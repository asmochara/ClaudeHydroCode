// ============================================================================
// simulation.cpp -- the coordinator: material-model construction, time-step
// control, the shared Coulomb logarithm, and the main run() loop that
// sequences the operator-split physics stages. The stages themselves live
// in per-physics files -- see the source-layout map in simulation.hpp.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>

// ============================================================================
// Construction: instantiate the material models requested by the deck
// (EOS, ionization, opacity) and build the initial mesh.
// ============================================================================
Simulation::Simulation(const InputDeck& deck) : deck_(deck) {
    twoTemperature_ = (deck_.control.temperatures == 2);
    radiationOn_ = deck_.radiation.enabled;
    laserOn_ = deck_.laser.enabled;
    if (laserOn_) {
        // Critical electron density, where the plasma frequency equals the
        // laser frequency and light can no longer propagate:
        //   n_crit = pi me c^2 / (e^2 lambda^2) = 1.11485e21 / lambda_um^2
        const double wavelength_um = deck_.laser.wavelength_um;
        criticalElectronDensity_percc =
            1.11485e21 / (wavelength_um * wavelength_um);
        buildRaySet();
    }

    // Instantiate materials. deck_.materials is a std::map, so iteration
    // order (and therefore material indexing) is deterministic: sorted by
    // material name.
    for (const auto& [name, spec] : deck_.materials) {
        Material material;
        material.name = name;
        material.A = spec.A;
        material.Z = spec.Z;
        material.fuel = spec.fuel;
        material.xD = spec.xD;
        material.xT = spec.xT;
        // Burn diagnostics run only if some material actually contains
        // deuterium or tritium (and the [burn] section didn't disable them).
        if (deck_.burn.enabled && (spec.xD > 0.0 || spec.xT > 0.0))
            burnOn_ = true;
        material.zbar =
            makeZbarModel(spec.ionization, spec.Z, spec.A, spec.zbar_table);
        if (twoTemperature_) {
            // 2T mode: separate ion and electron partial EOS.
            if (spec.eos == "ideal") {
                material.ion = std::make_shared<IdealIonEOS>(spec.gamma, spec.A);
                material.ele = std::make_shared<IdealElectronEOS>(
                    spec.gamma, spec.A, material.zbar, spec.degeneracy);
            } else {
                material.ion = std::make_shared<TableSpeciesEOS>(spec.table_ion);
                material.ele = std::make_shared<TableSpeciesEOS>(spec.table_electron);
            }
        } else {
            // 1T mode: one total-matter EOS.
            if (spec.eos == "ideal")
                material.eos = std::make_shared<IdealGasEOS>(spec.gamma, spec.A,
                                                             material.zbar);
            else
                material.eos = std::make_shared<TabulatedEOS>(spec.table_file);
        }
        if (radiationOn_) {
            if (!spec.opacity_table.empty())
                material.opacity = std::make_shared<TableOpacity>(spec.opacity_table);
            else
                material.opacity =
                    std::make_shared<ConstOpacity>(spec.kappa_R, spec.kappa_P);
        }
        materials_.push_back(std::move(material));
    }
    setupMesh();
}

// ============================================================================
// Coulomb logarithm (NRL Plasma Formulary, electron-ion). The two branches
// cover the classical and quantum-dominated impact-parameter regimes; the
// floor of 2 keeps cold/dense zones (where the formulary expressions go
// negative and lose meaning) at a sane strongly-coupled value. Shared by
// conduction, e-i coupling, and the laser absorption coefficient.
// ============================================================================
double Simulation::coulombLog(double electronDensity_percc, double temperature_eV,
                              double ionCharge) const {
    if (deck_.conduction.ln_lambda > 0.0) return deck_.conduction.ln_lambda;
    double lnLambda;
    if (temperature_eV > 10.0 * ionCharge * ionCharge)
        lnLambda = 24.0 - std::log(std::sqrt(electronDensity_percc) / temperature_eV);
    else
        lnLambda = 23.0 - std::log(std::sqrt(electronDensity_percc) * ionCharge *
                                   std::pow(temperature_eV, -1.5));
    return std::max(lnLambda, 2.0);
}

// ============================================================================
// Time-step control. Only the explicit hydro is stability-limited (the
// diffusion solves are implicit), so the constraint is the usual Courant
// condition on the sound crossing time of each zone, stiffened by the
// quadratic artificial viscosity's effective signal speed in compressing
// zones (the standard von Neumann-Richtmyer prescription). The step is also
// capped at dt_growth times the previous step so it recovers smoothly after
// a transient.
// ============================================================================
double Simulation::computeTimeStep() const {
    double dt_s = deck_.control.dt_max;
    for (int k = 0; k < numZones_; ++k) {
        const double zoneWidth_cm = nodeRadius_cm[k + 1] - nodeRadius_cm[k];
        // Compression speed: positive when the zone is being squeezed.
        const double compressionSpeed_cmps = std::max(
            0.0, -(nodeVelocity_cmps[k + 1] - nodeVelocity_cmps[k]));
        const double signalSpeed_cmps =
            zoneSoundSpeed_cmps[k] +
            4.0 * deck_.control.c_quad * compressionSpeed_cmps;
        dt_s = std::min(dt_s, deck_.control.cfl * zoneWidth_cm / signalSpeed_cmps);
    }
    if (timeStep_s > 0.0) dt_s = std::min(dt_s, timeStep_s * deck_.control.dt_growth);
    return dt_s;
}

// ============================================================================
// Main time loop. Each pass advances one time step through every enabled
// physics stage (see the class comment for why this order), then handles
// output. Only the hydro is explicit; everything else is implicit, so
// computeTimeStep() need only satisfy the acoustic CFL condition.
// ============================================================================
void Simulation::run() {
    std::filesystem::create_directories(deck_.output.directory);
    writeHistoryHeader();
    writeHistoryRow();
    int snapshotIndex = 0;
    writeSnapshot(snapshotIndex++);
    double nextDumpTime_s = (deck_.output.dt_dump > 0.0)
                                ? deck_.output.dt_dump
                                : 2.0 * deck_.control.t_end;  // i.e. never

    std::cout << "hydro1d: " << numZones_ << " zones, geometry d="
              << deck_.control.geometry << ", " << (twoTemperature_ ? "2T" : "1T")
              << ", conduction " << (deck_.conduction.enabled ? "on" : "off")
              << ", radiation " << (radiationOn_ ? "on" : "off") << ", laser "
              << (laserOn_ ? "on" : "off") << ", outer BC "
              << deck_.control.bc_outer << "\n";

    // Stages beyond the hydro leave (P, cs, Zbar) stale after they move
    // energy around; refresh once per step when any of them is active.
    const bool needsEndOfStepRefresh =
        deck_.conduction.enabled || twoTemperature_ || radiationOn_ || laserOn_;
    timeStep_s = deck_.control.dt_init;
    while (time_s < deck_.control.t_end && stepCount_ < deck_.control.max_steps) {
        timeStep_s = std::min(computeTimeStep(), deck_.control.t_end - time_s);
        if (stepCount_ == 0) timeStep_s = std::min(timeStep_s, deck_.control.dt_init);

        hydroStep(timeStep_s);
        laserStep(timeStep_s);
        couplingStep(timeStep_s);
        conductionStep(timeStep_s);
        radiationStep(timeStep_s);
        burnStep(timeStep_s);
        if (needsEndOfStepRefresh) {
            updateThermodynamics();
            applyTemperatureFloors();
        }

        time_s += timeStep_s;
        ++stepCount_;

        updateReport();
        if (stepCount_ % deck_.output.history_stride == 0) writeHistoryRow();
        if (time_s >= nextDumpTime_s - 1e-30) {
            writeSnapshot(snapshotIndex++);
            nextDumpTime_s += deck_.output.dt_dump;
        }
        if (stepCount_ % 20000 == 0)
            std::cout << "  step " << stepCount_ << "  t = " << time_s
                      << " s  dt = " << timeStep_s << " s\n";
    }

    writeHistoryRow();
    writeSnapshot(snapshotIndex);
    writeReport();
    std::cout << "hydro1d: done. t = " << time_s << " s in " << stepCount_
              << " steps. Output in " << deck_.output.directory << "/\n";
}
