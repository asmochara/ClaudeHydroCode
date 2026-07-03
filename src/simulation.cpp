// ============================================================================
// simulation.cpp -- the Lagrangian hydro engine and all operator-split
// physics stages. See simulation.hpp for the mesh layout, stage ordering,
// and the unit-suffix naming convention (_cm, _s, _gcc, _eV, _dyncm2, ...).
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"
#include "fusion.hpp"
#include "numerics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>

namespace {

// Braginskii electron-conduction coefficient gamma0 as a function of ion
// charge: kappa_e = gamma0 * ne kB Te tau_e / me. This rational fit passes
// through the tabulated Braginskii values (3.16 at Z=1) and approaches the
// Lorentz-gas limit at high Z.
double braginskiiGamma0(double ionCharge) {
    return 13.58 * (ionCharge + 0.24) / (ionCharge + 4.24);
}

// Levermore-Pomraning flux limiter for radiation diffusion. The argument
// R = |grad Er| / (kappa_R rho Er) measures how steep the radiation field is
// in units of the photon mean free path: lambda(R->0) = 1/3 recovers classic
// diffusion, lambda(R->inf) ~ 1/R caps the flux at the free-streaming value
// c*Er.
double levermorePomraningLambda(double R) {
    return (2.0 + R) / (6.0 + 3.0 * R + R * R);
}

// Solve the tridiagonal linear system
//   lower[k]*x[k-1] + diag[k]*x[k] + upper[k]*x[k+1] = rhs[k]
// by the Thomas algorithm (forward elimination + back substitution).
// diag and rhs are modified in place. O(n), no pivoting -- all our systems
// are diagonally dominant (backward-Euler diffusion matrices), so this is
// unconditionally safe here.
void thomasSolve(std::vector<double>& lower, std::vector<double>& diag,
                 std::vector<double>& upper, std::vector<double>& rhs,
                 std::vector<double>& solution) {
    const int n = static_cast<int>(diag.size());
    for (int k = 1; k < n; ++k) {
        const double w = lower[k] / diag[k - 1];
        diag[k] -= w * upper[k - 1];
        rhs[k] -= w * rhs[k - 1];
    }
    solution[n - 1] = rhs[n - 1] / diag[n - 1];
    for (int k = n - 2; k >= 0; --k)
        solution[k] = (rhs[k] - upper[k] * solution[k + 1]) / diag[k];
}

}  // namespace

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
// Geometry: face area and shell volume for planar (d=1), cylindrical (d=2),
// and spherical (d=3) symmetry. Planar quantities are per unit area,
// cylindrical per unit length, so "volume" has units cm, cm^2, cm^3
// respectively -- but every use pairs volume with density consistently, so
// zone masses are per-unit-area / per-unit-length masses in the reduced
// geometries and everything stays dimensionally coherent.
// ============================================================================
double Simulation::faceArea(double radius_cm) const {
    switch (deck_.control.geometry) {
        case 1: return 1.0;
        case 2: return 2.0 * phys::pi * radius_cm;
        default: return 4.0 * phys::pi * radius_cm * radius_cm;
    }
}

double Simulation::shellVolume(double innerRadius_cm, double outerRadius_cm) const {
    switch (deck_.control.geometry) {
        case 1:
            return outerRadius_cm - innerRadius_cm;
        case 2:
            return phys::pi * (outerRadius_cm * outerRadius_cm -
                               innerRadius_cm * innerRadius_cm);
        default:
            return 4.0 / 3.0 * phys::pi *
                   (outerRadius_cm * outerRadius_cm * outerRadius_cm -
                    innerRadius_cm * innerRadius_cm * innerRadius_cm);
    }
}

// ============================================================================
// Mesh setup: lay down the layers innermost-first, assign initial
// thermodynamic state, and freeze the Lagrangian zone masses.
// ============================================================================
void Simulation::setupMesh() {
    // Map a material name from a layer spec to its index in materials_.
    auto materialIndexByName = [&](const std::string& name) {
        for (size_t k = 0; k < materials_.size(); ++k)
            if (materials_[k].name == name) return static_cast<int>(k);
        throw std::runtime_error("unknown material " + name);
    };

    numZones_ = 0;
    for (const auto& layer : deck_.layers) numZones_ += layer.zones;
    nodeRadius_cm.assign(numZones_ + 1, 0.0);
    nodeVelocity_cmps.assign(numZones_ + 1, 0.0);
    nodeMass_g.assign(numZones_ + 1, 0.0);
    zoneMass_g.assign(numZones_, 0.0);
    zoneDensity_gcc.assign(numZones_, 0.0);
    zonePressure_dyncm2.assign(numZones_, 0.0);
    zoneSoundSpeed_cmps.assign(numZones_, 0.0);
    zoneViscousPressure_dyncm2.assign(numZones_, 0.0);
    zoneMeanIonization.assign(numZones_, 0.0);
    zoneMaterialIndex.assign(numZones_, 0);
    if (twoTemperature_) {
        ionSpecificEnergy_ergg.assign(numZones_, 0.0);
        electronSpecificEnergy_ergg.assign(numZones_, 0.0);
        ionTemperature_eV.assign(numZones_, 0.0);
        electronTemperature_eV.assign(numZones_, 0.0);
        ionPressure_dyncm2.assign(numZones_, 0.0);
        electronPressure_dyncm2.assign(numZones_, 0.0);
    } else {
        zoneSpecificEnergy_ergg.assign(numZones_, 0.0);
        zoneTemperature_eV.assign(numZones_, 0.0);
    }
    if (radiationOn_) radiationEnergyDensity_ergcc.assign(numZones_, 0.0);

    // ---- node positions ----------------------------------------------------
    // Within each layer the zone widths follow a geometric progression with
    // (outermost width)/(innermost width) = layer.ratio, so the user can
    // "feather" the mesh finer toward one side (e.g. finer zones at the
    // ablation surface). ratio = 1 gives uniform zoning.
    int zone = 0;
    nodeRadius_cm[0] = deck_.control.r_min;
    for (const auto& layer : deck_.layers) {
        const int layerZoneCount = layer.zones;
        // Common ratio g between adjacent zone widths, and the first width
        // w1 chosen so the widths sum exactly to the layer thickness.
        const double g = (layerZoneCount > 1)
                             ? std::pow(layer.ratio, 1.0 / (layerZoneCount - 1))
                             : 1.0;
        const double firstWidth_cm =
            (std::abs(g - 1.0) < 1e-12)
                ? layer.thickness / layerZoneCount
                : layer.thickness * (1.0 - g) / (1.0 - std::pow(g, layerZoneCount));
        double width_cm = firstWidth_cm;
        const int materialIndex = materialIndexByName(layer.material);
        const auto& material = materials_[materialIndex];

        // ---- initial temperature ---------------------------------------
        // The layer gives either T0 directly or P0, in which case we invert
        // the (total) pressure at the layer density for the temperature.
        double baseTemperature_eV = layer.T0;
        if (baseTemperature_eV <= 0.0 && layer.P0 > 0.0) {
            if (twoTemperature_) {
                baseTemperature_eV = invertMonotone(
                    [&](double trialT_eV) {
                        return material.ion->pressure(layer.rho0, trialT_eV) +
                               material.ele->pressure(layer.rho0, trialT_eV);
                    },
                    layer.P0, 1e-12, 1e9, 1.0);
            } else {
                baseTemperature_eV =
                    material.eos->temperatureFromPressure(layer.rho0, layer.P0);
            }
        }

        for (int k = 0; k < layerZoneCount; ++k, ++zone) {
            nodeRadius_cm[zone + 1] = nodeRadius_cm[zone] + width_cm;
            width_cm *= g;
            zoneMaterialIndex[zone] = materialIndex;
            zoneDensity_gcc[zone] = layer.rho0;
            const double floorT_eV = deck_.control.T_floor;
            if (twoTemperature_) {
                // Optional per-species overrides Ti0/Te0; otherwise both
                // species start at the common base temperature.
                ionTemperature_eV[zone] =
                    std::max((layer.Ti0 > 0.0) ? layer.Ti0 : baseTemperature_eV,
                             floorT_eV);
                electronTemperature_eV[zone] =
                    std::max((layer.Te0 > 0.0) ? layer.Te0 : baseTemperature_eV,
                             floorT_eV);
                ionSpecificEnergy_ergg[zone] =
                    material.ion->energy(zoneDensity_gcc[zone],
                                         ionTemperature_eV[zone]);
                electronSpecificEnergy_ergg[zone] =
                    material.ele->energy(zoneDensity_gcc[zone],
                                         electronTemperature_eV[zone]);
            } else {
                zoneTemperature_eV[zone] = std::max(baseTemperature_eV, floorT_eV);
                zoneSpecificEnergy_ergg[zone] =
                    material.eos->energy(zoneDensity_gcc[zone],
                                         zoneTemperature_eV[zone]);
            }
            if (radiationOn_) {
                // Radiation starts in equilibrium with the (electron)
                // temperature unless the layer sets Tr0 explicitly.
                const double radiationTemperature_eV =
                    (layer.Tr0 > 0.0)
                        ? layer.Tr0
                        : (twoTemperature_ ? electronTemperature_eV[zone]
                                           : zoneTemperature_eV[zone]);
                radiationEnergyDensity_ergcc[zone] =
                    phys::a_rad * radiationTemperature_eV * radiationTemperature_eV *
                    radiationTemperature_eV * radiationTemperature_eV;
            }
        }
    }
    // The geometric widths accumulate floating-point round-off; pin the
    // layer boundaries back to their exact positions so material interfaces
    // land where the deck says.
    {
        int boundaryNode = 0;
        double layerEdge_cm = deck_.control.r_min;
        for (const auto& layer : deck_.layers) {
            layerEdge_cm += layer.thickness;
            boundaryNode += layer.zones;
            nodeRadius_cm[boundaryNode] = layerEdge_cm;
        }
    }

    // ---- Lagrangian masses ---------------------------------------------------
    // Zone masses are fixed forever after this point. Node masses (used in
    // the momentum equation) are the half-masses of the two adjacent zones;
    // the boundary nodes carry half of their single neighboring zone.
    for (int k = 0; k < numZones_; ++k)
        zoneMass_g[k] = zoneDensity_gcc[k] *
                        shellVolume(nodeRadius_cm[k], nodeRadius_cm[k + 1]);
    for (int i = 1; i < numZones_; ++i)
        nodeMass_g[i] = 0.5 * (zoneMass_g[i - 1] + zoneMass_g[i]);
    nodeMass_g[0] = 0.5 * zoneMass_g[0];
    nodeMass_g[numZones_] = 0.5 * zoneMass_g[numZones_ - 1];

    // ---- shot-report bookkeeping ----------------------------------------------
    // Fuel zones drive the report's implosion metrics; the "hot spot" is the
    // innermost layer, and its outer boundary node (a fixed Lagrangian node,
    // so it tracks the gas/shell interface for all time) defines the
    // convergence ratio.
    zoneIsFuel_.assign(numZones_, 0);
    noFuelFlag_ = true;
    for (int k = 0; k < numZones_; ++k)
        if (materials_[zoneMaterialIndex[k]].fuel) {
            zoneIsFuel_[k] = 1;
            noFuelFlag_ = false;
        }
    if (noFuelFlag_)
        for (int k = 0; k < numZones_; ++k) zoneIsFuel_[k] = 1;  // fall back: all
    hotSpotNode_ = deck_.layers.front().zones;
    hotSpotRadius0_cm = nodeRadius_cm[hotSpotNode_];

    updateThermodynamics();

    // Total energy at t = 0 (the fluid starts at rest, so no kinetic term).
    // This anchors the energy-conservation diagnostic in the history file.
    initialTotalEnergy_erg = 0.0;
    for (int k = 0; k < numZones_; ++k) {
        initialTotalEnergy_erg +=
            zoneMass_g[k] * (twoTemperature_
                                 ? ionSpecificEnergy_ergg[k] +
                                       electronSpecificEnergy_ergg[k]
                                 : zoneSpecificEnergy_ergg[k]);
        if (radiationOn_)
            initialTotalEnergy_erg +=
                radiationEnergyDensity_ergcc[k] *
                shellVolume(nodeRadius_cm[k], nodeRadius_cm[k + 1]);
    }
}

// ============================================================================
// Thermodynamic refresh: invert temperature from (density, specific energy)
// -- energy is the conserved quantity the physics stages update -- then
// evaluate pressure, mean ionization, and the adiabatic sound speed used by
// the CFL condition and the artificial viscosity.
// ============================================================================
void Simulation::updateThermodynamics() {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        double soundSpeedSq_cm2s2;
        if (twoTemperature_) {
            // The previous temperature is passed as the inversion's initial
            // guess -- state changes little per step, so this converges in a
            // few iterations.
            ionTemperature_eV[k] = material.ion->temperature(
                zoneDensity_gcc[k], ionSpecificEnergy_ergg[k], ionTemperature_eV[k]);
            electronTemperature_eV[k] = material.ele->temperature(
                zoneDensity_gcc[k], electronSpecificEnergy_ergg[k],
                electronTemperature_eV[k]);
            ionPressure_dyncm2[k] =
                material.ion->pressure(zoneDensity_gcc[k], ionTemperature_eV[k]);
            electronPressure_dyncm2[k] =
                material.ele->pressure(zoneDensity_gcc[k], electronTemperature_eV[k]);
            zonePressure_dyncm2[k] = ionPressure_dyncm2[k] + electronPressure_dyncm2[k];
            // Mean ionization follows the electron temperature.
            zoneMeanIonization[k] =
                material.zbar->zbar(zoneDensity_gcc[k], electronTemperature_eV[k]);
            soundSpeedSq_cm2s2 =
                material.ion->cs2Contribution(zoneDensity_gcc[k],
                                              ionTemperature_eV[k]) +
                material.ele->cs2Contribution(zoneDensity_gcc[k],
                                              electronTemperature_eV[k]);
        } else {
            zoneTemperature_eV[k] = material.eos->temperature(
                zoneDensity_gcc[k], zoneSpecificEnergy_ergg[k], zoneTemperature_eV[k]);
            zonePressure_dyncm2[k] =
                material.eos->pressure(zoneDensity_gcc[k], zoneTemperature_eV[k]);
            zoneMeanIonization[k] =
                material.zbar->zbar(zoneDensity_gcc[k], zoneTemperature_eV[k]);
            soundSpeedSq_cm2s2 =
                material.eos->soundSpeed2(zoneDensity_gcc[k], zoneTemperature_eV[k]);
        }
        // Radiation adds gamma_rad * P_rad / rho = (4/3)(Er/3)/rho to the
        // sound speed (Eddington closure), which matters once aT^4 rivals
        // the matter pressure.
        if (radiationOn_)
            soundSpeedSq_cm2s2 += (4.0 / 9.0) * radiationEnergyDensity_ergcc[k] /
                                  zoneDensity_gcc[k];
        zoneSoundSpeed_cmps[k] = std::sqrt(soundSpeedSq_cm2s2);
    }
}

// ============================================================================
// Temperature floors: keep every temperature at or above control.T_floor so
// transport coefficients (which scale like powers of T) never see zero or
// negative temperatures. Flooring ADDS energy that did not come from any
// physical source; the total injected is accumulated so the energy-
// conservation diagnostic can account for it explicitly rather than
// misreporting it as a solver error.
// ============================================================================
void Simulation::applyTemperatureFloors() {
    const double floorT_eV = deck_.control.T_floor;
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        if (twoTemperature_) {
            if (ionTemperature_eV[k] < floorT_eV) {
                const double flooredEnergy_ergg =
                    material.ion->energy(zoneDensity_gcc[k], floorT_eV);
                floorEnergyInjected_erg +=
                    zoneMass_g[k] * (flooredEnergy_ergg - ionSpecificEnergy_ergg[k]);
                ionTemperature_eV[k] = floorT_eV;
                ionSpecificEnergy_ergg[k] = flooredEnergy_ergg;
            }
            if (electronTemperature_eV[k] < floorT_eV) {
                const double flooredEnergy_ergg =
                    material.ele->energy(zoneDensity_gcc[k], floorT_eV);
                floorEnergyInjected_erg +=
                    zoneMass_g[k] *
                    (flooredEnergy_ergg - electronSpecificEnergy_ergg[k]);
                electronTemperature_eV[k] = floorT_eV;
                electronSpecificEnergy_ergg[k] = flooredEnergy_ergg;
            }
        } else {
            if (zoneTemperature_eV[k] < floorT_eV) {
                const double flooredEnergy_ergg =
                    material.eos->energy(zoneDensity_gcc[k], floorT_eV);
                floorEnergyInjected_erg +=
                    zoneMass_g[k] * (flooredEnergy_ergg - zoneSpecificEnergy_ergg[k]);
                zoneTemperature_eV[k] = floorT_eV;
                zoneSpecificEnergy_ergg[k] = flooredEnergy_ergg;
            }
        }
        if (radiationOn_) {
            // A very small radiation floor (equilibrium with T_floor scaled
            // down by 1e-6) keeps Er positive without injecting meaningful
            // energy.
            const double minEr_ergcc =
                phys::a_rad * floorT_eV * floorT_eV * floorT_eV * floorT_eV * 1e-6;
            if (radiationEnergyDensity_ergcc[k] < minEr_ergcc) {
                floorEnergyInjected_erg +=
                    (minEr_ergcc - radiationEnergyDensity_ergcc[k]) *
                    shellVolume(nodeRadius_cm[k], nodeRadius_cm[k + 1]);
                radiationEnergyDensity_ergcc[k] = minEr_ergcc;
            }
        }
    }
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
// Hydro step (von Neumann-Richtmyer leapfrog):
//   1. accelerate nodes from the pressure + viscosity gradient  (u at n+1/2)
//   2. move nodes and recompute densities                       (r at n+1)
//   3. evaluate the artificial viscosity from the new velocities
//   4. update specific internal energies from PdV work, with a predictor-
//      corrector to time-center the pressure (2nd order on smooth flow)
// ============================================================================
void Simulation::hydroStep(double dt_s) {
    const auto& control = deck_.control;

    // ---- (1) momentum: u^{n+1/2} = u^{n-1/2} + dt * a^n ----------------------
    // The force on a node is its face area times the jump in total stress
    // (matter pressure + radiation pressure + artificial viscosity) across
    // it. Node 0 (the center, or the inner wall if r_min > 0) is held fixed.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 1; i < numZones_; ++i) {
        const double stressJump_dyncm2 =
            -(zonePressure_dyncm2[i] + radiationPressure_dyncm2(i) +
              zoneViscousPressure_dyncm2[i] - zonePressure_dyncm2[i - 1] -
              radiationPressure_dyncm2(i - 1) - zoneViscousPressure_dyncm2[i - 1]);
        nodeVelocity_cmps[i] +=
            dt_s * faceArea(nodeRadius_cm[i]) * stressJump_dyncm2 / nodeMass_g[i];
    }
    // Outer boundary: "wall" pins the node; "pressure" applies the drive
    // table; "free" is a vacuum boundary (used with the laser so the corona
    // can blow off and generate ablation pressure self-consistently).
    if (control.bc_outer == "pressure" || control.bc_outer == "free") {
        const double appliedPressure_dyncm2 =
            (control.bc_outer == "pressure") ? deck_.drive.pressure(time_s) : 0.0;
        nodeVelocity_cmps[numZones_] +=
            dt_s * faceArea(nodeRadius_cm[numZones_]) *
            (zonePressure_dyncm2[numZones_ - 1] +
             radiationPressure_dyncm2(numZones_ - 1) +
             zoneViscousPressure_dyncm2[numZones_ - 1] - appliedPressure_dyncm2) /
            nodeMass_g[numZones_];
        // Work done ON the system by the applied pressure (an inward-moving
        // boundary against an external pressure gains energy), tracked for
        // the energy-conservation diagnostic.
        if (appliedPressure_dyncm2 != 0.0)
            driveWorkDone_erg += -appliedPressure_dyncm2 *
                                 faceArea(nodeRadius_cm[numZones_]) *
                                 nodeVelocity_cmps[numZones_] * dt_s;
    }  // else wall: nodeVelocity_cmps[numZones_] stays 0

    // ---- (2) move nodes and recompute densities -------------------------------
    // Zone masses are constant (Lagrangian), so density follows from the new
    // volume alone. Node crossings ("mesh tangling") indicate the time step
    // was too large for the flow -- fatal, so fail loudly with advice.
    std::vector<double> oldDensity_gcc = zoneDensity_gcc;
    for (int i = 0; i <= numZones_; ++i)
        nodeRadius_cm[i] += dt_s * nodeVelocity_cmps[i];
    for (int k = 0; k < numZones_; ++k) {
        if (nodeRadius_cm[k + 1] <= nodeRadius_cm[k])
            throw std::runtime_error("mesh tangled at zone " + std::to_string(k) +
                                     ", t = " + std::to_string(time_s) +
                                     " s. Reduce cfl or increase viscosity.");
        zoneDensity_gcc[k] =
            zoneMass_g[k] / shellVolume(nodeRadius_cm[k], nodeRadius_cm[k + 1]);
    }

    // ---- (3) artificial viscosity at n+1/2 -------------------------------------
    // Combined quadratic (von Neumann-Richtmyer) + linear (Landshoff)
    // viscous pressure, active only in compression. It spreads shocks over a
    // few zones so the difference equations can integrate through them; the
    // linear term damps the residual post-shock ringing.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const double velocityJump_cmps =
            nodeVelocity_cmps[k + 1] - nodeVelocity_cmps[k];
        if (velocityJump_cmps < 0.0) {  // compressing
            const double midDensity_gcc =
                0.5 * (zoneDensity_gcc[k] + oldDensity_gcc[k]);
            zoneViscousPressure_dyncm2[k] =
                midDensity_gcc * (-velocityJump_cmps) *
                (control.c_quad * (-velocityJump_cmps) +
                 control.c_lin * zoneSoundSpeed_cmps[k]);
        } else {
            zoneViscousPressure_dyncm2[k] = 0.0;
        }
    }

    // ---- (4) internal energy: de = -(P + q) dV per unit mass -------------------
    // dV is the change in SPECIFIC volume (1/rho). Using the beginning-of-
    // step pressure alone would be only first-order accurate, so we take a
    // predictor step to estimate the end-of-step pressure and use the
    // average (time-centered pressure -> 2nd order on smooth flow). In 2T
    // mode each species does its own PdV work with its own partial pressure,
    // and the viscous (shock) heating goes entirely to the IONS -- physically,
    // a shock thermalizes the ion flow first and electrons heat later
    // through collisional coupling.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        const double specificVolumeChange_ccg =
            1.0 / zoneDensity_gcc[k] - 1.0 / oldDensity_gcc[k];
        if (twoTemperature_) {
            // Ions: predictor with old pressure (+ viscosity), then corrector
            // with the average of old and predicted pressures.
            const double ionEnergyPredicted_ergg =
                std::max(ionSpecificEnergy_ergg[k] -
                             (ionPressure_dyncm2[k] + zoneViscousPressure_dyncm2[k]) *
                                 specificVolumeChange_ccg,
                         0.0);
            const double ionTempPredicted_eV = material.ion->temperature(
                zoneDensity_gcc[k], ionEnergyPredicted_ergg, ionTemperature_eV[k]);
            const double ionPressurePredicted_dyncm2 = material.ion->pressure(
                zoneDensity_gcc[k], std::max(ionTempPredicted_eV, control.T_floor));
            ionSpecificEnergy_ergg[k] -=
                (0.5 * (ionPressure_dyncm2[k] + ionPressurePredicted_dyncm2) +
                 zoneViscousPressure_dyncm2[k]) *
                specificVolumeChange_ccg;

            // Electrons: same predictor-corrector, no viscous heating.
            const double eleEnergyPredicted_ergg =
                std::max(electronSpecificEnergy_ergg[k] -
                             electronPressure_dyncm2[k] * specificVolumeChange_ccg,
                         0.0);
            const double eleTempPredicted_eV = material.ele->temperature(
                zoneDensity_gcc[k], eleEnergyPredicted_ergg,
                electronTemperature_eV[k]);
            const double elePressurePredicted_dyncm2 = material.ele->pressure(
                zoneDensity_gcc[k], std::max(eleTempPredicted_eV, control.T_floor));
            electronSpecificEnergy_ergg[k] -=
                0.5 *
                (electronPressure_dyncm2[k] + elePressurePredicted_dyncm2) *
                specificVolumeChange_ccg;
        } else {
            const double energyPredicted_ergg =
                std::max(zoneSpecificEnergy_ergg[k] -
                             (zonePressure_dyncm2[k] + zoneViscousPressure_dyncm2[k]) *
                                 specificVolumeChange_ccg,
                         0.0);
            const double tempPredicted_eV = material.eos->temperature(
                zoneDensity_gcc[k], energyPredicted_ergg, zoneTemperature_eV[k]);
            const double pressurePredicted_dyncm2 = material.eos->pressure(
                zoneDensity_gcc[k], std::max(tempPredicted_eV, control.T_floor));
            zoneSpecificEnergy_ergg[k] -=
                (0.5 * (zonePressure_dyncm2[k] + pressurePredicted_dyncm2) +
                 zoneViscousPressure_dyncm2[k]) *
                specificVolumeChange_ccg;
        }
        // Radiation compresses adiabatically as a gamma = 4/3 gas:
        // Er V^{4/3} = const  =>  Er scales as rho^{4/3}. This is exact for
        // an isotropic photon gas and needs no predictor-corrector.
        if (radiationOn_)
            radiationEnergyDensity_ergcc[k] *=
                std::pow(zoneDensity_gcc[k] / oldDensity_gcc[k], 4.0 / 3.0);
    }

    updateThermodynamics();
    applyTemperatureFloors();
}

// ============================================================================
// Coulomb logarithm (NRL Plasma Formulary, electron-ion). The two branches
// cover the classical and quantum-dominated impact-parameter regimes; the
// floor of 2 keeps cold/dense zones (where the formulary expressions go
// negative and lose meaning) at a sane strongly-coupled value.
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
// Ray set: choose the impact parameters b_k so that every ray carries the
// same power P(t)/N. That holds when the b_k sit at the quantiles of the
// focal spot's cumulative power distribution
//   C(b) = integral_0^b 2 pi b' I(b') db'.
// Equal-power rays waste no resolution on dim parts of the spot and let the
// trace treat every ray identically.
// ============================================================================
void Simulation::buildRaySet() {
    const auto& laser = deck_.laser;
    const int rayCount = laser.rays;
    rayImpactParameter_cm.resize(rayCount);

    if (laser.profile == "flattop") {
        // Uniform intensity: C(b) ~ b^2, so the quantiles are closed-form.
        for (int k = 0; k < rayCount; ++k)
            rayImpactParameter_cm[k] =
                laser.beam_radius * std::sqrt((k + 0.5) / rayCount);
        return;
    }

    // Radial intensity profile I(b) in relative units, and the outer radius
    // beyond which it is treated as zero.
    double maxImpact_cm;
    std::function<double(double)> intensityProfile;
    if (laser.profile == "table") {
        const auto& table = laser.profile_table;
        maxImpact_cm = table.back().first;
        intensityProfile = [&table](double b_cm) {
            // Zero outside the tabulated range: a first radius > 0 therefore
            // produces an annular beam.
            if (b_cm < table.front().first || b_cm > table.back().first) return 0.0;
            return interpTimeTable(table, b_cm);
        };
    } else {  // gaussian | supergaussian, truncated where I/I0 = 1e-4
        const double order = (laser.profile == "gaussian") ? 2.0 : laser.sg_order;
        const double eFoldRadius_cm = laser.beam_radius;
        maxImpact_cm = eFoldRadius_cm * std::pow(std::log(1e4), 1.0 / order);
        intensityProfile = [order, eFoldRadius_cm](double b_cm) {
            return std::exp(-std::pow(b_cm / eFoldRadius_cm, order));
        };
    }

    // Build C(b) on a fine grid by the trapezoid rule, then invert each
    // quantile by linear interpolation.
    const int gridCount = 8192;
    std::vector<double> bGrid_cm(gridCount + 1), cumulativePower(gridCount + 1);
    for (int i = 0; i <= gridCount; ++i)
        bGrid_cm[i] = maxImpact_cm * i / gridCount;
    cumulativePower[0] = 0.0;
    for (int i = 1; i <= gridCount; ++i) {
        const double f0 = bGrid_cm[i - 1] * intensityProfile(bGrid_cm[i - 1]);
        const double f1 = bGrid_cm[i] * intensityProfile(bGrid_cm[i]);
        cumulativePower[i] =
            cumulativePower[i - 1] + 0.5 * (f0 + f1) * (bGrid_cm[i] - bGrid_cm[i - 1]);
    }
    if (cumulativePower[gridCount] <= 0.0)
        throw std::runtime_error("laser: focal-spot profile carries no power");
    for (int k = 0; k < rayCount; ++k) {
        const double target = (k + 0.5) / rayCount * cumulativePower[gridCount];
        const auto it =
            std::lower_bound(cumulativePower.begin(), cumulativePower.end(), target);
        const size_t i = std::max<size_t>(1, it - cumulativePower.begin());
        const double weight = (target - cumulativePower[i - 1]) /
                              std::max(cumulativePower[i] - cumulativePower[i - 1],
                                       1e-300);
        rayImpactParameter_cm[k] =
            bGrid_cm[i - 1] + weight * (bGrid_cm[i] - bGrid_cm[i - 1]);
    }
}

// ============================================================================
// Laser ray trace with refraction and inverse-bremsstrahlung absorption.
//
// GEOMETRY: uniform illumination from "infinitely many beams" reduces, in
// spherical symmetry, to a bundle of rays labeled by impact parameter b.
// Each ray obeys Bouguer's law (the spherical Snell's law)
//     mu(r) * r * sin(theta) = b        (b defined in vacuum where mu = 1)
// with refractive index mu = sqrt(1 - ne/n_crit). Within a zone of constant
// mu the ray is a straight chord whose distance of closest approach to the
// origin is d = b/mu. Walking inward zone by zone, a ray either:
//   - crosses the zone (d < inner radius): chord length from geometry;
//   - turns inside the zone (d >= inner radius): reaches depth d and comes
//     back out through the same zone;
//   - reflects at a zone face where mu drops enough that b/mu >= face radius
//     (total internal reflection), or where the zone is overdense
//     (ne >= n_crit, the critical surface).
// The outward path is the exact mirror of the inward path, so the ray
// retraces its chord list in reverse.
//
// ABSORPTION: along each chord the power decays as exp(-kappa_IB * length)
// with the inverse-bremsstrahlung coefficient
//     kappa_IB = nu_ei * (ne/n_crit) / (c * mu)          [1/cm]
// (nu_ei = NRL electron-ion collision frequency). Absorbed power is
// deposited in the traversed zone's electrons. A user-set fraction of any
// power reaching the critical surface is dumped there -- a stand-in for
// resonance absorption that also bootstraps coupling on a cold solid target
// before an underdense corona exists.
//
// LANGDON EFFECT: strong IB heating distorts the electron distribution away
// from Maxwellian (toward a super-Gaussian), which reduces the absorption.
// We apply the standard fit to Langdon's result (PRL 44, 575 (1980)):
//     kappa *= 1 - 0.553 / (1 + (0.27/alpha)^0.75),
//     alpha = Zbar * v_osc^2 / v_te^2,
// where v_osc is the electron quiver velocity in the laser field. The local
// intensity needed for v_osc is gathered from the trace itself (each chord
// contributes remaining ray power over its oblique tube cross-section) and
// used one step lagged -- the intensity field evolves far more slowly than
// the hydro time step, so the lag is harmless.
// ============================================================================
void Simulation::laserStep(double dt_s) {
    if (!laserOn_) return;
    const auto& laser = deck_.laser;
    const double totalPower_ergs = laser.powerAt(time_s);
    laserAbsorbedFraction_ = 0.0;
    rayDiag_.clear();
    if (totalPower_ergs <= 0.0) {
        laserIntensity_ergcm2s.assign(numZones_, 0.0);
        return;
    }
    laserIncident_erg += totalPower_ergs * dt_s;

    // ---- per-zone optics (refractive index, absorption coefficient) ----------
    if (laserIntensity_ergcm2s.empty()) laserIntensity_ergcm2s.assign(numZones_, 0.0);
    const double wavelength_cm = laser.wavelength_um * 1e-4;
    // Quiver-velocity coefficient: v_osc^2 = voscCoef * I  with I in
    // erg/cm^2/s. Derivation (linear polarization, time-averaged):
    //   I = (c/8pi) E^2,  v_osc = eE/(me omega),  omega = 2 pi c / lambda
    //   =>  v_osc^2 = 2 e^2 lambda^2 I / (pi me^2 c^3).
    const double voscCoef = 2.0 * phys::e_esu * phys::e_esu * wavelength_cm *
                            wavelength_cm /
                            (phys::pi * phys::m_e * phys::m_e * phys::c_light *
                             phys::c_light * phys::c_light);
    std::vector<double> refractiveIndex(numZones_), ibCoefficient_percm(numZones_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        const double electronTemp_eV =
            twoTemperature_ ? electronTemperature_eV[k] : zoneTemperature_eV[k];
        const double electronDensity_percc =
            zoneDensity_gcc[k] * zoneMeanIonization[k] / (material.A * phys::m_p);
        const double densityRatio = electronDensity_percc / criticalElectronDensity_percc;
        if (densityRatio < 1.0) {
            refractiveIndex[k] = std::sqrt(1.0 - densityRatio);
            const double effectiveCharge = std::max(zoneMeanIonization[k], 1.0);
            const double lnLambda =
                coulombLog(electronDensity_percc, electronTemp_eV, effectiveCharge);
            // NRL electron-ion collision frequency [1/s].
            const double collisionFrequency_pers =
                2.91e-6 * effectiveCharge * electronDensity_percc * lnLambda /
                std::pow(electronTemp_eV, 1.5);
            ibCoefficient_percm[k] = collisionFrequency_pers * densityRatio /
                                     (phys::c_light * refractiveIndex[k]);
            // Langdon reduction (see header comment), using the lagged local
            // intensity. Skipped when alpha is negligible.
            if (laser.langdon && laserIntensity_ergcm2s[k] > 0.0) {
                const double thermalSpeedSq_cm2s2 =
                    electronTemp_eV * phys::eV / phys::m_e;
                const double langdonAlpha = effectiveCharge * voscCoef *
                                            laserIntensity_ergcm2s[k] /
                                            thermalSpeedSq_cm2s2;
                if (langdonAlpha > 1e-4)
                    ibCoefficient_percm[k] *=
                        1.0 - 0.553 / (1.0 + std::pow(0.27 / langdonAlpha, 0.75));
            }
        } else {
            refractiveIndex[k] = 0.0;  // overdense: reflects at outer face
            ibCoefficient_percm[k] = 0.0;
        }
    }
    // Local intensity gathered during THIS trace, promoted to
    // laserIntensity_ergcm2s at the end for next step's Langdon factor.
    std::vector<double> gatheredIntensity_ergcm2s(numZones_, 0.0);

    std::vector<double> depositedPower_ergs(numZones_, 0.0);
    double absorbedPower_ergs = 0.0;
    const int rayCount = laser.rays;
    const double rayPower_ergs = totalPower_ergs / rayCount;

    // One chord of a ray's path: which zone it crosses and how long it is.
    struct Chord {
        int zone;
        double length_cm;
    };
    std::vector<Chord> inwardPath;
    for (int rayIdx = 0; rayIdx < rayCount; ++rayIdx) {
        const double impactParameter_cm = rayImpactParameter_cm[rayIdx];
        if (impactParameter_cm >= nodeRadius_cm[numZones_]) {
            // The ray misses the plasma entirely (spot larger than target).
            rayDiag_.push_back({impactParameter_cm, impactParameter_cm, 0.0});
            continue;
        }

        // ---- inward walk from the outer boundary ----------------------------
        inwardPath.clear();
        bool turnedInsideZone = false;
        int criticalSurfaceZone = -1;  // zone whose outer face reflected us
        double turningRadius_cm = nodeRadius_cm[0];
        for (int j = numZones_ - 1; j >= 0; --j) {
            const double innerRadius_cm = nodeRadius_cm[j];
            const double outerRadius_cm = nodeRadius_cm[j + 1];
            if (refractiveIndex[j] <= 0.0) {
                // Overdense zone: the ray reflects at its outer face (the
                // discrete critical surface).
                criticalSurfaceZone = j;
                turningRadius_cm = outerRadius_cm;
                break;
            }
            // Bouguer: the straight chord in this constant-mu zone passes
            // the origin at distance d = b/mu.
            const double closestApproach_cm =
                impactParameter_cm / refractiveIndex[j];
            if (closestApproach_cm >= outerRadius_cm) {
                // sin(theta) would exceed 1 at the face: total internal
                // reflection where the index drops across the interface.
                turningRadius_cm = outerRadius_cm;
                break;
            }
            // Half-chord length from the outer face to closest approach.
            const double halfChordAtOuter_cm =
                std::sqrt(outerRadius_cm * outerRadius_cm -
                          closestApproach_cm * closestApproach_cm);
            if (closestApproach_cm >= innerRadius_cm) {
                // The ray turns inside this zone: it travels to depth d and
                // back out, so the full in-zone path is twice the half-chord.
                inwardPath.push_back({j, 2.0 * halfChordAtOuter_cm});
                turningRadius_cm = closestApproach_cm;
                turnedInsideZone = true;
                break;
            }
            // The ray crosses the zone: chord length is the difference of
            // half-chords at the two faces.
            inwardPath.push_back(
                {j, halfChordAtOuter_cm -
                        std::sqrt(innerRadius_cm * innerRadius_cm -
                                  closestApproach_cm * closestApproach_cm)});
            // If j reaches 0 without turning (possible only when r_min > 0),
            // the ray hits the inner wall and reflects;
            // turningRadius_cm = nodeRadius_cm[0] was preset above.
        }

        // ---- attenuate along the path ---------------------------------------
        // Inward chords, then the optional critical-surface dump, then the
        // mirrored outward chords (the turning chord, if any, already covers
        // both directions so it is excluded from the reversed pass). Each
        // traversal also gathers local intensity: remaining ray power over
        // the ray tube's oblique cross-section area*cos(theta), with
        // cos(theta) floored at 0.1 near turning points where the true
        // intensity swelling is bounded by diffraction, not geometry.
        double remainingPower_ergs = rayPower_ergs;
        auto traverseChord = [&](const Chord& chord) {
            const double zoneCenter_cm =
                0.5 * (nodeRadius_cm[chord.zone] + nodeRadius_cm[chord.zone + 1]);
            const double sinTheta =
                (refractiveIndex[chord.zone] > 0.0)
                    ? impactParameter_cm /
                          (refractiveIndex[chord.zone] * zoneCenter_cm)
                    : 1.0;
            const double cosTheta = std::max(
                std::sqrt(std::max(1.0 - sinTheta * sinTheta, 0.0)), 0.1);
            gatheredIntensity_ergcm2s[chord.zone] +=
                remainingPower_ergs / (faceArea(zoneCenter_cm) * cosTheta);
            // expm1 keeps precision for optically thin chords (small tau).
            const double absorbed_ergs =
                remainingPower_ergs *
                (-std::expm1(-ibCoefficient_percm[chord.zone] * chord.length_cm));
            depositedPower_ergs[chord.zone] += absorbed_ergs;
            remainingPower_ergs -= absorbed_ergs;
        };
        for (const auto& chord : inwardPath) traverseChord(chord);
        if (criticalSurfaceZone >= 0 && laser.absorb_at_critical > 0.0) {
            const double dumped_ergs = remainingPower_ergs * laser.absorb_at_critical;
            depositedPower_ergs[criticalSurfaceZone] += dumped_ergs;
            remainingPower_ergs -= dumped_ergs;
        }
        const int outwardChordCount =
            static_cast<int>(inwardPath.size()) - (turnedInsideZone ? 1 : 0);
        for (int c = outwardChordCount - 1; c >= 0; --c) traverseChord(inwardPath[c]);

        // Whatever survives the round trip escapes back out of the plasma.
        absorbedPower_ergs += rayPower_ergs - remainingPower_ergs;
        rayDiag_.push_back({impactParameter_cm, turningRadius_cm,
                            (rayPower_ergs - remainingPower_ergs) / rayPower_ergs});
    }

    laserAbsorbedFraction_ = absorbedPower_ergs / totalPower_ergs;
    laserIntensity_ergcm2s = gatheredIntensity_ergcm2s;
    if (dt_s <= 0.0) return;  // trace-only mode (rayTraceReport)
    laserAbsorbed_erg += absorbedPower_ergs * dt_s;

    // ---- deposit into the electron (or 1T total-matter) energy ----------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        if (depositedPower_ergs[k] <= 0.0) continue;
        const auto& material = materials_[zoneMaterialIndex[k]];
        const double specificEnergyAdded_ergg =
            depositedPower_ergs[k] * dt_s / zoneMass_g[k];
        if (twoTemperature_) {
            electronSpecificEnergy_ergg[k] += specificEnergyAdded_ergg;
            electronTemperature_eV[k] = material.ele->temperature(
                zoneDensity_gcc[k], electronSpecificEnergy_ergg[k],
                electronTemperature_eV[k]);
        } else {
            zoneSpecificEnergy_ergg[k] += specificEnergyAdded_ergg;
            zoneTemperature_eV[k] = material.eos->temperature(
                zoneDensity_gcc[k], zoneSpecificEnergy_ergg[k],
                zoneTemperature_eV[k]);
        }
    }
}

void Simulation::rayTraceReport() {
    laserStep(0.0);  // dt = 0: trace and record diagnostics, deposit nothing
    std::printf("# ray trace at t = %g s, P = %g erg/s, n_crit = %g cm^-3\n",
                time_s, deck_.laser.powerAt(time_s), criticalElectronDensity_percc);
    std::printf("# %12s %14s %14s\n", "b [cm]", "r_turn [cm]", "f_abs");
    for (const auto& ray : rayDiag_)
        std::printf("  %12.6e %14.6e %14.6e\n", ray.impactParameter_cm,
                    ray.turningRadius_cm, ray.absorbedFraction);
    std::printf("# total absorbed fraction = %.6f\n", laserAbsorbedFraction_);
}

// ============================================================================
// Burn-off fusion diagnostics. Reaction rates use the Bosch-Hale Maxwellian
// reactivities evaluated at the ION temperature:
//   DT:   rate = nD * nT * <sv>_DT            [reactions / cm^3 / s]
//   DD:   rate = (1/2) nD^2 * <sv>_branch     (1/2 avoids double-counting
//                                              identical reactant pairs)
// Yields, fusion energy, and the burn-weighted <Ti> are accumulated, and the
// time of peak fusion power defines the burn bang time. NOTHING IS FED BACK:
// no charged-particle heating, no reactant depletion -- so running with burn
// on or off gives bit-identical hydrodynamics.
// ============================================================================
void Simulation::burnStep(double dt_s) {
    if (!burnOn_) return;
    double stepFusionPower_ergs = 0.0, stepYieldDT = 0.0, stepYieldDDn = 0.0;
    double stepWeightedTi_keV = 0.0, stepWeight = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
    reduction(+ : stepFusionPower_ergs, stepYieldDT, stepYieldDDn, \
                  stepWeightedTi_keV, stepWeight)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        if (material.xD <= 0.0) continue;
        const double ionTemp_keV =
            (twoTemperature_ ? ionTemperature_eV[k] : zoneTemperature_eV[k]) * 1e-3;
        // Below 0.2 keV the reactivity is beyond the fit's range and utterly
        // negligible (exponentially suppressed) -- skip the pow/exp work.
        if (ionTemp_keV < 0.2) continue;
        const double ionDensity_percc =
            zoneDensity_gcc[k] / (material.A * phys::m_p);
        const double deuteronDensity_percc = material.xD * ionDensity_percc;
        const double tritonDensity_percc = material.xT * ionDensity_percc;
        const double zoneVolume_cc =
            shellVolume(nodeRadius_cm[k], nodeRadius_cm[k + 1]);
        const double rateDT_perccps =
            deuteronDensity_percc * tritonDensity_percc * sigmavDT(ionTemp_keV);
        const double rateDDn_perccps = 0.5 * deuteronDensity_percc *
                                       deuteronDensity_percc * sigmavDDn(ionTemp_keV);
        const double rateDDp_perccps = 0.5 * deuteronDensity_percc *
                                       deuteronDensity_percc * sigmavDDp(ionTemp_keV);
        stepFusionPower_ergs +=
            (rateDT_perccps * fusion::Q_DT + rateDDn_perccps * fusion::Q_DDn +
             rateDDp_perccps * fusion::Q_DDp) *
            zoneVolume_cc;
        stepYieldDT += rateDT_perccps * zoneVolume_cc * dt_s;
        stepYieldDDn += rateDDn_perccps * zoneVolume_cc * dt_s;
        // Burn-weighted ion temperature: weight each zone's Ti by its DT
        // reactions this step (that is what an activation diagnostic sees).
        const double weight = rateDT_perccps * zoneVolume_cc * dt_s;
        stepWeightedTi_keV += weight * ionTemp_keV;
        stepWeight += weight;
    }
    fusionPower_ergs = stepFusionPower_ergs;
    neutronYieldDT += stepYieldDT;
    neutronYieldDDn += stepYieldDDn;
    fusionEnergy_erg += stepFusionPower_ergs * dt_s;
    burnWeightedTiSum_keV += stepWeightedTi_keV;
    burnWeightSum += stepWeight;
    if (stepFusionPower_ergs > peakFusionPower_ergs) {
        peakFusionPower_ergs = stepFusionPower_ergs;
        fusionBangTime_s = time_s;
    }
}

// ============================================================================
// Electron-ion temperature relaxation (2T only). The NRL-formulary
// equilibration rate
//   dTe/dt = nu_eq (Ti - Te),
//   nu_eq = 1.8e-19 sqrt(me mi) Z^2 ni lnLambda / (me Ti + mi Te)^{3/2}
// (masses in g, T in eV, n in cm^-3, nu in 1/s) is integrated pointwise with
// backward Euler: solving the coupled two-temperature relaxation implicitly
// gives the exact bounded update
//   (Ti - Te)^{new} = (Ti - Te) / (1 + q (1/cv_i + 1/cv_e)),
// which is unconditionally stable no matter how stiff the coupling (e.g. in
// cold dense matter where nu_eq is enormous). The energy exchanged is applied
// antisymmetrically, so the pair conserves energy to machine precision.
// ============================================================================
void Simulation::couplingStep(double dt_s) {
    if (!twoTemperature_) return;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        const double ionMass_g = material.A * phys::m_p;
        const double ionDensity_percc = zoneDensity_gcc[k] / ionMass_g;
        const double electronDensity_percc = ionDensity_percc * zoneMeanIonization[k];
        const double effectiveCharge = std::max(zoneMeanIonization[k], 1.0);
        const double lnLambda = coulombLog(electronDensity_percc,
                                           electronTemperature_eV[k], effectiveCharge);
        const double equilibrationRate_pers =
            1.8e-19 * std::sqrt(phys::m_e * ionMass_g) * zoneMeanIonization[k] *
            zoneMeanIonization[k] * ionDensity_percc * lnLambda /
            std::pow(phys::m_e * ionTemperature_eV[k] +
                         ionMass_g * electronTemperature_eV[k],
                     1.5);
        // q below is the coupling strength integrated over the step, per
        // unit mass: [erg/g/eV]. (3/2) ne kB per unit volume is the electron
        // heat capacity the formulary rate is defined against.
        const double couplingStrength_erggeV = dt_s * 1.5 * electronDensity_percc *
                                               phys::eV * equilibrationRate_pers /
                                               zoneDensity_gcc[k];
        const double ionHeatCapacity_erggeV =
            std::max(material.ion->cv(zoneDensity_gcc[k], ionTemperature_eV[k]),
                     1e-30);
        const double electronHeatCapacity_erggeV =
            std::max(material.ele->cv(zoneDensity_gcc[k], electronTemperature_eV[k]),
                     1e-30);
        const double newTemperatureGap_eV =
            (ionTemperature_eV[k] - electronTemperature_eV[k]) /
            (1.0 + couplingStrength_erggeV * (1.0 / ionHeatCapacity_erggeV +
                                              1.0 / electronHeatCapacity_erggeV));
        // Specific energy moved from ions to electrons this step.
        const double exchangedEnergy_ergg =
            couplingStrength_erggeV * newTemperatureGap_eV;
        electronSpecificEnergy_ergg[k] += exchangedEnergy_ergg;
        ionSpecificEnergy_ergg[k] -= exchangedEnergy_ergg;
        electronTemperature_eV[k] = material.ele->temperature(
            zoneDensity_gcc[k], electronSpecificEnergy_ergg[k],
            electronTemperature_eV[k]);
        ionTemperature_eV[k] = material.ion->temperature(
            zoneDensity_gcc[k], ionSpecificEnergy_ergg[k], ionTemperature_eV[k]);
    }
}

// ============================================================================
// Thermal conductivities. Both are expressed in units such that
//   heat flux = -kappa * d(kB*T)/dr    [erg/cm^2/s],
// i.e. kappa itself carries 1/(cm*s).
// ============================================================================
double Simulation::electronConductivity(int zone) const {
    // Spitzer-Harm: kappa_e = gamma0(Z) * ne * kB*Te * tau_e / me with the
    // NRL electron collision time tau_e = 3.44e5 Te^{3/2} / (ne lnLambda) [s]
    // (Te in eV, ne in cm^-3).
    const auto& material = materials_[zoneMaterialIndex[zone]];
    const double electronTemp_eV =
        twoTemperature_ ? electronTemperature_eV[zone] : zoneTemperature_eV[zone];
    const double electronDensity_percc =
        zoneDensity_gcc[zone] * zoneMeanIonization[zone] / (material.A * phys::m_p);
    const double effectiveCharge = std::max(zoneMeanIonization[zone], 1.0);
    const double lnLambda =
        coulombLog(electronDensity_percc, electronTemp_eV, effectiveCharge);
    const double collisionTime_s =
        3.44e5 * std::pow(electronTemp_eV, 1.5) / (electronDensity_percc * lnLambda);
    return braginskiiGamma0(effectiveCharge) * electronDensity_percc *
           (electronTemp_eV * phys::eV) * collisionTime_s / phys::m_e;
}

double Simulation::ionConductivity(int zone) const {
    // Braginskii: kappa_i = 3.9 * ni * kB*Ti * tau_i / mi with the NRL ion
    // collision time tau_i = 2.09e7 sqrt(A) Ti^{3/2} / (Z^4 ni lnLambda) [s].
    // Small next to the electron conductivity in equilibrium (sqrt(me/mi)),
    // but essential where Ti >> Te -- e.g. smoothing the converging-shock
    // ion-temperature spike at void closure.
    const auto& material = materials_[zoneMaterialIndex[zone]];
    const double ionMass_g = material.A * phys::m_p;
    const double ionDensity_percc = zoneDensity_gcc[zone] / ionMass_g;
    const double electronDensity_percc = ionDensity_percc * zoneMeanIonization[zone];
    const double effectiveCharge = std::max(zoneMeanIonization[zone], 1.0);
    const double lnLambda =
        coulombLog(electronDensity_percc, ionTemperature_eV[zone], effectiveCharge);
    const double collisionTime_s =
        2.09e7 * std::sqrt(material.A) * std::pow(ionTemperature_eV[zone], 1.5) /
        (effectiveCharge * effectiveCharge * effectiveCharge * effectiveCharge *
         ionDensity_percc * lnLambda);
    return 3.9 * ionDensity_percc * (ionTemperature_eV[zone] * phys::eV) *
           collisionTime_s / ionMass_g;
}

void Simulation::conductionStep(double dt_s) {
    if (!deck_.conduction.enabled || numZones_ < 2) return;
    solveConduction(dt_s, /*ionSpecies=*/false);  // electrons (or 1T matter)
    if (twoTemperature_ && deck_.conduction.ion_conduction)
        solveConduction(dt_s, /*ionSpecies=*/true);
}

// ============================================================================
// Flux-limited thermal conduction on one temperature field, backward Euler
// in time (so arbitrarily large diffusion numbers are stable -- Spitzer
// conduction in a keV hot spot is far too stiff for explicit stepping).
//
// Spatial discretization: finite-volume on the zones. The conductance of the
// face at node i between zones i-1 and i is
//   G_i = kappa_face * A_i * kB / (rc_i - rc_{i-1})   [erg/s per eV of dT]
// with kappa_face the HARMONIC mean of the zone conductivities (the correct
// series-resistance average, which keeps fluxes sane across sharp material
// interfaces where kappa jumps by orders of magnitude).
//
// Flux limiting: Spitzer-Harm is a small-gradient expansion; where the
// temperature scale length approaches the collision mean free path it wildly
// overpredicts the flux. The standard sharp limiter caps the flux against a
// user-set fraction f of the free-streaming value:
//   q = q_SH / (1 + |q_SH| / (f * n * kB*T * v_thermal)),
// applied by scaling the face conductivity, evaluated with beginning-of-step
// temperatures (frozen-coefficient linearization).
//
// Energy update: the SOLVED temperature field defines the face fluxes, and
// zone energies change by the flux divergence. Because each face's flux
// enters its two neighbors with opposite signs, total energy is conserved to
// round-off even when e(T) is nonlinear (TF ionization, degeneracy) -- which
// a naive per-zone e(T_new) - e(T_old) update would not guarantee.
// Temperatures are then re-inverted from the updated energies.
// ============================================================================
void Simulation::solveConduction(double dt_s, bool ionSpecies) {
    const double fluxLimiter = ionSpecies ? deck_.conduction.ion_flux_limiter
                                          : deck_.conduction.flux_limiter;
    std::vector<double>& temperature_eV =
        ionSpecies ? ionTemperature_eV
                   : (twoTemperature_ ? electronTemperature_eV : zoneTemperature_eV);
    std::vector<double>& specificEnergy_ergg =
        ionSpecies ? ionSpecificEnergy_ergg
                   : (twoTemperature_ ? electronSpecificEnergy_ergg
                                      : zoneSpecificEnergy_ergg);

    // ---- per-zone coefficients -------------------------------------------------
    std::vector<double> zoneCenter_cm(numZones_), heatCapacity_erggeV(numZones_);
    std::vector<double> conductivity_percms(numZones_);
    std::vector<double> ionFreeStreamFlux_ergcm2s(numZones_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        zoneCenter_cm[k] = 0.5 * (nodeRadius_cm[k] + nodeRadius_cm[k + 1]);
        heatCapacity_erggeV[k] = std::max(
            ionSpecies
                ? material.ion->cv(zoneDensity_gcc[k], temperature_eV[k])
                : (twoTemperature_
                       ? material.ele->cv(zoneDensity_gcc[k], temperature_eV[k])
                       : material.eos->cv(zoneDensity_gcc[k], temperature_eV[k])),
            1e-30);
        conductivity_percms[k] =
            ionSpecies ? ionConductivity(k) : electronConductivity(k);
        if (ionSpecies) {
            // Ion free-streaming flux is evaluated per zone because the ion
            // mass differs between materials (unlike electrons).
            const double ionMass_g = material.A * phys::m_p;
            const double ionDensity_percc = zoneDensity_gcc[k] / ionMass_g;
            const double thermalEnergy_erg = temperature_eV[k] * phys::eV;
            ionFreeStreamFlux_ergcm2s[k] =
                fluxLimiter * ionDensity_percc * thermalEnergy_erg *
                std::sqrt(thermalEnergy_erg / ionMass_g);
        }
    }

    // ---- face conductances G_i [erg/s/eV] at interior nodes i = 1..M-1 ---------
    std::vector<double> faceConductance_ergseV(numZones_ + 1, 0.0);
    for (int i = 1; i < numZones_; ++i) {
        const int zoneLeft = i - 1, zoneRight = i;
        double harmonicKappa_percms = 0.0;
        if (conductivity_percms[zoneLeft] > 0.0 && conductivity_percms[zoneRight] > 0.0)
            harmonicKappa_percms =
                2.0 * conductivity_percms[zoneLeft] * conductivity_percms[zoneRight] /
                (conductivity_percms[zoneLeft] + conductivity_percms[zoneRight]);
        const double centerSpacing_cm =
            zoneCenter_cm[zoneRight] - zoneCenter_cm[zoneLeft];
        // Temperature gradient in energy units, d(kB*T)/dr [erg/cm].
        const double thermalGradient_ergcm =
            (temperature_eV[zoneRight] - temperature_eV[zoneLeft]) * phys::eV /
            centerSpacing_cm;
        double freeStreamFlux_ergcm2s;
        if (ionSpecies) {
            freeStreamFlux_ergcm2s = 0.5 * (ionFreeStreamFlux_ergcm2s[zoneLeft] +
                                            ionFreeStreamFlux_ergcm2s[zoneRight]);
        } else {
            const auto& materialLeft = materials_[zoneMaterialIndex[zoneLeft]];
            const auto& materialRight = materials_[zoneMaterialIndex[zoneRight]];
            const double faceElectronDensity_percc =
                0.5 * (zoneDensity_gcc[zoneLeft] * zoneMeanIonization[zoneLeft] /
                           (materialLeft.A * phys::m_p) +
                       zoneDensity_gcc[zoneRight] * zoneMeanIonization[zoneRight] /
                           (materialRight.A * phys::m_p));
            const double faceThermalEnergy_erg =
                0.5 * (temperature_eV[zoneLeft] + temperature_eV[zoneRight]) *
                phys::eV;
            freeStreamFlux_ergcm2s = fluxLimiter * faceElectronDensity_percc *
                                     faceThermalEnergy_erg *
                                     std::sqrt(faceThermalEnergy_erg / phys::m_e);
        }
        const double spitzerFlux_ergcm2s =
            harmonicKappa_percms * std::abs(thermalGradient_ergcm);
        const double limitedKappa_percms =
            (freeStreamFlux_ergcm2s > 0.0)
                ? harmonicKappa_percms / (1.0 + spitzerFlux_ergcm2s /
                                                    freeStreamFlux_ergcm2s)
                : harmonicKappa_percms;
        faceConductance_ergseV[i] = limitedKappa_percms *
                                    faceArea(nodeRadius_cm[i]) * phys::eV /
                                    centerSpacing_cm;
    }
    // faceConductance_ergseV[0] and [numZones_] stay 0: insulated boundaries.

    // ---- backward-Euler tridiagonal system for T^{n+1} --------------------------
    //   (m cv/dt + G_i + G_{i+1}) T_k - G_i T_{k-1} - G_{i+1} T_{k+1}
    //       = (m cv/dt) T_k^n
    std::vector<double> lower(numZones_), diag(numZones_), upper(numZones_),
        rhs(numZones_), solvedTemperature_eV(numZones_);
    for (int k = 0; k < numZones_; ++k) {
        const double thermalInertia_ergseV =
            zoneMass_g[k] * heatCapacity_erggeV[k] / dt_s;
        lower[k] = -faceConductance_ergseV[k];
        upper[k] = -faceConductance_ergseV[k + 1];
        diag[k] = thermalInertia_ergseV + faceConductance_ergseV[k] +
                  faceConductance_ergseV[k + 1];
        rhs[k] = thermalInertia_ergseV * temperature_eV[k];
    }
    thomasSolve(lower, diag, upper, rhs, solvedTemperature_eV);

    // ---- conservative (flux-form) energy update ---------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const double netHeatInflow_ergs =
            (k > 0 ? faceConductance_ergseV[k] *
                         (solvedTemperature_eV[k - 1] - solvedTemperature_eV[k])
                   : 0.0) +
            (k < numZones_ - 1
                 ? faceConductance_ergseV[k + 1] *
                       (solvedTemperature_eV[k + 1] - solvedTemperature_eV[k])
                 : 0.0);
        specificEnergy_ergg[k] += dt_s * netHeatInflow_ergs / zoneMass_g[k];
        temperature_eV[k] =
            ionSpecies
                ? materials_[zoneMaterialIndex[k]].ion->temperature(
                      zoneDensity_gcc[k], specificEnergy_ergg[k],
                      solvedTemperature_eV[k])
                : (twoTemperature_
                       ? materials_[zoneMaterialIndex[k]].ele->temperature(
                             zoneDensity_gcc[k], specificEnergy_ergg[k],
                             solvedTemperature_eV[k])
                       : materials_[zoneMaterialIndex[k]].eos->temperature(
                             zoneDensity_gcc[k], specificEnergy_ergg[k],
                             solvedTemperature_eV[k]));
    }
}

// Heat capacity of the matter field the radiation exchanges energy with:
// electrons in 2T mode, the whole matter in 1T mode. [erg/g/eV]
double Simulation::matterHeatCapacity(int zone) const {
    const auto& material = materials_[zoneMaterialIndex[zone]];
    return std::max(
        twoTemperature_
            ? material.ele->cv(zoneDensity_gcc[zone], electronTemperature_eV[zone])
            : material.eos->cv(zoneDensity_gcc[zone], zoneTemperature_eV[zone]),
        1e-30);
}

// ============================================================================
// Grey flux-limited radiation diffusion, operator-split and backward Euler.
//
// The radiation energy density Er [erg/cm^3] obeys
//   dEr/dt = div( D grad Er ) + c kappa_P rho (a T^4 - Er)
// with diffusion coefficient D = c lambda(R) / (kappa_R rho), where lambda
// is the Levermore-Pomraning flux limiter (-> 1/3 in thick material, caps
// the flux at c Er in transparent material).
//
// MATTER COUPLING: the emission term a T^4 is linearized about the current
// matter temperature (T^4 -> T0^4 + 4 T0^3 dT) and the resulting implicit
// matter response is folded into the Er equation as the factor
//   fc = cv / (cv + 4 c kappa_P dt a T^3),
// so a single tridiagonal solve in Er captures the stiff emission/absorption
// exchange stably. The matter energy is then updated with EXACTLY the
// linearized source used in the solve, so matter+radiation energy is
// conserved to round-off.
//
// BOUNDARIES: zero-flux at both ends unless radiation.bc_outer = "vacuum",
// which adds a Marshak leak F = chi Er at the outer face, with
//   chi = (c/2) / (1 + 0.75 kappa_R rho dr)
// interpolating between the transparent (c/2, isotropic escape) and
// optically thick (diffusion-limited) limits. Leaked energy is tracked.
// ============================================================================
void Simulation::radiationStep(double dt_s) {
    if (!radiationOn_) return;
    // The matter temperature/energy the radiation couples to.
    std::vector<double>& matterTemperature_eV =
        twoTemperature_ ? electronTemperature_eV : zoneTemperature_eV;
    std::vector<double>& matterEnergy_ergg =
        twoTemperature_ ? electronSpecificEnergy_ergg : zoneSpecificEnergy_ergg;

    // ---- per-zone coefficients --------------------------------------------------
    std::vector<double> zoneVolume_cc(numZones_), zoneCenter_cm(numZones_);
    std::vector<double> rosseland_cm2g(numZones_), planck_cm2g(numZones_);
    std::vector<double> couplingFactor(numZones_), heatCapacity_erggeV(numZones_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        zoneVolume_cc[k] = shellVolume(nodeRadius_cm[k], nodeRadius_cm[k + 1]);
        zoneCenter_cm[k] = 0.5 * (nodeRadius_cm[k] + nodeRadius_cm[k + 1]);
        // The Rosseland floor keeps the diffusion coefficient finite in
        // effectively transparent zones (the flux limiter does the real work
        // there anyway).
        rosseland_cm2g[k] = std::max(
            material.opacity->rosseland(zoneDensity_gcc[k], matterTemperature_eV[k]),
            1e-10);
        planck_cm2g[k] = std::max(
            material.opacity->planck(zoneDensity_gcc[k], matterTemperature_eV[k]),
            0.0);
        heatCapacity_erggeV[k] = matterHeatCapacity(k);
        // Linearized-emission reduction factor fc (see header comment):
        // beta is the matter temperature's stiffness against radiating,
        // 4 c kappa_P dt a T^3, in the same units as cv.
        const double beta_erggeV =
            4.0 * phys::c_light * planck_cm2g[k] * dt_s * phys::a_rad *
            matterTemperature_eV[k] * matterTemperature_eV[k] *
            matterTemperature_eV[k];
        couplingFactor[k] =
            heatCapacity_erggeV[k] / (heatCapacity_erggeV[k] + beta_erggeV);
    }

    // ---- face diffusion conductances [cm^3/s] -----------------------------------
    // G_i converts an Er difference [erg/cm^3] into an energy flow [erg/s].
    std::vector<double> faceConductance_ccs(numZones_ + 1, 0.0);
    for (int i = 1; i < numZones_; ++i) {
        const int zoneLeft = i - 1, zoneRight = i;
        const double centerSpacing_cm =
            zoneCenter_cm[zoneRight] - zoneCenter_cm[zoneLeft];
        // Face inverse mean free path kappa_R * rho [1/cm].
        const double faceInverseMfp_percm =
            0.5 * (rosseland_cm2g[zoneLeft] * zoneDensity_gcc[zoneLeft] +
                   rosseland_cm2g[zoneRight] * zoneDensity_gcc[zoneRight]);
        const double faceEr_ergcc =
            std::max(0.5 * (radiationEnergyDensity_ergcc[zoneLeft] +
                            radiationEnergyDensity_ergcc[zoneRight]),
                     1e-300);
        // Levermore-Pomraning knudsen-like parameter R = |grad Er|/(kr rho Er),
        // evaluated with the lagged (beginning-of-step) field.
        const double gradientParameter =
            std::abs(radiationEnergyDensity_ergcc[zoneRight] -
                     radiationEnergyDensity_ergcc[zoneLeft]) /
            (centerSpacing_cm * faceInverseMfp_percm * faceEr_ergcc);
        const double diffusionCoefficient_cm2s =
            phys::c_light * levermorePomraningLambda(gradientParameter) /
            faceInverseMfp_percm;
        faceConductance_ccs[i] =
            diffusionCoefficient_cm2s * faceArea(nodeRadius_cm[i]) / centerSpacing_cm;
    }

    // ---- backward-Euler tridiagonal solve for Er^{n+1} ---------------------------
    std::vector<double> lower(numZones_), diag(numZones_), upper(numZones_),
        rhs(numZones_), solvedEr_ergcc(numZones_);
    for (int k = 0; k < numZones_; ++k) {
        const double volumeOverDt_ccs = zoneVolume_cc[k] / dt_s;
        // Effective emission/absorption coupling strength [cm^3/s].
        const double sourceStrength_ccs = zoneVolume_cc[k] * phys::c_light *
                                          planck_cm2g[k] * zoneDensity_gcc[k] *
                                          couplingFactor[k];
        const double equilibriumEr_ergcc =
            phys::a_rad * matterTemperature_eV[k] * matterTemperature_eV[k] *
            matterTemperature_eV[k] * matterTemperature_eV[k];
        lower[k] = -faceConductance_ccs[k];
        upper[k] = -faceConductance_ccs[k + 1];
        diag[k] = volumeOverDt_ccs + faceConductance_ccs[k] +
                  faceConductance_ccs[k + 1] + sourceStrength_ccs;
        rhs[k] = volumeOverDt_ccs * radiationEnergyDensity_ergcc[k] +
                 sourceStrength_ccs * equilibriumEr_ergcc;
    }
    // Vacuum (Marshak) leakage through the outer boundary, added implicitly
    // to the last zone's diagonal so the leak is evaluated at Er^{n+1}.
    double leakConductance_ccs = 0.0;
    if (deck_.radiation.bc_outer == "vacuum") {
        const double lastZoneWidth_cm =
            nodeRadius_cm[numZones_] - nodeRadius_cm[numZones_ - 1];
        const double leakSpeed_cmps =
            0.5 * phys::c_light /
            (1.0 + 0.75 * rosseland_cm2g[numZones_ - 1] *
                       zoneDensity_gcc[numZones_ - 1] * lastZoneWidth_cm);
        leakConductance_ccs = leakSpeed_cmps * faceArea(nodeRadius_cm[numZones_]);
        diag[numZones_ - 1] += leakConductance_ccs;
    }
    thomasSolve(lower, diag, upper, rhs, solvedEr_ergcc);

    if (leakConductance_ccs > 0.0)
        radiationLeaked_erg += leakConductance_ccs * solvedEr_ergcc[numZones_ - 1] *
                               dt_s;

    // ---- couple the exchanged energy back to the matter --------------------------
    // dTm below is EXACTLY the linearized response assumed in the solve, so
    // the matter gains precisely what the radiation lost (and vice versa).
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < numZones_; ++k) {
        const double equilibriumEr_ergcc =
            phys::a_rad * matterTemperature_eV[k] * matterTemperature_eV[k] *
            matterTemperature_eV[k] * matterTemperature_eV[k];
        const double matterTemperatureChange_eV =
            phys::c_light * planck_cm2g[k] * dt_s * couplingFactor[k] *
            (solvedEr_ergcc[k] - equilibriumEr_ergcc) / heatCapacity_erggeV[k];
        matterEnergy_ergg[k] += heatCapacity_erggeV[k] * matterTemperatureChange_eV;
        matterTemperature_eV[k] += matterTemperatureChange_eV;
        radiationEnergyDensity_ergcc[k] = solvedEr_ergcc[k];
    }
}

// ============================================================================
// Output: zone-by-zone snapshots, the per-step history file, and the
// end-of-run shot report.
// ============================================================================
void Simulation::writeSnapshot(int index) const {
    char filename[64];
    std::snprintf(filename, sizeof(filename), "snap_%05d.csv", index);
    std::ofstream out(std::filesystem::path(deck_.output.directory) / filename);
    out << "# t = " << time_s << " s, step = " << stepCount_ << "\n";
    out << "zone,material,r_left,r_right,r_center,u_left,u_right,rho,"
           "Ti_eV,Te_eV,Tr_eV,zbar,P,e,cs,q,alpha\n";
    out.precision(9);
    for (int k = 0; k < numZones_; ++k) {
        const auto& material = materials_[zoneMaterialIndex[k]];
        const double ionTemp_eV =
            twoTemperature_ ? ionTemperature_eV[k] : zoneTemperature_eV[k];
        const double electronTemp_eV =
            twoTemperature_ ? electronTemperature_eV[k] : zoneTemperature_eV[k];
        // Radiation temperature Tr = (Er/a)^{1/4}, zero when radiation off.
        const double radiationTemp_eV =
            radiationOn_
                ? std::pow(radiationEnergyDensity_ergcc[k] / phys::a_rad, 0.25)
                : 0.0;
        const double totalSpecificEnergy_ergg =
            twoTemperature_
                ? ionSpecificEnergy_ergg[k] + electronSpecificEnergy_ergg[k]
                : zoneSpecificEnergy_ergg[k];
        // Per-zone adiabat alpha = P / P_Fermi(rho): the compression metric.
        const double adiabat = zonePressure_dyncm2[k] /
                               fermiPressure0(zoneDensity_gcc[k], material.Z,
                                              material.A);
        out << k << ',' << material.name << ',' << nodeRadius_cm[k] << ','
            << nodeRadius_cm[k + 1] << ','
            << 0.5 * (nodeRadius_cm[k] + nodeRadius_cm[k + 1]) << ','
            << nodeVelocity_cmps[k] << ',' << nodeVelocity_cmps[k + 1] << ','
            << zoneDensity_gcc[k] << ',' << ionTemp_eV << ',' << electronTemp_eV
            << ',' << radiationTemp_eV << ',' << zoneMeanIonization[k] << ','
            << zonePressure_dyncm2[k] << ',' << totalSpecificEnergy_ergg << ','
            << zoneSoundSpeed_cmps[k] << ',' << zoneViscousPressure_dyncm2[k] << ','
            << adiabat << '\n';
    }
}

void Simulation::writeHistoryHeader() {
    historyFilePath_ =
        (std::filesystem::path(deck_.output.directory) / "history.csv").string();
    std::ofstream out(historyFilePath_);
    out << "step,t,dt,r_outer,u_outer,p_drive,P_laser,f_abs,rho_max,Ti_max,"
           "Te_max,Te_center,rhoR,E_int,E_kin,E_rad,W_drive,E_laser,E_floor,"
           "E_leak,E_err,P_fus,Y_n\n";
}

void Simulation::writeHistoryRow() {
    // Global sums for the energy-conservation diagnostic and quick-look
    // implosion metrics.
    double internalEnergy_erg = 0.0, kineticEnergy_erg = 0.0;
    double radiationEnergy_erg = 0.0, arealDensity_gcm2 = 0.0;
    double maxDensity_gcc = 0.0, maxIonTemp_eV = 0.0, maxElectronTemp_eV = 0.0;
    for (int k = 0; k < numZones_; ++k) {
        internalEnergy_erg +=
            zoneMass_g[k] * (twoTemperature_
                                 ? ionSpecificEnergy_ergg[k] +
                                       electronSpecificEnergy_ergg[k]
                                 : zoneSpecificEnergy_ergg[k]);
        if (radiationOn_)
            radiationEnergy_erg +=
                radiationEnergyDensity_ergcc[k] *
                shellVolume(nodeRadius_cm[k], nodeRadius_cm[k + 1]);
        arealDensity_gcm2 +=
            zoneDensity_gcc[k] * (nodeRadius_cm[k + 1] - nodeRadius_cm[k]);
        maxDensity_gcc = std::max(maxDensity_gcc, zoneDensity_gcc[k]);
        maxIonTemp_eV = std::max(
            maxIonTemp_eV, twoTemperature_ ? ionTemperature_eV[k] : zoneTemperature_eV[k]);
        maxElectronTemp_eV =
            std::max(maxElectronTemp_eV, twoTemperature_ ? electronTemperature_eV[k]
                                                         : zoneTemperature_eV[k]);
    }
    for (int i = 0; i <= numZones_; ++i)
        kineticEnergy_erg +=
            0.5 * nodeMass_g[i] * nodeVelocity_cmps[i] * nodeVelocity_cmps[i];
    const double totalEnergy_erg =
        internalEnergy_erg + kineticEnergy_erg + radiationEnergy_erg;
    // Relative energy-conservation error: everything the system holds now,
    // minus what it started with and every tracked source, plus every
    // tracked sink. Zero for a perfect scheme; the residual measures the
    // operator-splitting and hydro truncation error.
    const double scale_erg =
        std::max({std::abs(totalEnergy_erg), std::abs(driveWorkDone_erg),
                  std::abs(laserAbsorbed_erg), std::abs(initialTotalEnergy_erg),
                  1e-300});
    const double energyError = (totalEnergy_erg - initialTotalEnergy_erg -
                                driveWorkDone_erg - laserAbsorbed_erg -
                                floorEnergyInjected_erg + radiationLeaked_erg) /
                               scale_erg;

    std::ofstream out(historyFilePath_, std::ios::app);
    out.precision(9);
    out << stepCount_ << ',' << time_s << ',' << timeStep_s << ','
        << nodeRadius_cm[numZones_] << ',' << nodeVelocity_cmps[numZones_] << ','
        << deck_.drive.pressure(time_s) << ','
        << (laserOn_ ? deck_.laser.powerAt(time_s) : 0.0) << ','
        << laserAbsorbedFraction_ << ',' << maxDensity_gcc << ',' << maxIonTemp_eV
        << ',' << maxElectronTemp_eV << ','
        << (twoTemperature_ ? electronTemperature_eV[0] : zoneTemperature_eV[0])
        << ',' << arealDensity_gcm2 << ',' << internalEnergy_erg << ','
        << kineticEnergy_erg << ',' << radiationEnergy_erg << ','
        << driveWorkDone_erg << ',' << laserAbsorbed_erg << ','
        << floorEnergyInjected_erg << ',' << radiationLeaked_erg << ','
        << energyError << ',' << fusionPower_ergs << ','
        << neutronYieldDT + neutronYieldDDn << '\n';
}

// ============================================================================
// Fuel adiabat: mass-weighted alpha = P / P_Fermi(rho) over the DENSE part
// of the fuel (zones within 1/e of the current peak fuel density -- i.e. the
// compressed shell, excluding blown-off or unshocked material). P_Fermi is
// the T=0 fully-ionized electron Fermi pressure, the standard ICF reference
// (~2.2 rho^{5/3} Mbar for DT). alpha ~ 1 is fully degenerate fuel; ICF
// designs aim for alpha ~ 1-3 in flight.
// ============================================================================
double Simulation::fuelAdiabat() const {
    double peakFuelDensity_gcc = 0.0;
    for (int k = 0; k < numZones_; ++k)
        if (zoneIsFuel_[k])
            peakFuelDensity_gcc = std::max(peakFuelDensity_gcc, zoneDensity_gcc[k]);
    if (peakFuelDensity_gcc <= 0.0) return -1.0;
    const double denseThreshold_gcc = peakFuelDensity_gcc / 2.718281828;
    double massSum_g = 0.0, weightedAdiabatSum_g = 0.0;
    for (int k = 0; k < numZones_; ++k) {
        if (!zoneIsFuel_[k] || zoneDensity_gcc[k] < denseThreshold_gcc) continue;
        const auto& material = materials_[zoneMaterialIndex[k]];
        weightedAdiabatSum_g +=
            zoneMass_g[k] * zonePressure_dyncm2[k] /
            fermiPressure0(zoneDensity_gcc[k], material.Z, material.A);
        massSum_g += zoneMass_g[k];
    }
    return (massSum_g > 0.0) ? weightedAdiabatSum_g / massSum_g : -1.0;
}

// ============================================================================
// Shot-report accumulation, called once per step. Everything here is
// diagnostic-only; nothing feeds back into the physics.
// ============================================================================
void Simulation::updateReport() {
    auto& report = report_;

    // ---- global extrema over the whole run ------------------------------------
    double kineticEnergy_erg = 0.0;
    for (int i = 0; i <= numZones_; ++i)
        kineticEnergy_erg +=
            0.5 * nodeMass_g[i] * nodeVelocity_cmps[i] * nodeVelocity_cmps[i];
    report.maxKineticEnergy_erg =
        std::max(report.maxKineticEnergy_erg, kineticEnergy_erg);
    for (int k = 0; k < numZones_; ++k) {
        report.maxDensity_gcc = std::max(report.maxDensity_gcc, zoneDensity_gcc[k]);
        report.maxIonTemperature_eV =
            std::max(report.maxIonTemperature_eV,
                     twoTemperature_ ? ionTemperature_eV[k] : zoneTemperature_eV[k]);
        report.maxElectronTemperature_eV =
            std::max(report.maxElectronTemperature_eV,
                     twoTemperature_ ? electronTemperature_eV[k]
                                     : zoneTemperature_eV[k]);
    }
    report.peakDrivePressure_dyncm2 =
        std::max(report.peakDrivePressure_dyncm2, deck_.drive.pressure(time_s));
    if (laserOn_)
        report.peakLaserPower_ergs =
            std::max(report.peakLaserPower_ergs, deck_.laser.powerAt(time_s));

    // ---- fuel implosion speed (mass-averaged, inward positive) ----------------
    // Also accumulates the areal densities used just below.
    double fuelMass_g = 0.0, fuelMomentum_gcmps = 0.0;
    double fuelRhoR_gcm2 = 0.0, totalRhoR_gcm2 = 0.0;
    for (int k = 0; k < numZones_; ++k) {
        const double zoneWidth_cm = nodeRadius_cm[k + 1] - nodeRadius_cm[k];
        totalRhoR_gcm2 += zoneDensity_gcc[k] * zoneWidth_cm;
        if (!zoneIsFuel_[k]) continue;
        fuelMass_g += zoneMass_g[k];
        fuelMomentum_gcmps += zoneMass_g[k] * 0.5 *
                              (nodeVelocity_cmps[k] + nodeVelocity_cmps[k + 1]);
        fuelRhoR_gcm2 += zoneDensity_gcc[k] * zoneWidth_cm;
    }
    const double inwardSpeed_cmps =
        (fuelMass_g > 0.0) ? -fuelMomentum_gcmps / fuelMass_g : 0.0;
    if (inwardSpeed_cmps > report.peakImplosionSpeed_cmps) {
        report.peakImplosionSpeed_cmps = inwardSpeed_cmps;
        report.peakImplosionTime_s = time_s;
        // The design-relevant adiabat is the one the shell carries in
        // flight, so record it at the moment of peak velocity.
        report.fuelAdiabatAtPeakSpeed = fuelAdiabat();
    }

    // ---- peak fuel compression ("bang" proxy) and hot-spot state ---------------
    if (fuelRhoR_gcm2 > report.peakFuelRhoR_gcm2) {
        report.peakFuelRhoR_gcm2 = fuelRhoR_gcm2;
        report.totalRhoRAtBang_gcm2 = totalRhoR_gcm2;
        report.bangTime_s = time_s;
        report.hotSpotRadius_cm = nodeRadius_cm[hotSpotNode_];
        // Hot spot = all zones of the innermost layer (inside hotSpotNode_).
        double hotSpotMass_g = 0.0, weightedTi_geV = 0.0, weightedTe_geV = 0.0;
        double weightedP_gdyncm2 = 0.0, hotSpotRhoR_gcm2 = 0.0;
        for (int k = 0; k < hotSpotNode_ && k < numZones_; ++k) {
            hotSpotMass_g += zoneMass_g[k];
            weightedTi_geV +=
                zoneMass_g[k] *
                (twoTemperature_ ? ionTemperature_eV[k] : zoneTemperature_eV[k]);
            weightedTe_geV +=
                zoneMass_g[k] * (twoTemperature_ ? electronTemperature_eV[k]
                                                 : zoneTemperature_eV[k]);
            weightedP_gdyncm2 += zoneMass_g[k] * zonePressure_dyncm2[k];
            hotSpotRhoR_gcm2 +=
                zoneDensity_gcc[k] * (nodeRadius_cm[k + 1] - nodeRadius_cm[k]);
        }
        if (hotSpotMass_g > 0.0) {
            report.hotSpotTi_eV = weightedTi_geV / hotSpotMass_g;
            report.hotSpotTe_eV = weightedTe_geV / hotSpotMass_g;
            report.hotSpotPressure_dyncm2 = weightedP_gdyncm2 / hotSpotMass_g;
            report.hotSpotRhoR_gcm2 = hotSpotRhoR_gcm2;
        }
    }

    // ---- convergence ratio and in-flight aspect ratio ---------------------------
    report.minHotSpotRadius_cm =
        std::min(report.minHotSpotRadius_cm, nodeRadius_cm[hotSpotNode_]);
    // IFAR is measured once, the first time the hot-spot boundary passes 2/3
    // of its initial radius (the conventional in-flight sampling point).
    if (report.inFlightAspectRatio < 0.0 && hotSpotRadius0_cm > 0.0 &&
        nodeRadius_cm[hotSpotNode_] < (2.0 / 3.0) * hotSpotRadius0_cm) {
        double peakFuelDensity_gcc = 0.0;
        for (int k = 0; k < numZones_; ++k)
            if (zoneIsFuel_[k])
                peakFuelDensity_gcc =
                    std::max(peakFuelDensity_gcc, zoneDensity_gcc[k]);
        const double denseThreshold_gcc = peakFuelDensity_gcc / 2.718281828;
        // Shell = dense fuel (within 1/e of the peak); IFAR = its outer
        // radius over its total thickness.
        double shellOuterRadius_cm = 0.0, shellThickness_cm = 0.0;
        for (int k = 0; k < numZones_; ++k) {
            if (!zoneIsFuel_[k] || zoneDensity_gcc[k] < denseThreshold_gcc) continue;
            shellOuterRadius_cm = std::max(shellOuterRadius_cm, nodeRadius_cm[k + 1]);
            shellThickness_cm += nodeRadius_cm[k + 1] - nodeRadius_cm[k];
        }
        if (shellThickness_cm > 0.0) {
            report.inFlightAspectRatio = shellOuterRadius_cm / shellThickness_cm;
            report.ifarTime_s = time_s;
        }
    }
}

void Simulation::writeReport() const {
    const auto& report = report_;
    std::string text;
    char line[256];
    auto add = [&](const char* format, auto... args) {
        if constexpr (sizeof...(args) == 0) {
            text += format;
        } else {
            std::snprintf(line, sizeof(line), format, args...);
            text += line;
        }
        text += '\n';
    };

    add("=============== hydro1d shot report ===============");
    add("run: %d zones, %s, conduction %s, radiation %s, laser %s",
        numZones_, twoTemperature_ ? "2T" : "1T",
        deck_.conduction.enabled ? "on" : "off", radiationOn_ ? "on" : "off",
        laserOn_ ? "on" : "off");
    add("end: t = %.4g s in %ld steps", time_s, stepCount_);
    if (noFuelFlag_)
        add("NOTE: no material has 'fuel = true'; fuel metrics use all zones");

    if (laserOn_) {
        add("laser: E_inc = %.4g erg (%.1f kJ), absorbed = %.4g erg (%.1f kJ), "
            "coupling = %.1f%%, peak power = %.3g erg/s",
            laserIncident_erg, laserIncident_erg / 1e10, laserAbsorbed_erg,
            laserAbsorbed_erg / 1e10,
            (laserIncident_erg > 0.0) ? 100.0 * laserAbsorbed_erg / laserIncident_erg
                                      : 0.0,
            report.peakLaserPower_ergs);
    }
    if (report.peakDrivePressure_dyncm2 > 0.0)
        add("drive: peak pressure = %.3g dyn/cm^2 (%.1f Mbar), work = %.4g erg",
            report.peakDrivePressure_dyncm2, report.peakDrivePressure_dyncm2 / 1e12,
            driveWorkDone_erg);

    if (deck_.control.geometry == 3) {
        add("implosion:");
        add("  peak implosion speed  = %.1f km/s at %.4g s (fuel mass-avg)",
            report.peakImplosionSpeed_cmps / 1e5, report.peakImplosionTime_s);
        if (report.fuelAdiabatAtPeakSpeed > 0.0)
            add("  fuel adiabat then     = %.2f (mass-avg P/P_Fermi, dense shell)",
                report.fuelAdiabatAtPeakSpeed);
        if (report.inFlightAspectRatio > 0.0)
            add("  IFAR (at 2/3 R0)      = %.1f at %.4g s",
                report.inFlightAspectRatio, report.ifarTime_s);
        add("  convergence ratio     = %.1f (hot-spot boundary %.4g -> %.4g cm)",
            (report.minHotSpotRadius_cm > 0.0)
                ? hotSpotRadius0_cm / report.minHotSpotRadius_cm
                : -1.0,
            hotSpotRadius0_cm, report.minHotSpotRadius_cm);
        add("  bang time (peak fuel rhoR) = %.4g s", report.bangTime_s);
        add("  peak fuel rhoR        = %.3g g/cm^2 (total rhoR then: %.3g)",
            report.peakFuelRhoR_gcm2, report.totalRhoRAtBang_gcm2);
        add("  hot spot at bang: R = %.1f um, <Ti> = %.3g eV, <Te> = %.3g eV",
            report.hotSpotRadius_cm * 1e4, report.hotSpotTi_eV, report.hotSpotTe_eV);
        add("                    <P> = %.3g dyn/cm^2 (%.2f Gbar), rhoR = %.3g g/cm^2",
            report.hotSpotPressure_dyncm2, report.hotSpotPressure_dyncm2 / 1e15,
            report.hotSpotRhoR_gcm2);
    }
    if (burnOn_) {
        add("burn (diagnostic only, no self-heating):");
        add("  DT neutron yield      = %.4g  (DD-n yield: %.3g)", neutronYieldDT,
            neutronYieldDDn);
        add("  fusion energy         = %.4g erg (%.3g kJ)", fusionEnergy_erg,
            fusionEnergy_erg / 1e10);
        const double inputEnergy_erg = laserOn_ ? laserIncident_erg : driveWorkDone_erg;
        if (inputEnergy_erg > 0.0)
            add("  target gain           = %.3g (vs %s energy)",
                fusionEnergy_erg / inputEnergy_erg,
                laserOn_ ? "incident laser" : "drive work");
        add("  bang time (peak fusion power) = %.4g s, peak P_fus = %.4g erg/s",
            fusionBangTime_s, peakFusionPower_ergs);
        if (peakFusionPower_ergs > 0.0)
            add("  burn width (E_fus/P_fus,peak) = %.3g s",
                fusionEnergy_erg / peakFusionPower_ergs);
        if (burnWeightSum > 0.0)
            add("  burn-averaged Ti      = %.3g keV",
                burnWeightedTiSum_keV / burnWeightSum);
    }
    add("extrema: rho_max = %.4g g/cc, Ti_max = %.4g eV, Te_max = %.4g eV, "
        "E_kin_max = %.4g erg",
        report.maxDensity_gcc, report.maxIonTemperature_eV,
        report.maxElectronTemperature_eV, report.maxKineticEnergy_erg);
    add("energy bookkeeping: floors injected %.3g erg, radiation leaked %.3g erg",
        floorEnergyInjected_erg, radiationLeaked_erg);
    add("===================================================");

    std::cout << text;
    std::ofstream out(std::filesystem::path(deck_.output.directory) / "report.txt");
    out << text;
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
