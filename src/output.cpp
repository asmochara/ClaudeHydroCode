// ============================================================================
// output.cpp -- everything the code writes: zone-by-zone snapshots, the
// per-step history file, and the end-of-run shot report with its running
// design-metric accumulation. Part of the Simulation class; see
// simulation.hpp for the source layout and the unit-suffix convention.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

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
