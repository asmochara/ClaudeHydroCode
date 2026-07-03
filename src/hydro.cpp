// ============================================================================
// hydro.cpp -- the explicit Lagrangian hydrodynamics step (von Neumann-
// Richtmyer leapfrog): node acceleration from the stress gradient, node
// motion, artificial viscosity, and the predictor-corrector PdV energy
// update. Part of the Simulation class; see simulation.hpp for the source
// layout and the unit-suffix naming convention.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

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
