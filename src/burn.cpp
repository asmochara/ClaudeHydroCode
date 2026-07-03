// ============================================================================
// burn.cpp -- burn-off fusion diagnostics (no feedback on the hydro).
// Part of the Simulation class; see simulation.hpp for the source layout
// and the unit-suffix naming convention.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"
#include "fusion.hpp"

#include <cmath>

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
