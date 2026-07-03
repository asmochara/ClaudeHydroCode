// ============================================================================
// coupling.cpp -- electron-ion temperature relaxation (2T mode only).
// Part of the Simulation class; see simulation.hpp for the source layout
// and the unit-suffix naming convention.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>

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
