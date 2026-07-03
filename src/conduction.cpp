// ============================================================================
// conduction.cpp -- flux-limited thermal conduction: Spitzer-Harm electron
// conduction (1T and 2T) and Braginskii ion conduction (2T), each solved
// backward-Euler with a tridiagonal system and a conservative flux-form
// energy update. Part of the Simulation class; see simulation.hpp for the
// source layout and the unit-suffix naming convention.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"
#include "solvers.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
// Braginskii electron-conduction coefficient gamma0 as a function of ion
// charge: kappa_e = gamma0 * ne kB Te tau_e / me. This rational fit passes
// through the tabulated Braginskii values (3.16 at Z=1) and approaches the
// Lorentz-gas limit at high Z.
double braginskiiGamma0(double ionCharge) {
    return 13.58 * (ionCharge + 0.24) / (ionCharge + 4.24);
}
}  // namespace

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
