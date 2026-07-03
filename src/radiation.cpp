// ============================================================================
// radiation.cpp -- grey flux-limited radiation diffusion with linearized-
// Planck matter coupling, backward Euler with a tridiagonal solve. Part of
// the Simulation class; see simulation.hpp for the source layout and the
// unit-suffix naming convention.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"
#include "solvers.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
// Levermore-Pomraning flux limiter for radiation diffusion. The argument
// R = |grad Er| / (kappa_R rho Er) measures how steep the radiation field is
// in units of the photon mean free path: lambda(R->0) = 1/3 recovers classic
// diffusion, lambda(R->inf) ~ 1/R caps the flux at the free-streaming value
// c*Er.
double levermorePomraningLambda(double R) {
    return (2.0 + R) / (6.0 + 3.0 * R + R * R);
}
}  // namespace

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
