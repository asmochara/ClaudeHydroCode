// ============================================================================
// laser.cpp -- laser illumination: the equal-power ray set built from the
// focal-spot intensity profile, and the refracting ray trace with inverse-
// bremsstrahlung absorption, the critical-surface dump, and the Langdon
// absorption correction. Part of the Simulation class; see simulation.hpp
// for the source layout and the unit-suffix naming convention.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>
#include <vector>

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
