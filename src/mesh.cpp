// ============================================================================
// mesh.cpp -- mesh construction and thermodynamic state maintenance:
//   setupMesh()              lay down the layers, set initial conditions,
//                            freeze the Lagrangian zone masses
//   updateThermodynamics()   re-derive (T, P, Zbar, cs) from (rho, e)
//   applyTemperatureFloors() clamp temperatures at the user floor, tracking
//                            the injected energy
// Part of the Simulation class; see simulation.hpp for the source layout
// and the unit-suffix naming convention.
// ============================================================================

#include "simulation.hpp"
#include "constants.hpp"
#include "numerics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

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
