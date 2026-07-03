// ============================================================================
// eos.cpp -- equation-of-state implementations. Units throughout:
// density [g/cm^3], temperature [eV], pressure [dyn/cm^2], specific energy
// [erg/g], cv [erg/g/eV], sound speed squared [cm^2/s^2].
// ============================================================================

#include "eos.hpp"
#include "constants.hpp"
#include "numerics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// ============================================================================
// IdealGasEOS (1T total matter): P = rho (1+Zbar) kB T/(A m_p),
// e = P/((gamma-1) rho). With a T-dependent Zbar (Thomas-Fermi) the energy
// is no longer linear in T, so cv and the inversions are done numerically;
// with a fixed Zbar the generic routines reproduce the analytic results
// exactly (the secant inversion converges in one step on a linear function).
// ============================================================================

IdealGasEOS::IdealGasEOS(double gamma, double atomicMass_amu,
                         std::shared_ptr<const ZbarModel> ionizationModel)
    : gamma_(gamma),
      gasConstant_erggeV(phys::eV / (atomicMass_amu * phys::m_p)),
      ionizationModel_(std::move(ionizationModel)) {
    if (gamma <= 1.0) throw std::runtime_error("IdealGasEOS: gamma must be > 1");
}

double IdealGasEOS::pressure(double density_gcc, double temperature_eV) const {
    return density_gcc * gasConstant_erggeV *
           (1.0 + ionizationModel_->zbar(density_gcc, temperature_eV)) *
           temperature_eV;
}

double IdealGasEOS::energy(double density_gcc, double temperature_eV) const {
    return gasConstant_erggeV *
           (1.0 + ionizationModel_->zbar(density_gcc, temperature_eV)) *
           temperature_eV / (gamma_ - 1.0);
}

double IdealGasEOS::cv(double density_gcc, double temperature_eV) const {
    // Central finite difference of e(T); the interval is clipped at T = 0.
    const double dT_eV = 1e-3 * std::max(temperature_eV, 1e-12);
    return (energy(density_gcc, temperature_eV + dT_eV) -
            energy(density_gcc, std::max(temperature_eV - dT_eV, 0.0))) /
           (dT_eV + std::min(dT_eV, temperature_eV));
}

double IdealGasEOS::temperature(double density_gcc, double specificEnergy_ergg,
                                double guessTemperature_eV) const {
    return invertMonotone(
        [&](double trialT_eV) { return energy(density_gcc, trialT_eV); },
        specificEnergy_ergg, 0.0, 1e9, guessTemperature_eV);
}

double IdealGasEOS::temperatureFromPressure(double density_gcc,
                                            double pressure_dyncm2) const {
    return invertMonotone(
        [&](double trialT_eV) { return pressure(density_gcc, trialT_eV); },
        pressure_dyncm2, 0.0, 1e9,
        pressure_dyncm2 / (density_gcc * 2.0 * gasConstant_erggeV));
}

double IdealGasEOS::soundSpeed2(double density_gcc, double temperature_eV) const {
    // gamma P / rho; ionization derivatives are neglected, which is adequate
    // for the CFL condition and artificial viscosity this feeds.
    return gamma_ * pressure(density_gcc, temperature_eV) / density_gcc;
}

// ============================================================================
// TabulatedEOS (1T total matter from a two-block table file).
// ============================================================================

TabulatedEOS::TabulatedEOS(const std::string& path) : table_(path, 2) {}

double TabulatedEOS::pressure(double density_gcc, double temperature_eV) const {
    return table_.interp(0, density_gcc, temperature_eV);
}
double TabulatedEOS::energy(double density_gcc, double temperature_eV) const {
    return table_.interp(1, density_gcc, temperature_eV);
}

double TabulatedEOS::cv(double density_gcc, double temperature_eV) const {
    // Finite difference with the interval clamped inside the table's
    // temperature range.
    const double dT_eV = 1e-3 * std::max(temperature_eV, table_.Tmin());
    const double Tlow_eV = std::max(temperature_eV - dT_eV, table_.Tmin());
    const double Thigh_eV = std::min(temperature_eV + dT_eV, table_.Tmax());
    return (energy(density_gcc, Thigh_eV) - energy(density_gcc, Tlow_eV)) /
           (Thigh_eV - Tlow_eV);
}

double TabulatedEOS::temperature(double density_gcc, double specificEnergy_ergg,
                                 double guessTemperature_eV) const {
    return invertMonotone(
        [&](double trialT_eV) { return energy(density_gcc, trialT_eV); },
        specificEnergy_ergg, table_.Tmin(), table_.Tmax(), guessTemperature_eV);
}

double TabulatedEOS::temperatureFromPressure(double density_gcc,
                                             double pressure_dyncm2) const {
    return invertMonotone(
        [&](double trialT_eV) { return pressure(density_gcc, trialT_eV); },
        pressure_dyncm2, table_.Tmin(), table_.Tmax(),
        std::sqrt(table_.Tmin() * table_.Tmax()));
}

double TabulatedEOS::soundSpeed2(double density_gcc, double temperature_eV) const {
    // Thermodynamic identity, evaluated by finite differences of the table:
    //   cs^2 = (dP/drho)_T + T (dP/dT)_rho^2 / (rho^2 cv)
    const double dRho_gcc = 1e-3 * density_gcc;
    const double dT_eV = 1e-3 * std::max(temperature_eV, table_.Tmin());
    const double rhoLow_gcc = std::max(density_gcc - dRho_gcc, 1e-300);
    const double dPdRho =
        (pressure(density_gcc + dRho_gcc, temperature_eV) -
         pressure(rhoLow_gcc, temperature_eV)) /
        (density_gcc + dRho_gcc - rhoLow_gcc);
    const double Tlow_eV = std::max(temperature_eV - dT_eV, table_.Tmin());
    const double Thigh_eV = std::min(temperature_eV + dT_eV, table_.Tmax());
    const double dPdT = (pressure(density_gcc, Thigh_eV) -
                         pressure(density_gcc, Tlow_eV)) /
                        (Thigh_eV - Tlow_eV);
    const double heatCapacity_erggeV =
        std::max(cv(density_gcc, temperature_eV), 1e-300);
    const double soundSpeedSq =
        dPdRho + temperature_eV * dPdT * dPdT /
                     (density_gcc * density_gcc * heatCapacity_erggeV);
    // Floor keeps the CFL estimate sane on rough or noisy tables.
    return std::max(soundSpeedSq,
                    1e-6 * pressure(density_gcc, temperature_eV) / density_gcc);
}

// ============================================================================
// SpeciesEOS base class: generic numeric implementations that only need
// pressure() and energy() from the concrete species.
// ============================================================================

double SpeciesEOS::cv(double density_gcc, double temperature_eV) const {
    const double dT_eV = 1e-3 * std::max(temperature_eV, 1e-12);
    const double Tlow_eV = std::max(temperature_eV - dT_eV, Tlo());
    const double Thigh_eV = std::min(temperature_eV + dT_eV, Thi());
    // Floored positive so implicit solvers can always divide by cv.
    return std::max((energy(density_gcc, Thigh_eV) -
                     energy(density_gcc, Tlow_eV)) /
                        (Thigh_eV - Tlow_eV),
                    1e-30);
}

double SpeciesEOS::temperature(double density_gcc, double specificEnergy_ergg,
                               double guessTemperature_eV) const {
    return invertMonotone(
        [&](double trialT_eV) { return energy(density_gcc, trialT_eV); },
        specificEnergy_ergg, Tlo(), Thi(), guessTemperature_eV);
}

double SpeciesEOS::cs2Contribution(double density_gcc, double temperature_eV) const {
    // Same thermodynamic identity as TabulatedEOS::soundSpeed2, by finite
    // differences of this species' pressure.
    const double dRho_gcc = 1e-3 * density_gcc;
    const double dT_eV = 1e-3 * std::max(temperature_eV, 1e-12);
    const double rhoLow_gcc = std::max(density_gcc - dRho_gcc, 1e-300);
    const double dPdRho =
        (pressure(density_gcc + dRho_gcc, temperature_eV) -
         pressure(rhoLow_gcc, temperature_eV)) /
        (density_gcc + dRho_gcc - rhoLow_gcc);
    const double Tlow_eV = std::max(temperature_eV - dT_eV, Tlo());
    const double Thigh_eV = std::min(temperature_eV + dT_eV, Thi());
    const double dPdT = (pressure(density_gcc, Thigh_eV) -
                         pressure(density_gcc, Tlow_eV)) /
                        (Thigh_eV - Tlow_eV);
    const double soundSpeedSq =
        dPdRho + temperature_eV * dPdT * dPdT /
                     (density_gcc * density_gcc * cv(density_gcc, temperature_eV));
    return std::max(soundSpeedSq, 0.0);
}

// ============================================================================
// IdealIonEOS: classical ideal gas of the ions alone. Analytic overrides
// (everything is linear in T) avoid the generic numeric machinery.
// ============================================================================

IdealIonEOS::IdealIonEOS(double gamma, double atomicMass_amu)
    : gamma_(gamma),
      gasConstant_erggeV(phys::eV / (atomicMass_amu * phys::m_p)) {
    if (gamma <= 1.0) throw std::runtime_error("IdealIonEOS: gamma must be > 1");
}

double IdealIonEOS::pressure(double density_gcc, double temperature_eV) const {
    return density_gcc * gasConstant_erggeV * temperature_eV;
}
double IdealIonEOS::energy(double /*density_gcc*/, double temperature_eV) const {
    return gasConstant_erggeV * temperature_eV / (gamma_ - 1.0);
}
double IdealIonEOS::cv(double, double) const {
    return gasConstant_erggeV / (gamma_ - 1.0);
}
double IdealIonEOS::temperature(double, double specificEnergy_ergg, double) const {
    return std::max(specificEnergy_ergg, 0.0) * (gamma_ - 1.0) / gasConstant_erggeV;
}
double IdealIonEOS::cs2Contribution(double density_gcc, double temperature_eV) const {
    return gamma_ * pressure(density_gcc, temperature_eV) / density_gcc;
}

// ============================================================================
// IdealElectronEOS: classical electron gas, optionally with Fermi
// degeneracy (see the header for the physics discussion).
// ============================================================================

namespace {
// T=0 Fermi pressure of an electron gas of number density ne [cm^-3]:
//   E_F = hbar^2 (3 pi^2 ne)^{2/3} / (2 me),  P_F0 = (2/5) ne E_F.
inline double fermiPressureOfNe(double electronDensity_percc) {
    constexpr double hbarSqOver2me = phys::hbar * phys::hbar / (2.0 * phys::m_e);
    const double fermiEnergy_erg =
        hbarSqOver2me *
        std::pow(3.0 * phys::pi * phys::pi * electronDensity_percc, 2.0 / 3.0);
    return 0.4 * electronDensity_percc * fermiEnergy_erg;
}
}  // namespace

double fermiPressure0(double density_gcc, double nuclearCharge,
                      double atomicMass_amu) {
    return fermiPressureOfNe(density_gcc * nuclearCharge /
                             (atomicMass_amu * phys::m_p));
}

IdealElectronEOS::IdealElectronEOS(double gamma, double atomicMass_amu,
                                   std::shared_ptr<const ZbarModel> ionizationModel,
                                   bool includeDegeneracy)
    : gamma_(gamma),
      gasConstant_erggeV(phys::eV / (atomicMass_amu * phys::m_p)),
      includeDegeneracy_(includeDegeneracy),
      ionizationModel_(std::move(ionizationModel)) {
    if (gamma <= 1.0)
        throw std::runtime_error("IdealElectronEOS: gamma must be > 1");
}

double IdealElectronEOS::pressure(double density_gcc, double temperature_eV) const {
    const double zbar = ionizationModel_->zbar(density_gcc, temperature_eV);
    const double classicalPressure_dyncm2 =
        density_gcc * gasConstant_erggeV * zbar * temperature_eV;
    if (!includeDegeneracy_) return classicalPressure_dyncm2;
    // Quadrature interpolation between the classical and fully degenerate
    // limits. ne = rho*zbar/(A m_p) = rho*zbar*gasConstant/eV.
    const double fermiPressure_dyncm2 =
        fermiPressureOfNe(density_gcc * zbar * gasConstant_erggeV / phys::eV);
    return std::sqrt(classicalPressure_dyncm2 * classicalPressure_dyncm2 +
                     fermiPressure_dyncm2 * fermiPressure_dyncm2);
}

double IdealElectronEOS::energy(double density_gcc, double temperature_eV) const {
    // e = (3/2) P/rho is EXACT for a nonrelativistic Fermi gas at any
    // degeneracy (P = 2/3 of the energy density), so using it keeps P and e
    // mutually consistent. It implicitly assumes gamma = 5/3 for the
    // degenerate branch. Note the T -> 0 limit is the finite zero-point
    // energy (3/5) E_F per electron, not zero.
    if (includeDegeneracy_)
        return 1.5 * pressure(density_gcc, temperature_eV) / density_gcc;
    return gasConstant_erggeV *
           ionizationModel_->zbar(density_gcc, temperature_eV) * temperature_eV /
           (gamma_ - 1.0);
}

double IdealElectronEOS::cs2Contribution(double density_gcc,
                                         double temperature_eV) const {
    // gamma P / rho with the (possibly degeneracy-stiffened) pressure;
    // ionization derivatives neglected, as in the 1T ideal EOS.
    return gamma_ * pressure(density_gcc, temperature_eV) / density_gcc;
}

// ============================================================================
// TableSpeciesEOS: species partial EOS from a two-block table file.
// ============================================================================

TableSpeciesEOS::TableSpeciesEOS(const std::string& path) : table_(path, 2) {}

double TableSpeciesEOS::pressure(double density_gcc, double temperature_eV) const {
    return table_.interp(0, density_gcc, temperature_eV);
}
double TableSpeciesEOS::energy(double density_gcc, double temperature_eV) const {
    return table_.interp(1, density_gcc, temperature_eV);
}
