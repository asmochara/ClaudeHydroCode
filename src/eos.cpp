#include "eos.hpp"
#include "constants.hpp"
#include "numerics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// ---------------------------------------------------------------- IdealGasEOS

IdealGasEOS::IdealGasEOS(double gamma, double A, std::shared_ptr<const ZbarModel> zb)
    : gamma_(gamma), R0_(phys::eV / (A * phys::m_p)), zb_(std::move(zb)) {
    if (gamma <= 1.0) throw std::runtime_error("IdealGasEOS: gamma must be > 1");
}

double IdealGasEOS::pressure(double rho, double T) const {
    return rho * R0_ * (1.0 + zb_->zbar(rho, T)) * T;
}

double IdealGasEOS::energy(double rho, double T) const {
    return R0_ * (1.0 + zb_->zbar(rho, T)) * T / (gamma_ - 1.0);
}

double IdealGasEOS::cv(double rho, double T) const {
    const double dT = 1e-3 * std::max(T, 1e-12);
    return (energy(rho, T + dT) - energy(rho, std::max(T - dT, 0.0))) /
           (dT + std::min(dT, T));
}

double IdealGasEOS::temperature(double rho, double e, double Tguess) const {
    return invertMonotone([&](double T) { return energy(rho, T); },
                          e, 0.0, 1e9, Tguess);
}

double IdealGasEOS::temperatureFromPressure(double rho, double P) const {
    return invertMonotone([&](double T) { return pressure(rho, T); },
                          P, 0.0, 1e9, P / (rho * 2.0 * R0_));
}

double IdealGasEOS::soundSpeed2(double rho, double T) const {
    // Ionization derivatives neglected here; adequate for CFL/viscosity use.
    return gamma_ * pressure(rho, T) / rho;
}

// --------------------------------------------------------------- TabulatedEOS

TabulatedEOS::TabulatedEOS(const std::string& path) : tab_(path, 2) {}

double TabulatedEOS::pressure(double rho, double T) const { return tab_.interp(0, rho, T); }
double TabulatedEOS::energy(double rho, double T) const   { return tab_.interp(1, rho, T); }

double TabulatedEOS::cv(double rho, double T) const {
    const double dT = 1e-3 * std::max(T, tab_.Tmin());
    const double Tl = std::max(T - dT, tab_.Tmin());
    const double Th = std::min(T + dT, tab_.Tmax());
    return (energy(rho, Th) - energy(rho, Tl)) / (Th - Tl);
}

double TabulatedEOS::temperature(double rho, double e, double Tguess) const {
    return invertMonotone([&](double T) { return energy(rho, T); },
                          e, tab_.Tmin(), tab_.Tmax(), Tguess);
}

double TabulatedEOS::temperatureFromPressure(double rho, double P) const {
    return invertMonotone([&](double T) { return pressure(rho, T); },
                          P, tab_.Tmin(), tab_.Tmax(), std::sqrt(tab_.Tmin() * tab_.Tmax()));
}

double TabulatedEOS::soundSpeed2(double rho, double T) const {
    // cs^2 = (dP/drho)_T + T (dP/dT)_rho^2 / (rho^2 cv)
    const double drho = 1e-3 * rho;
    const double dT = 1e-3 * std::max(T, tab_.Tmin());
    const double rl = std::max(rho - drho, 1e-300);
    const double dPdrho = (pressure(rho + drho, T) - pressure(rl, T)) / (rho + drho - rl);
    const double Tl = std::max(T - dT, tab_.Tmin()), Th = std::min(T + dT, tab_.Tmax());
    const double dPdT = (pressure(rho, Th) - pressure(rho, Tl)) / (Th - Tl);
    const double c = std::max(cv(rho, T), 1e-300);
    const double cs2 = dPdrho + T * dPdT * dPdT / (rho * rho * c);
    return std::max(cs2, 1e-6 * pressure(rho, T) / rho);
}

// ----------------------------------------------------------------- SpeciesEOS

double SpeciesEOS::cv(double rho, double T) const {
    const double dT = 1e-3 * std::max(T, 1e-12);
    const double Tl = std::max(T - dT, Tlo());
    const double Th = std::min(T + dT, Thi());
    return std::max((energy(rho, Th) - energy(rho, Tl)) / (Th - Tl), 1e-30);
}

double SpeciesEOS::temperature(double rho, double e, double Tguess) const {
    return invertMonotone([&](double T) { return energy(rho, T); },
                          e, Tlo(), Thi(), Tguess);
}

double SpeciesEOS::cs2Contribution(double rho, double T) const {
    const double drho = 1e-3 * rho;
    const double dT = 1e-3 * std::max(T, 1e-12);
    const double rl = std::max(rho - drho, 1e-300);
    const double dPdrho = (pressure(rho + drho, T) - pressure(rl, T)) / (rho + drho - rl);
    const double Tl = std::max(T - dT, Tlo()), Th = std::min(T + dT, Thi());
    const double dPdT = (pressure(rho, Th) - pressure(rho, Tl)) / (Th - Tl);
    const double cs2 = dPdrho + T * dPdT * dPdT / (rho * rho * cv(rho, T));
    return std::max(cs2, 0.0);
}

// ----------------------------------------------------------------- IdealIonEOS

IdealIonEOS::IdealIonEOS(double gamma, double A)
    : gamma_(gamma), R0_(phys::eV / (A * phys::m_p)) {
    if (gamma <= 1.0) throw std::runtime_error("IdealIonEOS: gamma must be > 1");
}

double IdealIonEOS::pressure(double rho, double T) const { return rho * R0_ * T; }
double IdealIonEOS::energy(double /*rho*/, double T) const { return R0_ * T / (gamma_ - 1.0); }
double IdealIonEOS::cv(double, double) const { return R0_ / (gamma_ - 1.0); }
double IdealIonEOS::temperature(double, double e, double) const {
    return std::max(e, 0.0) * (gamma_ - 1.0) / R0_;
}
double IdealIonEOS::cs2Contribution(double rho, double T) const {
    return gamma_ * pressure(rho, T) / rho;
}

// ------------------------------------------------------------ IdealElectronEOS

IdealElectronEOS::IdealElectronEOS(double gamma, double A,
                                   std::shared_ptr<const ZbarModel> zb)
    : gamma_(gamma), R0_(phys::eV / (A * phys::m_p)), zb_(std::move(zb)) {
    if (gamma <= 1.0) throw std::runtime_error("IdealElectronEOS: gamma must be > 1");
}

double IdealElectronEOS::pressure(double rho, double T) const {
    return rho * R0_ * zb_->zbar(rho, T) * T;
}

double IdealElectronEOS::energy(double rho, double T) const {
    return R0_ * zb_->zbar(rho, T) * T / (gamma_ - 1.0);
}

double IdealElectronEOS::cs2Contribution(double rho, double T) const {
    // Ionization derivatives neglected, as in the 1T ideal EOS.
    return gamma_ * pressure(rho, T) / rho;
}

// -------------------------------------------------------------- TableSpeciesEOS

TableSpeciesEOS::TableSpeciesEOS(const std::string& path) : tab_(path, 2) {}

double TableSpeciesEOS::pressure(double rho, double T) const { return tab_.interp(0, rho, T); }
double TableSpeciesEOS::energy(double rho, double T) const   { return tab_.interp(1, rho, T); }
