// ============================================================================
// zbar.cpp -- mean-ionization models Zbar(density, electron temperature).
// Zbar feeds the ideal EOS (electron count), conduction and coupling rates,
// Coulomb logarithms, and the laser optics (electron density).
// ============================================================================

#include "zbar.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// R. M. More's fit to Thomas-Fermi average-atom ionization (in "Atomic and
// Molecular Physics of Controlled Thermonuclear Fusion", 1985; also
// reproduced in Atzeni & Meyer-ter-Vehn, "The Physics of Inertial Fusion",
// ch. 10). The intermediate variable names (T0, Tf, R, A1, B, C, Q1, Q, x)
// deliberately FOLLOW THE PUBLISHED FIT so the implementation can be checked
// line by line against the reference:
//   T0 = Te / Z^{4/3}            (scaled temperature, eV)
//   R  = rho / (Z A)             (scaled density, g/cc)
//   Zbar/Z = x / (1 + x + sqrt(1 + 2x)),  x = alpha * Q^beta
// with Q built from R and T0 through the fitted coefficients below.
// Known TF caveats: it over-ionizes cold low-density matter and pressure-
// ionizes solids (no shell structure) -- both inherent to the average-atom
// Thomas-Fermi model, not to this fit.
double ThomasFermiZbar::zbar(double density_gcc, double electronTemperature_eV) const {
    constexpr double alpha = 14.3139, beta = 0.6624;
    constexpr double a1 = 0.003323, a2 = 0.9718, a3 = 9.26148e-5, a4 = 3.10165;
    constexpr double b0 = -1.7630, b1 = 1.43175, b2 = 0.31546;
    constexpr double c1 = -0.366667, c2 = 0.983333;

    const double T0 = std::max(electronTemperature_eV, 1e-10) /
                      std::pow(Z_, 4.0 / 3.0);
    const double Tf = T0 / (1.0 + T0);
    const double R = std::max(density_gcc, 1e-30) / (Z_ * A_);
    const double Afac = a1 * std::pow(T0, a2) + a3 * std::pow(T0, a4);
    const double B = -std::exp(b0 + b1 * Tf + b2 * std::pow(Tf, 7));
    const double C = c1 * Tf + c2;
    const double Q1 = Afac * std::pow(R, B);
    const double Q = std::pow(std::pow(R, C) + std::pow(Q1, C), 1.0 / C);
    const double x = alpha * std::pow(Q, beta);
    const double zbarFit = Z_ * x / (1.0 + x + std::sqrt(1.0 + 2.0 * x));
    // The small floor keeps electron-derived quantities (ne, cv_e, coupling
    // rates) well-defined in cold neutral matter; the ceiling is full
    // ionization.
    return std::clamp(zbarFit, 1e-4, Z_);
}

TableZbar::TableZbar(const std::string& path, double maxCharge)
    : table_(path, 1), maxCharge_(maxCharge) {}

double TableZbar::zbar(double density_gcc, double electronTemperature_eV) const {
    return std::clamp(table_.interp(0, density_gcc, electronTemperature_eV), 1e-4,
                      maxCharge_);
}

std::shared_ptr<ZbarModel> makeZbarModel(const std::string& kind,
                                         double nuclearCharge,
                                         double atomicMass_amu,
                                         const std::string& tablePath) {
    if (kind == "fixed") return std::make_shared<FixedZbar>(nuclearCharge);
    if (kind == "tf")
        return std::make_shared<ThomasFermiZbar>(nuclearCharge, atomicMass_amu);
    if (kind == "table") return std::make_shared<TableZbar>(tablePath, nuclearCharge);
    throw std::runtime_error("unknown ionization model '" + kind + "'");
}
