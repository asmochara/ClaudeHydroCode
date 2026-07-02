#include "zbar.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// R. M. More's fit to Thomas-Fermi average-atom ionization (in "Atomic and
// Molecular Physics of Controlled Thermonuclear Fusion", 1985; also
// Atzeni & Meyer-ter-Vehn, "The Physics of Inertial Fusion", sec. 10).
double ThomasFermiZbar::zbar(double rho, double Te) const {
    constexpr double alpha = 14.3139, beta = 0.6624;
    constexpr double a1 = 0.003323, a2 = 0.9718, a3 = 9.26148e-5, a4 = 3.10165;
    constexpr double b0 = -1.7630, b1 = 1.43175, b2 = 0.31546;
    constexpr double c1 = -0.366667, c2 = 0.983333;

    const double T0 = std::max(Te, 1e-10) / std::pow(Z_, 4.0 / 3.0);
    const double Tf = T0 / (1.0 + T0);
    const double R = std::max(rho, 1e-30) / (Z_ * A_);
    const double Afac = a1 * std::pow(T0, a2) + a3 * std::pow(T0, a4);
    const double B = -std::exp(b0 + b1 * Tf + b2 * std::pow(Tf, 7));
    const double C = c1 * Tf + c2;
    const double Q1 = Afac * std::pow(R, B);
    const double Q = std::pow(std::pow(R, C) + std::pow(Q1, C), 1.0 / C);
    const double x = alpha * std::pow(Q, beta);
    const double zb = Z_ * x / (1.0 + x + std::sqrt(1.0 + 2.0 * x));
    // Small floor keeps electron quantities well-defined in cold matter.
    return std::clamp(zb, 1e-4, Z_);
}

TableZbar::TableZbar(const std::string& path, double Zmax)
    : tab_(path, 1), Zmax_(Zmax) {}

double TableZbar::zbar(double rho, double Te) const {
    return std::clamp(tab_.interp(0, rho, Te), 1e-4, Zmax_);
}

std::shared_ptr<ZbarModel> makeZbarModel(const std::string& kind, double Z, double A,
                                         const std::string& tablePath) {
    if (kind == "fixed") return std::make_shared<FixedZbar>(Z);
    if (kind == "tf") return std::make_shared<ThomasFermiZbar>(Z, A);
    if (kind == "table") return std::make_shared<TableZbar>(tablePath, Z);
    throw std::runtime_error("unknown ionization model '" + kind + "'");
}
