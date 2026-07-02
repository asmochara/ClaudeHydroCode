#include "eos.hpp"
#include "constants.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

// ---------------------------------------------------------------- IdealGasEOS

IdealGasEOS::IdealGasEOS(double gamma, double A, double Zbar)
    : gamma_(gamma), Rspec_((1.0 + Zbar) * phys::eV / (A * phys::m_p)) {
    if (gamma <= 1.0) throw std::runtime_error("IdealGasEOS: gamma must be > 1");
}

double IdealGasEOS::pressure(double rho, double T) const {
    return rho * Rspec_ * T;
}

double IdealGasEOS::energy(double /*rho*/, double T) const {
    return Rspec_ * T / (gamma_ - 1.0);
}

double IdealGasEOS::cv(double /*rho*/, double /*T*/) const {
    return Rspec_ / (gamma_ - 1.0);
}

double IdealGasEOS::temperature(double /*rho*/, double e, double /*Tguess*/) const {
    return e * (gamma_ - 1.0) / Rspec_;
}

double IdealGasEOS::temperatureFromPressure(double rho, double P) const {
    return P / (rho * Rspec_);
}

double IdealGasEOS::soundSpeed2(double rho, double T) const {
    return gamma_ * pressure(rho, T) / rho;
}

// --------------------------------------------------------------- TabulatedEOS

// ASCII format (comment lines starting with '#' allowed anywhere):
//   NR NT
//   rho values, ascending           (NR entries, g/cm^3)
//   T values, ascending             (NT entries, eV)
//   P(rho_i, T_j) row-major in rho  (NR*NT entries, dyn/cm^2)
//   e(rho_i, T_j) row-major in rho  (NR*NT entries, erg/g)
TabulatedEOS::TabulatedEOS(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("TabulatedEOS: cannot open " + path);

    // Stream tokens, skipping comments.
    std::vector<double> vals;
    std::string line;
    while (std::getline(in, line)) {
        auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        std::istringstream ss(line);
        double v;
        while (ss >> v) vals.push_back(v);
    }
    if (vals.size() < 2) throw std::runtime_error("TabulatedEOS: truncated file " + path);
    const size_t NR = static_cast<size_t>(vals[0]);
    const size_t NT = static_cast<size_t>(vals[1]);
    const size_t need = 2 + NR + NT + 2 * NR * NT;
    if (vals.size() != need)
        throw std::runtime_error("TabulatedEOS: expected " + std::to_string(need) +
                                 " numbers in " + path + ", got " + std::to_string(vals.size()));

    size_t k = 2;
    lnRho_.resize(NR);
    for (size_t i = 0; i < NR; ++i) lnRho_[i] = std::log(vals[k++]);
    lnT_.resize(NT);
    for (size_t j = 0; j < NT; ++j) lnT_[j] = std::log(vals[k++]);
    P_.assign(vals.begin() + k, vals.begin() + k + NR * NT);
    k += NR * NT;
    e_.assign(vals.begin() + k, vals.begin() + k + NR * NT);

    for (size_t i = 1; i < NR; ++i)
        if (lnRho_[i] <= lnRho_[i - 1])
            throw std::runtime_error("TabulatedEOS: rho grid not ascending");
    for (size_t j = 1; j < NT; ++j)
        if (lnT_[j] <= lnT_[j - 1])
            throw std::runtime_error("TabulatedEOS: T grid not ascending");

    Tmin_ = std::exp(lnT_.front());
    Tmax_ = std::exp(lnT_.back());
}

namespace {
// Locate x in ascending grid g; return index i with g[i] <= x <= g[i+1]
// (clamped) and the interpolation weight w in [0,1].
inline void locate(const std::vector<double>& g, double x, size_t& i, double& w) {
    if (x <= g.front()) { i = 0; w = 0.0; return; }
    if (x >= g.back())  { i = g.size() - 2; w = 1.0; return; }
    i = static_cast<size_t>(std::upper_bound(g.begin(), g.end(), x) - g.begin()) - 1;
    w = (x - g[i]) / (g[i + 1] - g[i]);
}
}  // namespace

double TabulatedEOS::interp(const std::vector<double>& tab, double rho, double T) const {
    const size_t NT = lnT_.size();
    size_t i, j;
    double wi, wj;
    locate(lnRho_, std::log(std::max(rho, 1e-300)), i, wi);
    locate(lnT_,   std::log(std::max(T,   1e-300)), j, wj);
    const double f00 = tab[i * NT + j],       f01 = tab[i * NT + j + 1];
    const double f10 = tab[(i + 1) * NT + j], f11 = tab[(i + 1) * NT + j + 1];
    return (1 - wi) * ((1 - wj) * f00 + wj * f01) + wi * ((1 - wj) * f10 + wj * f11);
}

double TabulatedEOS::pressure(double rho, double T) const { return interp(P_, rho, T); }
double TabulatedEOS::energy(double rho, double T) const   { return interp(e_, rho, T); }

double TabulatedEOS::cv(double rho, double T) const {
    const double dT = 1e-3 * std::max(T, Tmin_);
    const double Tl = std::max(T - dT, Tmin_);
    const double Th = std::min(T + dT, Tmax_);
    return (energy(rho, Th) - energy(rho, Tl)) / (Th - Tl);
}

double TabulatedEOS::temperature(double rho, double e, double Tguess) const {
    // Bisection on the (assumed monotone) e(rho,T); Newton-polish with cv.
    double lo = Tmin_, hi = Tmax_;
    if (e <= energy(rho, lo)) return lo;
    if (e >= energy(rho, hi)) return hi;
    double T = std::clamp(Tguess, lo, hi);
    for (int it = 0; it < 100; ++it) {
        const double em = energy(rho, T);
        if (em > e) hi = T; else lo = T;
        const double c = cv(rho, T);
        double Tn = (c > 0.0) ? T + (e - em) / c : 0.5 * (lo + hi);
        if (!(Tn > lo && Tn < hi)) Tn = 0.5 * (lo + hi);
        if (std::abs(Tn - T) < 1e-10 * T) return Tn;
        T = Tn;
    }
    return T;
}

double TabulatedEOS::temperatureFromPressure(double rho, double P) const {
    double lo = Tmin_, hi = Tmax_;
    if (P <= pressure(rho, lo)) return lo;
    if (P >= pressure(rho, hi)) return hi;
    for (int it = 0; it < 200; ++it) {
        const double mid = std::sqrt(lo * hi);
        if (pressure(rho, mid) > P) hi = mid; else lo = mid;
        if (hi - lo < 1e-10 * hi) break;
    }
    return std::sqrt(lo * hi);
}

double TabulatedEOS::soundSpeed2(double rho, double T) const {
    // cs^2 = (dP/drho)_T + T/(rho^2 cv) * (dP/dT)_rho^2   (thermodynamic identity)
    const double drho = 1e-3 * rho;
    const double dT = 1e-3 * std::max(T, Tmin_);
    const double dPdrho = (pressure(rho + drho, T) - pressure(std::max(rho - drho, 1e-300), T)) /
                          (rho + drho - std::max(rho - drho, 1e-300));
    const double Tl = std::max(T - dT, Tmin_), Th = std::min(T + dT, Tmax_);
    const double dPdT = (pressure(rho, Th) - pressure(rho, Tl)) / (Th - Tl);
    const double c = std::max(cv(rho, T), 1e-300);
    const double cs2 = dPdrho + T * dPdT * dPdT / (rho * rho * c);
    // Floor to keep the CFL estimate sane on rough tables.
    return std::max(cs2, 1e-6 * pressure(rho, T) / rho);
}
