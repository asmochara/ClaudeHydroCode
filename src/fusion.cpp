#include "fusion.hpp"

#include <algorithm>
#include <cmath>

namespace {
// Bosch-Hale parameterization:
//   theta = T / (1 - T(C2 + T(C4 + T C6)) / (1 + T(C3 + T(C5 + T C7))))
//   xi = (BG^2 / (4 theta))^(1/3)
//   <sv> = C1 * theta * sqrt(xi / (mrc2 T^3)) * exp(-3 xi)
double boschHale(double T, double BG, double mrc2, double C1, double C2,
                 double C3, double C4, double C5, double C6, double C7) {
    if (T < 0.2) return 0.0;
    T = std::min(T, 100.0);
    const double num = T * (C2 + T * (C4 + T * C6));
    const double den = 1.0 + T * (C3 + T * (C5 + T * C7));
    const double theta = T / (1.0 - num / den);
    const double xi = std::cbrt(BG * BG / (4.0 * theta));
    return C1 * theta * std::sqrt(xi / (mrc2 * T * T * T)) * std::exp(-3.0 * xi);
}
}  // namespace

double sigmavDT(double TkeV) {
    return boschHale(TkeV, 34.3827, 1.124656e6, 1.17302e-9, 1.51361e-2,
                     7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5);
}

double sigmavDDn(double TkeV) {
    return boschHale(TkeV, 31.3970, 9.37814e5, 5.43360e-12, 5.85778e-3,
                     7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0);
}

double sigmavDDp(double TkeV) {
    return boschHale(TkeV, 31.3970, 9.37814e5, 5.65718e-12, 3.41267e-3,
                     1.99167e-3, 0.0, 1.05060e-5, 0.0, 0.0);
}
