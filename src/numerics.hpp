#pragma once

#include <algorithm>
#include <cmath>

// Invert a monotonically increasing scalar function f on [lo, hi]:
// returns x with f(x) = target. Safeguarded secant iteration (falls back to
// bisection), so it is exact in one step for linear f and robust otherwise.
template <class F>
double invertMonotone(F&& f, double target, double lo, double hi, double guess) {
    if (target <= f(lo)) return lo;
    if (target >= f(hi)) return hi;
    double x0 = std::clamp(guess, lo, hi);
    double f0 = f(x0) - target;
    if (f0 > 0.0) hi = x0; else lo = x0;
    double x1 = std::clamp(x0 * 1.01 + 1e-14, lo, hi);
    if (x1 == x0) x1 = 0.5 * (lo + hi);
    for (int it = 0; it < 100; ++it) {
        double f1 = f(x1) - target;
        if (f1 > 0.0) hi = x1; else lo = x1;
        double x2;
        if (f1 != f0 && std::isfinite(f1) && std::isfinite(f0))
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        else
            x2 = 0.5 * (lo + hi);
        if (!(x2 > lo && x2 < hi)) x2 = 0.5 * (lo + hi);
        if (std::abs(x2 - x1) <= 1e-12 * std::max(std::abs(x2), 1e-30)) return x2;
        x0 = x1; f0 = f1; x1 = x2;
    }
    return 0.5 * (lo + hi);
}
