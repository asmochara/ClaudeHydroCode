#pragma once

#include <algorithm>
#include <cmath>

// ============================================================================
// invertMonotone: solve f(x) = target for x, where f is monotonically
// INCREASING on [lo, hi]. This is the workhorse behind every temperature
// inversion in the code (T from specific energy, T from pressure), so it has
// to be both fast on easy cases and unconditionally robust on hard ones.
//
// Method: safeguarded secant iteration. Each iterate updates a hard bracket
// [lo, hi] around the root; the secant proposal is accepted only if it lands
// strictly inside the bracket, otherwise the step falls back to bisection.
// Consequences:
//   - linear f (e.g. ideal-gas e(T) with fixed Zbar): exact in one step;
//   - smooth nonlinear f (Thomas-Fermi Zbar, Fermi degeneracy, tables):
//     converges superlinearly from the supplied guess, typically 3-8 calls;
//   - pathological f: degrades gracefully to bisection, never diverges.
// Targets at or beyond the bracket ends clamp to the corresponding end.
// ============================================================================
template <class F>
double invertMonotone(F&& f, double target, double lo, double hi, double guess) {
    if (target <= f(lo)) return lo;
    if (target >= f(hi)) return hi;
    double x0 = std::clamp(guess, lo, hi);
    double f0 = f(x0) - target;
    if (f0 > 0.0) hi = x0; else lo = x0;   // tighten the bracket immediately
    // Second point for the first secant step: a small nudge off the guess.
    double x1 = std::clamp(x0 * 1.01 + 1e-14, lo, hi);
    if (x1 == x0) x1 = 0.5 * (lo + hi);
    for (int iteration = 0; iteration < 100; ++iteration) {
        double f1 = f(x1) - target;
        if (f1 > 0.0) hi = x1; else lo = x1;
        double x2;
        if (f1 != f0 && std::isfinite(f1) && std::isfinite(f0))
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0);   // secant proposal
        else
            x2 = 0.5 * (lo + hi);                    // degenerate: bisect
        if (!(x2 > lo && x2 < hi)) x2 = 0.5 * (lo + hi);  // safeguard
        if (std::abs(x2 - x1) <= 1e-12 * std::max(std::abs(x2), 1e-30)) return x2;
        x0 = x1; f0 = f1; x1 = x2;
    }
    return 0.5 * (lo + hi);   // iteration budget exhausted: return midpoint
}
