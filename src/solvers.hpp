#pragma once

#include <vector>

// ============================================================================
// Shared linear-algebra helpers for the implicit physics stages.
// Header-inline so each caller's translation unit can inline them (the
// conduction and radiation solves both use thomasSolve every step).
// ============================================================================

// Solve the tridiagonal linear system
//   lower[k]*x[k-1] + diag[k]*x[k] + upper[k]*x[k+1] = rhs[k]
// by the Thomas algorithm (forward elimination + back substitution).
// diag and rhs are modified in place. O(n), no pivoting -- all our systems
// are diagonally dominant (backward-Euler diffusion matrices), so this is
// unconditionally safe here.
inline void thomasSolve(std::vector<double>& lower, std::vector<double>& diag,
                        std::vector<double>& upper, std::vector<double>& rhs,
                        std::vector<double>& solution) {
    const int n = static_cast<int>(diag.size());
    for (int k = 1; k < n; ++k) {
        const double w = lower[k] / diag[k - 1];
        diag[k] -= w * upper[k - 1];
        rhs[k] -= w * rhs[k - 1];
    }
    solution[n - 1] = rhs[n - 1] / diag[n - 1];
    for (int k = n - 2; k >= 0; --k)
        solution[k] = (rhs[k] - upper[k] * solution[k + 1]) / diag[k];
}
