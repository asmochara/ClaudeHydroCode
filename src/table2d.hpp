#pragma once

#include <string>
#include <vector>

// ============================================================================
// Table2D: a rectangular (density, temperature) table with one or more data
// blocks sharing the same grid. This single reader backs the tabulated EOS
// (blocks: P, e), the opacity tables (blocks: kappa_Rosseland, kappa_Planck),
// and the ionization tables (block: Zbar).
//
// ASCII file format ('#' comments allowed anywhere):
//   NR NT
//   density grid, ascending      (NR values, g/cm^3)
//   temperature grid, ascending  (NT values, eV)
//   block 0                      (NR*NT values, row-major in density)
//   block 1 ...
//
// Interpolation is bilinear in (ln rho, ln T); queries outside the grid
// clamp to the edge (constant extrapolation). Call logBlock() to switch a
// block's VALUES to logarithmic storage (interpLog then returns
// exp(interpolated log) -- appropriate for strictly positive quantities
// spanning decades, like opacities).
// ============================================================================
class Table2D {
public:
    Table2D(const std::string& path, int blockCount);
    double interp(int block, double density_gcc, double temperature_eV) const;
    void logBlock(int block);
    double interpLog(int block, double density_gcc, double temperature_eV) const;
    // Temperature range of the grid [eV] -- used to bracket inversions.
    double Tmin() const { return Tmin_; }
    double Tmax() const { return Tmax_; }
private:
    std::vector<double> logDensityGrid_;      // ln of the density grid
    std::vector<double> logTemperatureGrid_;  // ln of the temperature grid
    std::vector<std::vector<double>> blocks_; // row-major in density
    double Tmin_ = 0.0, Tmax_ = 0.0;
    std::string path_;                         // retained for error messages
};
