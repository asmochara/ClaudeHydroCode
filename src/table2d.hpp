#pragma once

#include <string>
#include <vector>

// Rectangular (rho, T) table with one or more data blocks sharing the same
// grid, loaded from the ASCII format described in the README:
//   NR NT
//   rho grid ascending [g/cc]
//   T grid ascending [eV]
//   block 0 (NR*NT, row-major in rho)
//   block 1 ...
// Interpolation is bilinear in (ln rho, ln T), linear in the block value
// (use logBlock() to switch a block to log-value interpolation).
class Table2D {
public:
    Table2D(const std::string& path, int nblocks);
    double interp(int block, double rho, double T) const;
    // Convert a block to ln(value) storage so interp effectively becomes
    // log-log-log; interpLog returns exp(interp). Requires positive entries.
    void logBlock(int block);
    double interpLog(int block, double rho, double T) const;
    double Tmin() const { return Tmin_; }
    double Tmax() const { return Tmax_; }
private:
    std::vector<double> lnRho_, lnT_;
    std::vector<std::vector<double>> blocks_;
    double Tmin_ = 0.0, Tmax_ = 0.0;
    std::string path_;
};
