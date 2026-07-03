// ============================================================================
// table2d.cpp -- the shared (density, temperature) table reader used by the
// tabulated EOS, opacity, and ionization models. See the header for the file
// format. The LOOKUP is bilinear in (ln rho, ln T) -- physical tables span
// many decades in both variables, so log spacing of the lookup coordinates
// is the natural choice. The VALUE is interpolated linearly unless the owner
// converts a block to log storage via logBlock() (used for opacities, which
// themselves span decades and must stay positive).
// ============================================================================

#include "table2d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

Table2D::Table2D(const std::string& path, int blockCount) : path_(path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Table2D: cannot open " + path);

    // Stream every number in the file into one flat list, stripping '#'
    // comments; the format is positional, so order is all that matters.
    std::vector<double> values;
    std::string line;
    while (std::getline(in, line)) {
        auto hashPos = line.find('#');
        if (hashPos != std::string::npos) line.erase(hashPos);
        std::istringstream stream(line);
        double v;
        while (stream >> v) values.push_back(v);
    }
    if (values.size() < 2)
        throw std::runtime_error("Table2D: truncated file " + path);
    const size_t densityCount = static_cast<size_t>(values[0]);
    const size_t temperatureCount = static_cast<size_t>(values[1]);
    const size_t expectedCount =
        2 + densityCount + temperatureCount +
        static_cast<size_t>(blockCount) * densityCount * temperatureCount;
    if (values.size() != expectedCount)
        throw std::runtime_error("Table2D: expected " + std::to_string(expectedCount) +
                                 " numbers in " + path + " (" +
                                 std::to_string(blockCount) + " blocks), got " +
                                 std::to_string(values.size()));

    // Grids are stored as logs; interp() takes logs of its query point once.
    size_t cursor = 2;
    logDensityGrid_.resize(densityCount);
    for (size_t i = 0; i < densityCount; ++i)
        logDensityGrid_[i] = std::log(values[cursor++]);
    logTemperatureGrid_.resize(temperatureCount);
    for (size_t j = 0; j < temperatureCount; ++j)
        logTemperatureGrid_[j] = std::log(values[cursor++]);
    for (int b = 0; b < blockCount; ++b) {
        blocks_.emplace_back(values.begin() + cursor,
                             values.begin() + cursor + densityCount * temperatureCount);
        cursor += densityCount * temperatureCount;
    }

    for (size_t i = 1; i < densityCount; ++i)
        if (logDensityGrid_[i] <= logDensityGrid_[i - 1])
            throw std::runtime_error("Table2D: rho grid not ascending in " + path);
    for (size_t j = 1; j < temperatureCount; ++j)
        if (logTemperatureGrid_[j] <= logTemperatureGrid_[j - 1])
            throw std::runtime_error("Table2D: T grid not ascending in " + path);

    Tmin_ = std::exp(logTemperatureGrid_.front());
    Tmax_ = std::exp(logTemperatureGrid_.back());
}

namespace {
// Locate x in the ascending grid: returns the bracketing lower index and the
// linear interpolation weight in [0, 1]. Queries outside the grid clamp to
// the nearest edge (weight 0 or 1), i.e. constant extrapolation.
inline void locate(const std::vector<double>& grid, double x, size_t& index,
                   double& weight) {
    if (x <= grid.front()) { index = 0; weight = 0.0; return; }
    if (x >= grid.back())  { index = grid.size() - 2; weight = 1.0; return; }
    index = static_cast<size_t>(std::upper_bound(grid.begin(), grid.end(), x) -
                                grid.begin()) - 1;
    weight = (x - grid[index]) / (grid[index + 1] - grid[index]);
}
}  // namespace

double Table2D::interp(int block, double density_gcc, double temperature_eV) const {
    const auto& tableBlock = blocks_[block];
    const size_t temperatureCount = logTemperatureGrid_.size();
    size_t densityIndex, temperatureIndex;
    double densityWeight, temperatureWeight;
    // The tiny floors avoid log(0) on degenerate queries.
    locate(logDensityGrid_, std::log(std::max(density_gcc, 1e-300)), densityIndex,
           densityWeight);
    locate(logTemperatureGrid_, std::log(std::max(temperature_eV, 1e-300)),
           temperatureIndex, temperatureWeight);
    // Standard bilinear stencil on the four surrounding grid points
    // (row-major in density: entry [i * NT + j]).
    const double f00 = tableBlock[densityIndex * temperatureCount + temperatureIndex];
    const double f01 =
        tableBlock[densityIndex * temperatureCount + temperatureIndex + 1];
    const double f10 =
        tableBlock[(densityIndex + 1) * temperatureCount + temperatureIndex];
    const double f11 =
        tableBlock[(densityIndex + 1) * temperatureCount + temperatureIndex + 1];
    return (1 - densityWeight) *
               ((1 - temperatureWeight) * f00 + temperatureWeight * f01) +
           densityWeight * ((1 - temperatureWeight) * f10 + temperatureWeight * f11);
}

void Table2D::logBlock(int block) {
    for (auto& value : blocks_[block]) {
        if (value <= 0.0)
            throw std::runtime_error("Table2D: non-positive entry in log block of " +
                                     path_);
        value = std::log(value);
    }
}

double Table2D::interpLog(int block, double density_gcc, double temperature_eV) const {
    return std::exp(interp(block, density_gcc, temperature_eV));
}
