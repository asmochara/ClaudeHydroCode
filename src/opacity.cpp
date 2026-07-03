#include "opacity.hpp"

TableOpacity::TableOpacity(const std::string& path) : table_(path, 2) {
    // Store both blocks as logs so the interpolation is log-log-log:
    // opacities are strictly positive and span decades, and interpolating
    // the log avoids negative or wildly overshooting values between nodes.
    table_.logBlock(0);   // Rosseland mean
    table_.logBlock(1);   // Planck mean
}

double TableOpacity::rosseland(double density_gcc, double temperature_eV) const {
    return table_.interpLog(0, density_gcc, temperature_eV);
}

double TableOpacity::planck(double density_gcc, double temperature_eV) const {
    return table_.interpLog(1, density_gcc, temperature_eV);
}
