#include "opacity.hpp"

TableOpacity::TableOpacity(const std::string& path) : tab_(path, 2) {
    tab_.logBlock(0);
    tab_.logBlock(1);
}

double TableOpacity::rosseland(double rho, double T) const {
    return tab_.interpLog(0, rho, T);
}

double TableOpacity::planck(double rho, double T) const {
    return tab_.interpLog(1, rho, T);
}
