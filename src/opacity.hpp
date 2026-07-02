#pragma once

#include "table2d.hpp"

#include <memory>
#include <string>

// Grey opacities: Rosseland and Planck means in cm^2/g as functions of
// (rho [g/cc], T [eV]).
class Opacity {
public:
    virtual ~Opacity() = default;
    virtual double rosseland(double rho, double T) const = 0;
    virtual double planck(double rho, double T) const = 0;
};

class ConstOpacity : public Opacity {
public:
    ConstOpacity(double kR, double kP) : kR_(kR), kP_(kP) {}
    double rosseland(double, double) const override { return kR_; }
    double planck(double, double) const override { return kP_; }
private:
    double kR_, kP_;
};

// Tabulated grey opacities, e.g. reduced from the Los Alamos OPLIB/TOPS
// astrophysical opacity tables. Two blocks in the standard table format:
// kappa_Rosseland then kappa_Planck, interpolated log-log-log.
class TableOpacity : public Opacity {
public:
    explicit TableOpacity(const std::string& path);
    double rosseland(double rho, double T) const override;
    double planck(double rho, double T) const override;
private:
    Table2D tab_;
};
