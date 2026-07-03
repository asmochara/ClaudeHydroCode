#pragma once

#include "table2d.hpp"

#include <memory>
#include <string>

// ============================================================================
// Grey opacities for the radiation-diffusion package: Rosseland and Planck
// means, both in cm^2/g, as functions of density [g/cm^3] and matter
// temperature [eV]. The Rosseland mean sets the diffusion coefficient
// (transport); the Planck mean sets the emission/absorption coupling rate.
// ============================================================================
class Opacity {
public:
    virtual ~Opacity() = default;
    virtual double rosseland(double density_gcc, double temperature_eV) const = 0;
    virtual double planck(double density_gcc, double temperature_eV) const = 0;
};

// Constant (grey, state-independent) opacities -- mainly for verification
// problems where an analytic answer is wanted.
class ConstOpacity : public Opacity {
public:
    ConstOpacity(double rosseland_cm2g, double planck_cm2g)
        : rosseland_cm2g_(rosseland_cm2g), planck_cm2g_(planck_cm2g) {}
    double rosseland(double, double) const override { return rosseland_cm2g_; }
    double planck(double, double) const override { return planck_cm2g_; }
private:
    double rosseland_cm2g_, planck_cm2g_;
};

// Tabulated grey opacities, e.g. reduced from the Los Alamos OPLIB/TOPS
// astrophysical opacity tables. Two blocks in the shared Table2D format
// (Rosseland then Planck), interpolated in log of the value as well since
// opacities span many decades.
class TableOpacity : public Opacity {
public:
    explicit TableOpacity(const std::string& path);
    double rosseland(double density_gcc, double temperature_eV) const override;
    double planck(double density_gcc, double temperature_eV) const override;
private:
    Table2D table_;
};
