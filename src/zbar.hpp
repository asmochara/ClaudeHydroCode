#pragma once

#include "table2d.hpp"

#include <memory>
#include <string>

// ============================================================================
// Mean-ionization models: Zbar as a function of density [g/cm^3] and
// electron temperature [eV]. Selected per material in the deck via
// 'ionization = fixed | tf | table'.
// ============================================================================
class ZbarModel {
public:
    virtual ~ZbarModel() = default;
    virtual double zbar(double density_gcc, double electronTemperature_eV) const = 0;
};

// Constant ionization: Zbar = Z always ('ionization = fixed', the default).
// With this choice the deck's Z key is read as the mean ionization itself.
class FixedZbar : public ZbarModel {
public:
    explicit FixedZbar(double fixedCharge) : fixedCharge_(fixedCharge) {}
    double zbar(double, double) const override { return fixedCharge_; }
private:
    double fixedCharge_;
};

// Thomas-Fermi average-atom ionization via R. M. More's fit (1985) -- the
// standard hydrocode TF fit. Inputs are the true nuclear charge Z and atomic
// weight A (use average-atom values for mixtures, e.g. Z=3.5, A=6.5 for CH).
class ThomasFermiZbar : public ZbarModel {
public:
    ThomasFermiZbar(double nuclearCharge, double atomicMass_amu)
        : Z_(nuclearCharge), A_(atomicMass_amu) {}
    double zbar(double density_gcc, double electronTemperature_eV) const override;
private:
    double Z_;   // nuclear charge
    double A_;   // atomic mass [amu]
};

// Tabulated Zbar(rho, T), e.g. reduced from TOPS/OPLIB ionization output.
// Single-block file in the shared Table2D format; values are clamped to
// [1e-4, Z].
class TableZbar : public ZbarModel {
public:
    TableZbar(const std::string& path, double maxCharge);
    double zbar(double density_gcc, double electronTemperature_eV) const override;
private:
    Table2D table_;
    double maxCharge_;   // nuclear charge, the physical ceiling
};

// Factory keyed on the deck's 'ionization' string.
std::shared_ptr<ZbarModel> makeZbarModel(const std::string& kind,
                                         double nuclearCharge,
                                         double atomicMass_amu,
                                         const std::string& tablePath);
