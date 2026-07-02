#pragma once

#include "table2d.hpp"

#include <memory>
#include <string>

// Mean-ionization models: Zbar(rho [g/cc], Te [eV]).
class ZbarModel {
public:
    virtual ~ZbarModel() = default;
    virtual double zbar(double rho, double Te) const = 0;
};

class FixedZbar : public ZbarModel {
public:
    explicit FixedZbar(double Z0) : Z0_(Z0) {}
    double zbar(double, double) const override { return Z0_; }
private:
    double Z0_;
};

// Thomas-Fermi average-atom ionization, R. M. More's fit (1985), as widely
// used in ICF/HEDP hydrocodes. Inputs are the nuclear charge Z and atomic
// weight A (average-atom values for mixtures).
class ThomasFermiZbar : public ZbarModel {
public:
    ThomasFermiZbar(double Z, double A) : Z_(Z), A_(A) {}
    double zbar(double rho, double Te) const override;
private:
    double Z_, A_;
};

// Tabulated Zbar(rho, T), e.g. derived from opacity/ionization tables
// (TOPS/OPLIB output). Single-block file in the standard table format.
class TableZbar : public ZbarModel {
public:
    TableZbar(const std::string& path, double Zmax);
    double zbar(double rho, double Te) const override;
private:
    Table2D tab_;
    double Zmax_;
};

std::shared_ptr<ZbarModel> makeZbarModel(const std::string& kind, double Z, double A,
                                         const std::string& tablePath);
