#pragma once

#include "opacity.hpp"
#include "table2d.hpp"
#include "zbar.hpp"

#include <memory>
#include <string>

// All quantities CGS with temperature in eV:
//   rho [g/cm^3], T [eV], P [dyn/cm^2], e [erg/g], cv [erg/g/eV], cs2 [cm^2/s^2]

// ---------------------------------------------------------------------------
// Single-temperature (total) EOS, used in 1T mode.
class EOS {
public:
    virtual ~EOS() = default;
    virtual double pressure(double rho, double T) const = 0;
    virtual double energy(double rho, double T) const = 0;
    virtual double cv(double rho, double T) const = 0;
    virtual double temperature(double rho, double e, double Tguess) const = 0;
    virtual double temperatureFromPressure(double rho, double P) const = 0;
    virtual double soundSpeed2(double rho, double T) const = 0;
};

// Ideal gas of ions + electrons with a (possibly rho,T-dependent) mean
// ionization: P = rho (1+Zbar(rho,T)) kB T / (A m_p), e = P/((gamma-1) rho).
class IdealGasEOS : public EOS {
public:
    IdealGasEOS(double gamma, double A, std::shared_ptr<const ZbarModel> zb);
    double pressure(double rho, double T) const override;
    double energy(double rho, double T) const override;
    double cv(double rho, double T) const override;
    double temperature(double rho, double e, double Tguess) const override;
    double temperatureFromPressure(double rho, double P) const override;
    double soundSpeed2(double rho, double T) const override;
private:
    double gamma_, R0_;  // R0 = eV/(A m_p) [erg/g/eV]
    std::shared_ptr<const ZbarModel> zb_;
};

// Tabulated total EOS: P and e blocks on a (rho,T) grid (README format).
class TabulatedEOS : public EOS {
public:
    explicit TabulatedEOS(const std::string& path);
    double pressure(double rho, double T) const override;
    double energy(double rho, double T) const override;
    double cv(double rho, double T) const override;
    double temperature(double rho, double e, double Tguess) const override;
    double temperatureFromPressure(double rho, double P) const override;
    double soundSpeed2(double rho, double T) const override;
private:
    Table2D tab_;
};

// ---------------------------------------------------------------------------
// Per-species EOS (ion or electron partial pressure/energy), used in 2T mode.
class SpeciesEOS {
public:
    virtual ~SpeciesEOS() = default;
    virtual double pressure(double rho, double T) const = 0;
    virtual double energy(double rho, double T) const = 0;
    virtual double cv(double rho, double T) const;        // numeric default
    virtual double temperature(double rho, double e, double Tguess) const;
    // This species' contribution to the adiabatic sound speed squared,
    // (dP/drho)_T + T (dP/dT)^2 / (rho^2 cv) by default.
    virtual double cs2Contribution(double rho, double T) const;
protected:
    virtual double Tlo() const { return 1e-12; }
    virtual double Thi() const { return 1e9; }
};

class IdealIonEOS : public SpeciesEOS {
public:
    IdealIonEOS(double gamma, double A);
    double pressure(double rho, double T) const override;
    double energy(double rho, double T) const override;
    double cv(double rho, double T) const override;
    double temperature(double rho, double e, double Tguess) const override;
    double cs2Contribution(double rho, double T) const override;
private:
    double gamma_, R0_;
};

class IdealElectronEOS : public SpeciesEOS {
public:
    IdealElectronEOS(double gamma, double A, std::shared_ptr<const ZbarModel> zb);
    double pressure(double rho, double T) const override;
    double energy(double rho, double T) const override;
    double cs2Contribution(double rho, double T) const override;
private:
    double gamma_, R0_;
    std::shared_ptr<const ZbarModel> zb_;
};

// Tabulated species EOS (e.g. SESAME/LEOS electron or ion sub-tables),
// same two-block file format as the total tabulated EOS.
class TableSpeciesEOS : public SpeciesEOS {
public:
    explicit TableSpeciesEOS(const std::string& path);
    double pressure(double rho, double T) const override;
    double energy(double rho, double T) const override;
protected:
    double Tlo() const override { return tab_.Tmin(); }
    double Thi() const override { return tab_.Tmax(); }
private:
    Table2D tab_;
};

// ---------------------------------------------------------------------------
// A material ties EOS models to the atomic/transport data.
struct Material {
    std::string name;
    double A = 1.0;                       // mean atomic mass [amu]
    double Z = 1.0;                       // nuclear charge (or fixed Zbar)
    std::shared_ptr<ZbarModel> zbar;
    std::shared_ptr<EOS> eos;             // 1T mode
    std::shared_ptr<SpeciesEOS> ion, ele; // 2T mode
    std::shared_ptr<Opacity> opacity;     // radiation (may be null)
};
