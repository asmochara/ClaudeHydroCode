#pragma once

#include <memory>
#include <string>
#include <vector>

// Equation-of-state interface. All quantities in CGS with temperature in eV:
//   rho [g/cm^3], T [eV], P [dyn/cm^2], e [erg/g], cv [erg/g/eV], cs2 [cm^2/s^2]
class EOS {
public:
    virtual ~EOS() = default;
    virtual double pressure(double rho, double T) const = 0;
    virtual double energy(double rho, double T) const = 0;
    virtual double cv(double rho, double T) const = 0;
    // Invert e(rho,T) for T; Tguess accelerates the search.
    virtual double temperature(double rho, double e, double Tguess) const = 0;
    // Invert P(rho,T) for T (used for pressure-specified initial conditions).
    virtual double temperatureFromPressure(double rho, double P) const = 0;
    // Adiabatic sound speed squared.
    virtual double soundSpeed2(double rho, double T) const = 0;
};

// Fully ionized (or fixed-ionization) ideal gas:
//   P = rho * (1+Zbar) * kB*T / (A*m_p),  e = P / ((gamma-1)*rho)
class IdealGasEOS : public EOS {
public:
    IdealGasEOS(double gamma, double A, double Zbar);
    double pressure(double rho, double T) const override;
    double energy(double rho, double T) const override;
    double cv(double rho, double T) const override;
    double temperature(double rho, double e, double Tguess) const override;
    double temperatureFromPressure(double rho, double P) const override;
    double soundSpeed2(double rho, double T) const override;
private:
    double gamma_;
    double Rspec_;  // (1+Zbar)*eV/(A*m_p)  [erg/g/eV]
};

// Tabulated EOS on a rectangular (rho, T) grid, bilinear interpolation in
// (ln rho, ln T). This is the hook for SESAME/LEOS-style tables; the ASCII
// format is documented in the README. e(rho,T) must increase with T.
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
    double interp(const std::vector<double>& tab, double rho, double T) const;
    std::vector<double> lnRho_, lnT_;   // grid (ascending)
    std::vector<double> P_, e_;         // row-major [irho*NT + iT]
    double Tmin_, Tmax_;
};

// A material ties an EOS to the atomic data needed by transport models.
struct Material {
    std::string name;
    double A = 1.0;      // mean atomic mass [amu]
    double Zbar = 1.0;   // mean ionization state
    std::shared_ptr<EOS> eos;
};
