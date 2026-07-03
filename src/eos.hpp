#pragma once

#include "opacity.hpp"
#include "table2d.hpp"
#include "zbar.hpp"

#include <memory>
#include <string>

// ============================================================================
// Equation-of-state interfaces.
//
// Two families:
//   EOS         - a single TOTAL-matter EOS, used in single-temperature (1T)
//                 mode where ions, electrons, and their energies are lumped.
//   SpeciesEOS  - the partial EOS of ONE species (ions or electrons), used
//                 in two-temperature (2T) mode where each species carries
//                 its own temperature and energy equation.
//
// UNITS (CGS + eV) for every method in this header:
//   density        [g/cm^3]     temperature   [eV]
//   pressure       [dyn/cm^2]   specificEnergy [erg/g]
//   heat capacity  [erg/g/eV]   sound speed^2  [cm^2/s^2]
//
// The hydro tracks (density, specificEnergy) as its conserved variables and
// inverts temperature from them each step, so every EOS must provide a
// robust temperature(density, energy) inversion. The previous temperature is
// always passed in as the initial guess.
// ============================================================================

// ---------------------------------------------------------------------------
// Single-temperature (total-matter) EOS.
// ---------------------------------------------------------------------------
class EOS {
public:
    virtual ~EOS() = default;
    virtual double pressure(double density_gcc, double temperature_eV) const = 0;
    virtual double energy(double density_gcc, double temperature_eV) const = 0;
    virtual double cv(double density_gcc, double temperature_eV) const = 0;
    // Invert e(rho, T) for T; guessTemperature_eV accelerates convergence.
    virtual double temperature(double density_gcc, double specificEnergy_ergg,
                               double guessTemperature_eV) const = 0;
    // Invert P(rho, T) for T (used for pressure-specified initial conditions).
    virtual double temperatureFromPressure(double density_gcc,
                                           double pressure_dyncm2) const = 0;
    // Adiabatic sound speed squared, for the CFL condition and viscosity.
    virtual double soundSpeed2(double density_gcc, double temperature_eV) const = 0;
};

// Ideal gas of ions plus electrons with a (possibly density/temperature-
// dependent) mean ionization from the attached ZbarModel:
//   P = rho * (1 + Zbar(rho,T)) * kB*T / (A m_p)
//   e = P / ((gamma - 1) rho)
// With a fixed Zbar this is the familiar gamma-law gas; with Thomas-Fermi
// ionization the pressure picks up the T-dependence of Zbar as well.
class IdealGasEOS : public EOS {
public:
    IdealGasEOS(double gamma, double atomicMass_amu,
                std::shared_ptr<const ZbarModel> ionizationModel);
    double pressure(double density_gcc, double temperature_eV) const override;
    double energy(double density_gcc, double temperature_eV) const override;
    double cv(double density_gcc, double temperature_eV) const override;
    double temperature(double density_gcc, double specificEnergy_ergg,
                       double guessTemperature_eV) const override;
    double temperatureFromPressure(double density_gcc,
                                   double pressure_dyncm2) const override;
    double soundSpeed2(double density_gcc, double temperature_eV) const override;
private:
    double gamma_;                  // adiabatic index (> 1)
    double gasConstant_erggeV;      // kB per unit ion mass = eV/(A m_p)
    std::shared_ptr<const ZbarModel> ionizationModel_;
};

// Tabulated total EOS: bilinear lookup of P(rho,T) and e(rho,T) on a
// rectangular grid (two blocks in the shared Table2D format; see README).
// This is the attachment point for SESAME/LEOS-derived data. e(rho,T) must
// increase monotonically with T at fixed rho for the inversion to work.
class TabulatedEOS : public EOS {
public:
    explicit TabulatedEOS(const std::string& path);
    double pressure(double density_gcc, double temperature_eV) const override;
    double energy(double density_gcc, double temperature_eV) const override;
    double cv(double density_gcc, double temperature_eV) const override;
    double temperature(double density_gcc, double specificEnergy_ergg,
                       double guessTemperature_eV) const override;
    double temperatureFromPressure(double density_gcc,
                                   double pressure_dyncm2) const override;
    double soundSpeed2(double density_gcc, double temperature_eV) const override;
private:
    Table2D table_;   // block 0: P [dyn/cm^2], block 1: e [erg/g]
};

// ---------------------------------------------------------------------------
// Per-species EOS (2T mode). The base class supplies generic numeric
// implementations of cv (finite difference of energy), the temperature
// inversion (safeguarded secant), and the sound-speed contribution (via the
// thermodynamic identity cs^2 = (dP/drho)_T + T (dP/dT)^2 / (rho^2 cv)), so
// concrete species only need pressure() and energy().
// ---------------------------------------------------------------------------
class SpeciesEOS {
public:
    virtual ~SpeciesEOS() = default;
    virtual double pressure(double density_gcc, double temperature_eV) const = 0;
    virtual double energy(double density_gcc, double temperature_eV) const = 0;
    virtual double cv(double density_gcc, double temperature_eV) const;
    virtual double temperature(double density_gcc, double specificEnergy_ergg,
                               double guessTemperature_eV) const;
    // This species' additive contribution to the total sound speed squared.
    virtual double cs2Contribution(double density_gcc, double temperature_eV) const;
protected:
    // Temperature bracket for the generic inversion; table-backed species
    // override these with the table's range.
    virtual double Tlo() const { return 1e-12; }
    virtual double Thi() const { return 1e9; }
};

// Classical ideal ion gas: P_i = rho kB T_i / (A m_p).
class IdealIonEOS : public SpeciesEOS {
public:
    IdealIonEOS(double gamma, double atomicMass_amu);
    double pressure(double density_gcc, double temperature_eV) const override;
    double energy(double density_gcc, double temperature_eV) const override;
    double cv(double density_gcc, double temperature_eV) const override;
    double temperature(double density_gcc, double specificEnergy_ergg,
                       double guessTemperature_eV) const override;
    double cs2Contribution(double density_gcc, double temperature_eV) const override;
private:
    double gamma_;
    double gasConstant_erggeV;   // eV/(A m_p) [erg/g/eV]
};

// Electron gas with optional Fermi degeneracy. With degeneracy on, the
// pressure interpolates between the classical and T=0 Fermi limits as
//   P_e = sqrt(P_classical^2 + P_Fermi0^2)
// (exact in both limits, ~15% worst-case error near degeneracy parameter
// theta ~ 1), and the energy uses e_e = (3/2) P_e / rho -- the EXACT
// relation for a nonrelativistic Fermi gas at any degeneracy, which keeps P
// and e thermodynamically consistent without extra machinery. Note there is
// NO cold-curve (chemical bonding) term: uncompressed solids carry their
// sub-Mbar zero-point pressure unbalanced. That is harmless once the drive
// exceeds ~1 Mbar; use tabulated EOS for cold-matter fidelity.
class IdealElectronEOS : public SpeciesEOS {
public:
    IdealElectronEOS(double gamma, double atomicMass_amu,
                     std::shared_ptr<const ZbarModel> ionizationModel,
                     bool includeDegeneracy);
    double pressure(double density_gcc, double temperature_eV) const override;
    double energy(double density_gcc, double temperature_eV) const override;
    double cs2Contribution(double density_gcc, double temperature_eV) const override;
private:
    double gamma_;
    double gasConstant_erggeV;   // eV/(A m_p) [erg/g/eV]
    bool includeDegeneracy_;
    std::shared_ptr<const ZbarModel> ionizationModel_;
};

// Tabulated species EOS (e.g. SESAME/LEOS electron or ion sub-tables),
// same two-block file format as the total tabulated EOS.
class TableSpeciesEOS : public SpeciesEOS {
public:
    explicit TableSpeciesEOS(const std::string& path);
    double pressure(double density_gcc, double temperature_eV) const override;
    double energy(double density_gcc, double temperature_eV) const override;
protected:
    double Tlo() const override { return table_.Tmin(); }
    double Thi() const override { return table_.Tmax(); }
private:
    Table2D table_;   // block 0: P [dyn/cm^2], block 1: e [erg/g]
};

// T=0 electron Fermi pressure at FULL ionization [dyn/cm^2]:
//   P_F0 = (2/5) ne E_F,  ne = rho * Z / (A m_p)
// This is the reference pressure for the ICF fuel adiabat alpha = P/P_F0
// (~2.2 rho^{5/3} Mbar for DT).
double fermiPressure0(double density_gcc, double nuclearCharge, double atomicMass_amu);

// ---------------------------------------------------------------------------
// A material bundles the EOS with the atomic data and transport/ionization/
// opacity models the physics stages need. Built once per deck material.
// ---------------------------------------------------------------------------
struct Material {
    std::string name;
    double A = 1.0;      // mean atomic mass [amu] (average-atom for mixtures)
    double Z = 1.0;      // nuclear charge (equals the fixed Zbar when the
                         // ionization model is "fixed")
    bool fuel = false;   // counts toward the shot report's fuel metrics
    double xD = 0.0;     // deuterium fraction of all ions (burn diagnostics)
    double xT = 0.0;     // tritium fraction of all ions
    std::shared_ptr<ZbarModel> zbar;       // mean-ionization model
    std::shared_ptr<EOS> eos;              // 1T mode only
    std::shared_ptr<SpeciesEOS> ion, ele;  // 2T mode only
    std::shared_ptr<Opacity> opacity;      // radiation only (else null)
};
