#pragma once

#include "eos.hpp"
#include "input.hpp"

#include <string>
#include <vector>

// 1D Lagrangian hydrodynamics on a staggered (von Neumann-Richtmyer) mesh:
// velocities and positions live on nodes, thermodynamic state in zones.
// Geometry: planar, cylindrical, or spherical. Shocks are captured with
// combined linear + quadratic artificial viscosity.
//
// Optional physics, each operator-split from the hydro update:
//  - two-temperature (ion/electron) model with NRL-formulary e-i coupling,
//    solved pointwise-implicitly (shock heating goes to the ions);
//  - flux-limited Spitzer-Harm electron conduction, backward Euler with a
//    tridiagonal solve;
//  - grey flux-limited radiation diffusion (Levermore-Pomraning limiter)
//    with linearized Planck emission/absorption coupling to the electrons,
//    backward Euler tridiagonal solve; opacities per material.
class Simulation {
public:
    explicit Simulation(const InputDeck& deck);
    void run();
    // Trace the laser through the initial state and print per-ray
    // diagnostics (impact parameter, turning radius, absorbed fraction).
    void rayTraceReport();

private:
    // Geometry helpers (d = 1, 2, 3).
    double area(double r) const;
    double volume(double r0, double r1) const;

    void setupMesh();
    double computeDt() const;
    void hydroStep(double dt);
    void laserStep(double dt);       // ray-traced deposition (dt=0: trace only)
    void couplingStep(double dt);    // electron-ion temperature relaxation (2T)
    void conductionStep(double dt);  // electron (+ ion, in 2T) conduction
    void solveConduction(double dt, bool ion);
    void radiationStep(double dt);   // grey FLD + matter coupling
    void updateEosDerived();
    void applyFloors();

    void writeSnapshot(int index) const;
    void writeHistoryHeader();
    void writeHistoryRow();

    double zoneKappa(int j) const;      // Spitzer electron conductivity [1/(cm s)]
    double zoneKappaIon(int j) const;   // Braginskii ion conductivity [1/(cm s)]
    double coulombLog(double ne, double T, double Z) const;
    double matterCv(int k) const;       // heat capacity coupled to radiation
    double radP(int k) const { return rad_ ? Er[k] / 3.0 : 0.0; }

    InputDeck deck_;
    std::vector<Material> mats_;
    bool twoT_ = false;
    bool rad_ = false;
    bool laser_ = false;
    double ncrit_ = 0.0;          // critical electron density [1/cm^3]

    int M = 0;                    // number of zones (M+1 nodes)
    // Node-centered:
    std::vector<double> r, u;     // position [cm], velocity [cm/s]
    std::vector<double> dmNode;   // node-coupled mass
    // Zone-centered:
    std::vector<double> dm, rho, cs, q;
    std::vector<double> e, T;             // 1T mode: total matter energy/temp
    std::vector<double> ei, ee, Ti, Te;   // 2T mode: per-species
    std::vector<double> P;                // matter pressure (Pi+Pe in 2T)
    std::vector<double> Pi_, Pe_;         // species pressures (2T)
    std::vector<double> zb;               // mean ionization
    std::vector<double> Er;               // radiation energy density [erg/cm^3]
    std::vector<int> matid;

    double t = 0.0;
    double dt_ = 0.0;
    long step = 0;
    double driveWork = 0.0;       // cumulative boundary work done on the system
    double Eleak = 0.0;           // cumulative radiation lost through boundary
    double Elaser = 0.0;          // cumulative laser energy absorbed
    double fabs_ = 0.0;           // instantaneous absorbed laser fraction
    double E0 = 0.0;              // initial total energy
    std::string histPath_;

    struct RayInfo { double b, rmin, fabs; };
    std::vector<RayInfo> rayDiag_;  // per-ray diagnostics from last laserStep
};
