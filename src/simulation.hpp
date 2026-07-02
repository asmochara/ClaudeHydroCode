#pragma once

#include "eos.hpp"
#include "input.hpp"

#include <string>
#include <vector>

// 1D Lagrangian hydrodynamics on a staggered (von Neumann-Richtmyer) mesh:
// velocities and positions live on nodes, thermodynamic state in zones.
// Geometry: planar, cylindrical, or spherical. Shocks are captured with
// combined linear + quadratic artificial viscosity. Thermal conduction is
// flux-limited Spitzer-Harm, integrated implicitly (backward Euler,
// tridiagonal solve) and operator-split from the hydro update.
class Simulation {
public:
    explicit Simulation(const InputDeck& deck);
    void run();

private:
    // Geometry helpers (d = 1, 2, 3).
    double area(double r) const;
    double volume(double r0, double r1) const;

    void setupMesh();
    double computeDt() const;
    void hydroStep(double dt);
    void conductionStep(double dt);
    void updateEosDerived();
    void applyFloors();

    void writeSnapshot(int index) const;
    void writeHistoryHeader();
    void writeHistoryRow();

    double zoneKappa(int j) const;      // Spitzer conductivity [1/(cm s)]
    double coulombLog(double ne, double T, double Z) const;

    InputDeck deck_;
    std::vector<Material> mats_;

    int M = 0;                    // number of zones (M+1 nodes)
    // Node-centered:
    std::vector<double> r, u;     // position [cm], velocity [cm/s]
    std::vector<double> dmNode;   // node-coupled mass [g per unit solid measure]
    // Zone-centered:
    std::vector<double> dm, rho, e, T, P, cs, q;
    std::vector<int> matid;

    double t = 0.0;
    double dt_ = 0.0;
    long step = 0;
    double driveWork = 0.0;       // cumulative boundary work done on the system
    double E0 = 0.0;              // initial total energy
    std::string histPath_;
};
