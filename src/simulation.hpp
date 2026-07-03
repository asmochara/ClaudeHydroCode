#pragma once

#include "eos.hpp"
#include "input.hpp"

#include <string>
#include <vector>

// ============================================================================
// Simulation: the 1D Lagrangian hydrodynamics engine and all attached physics.
//
// MESH LAYOUT (staggered, von Neumann-Richtmyer):
//   The mesh has numZones_ zones and numZones_+1 nodes. Kinematic quantities
//   (position, velocity) live on NODES; thermodynamic quantities (density,
//   temperature, pressure, ...) live on ZONES. Zone k lies between nodes k
//   and k+1. Zone masses are fixed for all time (Lagrangian frame): the mesh
//   moves with the fluid, so there is no advection and material interfaces
//   stay exactly on node boundaries.
//
// TIME INTEGRATION (one call to each per step, in this order):
//   hydroStep       - explicit leapfrog momentum + energy update with
//                     artificial viscosity (the only explicitly time-step-
//                     limited physics; see computeTimeStep)
//   laserStep       - ray-traced laser deposition into the electrons
//   couplingStep    - electron-ion temperature relaxation (2T only),
//                     pointwise implicit
//   conductionStep  - flux-limited Spitzer-Harm electron conduction and
//                     (2T only) Braginskii ion conduction, each backward-
//                     Euler with a tridiagonal solve
//   radiationStep   - grey flux-limited radiation diffusion coupled to the
//                     electrons, backward-Euler tridiagonal solve
//   burnStep        - burn-off fusion diagnostics (no feedback on the hydro)
//
// UNITS: CGS everywhere, with temperatures in eV. Variable names carry the
// unit as a suffix:
//   _cm      centimeters              _s        seconds
//   _g       grams                    _erg      erg
//   _eV/_keV electron-volts           _cmps     cm/s
//   _gcc     g/cm^3                   _dyncm2   dyn/cm^2 (1 Mbar = 1e12)
//   _ergg    erg/g                    _ergcc    erg/cm^3
//   _ergs    erg/s                    _percc    1/cm^3
//   _gcm2    g/cm^2 (areal density)   _ergcm2s  erg/cm^2/s (intensity)
// ============================================================================
class Simulation {
public:
    explicit Simulation(const InputDeck& deck);

    // Advance from t = 0 to control.t_end, writing snapshots, the history
    // file, and the end-of-run shot report.
    void run();

    // Diagnostic mode (CLI --raytrace): trace the laser through the initial
    // state and print each ray's impact parameter, turning radius, and
    // absorbed fraction without advancing the hydro.
    void rayTraceReport();

private:
    // ---- geometry helpers --------------------------------------------------
    // The geometry index d = 1 (planar), 2 (cylindrical), 3 (spherical) sets
    // the area and volume elements. In planar geometry quantities are per
    // unit area; in cylindrical, per unit length.
    double faceArea(double radius_cm) const;                       // [cm^2]
    double shellVolume(double innerRadius_cm, double outerRadius_cm) const;  // [cm^3]

    // ---- setup -------------------------------------------------------------
    void setupMesh();          // build zones from the layer specs, set ICs
    void buildRaySet();        // equal-power ray impact parameters from the
                               // focal-spot intensity profile

    // ---- physics stages (see class comment for ordering) --------------------
    double computeTimeStep() const;
    void hydroStep(double dt_s);
    void laserStep(double dt_s);       // dt_s = 0: trace only, deposit nothing
    void couplingStep(double dt_s);
    void burnStep(double dt_s);
    void conductionStep(double dt_s);
    void solveConduction(double dt_s, bool ionSpecies);
    void radiationStep(double dt_s);

    // Re-derive temperature, pressure, mean ionization, and sound speed from
    // (density, specific energy) so every stage starts from a consistent
    // thermodynamic state.
    void updateThermodynamics();
    // Clamp temperatures (and the radiation field) at the user floor; the
    // energy this injects is tracked in floorEnergyInjected_erg so the
    // energy-conservation diagnostic stays honest.
    void applyTemperatureFloors();

    // ---- transport coefficients ---------------------------------------------
    // Thermal conductivities in units such that flux = -kappa * d(kB*T)/dr,
    // i.e. kappa in 1/(cm*s); see the implementations for the formulas.
    double electronConductivity(int zone) const;
    double ionConductivity(int zone) const;
    // NRL-formulary electron-ion Coulomb logarithm (or the fixed deck value).
    double coulombLog(double electronDensity_percc, double temperature_eV,
                      double ionCharge) const;
    // Heat capacity of whatever the radiation couples to (electrons in 2T,
    // total matter in 1T). [erg/g/eV]
    double matterHeatCapacity(int zone) const;
    // Radiation pressure Er/3 (Eddington closure); zero when radiation is off.
    double radiationPressure_dyncm2(int zone) const {
        return radiationOn_ ? radiationEnergyDensity_ergcc[zone] / 3.0 : 0.0;
    }

    // ---- output -------------------------------------------------------------
    void writeSnapshot(int index) const;
    void writeHistoryHeader();
    void writeHistoryRow();
    void updateReport();       // per-step running design metrics
    void writeReport() const;  // end-of-run shot report (stdout + report.txt)
    // Mass-weighted adiabat alpha = P/P_Fermi of the dense fuel shell.
    double fuelAdiabat() const;

    // ---- configuration ------------------------------------------------------
    InputDeck deck_;
    std::vector<Material> materials_;   // instantiated EOS/transport models
    bool twoTemperature_ = false;       // separate ion/electron temperatures?
    bool radiationOn_ = false;
    bool laserOn_ = false;
    bool burnOn_ = false;

    // ---- laser configuration -------------------------------------------------
    double criticalElectronDensity_percc = 0.0;  // ne where the laser reflects
    std::vector<double> rayImpactParameter_cm;   // one entry per ray
    // Local laser intensity per zone from the previous step's trace (lagged
    // one step), used only by the Langdon absorption correction.
    std::vector<double> laserIntensity_ergcm2s;

    // ---- mesh and state arrays ------------------------------------------------
    int numZones_ = 0;                            // zones; nodes = numZones_+1
    // Node-centered (size numZones_+1):
    std::vector<double> nodeRadius_cm;
    std::vector<double> nodeVelocity_cmps;
    std::vector<double> nodeMass_g;               // half-zone masses coupled to
                                                  // each node (momentum eq.)
    // Zone-centered (size numZones_):
    std::vector<double> zoneMass_g;               // fixed Lagrangian masses
    std::vector<double> zoneDensity_gcc;
    std::vector<double> zoneSoundSpeed_cmps;      // adiabatic, incl. radiation
    std::vector<double> zoneViscousPressure_dyncm2;  // artificial viscosity q
    std::vector<double> zoneSpecificEnergy_ergg;  // 1T mode: total matter
    std::vector<double> zoneTemperature_eV;       // 1T mode
    std::vector<double> ionSpecificEnergy_ergg;   // 2T mode
    std::vector<double> electronSpecificEnergy_ergg;
    std::vector<double> ionTemperature_eV;        // 2T mode
    std::vector<double> electronTemperature_eV;
    std::vector<double> zonePressure_dyncm2;      // total MATTER pressure
                                                  // (ion+electron in 2T)
    std::vector<double> ionPressure_dyncm2;       // 2T mode
    std::vector<double> electronPressure_dyncm2;  // 2T mode
    std::vector<double> zoneMeanIonization;       // Zbar from the material's
                                                  // ionization model
    std::vector<double> radiationEnergyDensity_ergcc;  // grey Er
    std::vector<int> zoneMaterialIndex;           // into materials_

    // ---- time and global energy bookkeeping -----------------------------------
    double time_s = 0.0;
    double timeStep_s = 0.0;
    long stepCount_ = 0;
    // Cumulative energy sources/sinks for the conservation diagnostic
    //   E_err = (E_now - E_initial - work - laser - floors + leak) / scale
    double driveWorkDone_erg = 0.0;      // boundary-pressure work on the system
    double radiationLeaked_erg = 0.0;    // lost through the vacuum radiation BC
    double laserAbsorbed_erg = 0.0;      // deposited by the ray trace
    double laserIncident_erg = 0.0;      // integrated incident laser power
    double floorEnergyInjected_erg = 0.0;  // added by temperature/Er floors
    double initialTotalEnergy_erg = 0.0;
    double laserAbsorbedFraction_ = 0.0;   // instantaneous, for history/report
    std::string historyFilePath_;

    // Per-ray diagnostics captured by the most recent laserStep call.
    struct RayInfo {
        double impactParameter_cm;
        double turningRadius_cm;   // deepest radius reached
        double absorbedFraction;
    };
    std::vector<RayInfo> rayDiag_;

    // ---- burn-off fusion accumulators (diagnostic only) -----------------------
    double fusionPower_ergs = 0.0;       // instantaneous, whole problem
    double peakFusionPower_ergs = 0.0;
    double fusionBangTime_s = -1.0;      // time of peak fusion power
    double neutronYieldDT = 0.0;         // cumulative T(d,n)4He reactions
    double neutronYieldDDn = 0.0;        // cumulative D(d,n)3He reactions
    double fusionEnergy_erg = 0.0;       // cumulative (would-be) fusion energy
    double burnWeightedTiSum_keV = 0.0;  // for burn-averaged <Ti>
    double burnWeightSum = 0.0;

    // ---- shot-report bookkeeping ----------------------------------------------
    std::vector<char> zoneIsFuel_;     // per-zone: material flagged fuel = true
    bool noFuelFlag_ = false;          // nothing flagged; fell back to all zones
    int hotSpotNode_ = 0;              // outer node of the innermost layer
    double hotSpotRadius0_cm = 0.0;    // its initial radius (convergence ratio)
    struct Report {
        double peakImplosionSpeed_cmps = 0.0;   // mass-averaged over fuel
        double peakImplosionTime_s = -1.0;
        double fuelAdiabatAtPeakSpeed = -1.0;   // alpha = P/P_Fermi then
        double peakFuelRhoR_gcm2 = 0.0;         // "bang" proxy without burn
        double bangTime_s = -1.0;               // time of peak fuel rhoR
        double totalRhoRAtBang_gcm2 = 0.0;
        double hotSpotRadius_cm = 0.0;          // all hot-spot values at bang
        double hotSpotTi_eV = 0.0;
        double hotSpotTe_eV = 0.0;
        double hotSpotPressure_dyncm2 = 0.0;
        double hotSpotRhoR_gcm2 = 0.0;
        double minHotSpotRadius_cm = 1e300;     // for the convergence ratio
        double inFlightAspectRatio = -1.0;      // IFAR at 2/3 initial radius
        double ifarTime_s = -1.0;
        double maxDensity_gcc = 0.0;            // global extrema over the run
        double maxIonTemperature_eV = 0.0;
        double maxElectronTemperature_eV = 0.0;
        double maxKineticEnergy_erg = 0.0;
        double peakDrivePressure_dyncm2 = 0.0;
        double peakLaserPower_ergs = 0.0;
    } report_;
};
