#pragma once

#include <map>
#include <string>
#include <utility>
#include <vector>

// Parsed representation of an input deck. Units: CGS + eV.

struct ControlSpec {
    double t_end = 0.0;          // [s]
    double dt_init = 1e-14;      // [s]
    double dt_max = 1e30;        // [s]
    double cfl = 0.4;            // Courant factor
    double dt_growth = 1.1;      // max dt growth per step
    long   max_steps = 100000000;
    int    geometry = 3;         // 1 planar, 2 cylindrical, 3 spherical
    double r_min = 0.0;          // inner radius [cm]
    std::string bc_outer = "wall";  // "wall" or "pressure" (drive table)
    double c_quad = 2.0;         // quadratic artificial viscosity coefficient
    double c_lin = 0.3;          // linear artificial viscosity coefficient
    double T_floor = 1e-4;       // [eV]
    int temperatures = 1;        // 1 = single-T, 2 = separate ion/electron T
};

struct ConductionSpec {
    bool enabled = false;
    double flux_limiter = 0.06;  // fraction of free-streaming flux (electrons)
    double ln_lambda = -1.0;     // fixed Coulomb log; <0 means compute (NRL)
    bool ion_conduction = true;  // Braginskii ion conduction (2T mode only)
    double ion_flux_limiter = 0.3;  // fraction of ion free-streaming flux
};

struct RadiationSpec {
    bool enabled = false;        // grey flux-limited radiation diffusion
    std::string bc_outer = "insulated";  // "insulated" or "vacuum" (Marshak leak)
};

struct OutputSpec {
    std::string directory = "output";
    double dt_dump = -1.0;       // snapshot interval [s]; <0 = start/end only
    int history_stride = 10;     // history row every N steps
};

struct MaterialSpec {
    std::string name;
    std::string eos = "ideal";   // "ideal" or "table"
    double gamma = 5.0 / 3.0;
    double A = 1.0;              // amu (average-atom for mixtures)
    double Z = 1.0;              // nuclear charge (fixed Zbar if ionization=fixed)
    std::string table_file;      // total EOS table (1T mode)
    std::string table_ion;       // ion EOS table (2T mode, eos=table)
    std::string table_electron;  // electron EOS table (2T mode, eos=table)
    std::string ionization = "fixed";  // "fixed" | "tf" (Thomas-Fermi) | "table"
    std::string zbar_table;      // Zbar(rho,T) table (ionization=table)
    std::string opacity_table;   // kappa_R,kappa_P table; else constants below
    double kappa_R = -1.0;       // [cm^2/g] constant Rosseland mean
    double kappa_P = -1.0;       // [cm^2/g] constant Planck mean
};

struct LayerSpec {
    std::string material;
    double thickness = 0.0;      // [cm]
    int zones = 0;
    double rho0 = 0.0;           // [g/cm^3]
    double T0 = -1.0;            // [eV]        (give T0 or P0)
    double P0 = -1.0;            // [dyn/cm^2]
    double Ti0 = -1.0;           // [eV] override ion T (2T mode)
    double Te0 = -1.0;           // [eV] override electron T (2T mode)
    double Tr0 = -1.0;           // [eV] initial radiation temperature
    double ratio = 1.0;          // outermost/innermost zone width (geometric)
};

struct DriveSpec {
    // Piecewise-linear applied pressure at the outer boundary:
    // (time [s], pressure [dyn/cm^2]) pairs, ascending in time.
    std::vector<std::pair<double, double>> table;
    double pressure(double t) const;
    bool empty() const { return table.empty(); }
};

struct InputDeck {
    ControlSpec control;
    ConductionSpec conduction;
    RadiationSpec radiation;
    OutputSpec output;
    std::map<std::string, MaterialSpec> materials;
    std::vector<LayerSpec> layers;
    DriveSpec drive;
};

// Parse the deck at `path`; throws std::runtime_error with a useful message
// on malformed input.
InputDeck parseDeck(const std::string& path);
