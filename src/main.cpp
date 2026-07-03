// ============================================================================
// main.cpp -- command-line entry point.
//   hydro1d INPUT_DECK              run a simulation
//   hydro1d --tf Z A rho T          print the Thomas-Fermi Zbar (spot checks
//                                   against ionization references)
//   hydro1d --raytrace INPUT_DECK   trace the laser through the initial state
//                                   and print per-ray turning radii and
//                                   absorbed fractions (refraction checks)
// ============================================================================

#include "input.hpp"
#include "simulation.hpp"
#include "zbar.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    // --- utility: Thomas-Fermi mean ionization at given conditions ---------
    if (argc == 6 && std::string(argv[1]) == "--tf") {
        const double nuclearCharge = std::atof(argv[2]);
        const double atomicMass_amu = std::atof(argv[3]);
        const double density_gcc = std::atof(argv[4]);
        const double temperature_eV = std::atof(argv[5]);
        std::cout << ThomasFermiZbar(nuclearCharge, atomicMass_amu)
                         .zbar(density_gcc, temperature_eV)
                  << "\n";
        return 0;
    }
    // --- utility: laser ray diagnostics on the initial state ---------------
    if (argc == 3 && std::string(argv[1]) == "--raytrace") {
        try {
            Simulation simulation(parseDeck(argv[2]));
            simulation.rayTraceReport();
        } catch (const std::exception& ex) {
            std::cerr << "hydro1d: error: " << ex.what() << "\n";
            return 1;
        }
        return 0;
    }
    if (argc != 2) {
        std::cerr << "usage: hydro1d INPUT_DECK\n"
                  << "       hydro1d --tf Z A rho[g/cc] T[eV]   (print TF Zbar)\n"
                  << "       hydro1d --raytrace INPUT_DECK      (laser ray diagnostics)\n";
        return 2;
    }
    // --- normal run ---------------------------------------------------------
    try {
#ifdef _OPENMP
        std::cout << "hydro1d: OpenMP enabled, max threads = "
                  << omp_get_max_threads() << "\n";
#endif
        Simulation simulation(parseDeck(argv[1]));
        simulation.run();
    } catch (const std::exception& ex) {
        std::cerr << "hydro1d: error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
