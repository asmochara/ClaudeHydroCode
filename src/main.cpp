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
    // Utility: print the Thomas-Fermi mean ionization for given conditions.
    if (argc == 6 && std::string(argv[1]) == "--tf") {
        const double Z = std::atof(argv[2]), A = std::atof(argv[3]);
        const double rho = std::atof(argv[4]), T = std::atof(argv[5]);
        std::cout << ThomasFermiZbar(Z, A).zbar(rho, T) << "\n";
        return 0;
    }
    if (argc != 2) {
        std::cerr << "usage: hydro1d INPUT_DECK\n"
                  << "       hydro1d --tf Z A rho[g/cc] T[eV]   (print TF Zbar)\n";
        return 2;
    }
    try {
#ifdef _OPENMP
        std::cout << "hydro1d: OpenMP enabled, max threads = "
                  << omp_get_max_threads() << "\n";
#endif
        Simulation sim(parseDeck(argv[1]));
        sim.run();
    } catch (const std::exception& ex) {
        std::cerr << "hydro1d: error: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
