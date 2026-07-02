#include "input.hpp"
#include "simulation.hpp"

#include <exception>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "usage: hydro1d INPUT_DECK\n";
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
