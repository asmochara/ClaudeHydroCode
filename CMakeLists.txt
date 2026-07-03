cmake_minimum_required(VERSION 3.14)
project(hydro1d CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(HYDRO1D_USE_OPENMP "Enable OpenMP threading" ON)

# The Simulation class is split across per-physics translation units; see
# the source-layout map at the top of src/simulation.hpp.
add_executable(hydro1d
  src/main.cpp
  src/input.cpp
  src/table2d.cpp
  src/eos.cpp
  src/fusion.cpp
  src/zbar.cpp
  src/opacity.cpp
  src/simulation.cpp
  src/mesh.cpp
  src/hydro.cpp
  src/laser.cpp
  src/coupling.cpp
  src/conduction.cpp
  src/radiation.cpp
  src/burn.cpp
  src/output.cpp
)

# Warning flags per compiler family (MSVC does not understand -Wall/-Wextra).
if(MSVC)
  target_compile_options(hydro1d PRIVATE /W3)
else()
  target_compile_options(hydro1d PRIVATE -Wall -Wextra)
endif()

if(HYDRO1D_USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(hydro1d PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()
