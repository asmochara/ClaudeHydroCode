cmake_minimum_required(VERSION 3.14)
project(hydro1d CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(HYDRO1D_USE_OPENMP "Enable OpenMP threading" ON)

add_executable(hydro1d
  src/main.cpp
  src/input.cpp
  src/eos.cpp
  src/simulation.cpp
)

target_compile_options(hydro1d PRIVATE -Wall -Wextra)

if(HYDRO1D_USE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(hydro1d PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()
