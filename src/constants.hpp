#pragma once

// Physical constants in CGS units. Temperatures are carried in eV
// throughout the code; multiply by phys::eV to get ergs (i.e. k_B*T).
namespace phys {
constexpr double eV   = 1.602176634e-12;    // erg per eV
constexpr double m_p  = 1.67262192369e-24;  // proton mass [g]
constexpr double m_e  = 9.1093837015e-28;   // electron mass [g]
constexpr double pi   = 3.14159265358979323846;
constexpr double c_light = 2.99792458e10;   // [cm/s]
constexpr double hbar = 1.054571817e-27;    // [erg s]
// Radiation constant a = 4 sigma_SB / c expressed per eV^4:
// E_rad(equilibrium) = a_rad * T[eV]^4  [erg/cm^3]
constexpr double a_rad = 137.20172;
}
