#pragma once

// ============================================================================
// Physical constants, CGS units (CODATA values). Temperatures are carried in
// eV throughout the code: multiply a temperature in eV by phys::eV to get
// the thermal energy kB*T in erg.
// ============================================================================
namespace phys {
constexpr double eV   = 1.602176634e-12;    // erg per eV (kB*T conversion)
constexpr double m_p  = 1.67262192369e-24;  // proton mass [g]
constexpr double m_e  = 9.1093837015e-28;   // electron mass [g]
constexpr double pi   = 3.14159265358979323846;
constexpr double c_light = 2.99792458e10;   // speed of light [cm/s]
constexpr double hbar = 1.054571817e-27;    // reduced Planck constant [erg s]
constexpr double e_esu = 4.80320471257e-10; // electron charge [esu]
// Radiation constant a = 4 sigma_SB / c, expressed per eV^4 so that the
// equilibrium radiation energy density is simply
//   E_rad = a_rad * T[eV]^4   [erg/cm^3]
constexpr double a_rad = 137.20172;
}  // namespace phys
