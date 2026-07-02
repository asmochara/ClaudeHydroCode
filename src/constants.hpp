#pragma once

// Physical constants in CGS units. Temperatures are carried in eV
// throughout the code; multiply by phys::eV to get ergs (i.e. k_B*T).
namespace phys {
constexpr double eV   = 1.602176634e-12;    // erg per eV
constexpr double m_p  = 1.67262192369e-24;  // proton mass [g]
constexpr double m_e  = 9.1093837015e-28;   // electron mass [g]
constexpr double pi   = 3.14159265358979323846;
}
