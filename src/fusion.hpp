#pragma once

// Bosch & Hale (Nucl. Fusion 32, 611 (1992)) Maxwellian reactivities
// <sigma*v> in cm^3/s, T in keV. Fits are valid for 0.2-100 keV; inputs are
// clamped to that range (reactivity below 0.2 keV is treated as zero).
double sigmavDT(double TkeV);    // T(d,n)4He
double sigmavDDn(double TkeV);   // D(d,n)3He
double sigmavDDp(double TkeV);   // D(d,p)T

// Total reaction Q-values [erg] (all products).
namespace fusion {
constexpr double MeV = 1.602176634e-6;      // erg
constexpr double Q_DT = 17.589 * MeV;       // alpha 3.52 + n 14.07 MeV
constexpr double Q_DDn = 3.268 * MeV;       // 3He 0.82 + n 2.45 MeV
constexpr double Q_DDp = 4.032 * MeV;       // T 1.01 + p 3.02 MeV
}
