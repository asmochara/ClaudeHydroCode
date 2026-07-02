#include "simulation.hpp"
#include "constants.hpp"
#include "numerics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
// Braginskii electron-conduction coefficient vs Z, fit through the tabulated
// values (3.16 at Z=1, approaching the Lorentz limit at high Z).
double braginskiiGamma0(double Z) {
    return 13.58 * (Z + 0.24) / (Z + 4.24);
}

// Levermore-Pomraning flux limiter for radiation diffusion.
double lpLambda(double R) {
    return (2.0 + R) / (6.0 + 3.0 * R + R * R);
}

// Solve the tridiagonal system a x_{k-1} + b x_k + c x_{k+1} = d in place.
void thomasSolve(std::vector<double>& a, std::vector<double>& b,
                 std::vector<double>& c, std::vector<double>& d,
                 std::vector<double>& x) {
    const int n = static_cast<int>(b.size());
    for (int k = 1; k < n; ++k) {
        const double w = a[k] / b[k - 1];
        b[k] -= w * c[k - 1];
        d[k] -= w * d[k - 1];
    }
    x[n - 1] = d[n - 1] / b[n - 1];
    for (int k = n - 2; k >= 0; --k) x[k] = (d[k] - c[k] * x[k + 1]) / b[k];
}
}  // namespace

Simulation::Simulation(const InputDeck& deck) : deck_(deck) {
    twoT_ = (deck_.control.temperatures == 2);
    rad_ = deck_.radiation.enabled;

    // Instantiate materials in a deterministic order and remember indices.
    for (const auto& [name, spec] : deck_.materials) {
        Material m;
        m.name = name;
        m.A = spec.A;
        m.Z = spec.Z;
        m.zbar = makeZbarModel(spec.ionization, spec.Z, spec.A, spec.zbar_table);
        if (twoT_) {
            if (spec.eos == "ideal") {
                m.ion = std::make_shared<IdealIonEOS>(spec.gamma, spec.A);
                m.ele = std::make_shared<IdealElectronEOS>(spec.gamma, spec.A, m.zbar);
            } else {
                m.ion = std::make_shared<TableSpeciesEOS>(spec.table_ion);
                m.ele = std::make_shared<TableSpeciesEOS>(spec.table_electron);
            }
        } else {
            if (spec.eos == "ideal")
                m.eos = std::make_shared<IdealGasEOS>(spec.gamma, spec.A, m.zbar);
            else
                m.eos = std::make_shared<TabulatedEOS>(spec.table_file);
        }
        if (rad_) {
            if (!spec.opacity_table.empty())
                m.opacity = std::make_shared<TableOpacity>(spec.opacity_table);
            else
                m.opacity = std::make_shared<ConstOpacity>(spec.kappa_R, spec.kappa_P);
        }
        mats_.push_back(std::move(m));
    }
    setupMesh();
}

double Simulation::area(double r_) const {
    switch (deck_.control.geometry) {
        case 1: return 1.0;
        case 2: return 2.0 * phys::pi * r_;
        default: return 4.0 * phys::pi * r_ * r_;
    }
}

double Simulation::volume(double r0, double r1) const {
    switch (deck_.control.geometry) {
        case 1: return r1 - r0;
        case 2: return phys::pi * (r1 * r1 - r0 * r0);
        default: return 4.0 / 3.0 * phys::pi * (r1 * r1 * r1 - r0 * r0 * r0);
    }
}

void Simulation::setupMesh() {
    auto matIndex = [&](const std::string& name) {
        for (size_t k = 0; k < mats_.size(); ++k)
            if (mats_[k].name == name) return static_cast<int>(k);
        throw std::runtime_error("unknown material " + name);
    };

    M = 0;
    for (const auto& l : deck_.layers) M += l.zones;
    r.assign(M + 1, 0.0);
    u.assign(M + 1, 0.0);
    dmNode.assign(M + 1, 0.0);
    dm.assign(M, 0.0);
    rho.assign(M, 0.0);
    P.assign(M, 0.0);
    cs.assign(M, 0.0);
    q.assign(M, 0.0);
    zb.assign(M, 0.0);
    matid.assign(M, 0);
    if (twoT_) {
        ei.assign(M, 0.0); ee.assign(M, 0.0);
        Ti.assign(M, 0.0); Te.assign(M, 0.0);
        Pi_.assign(M, 0.0); Pe_.assign(M, 0.0);
    } else {
        e.assign(M, 0.0);
        T.assign(M, 0.0);
    }
    if (rad_) Er.assign(M, 0.0);

    // Node positions: per layer, zone widths follow a geometric progression
    // with (outermost width)/(innermost width) = ratio.
    int j = 0;
    r[0] = deck_.control.r_min;
    for (const auto& l : deck_.layers) {
        const int N = l.zones;
        const double g = (N > 1) ? std::pow(l.ratio, 1.0 / (N - 1)) : 1.0;
        const double w1 = (std::abs(g - 1.0) < 1e-12)
                              ? l.thickness / N
                              : l.thickness * (1.0 - g) / (1.0 - std::pow(g, N));
        double w = w1;
        const int mi = matIndex(l.material);
        const auto& mat = mats_[mi];

        // Base temperature from T0 or from inverting the total pressure.
        double Tbase = l.T0;
        if (Tbase <= 0.0 && l.P0 > 0.0) {
            if (twoT_) {
                Tbase = invertMonotone(
                    [&](double Tx) {
                        return mat.ion->pressure(l.rho0, Tx) +
                               mat.ele->pressure(l.rho0, Tx);
                    },
                    l.P0, 1e-12, 1e9, 1.0);
            } else {
                Tbase = mat.eos->temperatureFromPressure(l.rho0, l.P0);
            }
        }

        for (int k = 0; k < N; ++k, ++j) {
            r[j + 1] = r[j] + w;
            w *= g;
            matid[j] = mi;
            rho[j] = l.rho0;
            const double Tf = deck_.control.T_floor;
            if (twoT_) {
                Ti[j] = std::max((l.Ti0 > 0.0) ? l.Ti0 : Tbase, Tf);
                Te[j] = std::max((l.Te0 > 0.0) ? l.Te0 : Tbase, Tf);
                ei[j] = mat.ion->energy(rho[j], Ti[j]);
                ee[j] = mat.ele->energy(rho[j], Te[j]);
            } else {
                T[j] = std::max(Tbase, Tf);
                e[j] = mat.eos->energy(rho[j], T[j]);
            }
            if (rad_) {
                const double Tr = (l.Tr0 > 0.0) ? l.Tr0
                                                : (twoT_ ? Te[j] : T[j]);
                Er[j] = phys::a_rad * Tr * Tr * Tr * Tr;
            }
        }
    }
    // Recompute exact layer boundaries to avoid width round-off drift.
    {
        int jj = 0;
        double edge = deck_.control.r_min;
        for (const auto& l : deck_.layers) {
            edge += l.thickness;
            jj += l.zones;
            r[jj] = edge;
        }
    }

    for (int k = 0; k < M; ++k) dm[k] = rho[k] * volume(r[k], r[k + 1]);
    for (int i = 1; i < M; ++i) dmNode[i] = 0.5 * (dm[i - 1] + dm[i]);
    dmNode[0] = 0.5 * dm[0];
    dmNode[M] = 0.5 * dm[M - 1];

    updateEosDerived();

    E0 = 0.0;  // starts at rest
    for (int k = 0; k < M; ++k) {
        E0 += dm[k] * (twoT_ ? ei[k] + ee[k] : e[k]);
        if (rad_) E0 += Er[k] * volume(r[k], r[k + 1]);
    }
}

void Simulation::updateEosDerived() {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const auto& m = mats_[matid[k]];
        double cs2;
        if (twoT_) {
            Ti[k] = m.ion->temperature(rho[k], ei[k], Ti[k]);
            Te[k] = m.ele->temperature(rho[k], ee[k], Te[k]);
            Pi_[k] = m.ion->pressure(rho[k], Ti[k]);
            Pe_[k] = m.ele->pressure(rho[k], Te[k]);
            P[k] = Pi_[k] + Pe_[k];
            zb[k] = m.zbar->zbar(rho[k], Te[k]);
            cs2 = m.ion->cs2Contribution(rho[k], Ti[k]) +
                  m.ele->cs2Contribution(rho[k], Te[k]);
        } else {
            T[k] = m.eos->temperature(rho[k], e[k], T[k]);
            P[k] = m.eos->pressure(rho[k], T[k]);
            zb[k] = m.zbar->zbar(rho[k], T[k]);
            cs2 = m.eos->soundSpeed2(rho[k], T[k]);
        }
        if (rad_) cs2 += (4.0 / 9.0) * Er[k] / rho[k];
        cs[k] = std::sqrt(cs2);
    }
}

void Simulation::applyFloors() {
    const double Tf = deck_.control.T_floor;
    for (int k = 0; k < M; ++k) {
        const auto& m = mats_[matid[k]];
        if (twoT_) {
            if (Ti[k] < Tf) { Ti[k] = Tf; ei[k] = m.ion->energy(rho[k], Tf); }
            if (Te[k] < Tf) { Te[k] = Tf; ee[k] = m.ele->energy(rho[k], Tf); }
        } else {
            if (T[k] < Tf) { T[k] = Tf; e[k] = m.eos->energy(rho[k], Tf); }
        }
        if (rad_) Er[k] = std::max(Er[k], phys::a_rad * Tf * Tf * Tf * Tf * 1e-6);
    }
}

double Simulation::computeDt() const {
    double dt = deck_.control.dt_max;
    for (int k = 0; k < M; ++k) {
        const double dr = r[k + 1] - r[k];
        const double du = std::max(0.0, -(u[k + 1] - u[k]));  // compression speed
        const double sig = cs[k] + 4.0 * deck_.control.c_quad * du;
        dt = std::min(dt, deck_.control.cfl * dr / sig);
    }
    if (dt_ > 0.0) dt = std::min(dt, dt_ * deck_.control.dt_growth);
    return dt;
}

void Simulation::hydroStep(double dt) {
    const auto& c = deck_.control;

    // --- momentum: u^{n+1/2} = u^{n-1/2} + dt * a^n -----------------------
    // Total stress = matter pressure + radiation pressure (Er/3) + viscosity.
    // Node 0 (center or inner wall) is fixed.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 1; i < M; ++i) {
        const double f = -(P[i] + radP(i) + q[i] - P[i - 1] - radP(i - 1) - q[i - 1]);
        u[i] += dt * area(r[i]) * f / dmNode[i];
    }
    if (c.bc_outer == "pressure") {
        const double pd = deck_.drive.pressure(t);
        u[M] += dt * area(r[M]) * (P[M - 1] + radP(M - 1) + q[M - 1] - pd) / dmNode[M];
        // Work done on the system by the applied pressure.
        driveWork += -pd * area(r[M]) * u[M] * dt;
    }  // else wall: u[M] stays 0

    // --- move nodes -------------------------------------------------------
    std::vector<double> rhoOld = rho;
    for (int i = 0; i <= M; ++i) r[i] += dt * u[i];
    for (int k = 0; k < M; ++k) {
        if (r[k + 1] <= r[k])
            throw std::runtime_error("mesh tangled at zone " + std::to_string(k) +
                                     ", t = " + std::to_string(t) +
                                     " s. Reduce cfl or increase viscosity.");
        rho[k] = dm[k] / volume(r[k], r[k + 1]);
    }

    // --- artificial viscosity at n+1/2 -------------------------------------
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const double du = u[k + 1] - u[k];
        if (du < 0.0) {
            const double rb = 0.5 * (rho[k] + rhoOld[k]);
            q[k] = rb * (-du) * (c.c_quad * (-du) + c.c_lin * cs[k]);
        } else {
            q[k] = 0.0;
        }
    }

    // --- energy: de = -(P+q) dV per species, predictor-corrector ----------
    // In 2T mode the artificial-viscosity (shock) heating goes to the ions.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const auto& m = mats_[matid[k]];
        const double dV = 1.0 / rho[k] - 1.0 / rhoOld[k];  // specific volume change
        if (twoT_) {
            const double eiP = std::max(ei[k] - (Pi_[k] + q[k]) * dV, 0.0);
            const double TiP = m.ion->temperature(rho[k], eiP, Ti[k]);
            const double PiP = m.ion->pressure(rho[k], std::max(TiP, c.T_floor));
            ei[k] -= (0.5 * (Pi_[k] + PiP) + q[k]) * dV;

            const double eeP = std::max(ee[k] - Pe_[k] * dV, 0.0);
            const double TeP = m.ele->temperature(rho[k], eeP, Te[k]);
            const double PeP = m.ele->pressure(rho[k], std::max(TeP, c.T_floor));
            ee[k] -= 0.5 * (Pe_[k] + PeP) * dV;
        } else {
            const double ePred = std::max(e[k] - (P[k] + q[k]) * dV, 0.0);
            const double Tp = m.eos->temperature(rho[k], ePred, T[k]);
            const double Pp = m.eos->pressure(rho[k], std::max(Tp, c.T_floor));
            e[k] -= (0.5 * (P[k] + Pp) + q[k]) * dV;
        }
        // Radiation compresses adiabatically as a gamma = 4/3 gas.
        if (rad_) Er[k] *= std::pow(rho[k] / rhoOld[k], 4.0 / 3.0);
    }

    updateEosDerived();
    applyFloors();
}

double Simulation::coulombLog(double ne, double Te_, double Z) const {
    if (deck_.conduction.ln_lambda > 0.0) return deck_.conduction.ln_lambda;
    // NRL formulary electron-ion Coulomb logarithm, floored for cold/dense zones.
    double ll;
    if (Te_ > 10.0 * Z * Z)
        ll = 24.0 - std::log(std::sqrt(ne) / Te_);
    else
        ll = 23.0 - std::log(std::sqrt(ne) * Z * std::pow(Te_, -1.5));
    return std::max(ll, 2.0);
}

void Simulation::couplingStep(double dt) {
    if (!twoT_) return;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const auto& m = mats_[matid[k]];
        const double mi = m.A * phys::m_p;
        const double ni = rho[k] / mi;
        const double ne = ni * zb[k];
        const double Zeff = std::max(zb[k], 1.0);
        const double lnL = coulombLog(ne, Te[k], Zeff);
        // NRL formulary electron-ion temperature equilibration rate:
        // dTe/dt = nu (Ti - Te), nu in s^-1 with m in g, T in eV, n in cm^-3.
        const double nu = 1.8e-19 * std::sqrt(phys::m_e * mi) * zb[k] * zb[k] * ni *
                          lnL / std::pow(phys::m_e * Ti[k] + mi * Te[k], 1.5);
        // Backward-Euler pointwise solve of the two-temperature relaxation.
        const double qm = dt * 1.5 * ne * phys::eV * nu / rho[k];  // erg/g/eV
        const double cvi = std::max(m.ion->cv(rho[k], Ti[k]), 1e-30);
        const double cve = std::max(m.ele->cv(rho[k], Te[k]), 1e-30);
        const double dT = (Ti[k] - Te[k]) / (1.0 + qm * (1.0 / cvi + 1.0 / cve));
        const double Q = qm * dT;  // specific energy moved ion -> electron
        ee[k] += Q;
        ei[k] -= Q;
        Te[k] = m.ele->temperature(rho[k], ee[k], Te[k]);
        Ti[k] = m.ion->temperature(rho[k], ei[k], Ti[k]);
    }
}

double Simulation::zoneKappa(int j) const {
    // Spitzer-Harm electron thermal conductivity, kappa in units such that
    // flux = -kappa * d(kB*Te)/dr  [erg/cm^2/s], i.e. kappa in 1/(cm s):
    //   kappa = gamma0(Z) * ne * kB*Te * tau_e / me
    //   tau_e = 3.44e5 * Te[eV]^{3/2} / (ne * lnLambda)   [s]  (NRL formulary)
    const auto& m = mats_[matid[j]];
    const double Tel = twoT_ ? Te[j] : T[j];
    const double ne = rho[j] * zb[j] / (m.A * phys::m_p);
    const double Zeff = std::max(zb[j], 1.0);
    const double lnL = coulombLog(ne, Tel, Zeff);
    const double tau = 3.44e5 * std::pow(Tel, 1.5) / (ne * lnL);
    return braginskiiGamma0(Zeff) * ne * (Tel * phys::eV) * tau / phys::m_e;
}

double Simulation::zoneKappaIon(int j) const {
    // Braginskii ion thermal conductivity (same flux convention as zoneKappa):
    //   kappa_i = 3.9 * ni * kB*Ti * tau_i / mi
    //   tau_i = 2.09e7 * sqrt(mu) * Ti[eV]^{3/2} / (Z^4 * ni * lnLambda)  [s]
    // (NRL formulary ion collision time, mu = mi/mp).
    const auto& m = mats_[matid[j]];
    const double mi = m.A * phys::m_p;
    const double ni = rho[j] / mi;
    const double ne = ni * zb[j];
    const double Zeff = std::max(zb[j], 1.0);
    const double lnL = coulombLog(ne, Ti[j], Zeff);
    const double tau = 2.09e7 * std::sqrt(m.A) * std::pow(Ti[j], 1.5) /
                       (Zeff * Zeff * Zeff * Zeff * ni * lnL);
    return 3.9 * ni * (Ti[j] * phys::eV) * tau / mi;
}

void Simulation::conductionStep(double dt) {
    if (!deck_.conduction.enabled || M < 2) return;
    solveConduction(dt, /*ion=*/false);
    if (twoT_ && deck_.conduction.ion_conduction) solveConduction(dt, /*ion=*/true);
}

// Flux-limited thermal conduction on one temperature field (electrons in 1T
// and 2T mode; ions in 2T mode), backward Euler with a tridiagonal solve.
void Simulation::solveConduction(double dt, bool ion) {
    const double f = ion ? deck_.conduction.ion_flux_limiter
                         : deck_.conduction.flux_limiter;
    std::vector<double>& Tc = ion ? Ti : (twoT_ ? Te : T);
    std::vector<double>& ec = ion ? ei : (twoT_ ? ee : e);
    auto energyOf = [&](int k, double Tx) {
        const auto& m = mats_[matid[k]];
        return ion ? m.ion->energy(rho[k], Tx)
                   : (twoT_ ? m.ele->energy(rho[k], Tx)
                            : m.eos->energy(rho[k], Tx));
    };

    // Zone centers, heat capacities, conductivities, free-streaming fluxes.
    std::vector<double> rc(M), cv(M), kap(M), qf(M);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const auto& m = mats_[matid[k]];
        rc[k] = 0.5 * (r[k] + r[k + 1]);
        cv[k] = std::max(ion ? m.ion->cv(rho[k], Tc[k])
                             : (twoT_ ? m.ele->cv(rho[k], Tc[k])
                                      : m.eos->cv(rho[k], Tc[k])), 1e-30);
        kap[k] = ion ? zoneKappaIon(k) : zoneKappa(k);
        if (ion) {
            const double mi = m.A * phys::m_p;
            const double ni = rho[k] / mi;
            const double kT = Tc[k] * phys::eV;
            qf[k] = f * ni * kT * std::sqrt(kT / mi);
        }
    }

    // Face conductances G_i (erg/s/eV) at interior nodes i = 1..M-1, with the
    // conductivity harmonically averaged across the face (correct behavior at
    // material interfaces) and a sharp flux limiter:
    //   q = q_SH / (1 + |q_SH| / (f * n * kB*T * v_th))
    // evaluated with the beginning-of-step temperature field.
    std::vector<double> G(M + 1, 0.0);
    for (int i = 1; i < M; ++i) {
        const int jl = i - 1, jr = i;
        double kf = 0.0;
        if (kap[jl] > 0.0 && kap[jr] > 0.0)
            kf = 2.0 * kap[jl] * kap[jr] / (kap[jl] + kap[jr]);
        const double drc = rc[jr] - rc[jl];
        const double gradKT = (Tc[jr] - Tc[jl]) * phys::eV / drc;  // d(kB T)/dr
        double qfs;
        if (ion) {
            qfs = 0.5 * (qf[jl] + qf[jr]);  // carrier mass varies by material
        } else {
            const auto& ml = mats_[matid[jl]];
            const auto& mr = mats_[matid[jr]];
            const double nef = 0.5 * (rho[jl] * zb[jl] / (ml.A * phys::m_p) +
                                      rho[jr] * zb[jr] / (mr.A * phys::m_p));
            const double kTf = 0.5 * (Tc[jl] + Tc[jr]) * phys::eV;
            qfs = f * nef * kTf * std::sqrt(kTf / phys::m_e);
        }
        const double qsh = kf * std::abs(gradKT);
        const double keff = (qfs > 0.0) ? kf / (1.0 + qsh / qfs) : kf;
        G[i] = keff * area(r[i]) * phys::eV / drc;
    }

    // Backward-Euler tridiagonal system for T^{n+1}.
    std::vector<double> a(M), b(M), cc(M), d(M), Tn(M);
    for (int k = 0; k < M; ++k) {
        const double diag0 = dm[k] * cv[k] / dt;
        a[k] = -G[k];
        cc[k] = -G[k + 1];
        b[k] = diag0 + G[k] + G[k + 1];
        d[k] = diag0 * Tc[k];
    }
    thomasSolve(a, b, cc, d, Tn);

    // Update internal energy consistently with the EOS at constant density.
    const double Tf = deck_.control.T_floor;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const double Tnew = std::max(Tn[k], Tf);
        ec[k] += energyOf(k, Tnew) - energyOf(k, Tc[k]);
        Tc[k] = Tnew;
    }
}

double Simulation::matterCv(int k) const {
    const auto& m = mats_[matid[k]];
    return std::max(twoT_ ? m.ele->cv(rho[k], Te[k])
                          : m.eos->cv(rho[k], T[k]), 1e-30);
}

void Simulation::radiationStep(double dt) {
    if (!rad_) return;
    std::vector<double>& Tm = twoT_ ? Te : T;   // matter temp coupled to radiation
    std::vector<double>& em = twoT_ ? ee : e;

    // Zone data: volumes, centers, opacities, coupling factors.
    std::vector<double> V(M), rc(M), kR(M), kP(M), fc(M), cvv(M);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const auto& m = mats_[matid[k]];
        V[k] = volume(r[k], r[k + 1]);
        rc[k] = 0.5 * (r[k] + r[k + 1]);
        kR[k] = std::max(m.opacity->rosseland(rho[k], Tm[k]), 1e-10);
        kP[k] = std::max(m.opacity->planck(rho[k], Tm[k]), 0.0);
        cvv[k] = matterCv(k);
        // Linearized-emission reduction factor: with backward-Euler coupling
        // to the matter, the effective source is fc * c*kP*rho*(a Tm^4 - Er).
        const double beta = 4.0 * phys::c_light * kP[k] * dt * phys::a_rad *
                            Tm[k] * Tm[k] * Tm[k];
        fc[k] = cvv[k] / (cvv[k] + beta);
    }

    // Face diffusion conductances (erg s^-1 per erg cm^-3 of Er difference)
    // with the Levermore-Pomraning flux limiter on the lagged field.
    std::vector<double> G(M + 1, 0.0);
    for (int i = 1; i < M; ++i) {
        const int jl = i - 1, jr = i;
        const double drc = rc[jr] - rc[jl];
        const double krf = 0.5 * (kR[jl] * rho[jl] + kR[jr] * rho[jr]);  // 1/cm
        const double Ef = std::max(0.5 * (Er[jl] + Er[jr]), 1e-300);
        const double R = std::abs(Er[jr] - Er[jl]) / (drc * krf * Ef);
        const double D = phys::c_light * lpLambda(R) / krf;  // cm^2/s
        G[i] = D * area(r[i]) / drc;
    }

    // Backward-Euler tridiagonal solve for Er^{n+1}.
    std::vector<double> a(M), b(M), cc(M), d(M), En(M);
    for (int k = 0; k < M; ++k) {
        const double diag0 = V[k] / dt;
        const double S = V[k] * phys::c_light * kP[k] * rho[k] * fc[k];
        const double aT4 = phys::a_rad * Tm[k] * Tm[k] * Tm[k] * Tm[k];
        a[k] = -G[k];
        cc[k] = -G[k + 1];
        b[k] = diag0 + G[k] + G[k + 1] + S;
        d[k] = diag0 * Er[k] + S * aT4;
    }
    // Vacuum (Marshak) leakage through the outer boundary: F = chi * Er with
    // chi interpolating the transparent (c/2) and optically thick limits.
    double chiA = 0.0;
    if (deck_.radiation.bc_outer == "vacuum") {
        const double drl = r[M] - r[M - 1];
        const double chi = 0.5 * phys::c_light /
                           (1.0 + 0.75 * kR[M - 1] * rho[M - 1] * drl);
        chiA = chi * area(r[M]);
        b[M - 1] += chiA;
    }
    thomasSolve(a, b, cc, d, En);

    if (chiA > 0.0) Eleak += chiA * En[M - 1] * dt;

    // Couple absorbed/emitted energy back to the matter (exactly the
    // linearized source used in the solve, so the exchange conserves energy).
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const double aT4 = phys::a_rad * Tm[k] * Tm[k] * Tm[k] * Tm[k];
        const double dTm = phys::c_light * kP[k] * dt * fc[k] * (En[k] - aT4) / cvv[k];
        em[k] += cvv[k] * dTm;
        Tm[k] += dTm;
        Er[k] = En[k];
    }
}

void Simulation::writeSnapshot(int index) const {
    char name[64];
    std::snprintf(name, sizeof(name), "snap_%05d.csv", index);
    std::ofstream out(std::filesystem::path(deck_.output.directory) / name);
    out << "# t = " << t << " s, step = " << step << "\n";
    out << "zone,material,r_left,r_right,r_center,u_left,u_right,rho,"
           "Ti_eV,Te_eV,Tr_eV,zbar,P,e,cs,q\n";
    out.precision(9);
    for (int k = 0; k < M; ++k) {
        const double Tik = twoT_ ? Ti[k] : T[k];
        const double Tek = twoT_ ? Te[k] : T[k];
        const double Trk = rad_ ? std::pow(Er[k] / phys::a_rad, 0.25) : 0.0;
        const double ek = twoT_ ? ei[k] + ee[k] : e[k];
        out << k << ',' << mats_[matid[k]].name << ','
            << r[k] << ',' << r[k + 1] << ',' << 0.5 * (r[k] + r[k + 1]) << ','
            << u[k] << ',' << u[k + 1] << ','
            << rho[k] << ',' << Tik << ',' << Tek << ',' << Trk << ','
            << zb[k] << ',' << P[k] << ',' << ek << ',' << cs[k] << ',' << q[k]
            << '\n';
    }
}

void Simulation::writeHistoryHeader() {
    histPath_ = (std::filesystem::path(deck_.output.directory) / "history.csv").string();
    std::ofstream out(histPath_);
    out << "step,t,dt,r_outer,u_outer,p_drive,rho_max,Ti_max,Te_max,Te_center,"
           "rhoR,E_int,E_kin,E_rad,W_drive,E_leak,E_err\n";
}

void Simulation::writeHistoryRow() {
    double Eint = 0.0, Ekin = 0.0, Erad = 0.0, rhoR = 0.0, rhomax = 0.0;
    double Timax = 0.0, Temax = 0.0;
    for (int k = 0; k < M; ++k) {
        Eint += dm[k] * (twoT_ ? ei[k] + ee[k] : e[k]);
        if (rad_) Erad += Er[k] * volume(r[k], r[k + 1]);
        rhoR += rho[k] * (r[k + 1] - r[k]);
        rhomax = std::max(rhomax, rho[k]);
        Timax = std::max(Timax, twoT_ ? Ti[k] : T[k]);
        Temax = std::max(Temax, twoT_ ? Te[k] : T[k]);
    }
    for (int i = 0; i <= M; ++i) Ekin += 0.5 * dmNode[i] * u[i] * u[i];
    const double Etot = Eint + Ekin + Erad;
    const double scale = std::max({std::abs(Etot), std::abs(driveWork), std::abs(E0), 1e-300});
    const double err = (Etot - E0 - driveWork + Eleak) / scale;

    std::ofstream out(histPath_, std::ios::app);
    out.precision(9);
    out << step << ',' << t << ',' << dt_ << ',' << r[M] << ',' << u[M] << ','
        << deck_.drive.pressure(t) << ',' << rhomax << ',' << Timax << ','
        << Temax << ',' << (twoT_ ? Te[0] : T[0]) << ',' << rhoR << ','
        << Eint << ',' << Ekin << ',' << Erad << ',' << driveWork << ','
        << Eleak << ',' << err << '\n';
}

void Simulation::run() {
    std::filesystem::create_directories(deck_.output.directory);
    writeHistoryHeader();
    writeHistoryRow();
    int snapIndex = 0;
    writeSnapshot(snapIndex++);
    double nextDump = (deck_.output.dt_dump > 0.0) ? deck_.output.dt_dump
                                                   : 2.0 * deck_.control.t_end;

    std::cout << "hydro1d: " << M << " zones, geometry d=" << deck_.control.geometry
              << ", " << (twoT_ ? "2T" : "1T")
              << ", conduction " << (deck_.conduction.enabled ? "on" : "off")
              << ", radiation " << (rad_ ? "on" : "off")
              << ", outer BC " << deck_.control.bc_outer << "\n";

    const bool refresh = deck_.conduction.enabled || twoT_ || rad_;
    dt_ = deck_.control.dt_init;
    while (t < deck_.control.t_end && step < deck_.control.max_steps) {
        dt_ = std::min(computeDt(), deck_.control.t_end - t);
        if (step == 0) dt_ = std::min(dt_, deck_.control.dt_init);

        hydroStep(dt_);
        couplingStep(dt_);
        conductionStep(dt_);
        radiationStep(dt_);
        if (refresh) {
            updateEosDerived();
            applyFloors();
        }

        t += dt_;
        ++step;

        if (step % deck_.output.history_stride == 0) writeHistoryRow();
        if (t >= nextDump - 1e-30) {
            writeSnapshot(snapIndex++);
            nextDump += deck_.output.dt_dump;
        }
        if (step % 20000 == 0)
            std::cout << "  step " << step << "  t = " << t << " s  dt = " << dt_
                      << " s\n";
    }

    writeHistoryRow();
    writeSnapshot(snapIndex);
    std::cout << "hydro1d: done. t = " << t << " s in " << step << " steps. Output in "
              << deck_.output.directory << "/\n";
}
