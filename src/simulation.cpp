#include "simulation.hpp"
#include "constants.hpp"

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
}  // namespace

Simulation::Simulation(const InputDeck& deck) : deck_(deck) {
    // Instantiate materials in a deterministic order and remember indices.
    for (const auto& [name, spec] : deck_.materials) {
        Material m;
        m.name = name;
        m.A = spec.A;
        m.Zbar = spec.Z;
        if (spec.eos == "ideal")
            m.eos = std::make_shared<IdealGasEOS>(spec.gamma, spec.A, spec.Z);
        else
            m.eos = std::make_shared<TabulatedEOS>(spec.table_file);
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
    e.assign(M, 0.0);
    T.assign(M, 0.0);
    P.assign(M, 0.0);
    cs.assign(M, 0.0);
    q.assign(M, 0.0);
    matid.assign(M, 0);

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
        for (int k = 0; k < N; ++k, ++j) {
            r[j + 1] = r[j] + w;
            w *= g;
            matid[j] = mi;
            rho[j] = l.rho0;
            const auto& eos = *mats_[mi].eos;
            const double T0 = (l.T0 > 0.0) ? l.T0
                                           : eos.temperatureFromPressure(l.rho0, l.P0);
            T[j] = std::max(T0, deck_.control.T_floor);
            e[j] = eos.energy(rho[j], T[j]);
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

    double Eint = 0.0;
    for (int k = 0; k < M; ++k) Eint += dm[k] * e[k];
    E0 = Eint;  // starts at rest
}

void Simulation::updateEosDerived() {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const auto& eos = *mats_[matid[k]].eos;
        T[k] = eos.temperature(rho[k], e[k], T[k]);
        P[k] = eos.pressure(rho[k], T[k]);
        cs[k] = std::sqrt(eos.soundSpeed2(rho[k], T[k]));
    }
}

void Simulation::applyFloors() {
    const double Tf = deck_.control.T_floor;
    for (int k = 0; k < M; ++k) {
        if (T[k] < Tf) {
            T[k] = Tf;
            e[k] = mats_[matid[k]].eos->energy(rho[k], Tf);
        }
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
    // Node 0 (center or inner wall) is fixed.
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 1; i < M; ++i) {
        const double f = -(P[i] + q[i] - P[i - 1] - q[i - 1]);
        u[i] += dt * area(r[i]) * f / dmNode[i];
    }
    if (c.bc_outer == "pressure") {
        const double pd = deck_.drive.pressure(t);
        u[M] += dt * area(r[M]) * (P[M - 1] + q[M - 1] - pd) / dmNode[M];
    }  // else wall: u[M] stays 0

    // Boundary work done on the system by the applied pressure
    // (force -A*p_drive acting through u[M]*dt).
    if (c.bc_outer == "pressure")
        driveWork += -deck_.drive.pressure(t) * area(r[M]) * u[M] * dt;

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

    // --- energy: de = -(P+q) dV, predictor-corrector on the pressure ------
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const auto& eos = *mats_[matid[k]].eos;
        const double dV = 1.0 / rho[k] - 1.0 / rhoOld[k];  // specific volume change
        const double ePred = e[k] - (P[k] + q[k]) * dV;
        const double Tp = eos.temperature(rho[k], std::max(ePred, 0.0), T[k]);
        const double Pp = eos.pressure(rho[k], std::max(Tp, c.T_floor));
        e[k] -= (0.5 * (P[k] + Pp) + q[k]) * dV;
    }

    // Node masses follow the (fixed) zone masses; recompute only for clarity
    // if zone masses ever become dynamic. dm is constant in Lagrangian frame.

    updateEosDerived();
    applyFloors();
}

double Simulation::coulombLog(double ne, double Te, double Z) const {
    if (deck_.conduction.ln_lambda > 0.0) return deck_.conduction.ln_lambda;
    // NRL formulary electron-ion Coulomb logarithm, floored for cold/dense zones.
    double ll;
    if (Te > 10.0 * Z * Z)
        ll = 24.0 - std::log(std::sqrt(ne) / Te);
    else
        ll = 23.0 - std::log(std::sqrt(ne) * Z * std::pow(Te, -1.5));
    return std::max(ll, 2.0);
}

double Simulation::zoneKappa(int j) const {
    // Spitzer-Harm electron thermal conductivity, kappa in units such that
    // flux = -kappa * d(kB*Te)/dr  [erg/cm^2/s], i.e. kappa in 1/(cm s):
    //   kappa = gamma0(Z) * ne * kB*Te * tau_e / me
    //   tau_e = 3.44e5 * Te[eV]^{3/2} / (ne * lnLambda)   [s]  (NRL formulary)
    const auto& m = mats_[matid[j]];
    const double ne = rho[j] * m.Zbar / (m.A * phys::m_p);
    const double Te = T[j];
    const double lnL = coulombLog(ne, Te, m.Zbar);
    const double tau = 3.44e5 * std::pow(Te, 1.5) / (ne * lnL);
    return braginskiiGamma0(m.Zbar) * ne * (Te * phys::eV) * tau / phys::m_e;
}

void Simulation::conductionStep(double dt) {
    if (!deck_.conduction.enabled || M < 2) return;
    const double f = deck_.conduction.flux_limiter;

    // Zone centers and heat capacities.
    std::vector<double> rc(M), cv(M), kap(M);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        rc[k] = 0.5 * (r[k] + r[k + 1]);
        cv[k] = mats_[matid[k]].eos->cv(rho[k], T[k]);
        kap[k] = zoneKappa(k);
    }

    // Face conductances G_i (erg/s/eV) at interior nodes i = 1..M-1, with
    // Spitzer conductivity harmonically averaged across the face (correct
    // behavior at material interfaces) and a sharp flux limiter:
    //   q = q_SH / (1 + |q_SH| / (f * ne * kB*Te * v_te))
    // evaluated with the beginning-of-step temperature field.
    std::vector<double> G(M + 1, 0.0);
    for (int i = 1; i < M; ++i) {
        const int jl = i - 1, jr = i;
        double kf = 0.0;
        if (kap[jl] > 0.0 && kap[jr] > 0.0)
            kf = 2.0 * kap[jl] * kap[jr] / (kap[jl] + kap[jr]);
        const double drc = rc[jr] - rc[jl];
        const double gradKT = (T[jr] - T[jl]) * phys::eV / drc;  // d(kB T)/dr
        const auto& ml = mats_[matid[jl]];
        const auto& mr = mats_[matid[jr]];
        const double nef = 0.5 * (rho[jl] * ml.Zbar / (ml.A * phys::m_p) +
                                  rho[jr] * mr.Zbar / (mr.A * phys::m_p));
        const double Tf = 0.5 * (T[jl] + T[jr]);
        const double kTf = Tf * phys::eV;
        const double qfs = f * nef * kTf * std::sqrt(kTf / phys::m_e);
        const double qsh = kf * std::abs(gradKT);
        const double keff = (qfs > 0.0) ? kf / (1.0 + qsh / qfs) : kf;
        G[i] = keff * area(r[i]) * phys::eV / drc;
    }

    // Backward-Euler tridiagonal system for T^{n+1} (Thomas algorithm):
    //   (dm*cv/dt + G_i + G_{i+1}) T_j - G_i T_{j-1} - G_{i+1} T_{j+1} = dm*cv/dt T_j^n
    std::vector<double> a(M), b(M), cc(M), d(M);
    for (int k = 0; k < M; ++k) {
        const double diag0 = dm[k] * cv[k] / dt;
        a[k] = -G[k];
        cc[k] = -G[k + 1];
        b[k] = diag0 + G[k] + G[k + 1];
        d[k] = diag0 * T[k];
    }
    for (int k = 1; k < M; ++k) {
        const double w = a[k] / b[k - 1];
        b[k] -= w * cc[k - 1];
        d[k] -= w * d[k - 1];
    }
    std::vector<double> Tn(M);
    Tn[M - 1] = d[M - 1] / b[M - 1];
    for (int k = M - 2; k >= 0; --k) Tn[k] = (d[k] - cc[k] * Tn[k + 1]) / b[k];

    // Update internal energy consistently with the EOS at constant density.
    const double Tf = deck_.control.T_floor;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int k = 0; k < M; ++k) {
        const auto& eos = *mats_[matid[k]].eos;
        const double Tnew = std::max(Tn[k], Tf);
        e[k] += eos.energy(rho[k], Tnew) - eos.energy(rho[k], T[k]);
        T[k] = Tnew;
        P[k] = eos.pressure(rho[k], Tnew);
        cs[k] = std::sqrt(eos.soundSpeed2(rho[k], Tnew));
    }
}

void Simulation::writeSnapshot(int index) const {
    char name[64];
    std::snprintf(name, sizeof(name), "snap_%05d.csv", index);
    std::ofstream out(std::filesystem::path(deck_.output.directory) / name);
    out << "# t = " << t << " s, step = " << step << "\n";
    out << "zone,material,r_left,r_right,r_center,u_left,u_right,rho,T_eV,P,e,cs,q\n";
    out.precision(9);
    for (int k = 0; k < M; ++k) {
        out << k << ',' << mats_[matid[k]].name << ','
            << r[k] << ',' << r[k + 1] << ',' << 0.5 * (r[k] + r[k + 1]) << ','
            << u[k] << ',' << u[k + 1] << ','
            << rho[k] << ',' << T[k] << ',' << P[k] << ','
            << e[k] << ',' << cs[k] << ',' << q[k] << '\n';
    }
}

void Simulation::writeHistoryHeader() {
    histPath_ = (std::filesystem::path(deck_.output.directory) / "history.csv").string();
    std::ofstream out(histPath_);
    out << "step,t,dt,r_outer,u_outer,p_drive,rho_max,T_max,T_center,rhoR,"
           "E_int,E_kin,W_drive,E_err\n";
}

void Simulation::writeHistoryRow() {
    double Eint = 0.0, Ekin = 0.0, rhoR = 0.0, rhomax = 0.0, Tmax = 0.0;
    for (int k = 0; k < M; ++k) {
        Eint += dm[k] * e[k];
        rhoR += rho[k] * (r[k + 1] - r[k]);
        rhomax = std::max(rhomax, rho[k]);
        Tmax = std::max(Tmax, T[k]);
    }
    for (int i = 0; i <= M; ++i) Ekin += 0.5 * dmNode[i] * u[i] * u[i];
    const double Etot = Eint + Ekin;
    const double scale = std::max({std::abs(Etot), std::abs(driveWork), std::abs(E0)});
    const double err = (Etot - E0 - driveWork) / scale;

    std::ofstream out(histPath_, std::ios::app);
    out.precision(9);
    out << step << ',' << t << ',' << dt_ << ',' << r[M] << ',' << u[M] << ','
        << deck_.drive.pressure(t) << ',' << rhomax << ',' << Tmax << ','
        << T[0] << ',' << rhoR << ',' << Eint << ',' << Ekin << ','
        << driveWork << ',' << err << '\n';
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
              << ", conduction " << (deck_.conduction.enabled ? "on" : "off")
              << ", outer BC " << deck_.control.bc_outer << "\n";

    dt_ = deck_.control.dt_init;
    while (t < deck_.control.t_end && step < deck_.control.max_steps) {
        dt_ = std::min(computeDt(), deck_.control.t_end - t);
        if (step == 0) dt_ = std::min(dt_, deck_.control.dt_init);

        hydroStep(dt_);
        conductionStep(dt_);

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
