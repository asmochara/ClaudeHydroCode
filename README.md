# hydro1d — 1D Lagrangian hydrodynamics for ICF-style implosions

A compact C++17 Lagrangian hydrodynamics code for shocks and compressible
flow in 1D planar, cylindrical, or spherical geometry, aimed at simulating
inertial-confinement-fusion capsule implosions driven by an applied
(ablation-like) pressure source. No radiation transport yet.

## Physics and numerics

- **Hydro scheme:** staggered-grid (von Neumann–Richtmyer) Lagrangian
  differencing. Velocities/positions live on nodes, thermodynamic state in
  zones; zone masses are fixed. Explicit leapfrog time integration with a
  predictor–corrector (time-centered pressure) internal-energy update, so
  the scheme is second-order-accurate on smooth flow for any EOS.
- **Shock capturing:** combined quadratic (von Neumann–Richtmyer) + linear
  (Landshoff) artificial viscosity, active only in compression
  (`c_quad`, `c_lin` in the deck).
- **Geometry:** planar / cylindrical / spherical via the `geometry` key;
  the planar option exists mainly so shock-tube validation problems can be
  run with the same executable.
- **Equation of state:** per-material, either
  - `ideal`: fully/partially ionized ideal gas
    `P = rho (1+Zbar) kB T / (A m_p)`, `e = P / ((gamma-1) rho)`, or
  - `table`: bilinear lookup of `P(rho,T)` and `e(rho,T)` on a rectangular
    `(rho, T)` grid (interpolated in log rho, log T), with Newton/bisection
    inversion for `T(rho,e)` and a thermodynamically consistent sound speed
    `cs^2 = (dP/drho)_T + T (dP/dT)_rho^2 / (rho^2 cv)`.
    This is the hook for SESAME/LEOS-derived tables (format below).
- **Thermal conduction:** single-temperature flux-limited Spitzer–Härm
  electron conduction,
  `kappa = gamma0(Z) ne kB Te tau_e / me` with the NRL-formulary collision
  time and Coulomb logarithm and `gamma0(Z) = 13.58 (Z+0.24)/(Z+4.24)`
  (a fit through the Braginskii coefficients, 3.2 at Z=1). The heat flux is
  limited against the free-streaming flux with the standard sharp limiter
  `q = q_SH / (1 + |q_SH| / (f ne kB Te v_te))`, `f = flux_limiter`
  (default 0.06). Conduction is operator-split from the hydro and
  integrated implicitly (backward Euler, tridiagonal Thomas solve), with
  conductivities harmonically averaged at zone faces so material
  interfaces behave correctly. Boundaries are insulated.
- **Pressure drive:** a piecewise-linear-in-time applied pressure at the
  outer boundary (`[drive]` table + `bc_outer = pressure`) stands in for
  the ablation pressure that starts the confinement. The cumulative drive
  work is tracked and an energy-conservation error is reported in the
  history file.
- **Units:** CGS everywhere, temperatures in eV
  (1 Mbar = 1e12 dyn/cm^2).

## Building

Requires CMake ≥ 3.14 and a C++17 compiler. OpenMP is used if found.

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Running

```sh
./build/hydro1d examples/icf_capsule.deck
```

Two example decks are provided:

- `examples/sod.deck` — planar Sod shock tube. At t = 0.2 the computed
  contact plateau (0.42633 vs 0.42632 exact), post-shock density
  (0.26552 vs 0.26557), and shock position (0.850 vs 0.8504) match the
  analytic solution to <0.1%; global energy error ~3e-5.
- `examples/icf_capsule.deck` — a simplified spherical DT-gas / DT-ice /
  CH capsule driven by a 100 Mbar pressure ramp, conduction on. It
  produces an implosion with convergence ratio ~35, a multi-keV hot spot,
  peak fuel density of several hundred g/cc, and peak ρR ≈ 1.7 g/cm².
- `examples/sod_table.deck` — the same Sod problem run through the
  tabulated-EOS path (generate the table first, see below); it reproduces
  the ideal-gas run to interpolation error (~0.06% L1 in density).

A quick-look plotting script (requires matplotlib + pandas) is included:

```sh
python3 tools/plot_snapshot.py out_icf/snap_00000.csv out_icf/snap_00021.csv
```

## Input deck format

Plain text, `key = value` inside `[section]` headers, `#` comments.
See the examples for complete decks.

| Section | Keys |
|---|---|
| `[control]` | `t_end`, `dt_init`, `dt_max`, `cfl`, `dt_growth`, `max_steps`, `geometry` (planar/cylindrical/spherical), `r_min`, `bc_outer` (wall/pressure), `c_quad`, `c_lin`, `T_floor` [eV] |
| `[conduction]` | `enabled`, `flux_limiter`, `ln_lambda` (number or `auto` for NRL formulary) |
| `[drive]` | `table = t0 p0 t1 p1 ...` (s, dyn/cm^2), linearly interpolated, end values held |
| `[output]` | `directory`, `dt_dump` (s), `history_stride` |
| `[material NAME]` | `eos` (ideal/table), `gamma`, `A` (amu), `Z` (mean ionization), `table` (file, for `eos = table`) |
| `[layer]` (repeatable, innermost first) | `material`, `thickness` (cm), `zones`, `rho0` (g/cc), `T0` (eV) **or** `P0` (dyn/cm^2), `ratio` (outer/inner zone-width ratio for geometric zoning) |

## Tabulated EOS format

ASCII, `#` comments allowed anywhere:

```
NR NT
rho grid, ascending           (NR values, g/cc)
T grid, ascending             (NT values, eV)
P(rho_i, T_j) row-major in rho  (NR*NT values, dyn/cm^2)
e(rho_i, T_j) row-major in rho  (NR*NT values, erg/g)
```

`e` must be monotonically increasing in T at fixed rho.
`tools/make_ideal_table.py` writes an ideal-gas table in this format, both
as a self-test of the lookup path and as a template for converting real
SESAME/LEOS data.

## Output

- `snap_NNNNN.csv` — zone-by-zone state (radii, velocities, rho, T, P, e,
  cs, artificial viscosity, material) at each dump time.
- `history.csv` — per-step time series: outer radius/velocity, drive
  pressure, max density, max/central temperature, ρR, internal and kinetic
  energy, cumulative drive work, and the relative energy-conservation
  error `E_err = (E_int + E_kin - E_0 - W_drive)/E`.

## Parallelism and design notes

The zone loops (EOS evaluation, viscosity, energy update, conduction
coefficients) are threaded with OpenMP; the code is memory-light and, as
intended, entirely comfortable on a single core (the ICF example runs in
under a second). The physics is deliberately factored — `EOS` is an
abstract interface, conduction is a self-contained operator-split stage —
so the natural upgrades slot in without restructuring:

- **LLNL libraries:** the current core is dependency-free on purpose (300
  zones on one core doesn't benefit from them), but the loop structure maps
  directly onto [RAJA](https://github.com/LLNL/RAJA) kernels if GPU/many-core
  portability is ever wanted, and stiffer physics (radiation diffusion,
  multi-group) would be a good fit for
  [SUNDIALS](https://github.com/LLNL/sundials) implicit integrators instead
  of the built-in tridiagonal solve. The `TabulatedEOS` reader is the
  attachment point for LEOS/SESAME table data.
- **Obvious next physics steps:** separate ion/electron temperatures,
  radiation diffusion, ionization models (Thomas–Fermi Zbar), and DT burn.
