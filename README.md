# hydro1d — 1D Lagrangian hydrodynamics for ICF-style implosions

A compact C++17 Lagrangian hydrodynamics code for shocks and compressible
flow in 1D planar, cylindrical, or spherical geometry, aimed at simulating
inertial-confinement-fusion capsule implosions driven by an applied
(ablation-like) pressure source.

## Physics and numerics

- **Hydro scheme:** staggered-grid (von Neumann–Richtmyer) Lagrangian
  differencing. Velocities/positions live on nodes, thermodynamic state in
  zones; zone masses are fixed. Explicit leapfrog time integration with a
  predictor–corrector (time-centered pressure) internal-energy update, so
  the scheme is second-order-accurate on smooth flow for any EOS.
- **Shock capturing:** combined quadratic (von Neumann–Richtmyer) + linear
  (Landshoff) artificial viscosity, active only in compression
  (`c_quad`, `c_lin` in the deck). In 2T mode the viscous (shock) heating
  goes to the ions, as it should.
- **Geometry:** planar / cylindrical / spherical via the `geometry` key;
  the planar option exists mainly so shock-tube validation problems can be
  run with the same executable.
- **Temperatures:** single-temperature (`temperatures = 1`) or separate
  ion/electron temperatures (`temperatures = 2`). In 2T mode each species
  carries its own energy equation, and the temperatures relax at the
  NRL-formulary electron–ion equilibration rate, integrated
  pointwise-implicitly (unconditionally stable).
- **Equation of state:** per-material, either
  - `ideal`: ideal ion + electron gases with mean ionization Zbar from the
    ionization model (below); in 1T mode
    `P = rho (1+Zbar) kB T / (A m_p)`, in 2T mode the ion and electron
    partial EOS are separate; or
  - `table`: bilinear lookup of `P(rho,T)` and `e(rho,T)` on a rectangular
    grid (interpolated in log rho, log T) with robust `T(rho,e)` inversion
    and a thermodynamically consistent sound speed. In 2T mode give
    separate `table_ion` and `table_electron` files (e.g. SESAME/LEOS
    sub-tables). This is the hook for real tabular EOS data.
- **Ionization:** per-material `ionization = fixed | tf | table`.
  `fixed` uses Zbar = Z; `tf` is the Thomas–Fermi average-atom fit of
  R. M. More (1985) — the standard hydrocode TF fit — using the material's
  nuclear charge Z and atomic weight A (average-atom values for mixtures);
  `table` reads Zbar(rho,T) from a single-block table file, e.g. reduced
  from TOPS/OPLIB ionization output. Zbar feeds the ideal EOS, conduction,
  coupling, and Coulomb logarithms. (`hydro1d --tf Z A rho T` prints the
  TF Zbar for quick checks.)
- **Thermal conduction:** flux-limited Spitzer–Härm electron conduction,
  `kappa_e = gamma0(Zbar) ne kB Te tau_e / me` with the NRL-formulary
  collision time and Coulomb logarithm and
  `gamma0(Z) = 13.58 (Z+0.24)/(Z+4.24)` (a fit through the Braginskii
  coefficients, 3.2 at Z=1). The heat flux is limited against the
  free-streaming flux with the standard sharp limiter
  `q = q_SH / (1 + |q_SH| / (f ne kB Te v_te))`, `f = flux_limiter`
  (default 0.06). In 2T mode, Braginskii **ion conduction**
  (`kappa_i = 3.9 ni kB Ti tau_i / mi`, `ion_flux_limiter` default 0.3) is
  also applied to the ion temperature — important for smoothing the
  converging-shock ion-temperature spike at void closure. Both are
  operator-split, integrated implicitly (backward Euler, tridiagonal
  Thomas solve), with face conductivities harmonically averaged so
  material interfaces behave correctly. Boundaries are insulated.
- **Radiation:** grey flux-limited radiation diffusion
  (Levermore–Pomraning limiter) for the radiation energy density, with
  linearized-Planck emission/absorption coupling to the electron (or 1T
  matter) temperature, solved backward-Euler with a tridiagonal solve; the
  discrete matter–radiation exchange is exactly energy-conserving.
  Radiation pressure (Er/3) enters the momentum equation and Er is
  compressed adiabatically (gamma = 4/3) during the hydro update.
  Opacities per material: constant `kappa_R`/`kappa_P` [cm^2/g] or an
  `opacity_table` of Rosseland and Planck means on the standard grid —
  the intended source is the Los Alamos OPLIB tables via TOPS
  (https://aphysics2.lanl.gov/apps/), reduced to grey means.
  `radiation.bc_outer` is `insulated` or `vacuum` (Marshak leakage,
  tracked in the energy budget).
- **Laser ray tracing:** spherically symmetric direct-drive illumination
  (`[laser]`). Uniform "infinitely many beams" illumination reduces in 1D
  to a bundle of rays sampling the focal spot's impact parameter b,
  distributed according to the spot's radial intensity profile: `profile =
  flattop` (sharp edge at `beam_radius`), `gaussian` / `supergaussian`
  (I ∝ exp(−(b/`beam_radius`)^`sg_order`) with `beam_radius` the
  1/e-intensity radius, truncated at 1e-4 of peak), or `table` (arbitrary
  piecewise-linear radius/intensity pairs — a first radius > 0 gives an
  annular beam). Rays carry equal power by sampling the quantiles of the
  profile's cumulative power distribution.
  Each ray obeys Bouguer's law `mu r sin(theta) = b` in the spherically
  stratified plasma with `mu = sqrt(1 - ne/n_crit)`,
  `n_crit = 1.115e21/lambda_um^2`, so rays refract through the corona,
  turn at their Bouguer radius (or reflect at the critical surface / a
  total-internal-reflection interface), and retrace the mirrored path
  outward. Inverse-bremsstrahlung absorption
  (`kappa = nu_ei (ne/nc) / (c mu)`, NRL collision frequency) attenuates
  each chord and deposits into the electrons; a user-set fraction of the
  power reaching the critical surface is dumped there
  (`absorb_at_critical`, a resonance-absorption stand-in that also
  bootstraps absorption before a corona exists). Inputs: total power vs
  time on target (all beams summed), wavelength, beam radius, ray count.
  Unabsorbed light escapes; absorbed/incident power bookkeeping is exact
  and reported in the history. Use `bc_outer = free` so the corona can
  blow off and the ablation pressure emerges self-consistently — this
  replaces the prescribed pressure drive.
  (`hydro1d --raytrace DECK` prints per-ray turning radii and absorbed
  fractions for the initial state.)
- **Pressure drive (alternative):** a piecewise-linear-in-time applied
  pressure at the outer boundary (`[drive]` table + `bc_outer = pressure`)
  stands in for the ablation pressure when you don't want to model the
  laser. The cumulative drive work is tracked and an energy-conservation
  error is reported in the history file.
- **Units:** CGS everywhere, temperatures in eV
  (1 Mbar = 1e12 dyn/cm^2).

Not included (yet): DT burn/alpha heating, multigroup radiation,
laser-plasma instabilities (the ray trace covers only refraction, inverse
bremsstrahlung, and a critical-surface dump), electron degeneracy in the
ideal EOS (use tables for degenerate/cold-matter regimes).

## Building

Requires CMake ≥ 3.14 and a C++17 compiler. OpenMP is used if found.

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Running

```sh
./build/hydro1d examples/icf_capsule_2t.deck
```

Example decks (all runtimes seconds on one core):

- `examples/sod.deck` — planar Sod shock tube. At t = 0.2 the computed
  contact plateau (0.42633 vs 0.42632 exact), post-shock density
  (0.26552 vs 0.26557), and shock position (0.850 vs 0.8504) match the
  analytic solution to <0.1%; global energy error ~3e-5.
- `examples/sod_table.deck` — the same problem through the tabulated-EOS
  path (generate `examples/ideal_g1.4.eos` first with
  `tools/make_ideal_table.py`); matches the ideal run to interpolation
  error (~0.06% L1 in density).
- `examples/relax_ei.deck` — 2T verification: a static uniform DT plasma
  with Ti = 300 eV, Te = 100 eV relaxes to 200 eV at the NRL rate with
  machine-precision energy conservation.
- `examples/rad_equil.deck` — radiation verification: an optically thick
  slab with Tr = 200 eV, Te = 100 eV relaxes to exact equilibrium
  (Er = a Te^4, Te → 100.18 eV as predicted by the energy budget) with
  machine-precision conservation.
- `examples/icf_capsule.deck` — 1T spherical DT-gas / DT-ice / CH capsule
  driven by a 100 Mbar pressure ramp, conduction on: convergence ratio
  ~35, multi-keV hot spot, peak fuel density ~420 g/cc, peak
  rhoR ≈ 1.7 g/cm².
- `examples/icf_capsule_2t.deck` — the same capsule with the full stack:
  2T, ion + electron conduction, Thomas–Fermi ionization, and radiation
  diffusion with Kramers-like demo opacity tables (generate
  `examples/dt.opac` / `examples/ch.opac` first with
  `tools/make_test_opacity.py`, see the deck header). Shows the expected
  2T signatures: shocked ions run hotter than electrons in flight, a
  brief multi-keV ion flash at void closure smoothed by ion conduction,
  and a hot spot cooled by radiating into the optically thick fuel.
- `examples/raytrace_test.deck` — ray-tracing verification (run with
  `--raytrace`): a uniform sphere with refractive index 0.6; the reported
  turning radii match Bouguer's law (r_turn = b/0.6, total external
  reflection for b > 0.6R) to machine precision.
- `examples/icf_direct_drive.deck` — the same capsule driven by the laser
  instead of an applied pressure (351 nm, ~60 TW peak, `bc_outer = free`).
  Absorption bootstraps from the critical-surface dump on cold solid to
  ~65% inverse bremsstrahlung once the corona forms, then falls as the
  expanding corona refracts rays away from the shrinking target — all
  resolved by the ray trace. The un-optimized pulse plus radiative preheat
  from the crude demo opacities give a modest but complete implosion
  (~155 km/s shell, convergence ratio ~26, ~50 g/cc fuel, ~1.5 keV hot
  spot); pulse-shape and opacity fidelity are left to the user.

Quick-look plotting (requires matplotlib + pandas):

```sh
python3 tools/plot_snapshot.py out_icf_2t/snap_00000.csv out_icf_2t/snap_00021.csv
```

## Input deck format

Plain text, `key = value` inside `[section]` headers, `#` comments.
See the examples for complete decks.

| Section | Keys |
|---|---|
| `[control]` | `t_end`, `dt_init`, `dt_max`, `cfl`, `dt_growth`, `max_steps`, `geometry` (planar/cylindrical/spherical), `temperatures` (1/2), `r_min`, `bc_outer` (wall/pressure/free), `c_quad`, `c_lin`, `T_floor` [eV] |
| `[conduction]` | `enabled`, `flux_limiter`, `ln_lambda` (number or `auto`), `ion_conduction` (2T, default true), `ion_flux_limiter` |
| `[radiation]` | `enabled`, `bc_outer` (insulated/vacuum) |
| `[drive]` | `table = t0 p0 t1 p1 ...` (s, dyn/cm^2), linearly interpolated, end values held |
| `[laser]` | `enabled`, `wavelength_um`, `profile` (flattop/gaussian/supergaussian/table), `beam_radius` (cm), `sg_order`, `profile_table = r0 I0 r1 I1 ...` (cm, relative intensity), `rays`, `absorb_at_critical` (0–1), `power = t0 P0 t1 P1 ...` (s, erg/s; total on target, 1 TW = 1e19 erg/s) |
| `[output]` | `directory`, `dt_dump` (s), `history_stride` |
| `[material NAME]` | `eos` (ideal/table), `gamma`, `A` (amu), `Z` (nuclear charge; the fixed Zbar when `ionization = fixed`), `table` (1T EOS file), `table_ion`/`table_electron` (2T EOS files), `ionization` (fixed/tf/table), `zbar_table`, `opacity_table` **or** `kappa_R` + `kappa_P` (cm^2/g) |
| `[layer]` (repeatable, innermost first) | `material`, `thickness` (cm), `zones`, `rho0` (g/cc), `T0` (eV) **or** `P0` (dyn/cm^2), optional `Ti0`/`Te0` (2T), `Tr0` (radiation), `ratio` (outer/inner zone-width ratio) |

## Table file format

All tables (EOS, opacity, ionization) share one ASCII layout with `#`
comments allowed anywhere:

```
NR NT
rho grid, ascending                 (NR values, g/cc)
T grid, ascending                   (NT values, eV)
block 1 (rho_i, T_j) row-major in rho   (NR*NT values)
block 2 ...
```

- **EOS** (total, ion, or electron): two blocks — P [dyn/cm^2] then
  e [erg/g]; e must increase with T at fixed rho.
- **Opacity**: two blocks — Rosseland then Planck mean [cm^2/g]
  (interpolated log-log-log).
- **Ionization**: one block — Zbar.

`tools/make_ideal_table.py` writes an ideal-gas EOS table (self-test of
the lookup path and a template for SESAME/LEOS conversions);
`tools/make_test_opacity.py` writes a Kramers-like demo opacity table and
documents the intended OPLIB/TOPS workflow.

## Output

- `snap_NNNNN.csv` — zone-by-zone state at each dump time: radii,
  velocities, rho, Ti, Te, Tr (radiation temperature), Zbar, matter
  pressure, specific internal energy, sound speed, artificial viscosity,
  material.
- `history.csv` — per-step time series: outer radius/velocity, drive
  pressure, laser power and instantaneous absorbed fraction, max density,
  max ion/electron temperature, central electron temperature, rhoR,
  internal/kinetic/radiation energy, cumulative drive work, cumulative
  absorbed laser energy, radiation leaked through the boundary, and the
  relative energy-conservation error
  `E_err = (E_int + E_kin + E_rad - E_0 - W_drive - E_laser + E_leak)/E`.

## Parallelism and design notes

The zone loops (EOS evaluation, viscosity, energy update, transport
coefficients) are threaded with OpenMP; the tridiagonal solves are serial
but trivially cheap in 1D. The code is memory-light and, as intended,
entirely comfortable on a single core (the full-physics capsule runs in
tens of seconds).

The physics is deliberately factored — `EOS`/`SpeciesEOS`, `ZbarModel`,
and `Opacity` are abstract interfaces sharing one `Table2D` reader, and
each transport process is a self-contained operator-split stage — so the
natural upgrades slot in without restructuring:

- **LLNL libraries:** the core is dependency-free on purpose (a few
  hundred zones on one core don't benefit), but the loop structure maps
  directly onto [RAJA](https://github.com/LLNL/RAJA) kernels if GPU/many-core
  portability is ever wanted, and stiffer coupled physics (multigroup
  radiation) would be a good fit for
  [SUNDIALS](https://github.com/LLNL/sundials) implicit integrators instead
  of the built-in tridiagonal solves. The table readers are the
  attachment points for LEOS/SESAME EOS data and OPLIB/TOPS opacities.
- **Obvious next physics steps:** DT burn with alpha heating, multigroup
  radiation diffusion, and refined laser coupling (Langdon effect,
  resonance-absorption models, cross-beam energy transfer proxies).
