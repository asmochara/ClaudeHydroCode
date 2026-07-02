#!/usr/bin/env python3
"""Generate a Kramers-like grey opacity table for TESTING/DEMO use.

    kappa_R = min(cap, k0 * rho * T^-3.5 + kappa_es)
    kappa_P = min(cap, 30 * k0 * rho * T^-3.5)

with k0 = 2.2e8 * Z^3 / A (free-free scaling) and a Thomson-scattering floor
kappa_es = 0.4*Z/A on the Rosseland mean. This captures the right orders of
magnitude and trends for hydrogenic plasmas but is NOT production data: for
real work, reduce Rosseland/Planck means from the Los Alamos OPLIB tables
(TOPS, https://aphysics2.lanl.gov/apps/) onto this grid format:

    NR NT
    rho grid ascending [g/cc]
    T grid ascending [eV]
    kappa_Rosseland(rho_i, T_j) row-major in rho [cm^2/g]
    kappa_Planck(rho_i, T_j) row-major in rho [cm^2/g]
"""
import argparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("output")
p.add_argument("--Z", type=float, default=1.0)
p.add_argument("--A", type=float, default=1.0)
p.add_argument("--cap", type=float, default=1e6, help="max opacity [cm^2/g]")
p.add_argument("--rho-min", type=float, default=1e-6)
p.add_argument("--rho-max", type=float, default=1e4)
p.add_argument("--T-min", type=float, default=0.01)
p.add_argument("--T-max", type=float, default=1e5)
p.add_argument("--nr", type=int, default=60)
p.add_argument("--nt", type=int, default=80)
args = p.parse_args()


def geomspace(lo, hi, n):
    return [lo * (hi / lo) ** (i / (n - 1)) for i in range(n)]


def writerow(f, vals):
    f.write(" ".join(f"{v:.6e}" for v in vals) + "\n")


rho = geomspace(args.rho_min, args.rho_max, args.nr)
T = geomspace(args.T_min, args.T_max, args.nt)
k0 = 2.2e8 * args.Z**3 / args.A
kes = 0.4 * args.Z / args.A

with open(args.output, "w") as f:
    f.write(f"# TEST Kramers-like grey opacity: Z={args.Z} A={args.A} "
            f"(not production data)\n")
    f.write(f"{args.nr} {args.nt}\n")
    f.write("# rho grid [g/cc]\n")
    writerow(f, rho)
    f.write("# T grid [eV]\n")
    writerow(f, T)
    f.write("# kappa_Rosseland [cm^2/g]\n")
    for r in rho:
        writerow(f, [min(args.cap, k0 * r * t**-3.5 + kes) for t in T])
    f.write("# kappa_Planck [cm^2/g]\n")
    for r in rho:
        writerow(f, [min(args.cap, 30.0 * k0 * r * t**-3.5) for t in T])

print(f"wrote {args.output}: {args.nr}x{args.nt} test opacity table")
