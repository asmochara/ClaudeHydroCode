#!/usr/bin/env python3
"""Generate a tabulated-EOS file that reproduces an ideal gas.

Useful for validating the table-lookup EOS path against the analytic ideal
EOS, and as a template for writing real (SESAME/LEOS-derived) tables.

Format (see README):
    NR NT
    rho grid (ascending, g/cc)
    T grid (ascending, eV)
    P(rho_i, T_j) row-major in rho  [dyn/cm^2]
    e(rho_i, T_j) row-major in rho  [erg/g]
"""
import argparse
import math

EV = 1.602176634e-12    # erg
MP = 1.67262192369e-24  # g

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("output")
p.add_argument("--gamma", type=float, default=1.4)
p.add_argument("--A", type=float, default=1.0)
p.add_argument("--Z", type=float, default=1.0)
p.add_argument("--rho-min", type=float, default=1e-2)
p.add_argument("--rho-max", type=float, default=1e2)
p.add_argument("--T-min", type=float, default=1e-15)
p.add_argument("--T-max", type=float, default=1e-9)
p.add_argument("--nr", type=int, default=120)
p.add_argument("--nt", type=int, default=120)
args = p.parse_args()


def geomspace(lo, hi, n):
    return [lo * (hi / lo) ** (i / (n - 1)) for i in range(n)]


def writerow(f, vals):
    f.write(" ".join(f"{v:.10e}" for v in vals) + "\n")


rho = geomspace(args.rho_min, args.rho_max, args.nr)
T = geomspace(args.T_min, args.T_max, args.nt)
Rspec = (1.0 + args.Z) * EV / (args.A * MP)

with open(args.output, "w") as f:
    f.write(f"# ideal-gas table: gamma={args.gamma} A={args.A} Z={args.Z}\n")
    f.write(f"{args.nr} {args.nt}\n")
    f.write("# rho grid [g/cc]\n")
    writerow(f, rho)
    f.write("# T grid [eV]\n")
    writerow(f, T)
    f.write("# P [dyn/cm^2]\n")
    for r in rho:
        writerow(f, [r * Rspec * t for t in T])
    f.write("# e [erg/g]\n")
    for _ in rho:
        writerow(f, [Rspec * t / (args.gamma - 1.0) for t in T])

print(f"wrote {args.output}: {args.nr}x{args.nt} ideal-gas table")
