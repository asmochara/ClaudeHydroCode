#!/usr/bin/env python3
"""Quick-look plot of hydro1d snapshot files.

  python3 tools/plot_snapshot.py out_icf/snap_00014.csv [more snapshots...]

Plots rho, T, P, and velocity vs radius; requires matplotlib + pandas.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

if len(sys.argv) < 2:
    sys.exit(__doc__)

fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for path in sys.argv[1:]:
    with open(path) as f:
        tline = f.readline()  # "# t = ... s, step = ..."
    label = tline.split(",")[0].lstrip("# ").strip()
    df = pd.read_csv(path, comment="#")
    uc = 0.5 * (df.u_left + df.u_right)
    axes[0, 0].plot(df.r_center, df.rho, label=label)
    axes[0, 1].plot(df.r_center, df.Ti_eV, label=f"Ti {label}")
    axes[0, 1].plot(df.r_center, df.Te_eV, "--", label=f"Te {label}")
    if df.Tr_eV.max() > 0:
        axes[0, 1].plot(df.r_center, df.Tr_eV, ":", label=f"Tr {label}")
    axes[1, 0].plot(df.r_center, df.P, label=label)
    axes[1, 1].plot(df.r_center, uc, label=label)

axes[0, 0].set_ylabel(r"$\rho$ [g/cm$^3$]")
axes[0, 1].set_ylabel(r"$T_i,\ T_e,\ T_r$ [eV]")
axes[1, 0].set_ylabel(r"$P$ [dyn/cm$^2$]")
axes[1, 1].set_ylabel(r"$u$ [cm/s]")
for ax in axes.flat:
    ax.set_xlabel("r [cm]")
    ax.grid(alpha=0.3)
axes[0, 0].legend(fontsize=7)
fig.tight_layout()
out = "snapshot.png"
fig.savefig(out, dpi=140)
print(f"wrote {out}")
