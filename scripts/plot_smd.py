#!/usr/bin/env python3
"""
plot_smd.py

Plots a steered QM/MM MD run of CP2K: the collective variable (CV) against the moving
restraint target, the work done by the restraint and, if the .ener file is given, the
potential energy (the restraint energy is written with --csv). The panels share the time axis (one quantity per panel).

Inputs
  - PROJECT-COLVAR.metadynLog : column 1 is the time (fs), column 2 the CV in atomic units
    (bohr for distances), written when &FREE_ENERGY/&METADYN/&PRINT/&COLVAR is present
  - PROJECT-1.ener (optional)   : Step, Time[fs], Kin[a.u.], Temp[K], Pot[a.u.], Cons Qty[a.u.], UsedTime[s]

The target follows CP2K's &COLLECTIVE section: TARGET + TARGET_GROWTH*t, stopped at TARGET_LIMIT.
The restraint energy is K*(CV - target)^2 (the CP2K form, without 1/2) and the work done by the
moving restraint is W(t) = sum of -2K*(CV - target)*d(target).

Usage examples:
    python3 plot_smd.py -c SMD-NVT-EE-GPW-COLVAR.metadynLog
    python3 plot_smd.py -c SMD-NVT-EE-GPW-COLVAR.metadynLog -e SMD-NVT-EE-GPW-1.ener \\
        --target -2.0 --growth 0.002 --limit 2.0 -k 100 -o smd.png
    python3 plot_smd.py -c SMD-NVT-EE-GPW-COLVAR.metadynLog --no-show --csv smd.csv
"""

import argparse
import sys

import numpy as np

BOHR = 0.529177210903          # angstrom
HARTREE = 627.509474           # kcal/mol
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3de"
SERIES_1, SERIES_2 = "#2a78d6", "#eb6834"


def parse_args():
    p = argparse.ArgumentParser(description="Plot the CV, restraint and energies of a CP2K steered MD run.",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-c", "--colvar", required=True, help="CP2K metadynLog file with the CV")
    p.add_argument("-e", "--ener", default=None, help="CP2K .ener file of the same run (adds the potential energy)")
    p.add_argument("--cv-col", type=int, default=2, help="Column of the CV in the metadynLog (1 = time)")
    p.add_argument("--atomic", action="store_true", help="The CV is not a distance: do not convert bohr to angstrom")
    p.add_argument("--target", type=float, default=-2.0, help="TARGET at t=0 (angstrom)")
    p.add_argument("--growth", type=float, default=0.002, help="TARGET_GROWTH (angstrom/fs)")
    p.add_argument("--limit", type=float, default=2.0, help="TARGET_LIMIT (angstrom)")
    p.add_argument("-k", "--kforce", type=float, default=100.0, help="Restraint constant K (kcal/mol/angstrom^2)")
    p.add_argument("--t0", type=float, default=0.0, help="Time (fs) at which the target starts to move")
    p.add_argument("--cv-label", default="CV", help="Label of the CV axis")
    p.add_argument("-o", "--output", default="smd.png", help="Output figure (saved at --dpi)")
    p.add_argument("--dpi", type=int, default=600, help="Resolution of the saved figure")
    p.add_argument("--csv", default=None, help="Also write time, CV, target, restraint energy and work to this file")
    p.add_argument("--no-show", action="store_true", help="Do not open a window with the plot (only save it)")
    return p.parse_args()


def load_table(path, min_cols):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            try:
                vals = [float(x) for x in line.split()]
            except ValueError:
                continue
            if len(vals) >= min_cols:
                rows.append(vals[:max(min_cols, 7)])
    if not rows:
        sys.exit(f"ERROR: no data read from {path}")
    width = min(len(r) for r in rows)
    return np.array([r[:width] for r in rows])


def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=9)
    ax.yaxis.label.set_color(INK)


def main():
    a = parse_args()
    import matplotlib
    if a.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dat = load_table(a.colvar, a.cv_col)
    t = dat[:, 0]
    cv = dat[:, a.cv_col - 1] * (1.0 if a.atomic else BOHR)
    sign = 1.0 if a.limit >= a.target else -1.0
    target = a.target + a.growth * np.clip(t - a.t0, 0, None)
    target = np.minimum(target, a.limit) if sign > 0 else np.maximum(target, a.limit)
    e_rest = a.kforce * (cv - target) ** 2
    dtarget = np.diff(target, prepend=target[0])
    work = np.cumsum(-2.0 * a.kforce * (cv - target) * dtarget)

    print(f"{a.colvar}: {len(t)} points, t = {t[0]:.1f} - {t[-1]:.1f} fs")
    print(f"  CV  start {cv[0]:8.3f}   end {cv[-1]:8.3f}   min {cv.min():8.3f}   max {cv.max():8.3f}  (angstrom)")
    print(f"  max |CV - target| {np.abs(cv - target).max():6.3f} A at t = {t[np.argmax(np.abs(cv - target))]:.1f} fs")
    print(f"  work at the end {work[-1]:8.2f} kcal/mol, maximum {work.max():8.2f} kcal/mol at t = {t[np.argmax(work)]:.1f} fs")

    panels = 3 if a.ener else 2
    fig, axes = plt.subplots(panels, 1, figsize=(7.0, 2.4 * panels), sharex=True, constrained_layout=True)

    ax = axes[0]
    ax.plot(t, target, color=MUTED, linewidth=1.5, linestyle=(0, (4, 3)), label="restraint target")
    ax.plot(t, cv, color=SERIES_1, linewidth=2.0, label="CV")
    ax.set_ylabel(f"{a.cv_label} (Å)" if not a.atomic else a.cv_label)
    ax.legend(frameon=False, fontsize=9, loc="best")
    style(ax)

    ax = axes[1]
    ax.plot(t, work, color=SERIES_1, linewidth=2.0)
    ax.set_ylabel("Work (kcal/mol)")
    style(ax)

    if a.ener:
        en = load_table(a.ener, 5)
        te, pot = en[:, 1], (en[:, 4] - en[0, 4]) * HARTREE
        ax = axes[2]
        ax.plot(te, pot, color=SERIES_1, linewidth=1.2)
        ax.set_ylabel("E$_{pot}$ - E$_{pot}$(0)\n(kcal/mol)")
        style(ax)
        print(f"{a.ener}: {len(te)} steps, potential energy range {pot.min():.1f} to {pot.max():.1f} kcal/mol")

    axes[-1].set_xlabel("Time (fs)", color=INK)
    fig.savefig(a.output, dpi=a.dpi)
    print(f"Figure saved to {a.output} ({a.dpi} dpi)")

    if a.csv:
        np.savetxt(a.csv, np.column_stack([t, cv, target, e_rest, work]), fmt="%12.4f",
                   header="time[fs] CV[A] target[A] E_restraint[kcal/mol] work[kcal/mol]")
        print(f"Data written to {a.csv}")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
