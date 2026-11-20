#!/usr/bin/env python3
"""
plot_cp2k.py

Plots of the CP2K runs of the tutorial, shown on screen and saved at 600 dpi.

  opt   energy against the iteration number of one or more geometry optimizations (GEO_OPT),
        one panel per output file. The energy of iteration 0 is the first energy of the run,
        the others are the "Total Energy" of the "Informations at step" blocks.
  temp  temperature of the system against the simulation time of an MD run, from the .ener file,
        with a moving average and the target temperature.

The figures use a dark background by default (--theme dark, matching the tutorial pages);
--theme light gives a white background for documents and slides.

Usage examples
  python3 plot_cp2k.py opt mm_em.out
  python3 plot_cp2k.py opt qmmm_em_me.out qmmm_em_ee.out --titles "Mechanical embedding,Electrostatic embedding"
  python3 plot_cp2k.py temp QMMM_MD_EQ-1.ener
  python3 plot_cp2k.py temp QMMM_MD_EQ-1.ener --window 200 --target 300 --no-show --prefix md_temp
  python3 plot_cp2k.py temp QMMM_MD_EQ-1.ener --theme light
"""

import argparse
import os
import re
import sys

import numpy as np

HARTREE_TO_KCAL = 627.509474
HARTREE_TO_KJ = 2625.499639

# background, text, axes/ticks, grid and data colours of each theme
# (dark: same surface and text colours as the tutorial pages)
THEMES = {
    "dark": {"bg": "#0e1014", "ink": "#f3f5f8", "muted": "#9ba3b1", "grid": "#262a32", "series": "#3987e5",
             "faint": 0.5},
    "light": {"bg": "#ffffff", "ink": "#1f1f1e", "muted": "#6b6a64", "grid": "#e4e3de", "series": "#2a78d6",
              "faint": 0.35},
}
# opacity of the instantaneous temperature trace drawn under its moving average: FAINT
BG, INK, MUTED, GRID, SERIES_1, FAINT = (THEMES["dark"][k] for k in ("bg", "ink", "muted", "grid", "series", "faint"))


def set_theme(name, matplotlib):
    """Colours of the chosen theme, also applied to the figure, axes, legend and saved file."""
    global BG, INK, MUTED, GRID, SERIES_1, FAINT
    th = THEMES[name]
    BG, INK, MUTED, GRID, SERIES_1, FAINT = th["bg"], th["ink"], th["muted"], th["grid"], th["series"], th["faint"]
    matplotlib.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "savefig.edgecolor": BG,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
        "xtick.color": MUTED, "ytick.color": MUTED, "grid.color": GRID,
    })


# ------------------------------------------------------------------------------------------ input
def parse_args():
    ap = argparse.ArgumentParser(description="Plots of CP2K geometry optimizations and MD runs.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--prefix", default=None,
                        help="Name of the output .png and .dat files (default: from the first input file)")
    common.add_argument("--dpi", type=int, default=600, help="Resolution of the saved figure")
    common.add_argument("--no-show", action="store_true", help="Do not open a window with the plot (only save it)")
    common.add_argument("--theme", choices=sorted(THEMES), default="dark",
                        help="Colour theme: dark background (tutorial pages) or white background")

    po = sub.add_parser("opt", parents=[common], formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                        help="Energy against the iteration number of GEO_OPT runs")
    po.add_argument("outputs", nargs="+", help="CP2K output files of GEO_OPT runs (one panel each)")
    po.add_argument("--titles", default=None,
                    help="Comma-separated panel titles (default: the file names)")
    po.add_argument("--unit", choices=["hartree", "kcal", "kj"], default="hartree", help="Energy unit")
    po.add_argument("--relative", action="store_true",
                    help="Plot the energy relative to the last iteration instead of the absolute energy")
    po.add_argument("--logy", action="store_true",
                    help="Logarithmic energy axis (implies --relative; the last point is left out)")

    pt = sub.add_parser("temp", parents=[common], formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                        help="Temperature against the simulation time of an MD run")
    pt.add_argument("ener", help="CP2K .ener file (Step, Time[fs], Kin, Temp[K], Pot, Cons Qty, UsedTime)")
    pt.add_argument("--window", type=float, default=100.0, help="Width of the moving average (fs)")
    pt.add_argument("--target", type=float, default=300.0, help="Target temperature (K) drawn as a reference line")
    pt.add_argument("--xunit", choices=["fs", "ps"], default="ps", help="Unit of the time axis")
    return ap.parse_args()


def read_opt(path):
    """Energies (hartree) of a GEO_OPT output: the first energy of the run, then one per iteration."""
    with open(path, "r", errors="replace") as fh:
        text = fh.read()
    first = re.search(r"ENERGY\| Total FORCE_EVAL \( \w+ \) energy \[a\.u\.\]:\s+(\S+)", text)
    steps = re.findall(r"Informations at step =\s+(\d+).*?Total Energy\s+=\s+(\S+)", text, re.S)
    if first is None or not steps:
        sys.exit(f"ERROR: no geometry optimization steps found in {path}")
    it = np.array([0] + [int(s) for s, _ in steps])
    e = np.array([float(first.group(1))] + [float(v) for _, v in steps])
    status = "converged" if ("GEOMETRY OPTIMIZATION COMPLETED" in text or "run CONVERGED" in text) else (
        "maximum number of steps reached" if "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in text
        else "not finished")
    return it, e, status


def read_ener(path):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    return data[:, 1], data[:, 3]


def moving_average(t, y, width):
    """Centred moving average over a time window (fs), shrunk at the ends."""
    out = np.empty_like(y)
    h = width / 2.0
    for i in range(len(y)):
        m = (t >= t[i] - h) & (t <= t[i] + h)
        out[i] = y[m].mean()
    return out


# ------------------------------------------------------------------------------------------- plot
def style(ax):
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=9)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def default_prefix(path, suffix):
    base = os.path.basename(path)
    for ext in (".out", ".ener", ".log"):
        if base.endswith(ext):
            base = base[: -len(ext)]
    return f"{base}_{suffix}"


def plot_opt(a, plt):
    titles = a.titles.split(",") if a.titles else [os.path.basename(f) for f in a.outputs]
    if len(titles) != len(a.outputs):
        sys.exit("ERROR: --titles must have one title per output file")
    conv = {"hartree": 1.0, "kcal": HARTREE_TO_KCAL, "kj": HARTREE_TO_KJ}[a.unit]
    lab = {"hartree": "hartree", "kcal": "kcal/mol", "kj": "kJ/mol"}[a.unit]
    relative = a.relative or a.logy

    n = len(a.outputs)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 3.6), constrained_layout=True, squeeze=False)
    rows = []
    for ax, f, title in zip(axes[0], a.outputs, titles):
        it, e, status = read_opt(f)
        y = (e - e[-1]) * conv if relative else e * conv
        print(f"{f}: {len(it) - 1} iterations ({status}), energy {e[0]:.6f} -> {e[-1]:.6f} hartree "
              f"({(e[-1] - e[0]) * HARTREE_TO_KCAL:.1f} kcal/mol)")
        if a.logy:
            ax.semilogy(it[:-1], y[:-1], color=SERIES_1, linewidth=2.0)
        else:
            ax.plot(it, y, color=SERIES_1, linewidth=2.0)
        ax.set_title(title, color=INK, fontsize=10)
        ax.set_xlabel("Iteration")
        ax.set_ylabel(f"E - E(last) ({lab})" if relative else f"Energy ({lab})")
        ax.margins(x=0.01)
        style(ax)
        rows.append((f, it, e))

    prefix = a.prefix or default_prefix(a.outputs[0], "energy")
    with open(f"{prefix}.dat", "w") as fh:
        for f, it, e in rows:
            fh.write(f"# {f}\n# iteration  energy[hartree]\n")
            for i, v in zip(it, e):
                fh.write(f"{i:8d} {v:20.10f}\n")
            fh.write("\n\n")
    return fig, prefix


def plot_temp(a, plt):
    t, temp = read_ener(a.ener)
    tconv, tlab = (1.0, "fs") if a.xunit == "fs" else (1e-3, "ps")
    avg = moving_average(t, temp, a.window)
    late = t >= t[-1] / 2.0
    print(f"{a.ener}: {len(t)} steps, t = {t[0]:.0f} to {t[-1]:.0f} fs, T = {temp.mean():.1f} +/- {temp.std():.1f} K "
          f"over the run, {temp[late].mean():.1f} +/- {temp[late].std():.1f} K over the second half")

    fig, ax = plt.subplots(figsize=(7.0, 3.4), constrained_layout=True)
    ax.plot(t * tconv, temp, color=SERIES_1, alpha=FAINT, linewidth=0.8, label="instantaneous")
    ax.plot(t * tconv, avg, color=SERIES_1, linewidth=2.0, label=f"{a.window:g} fs moving average")
    ax.axhline(a.target, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)), label=f"target {a.target:g} K")
    ax.set_xlabel(f"Time ({tlab})")
    ax.set_ylabel("Temperature (K)")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.margins(x=0.01)
    style(ax)

    prefix = a.prefix or default_prefix(a.ener, "temperature")
    with open(f"{prefix}.dat", "w") as fh:
        fh.write(f"# {a.ener}\n# time[fs]  T[K]  T_moving_average[K]\n")
        for x, y, m in zip(t, temp, avg):
            fh.write(f"{x:10.2f} {y:10.3f} {m:10.3f}\n")
    return fig, prefix


def main():
    a = parse_args()
    import matplotlib
    if a.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    set_theme(a.theme, matplotlib)

    fig, prefix = plot_opt(a, plt) if a.mode == "opt" else plot_temp(a, plt)
    fig.savefig(f"{prefix}.png", dpi=a.dpi)
    print(f"wrote {prefix}.dat")
    print(f"wrote {prefix}.png ({a.dpi} dpi)")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
