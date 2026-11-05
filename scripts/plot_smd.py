#!/usr/bin/env python3
"""
plot_smd.py

Time view of a CP2K steered MD run (moving harmonic restraint on a COLVAR, &CONSTRAINT/&COLLECTIVE
with TARGET_GROWTH), to be read side by side with the trajectory in VMD. The panels share the time
axis, one quantity per panel:

  1. the CV and the restraint centre lambda(t) = TARGET + TARGET_GROWTH*t (stopped at TARGET_LIMIT)
  2. the work done by the moving restraint, W(t) = int -2K (s - lambda) dlambda, with the same
     trapezoidal integration as smd_jarzynski.py (so both scripts give the same numbers)
  3. optionally, distances between pairs of atoms read from the DCD trajectory (--dcd --pairs)

The restraint parameters (TARGET, TARGET_GROWTH, TARGET_LIMIT, K), the time step and the stride of
the trajectory are read from the input echo of the CP2K output (-o), as in smd_jarzynski.py, or
given on the command line in the same internal CP2K units.

Inputs
  PROJECT-COLVAR.metadynLog : column 0 is the time (fs), column 1 the CV in internal units (bohr)
  PROJECT.out               : CP2K output with the input echo (ECHO_INPUT)
  PROJECT-pos-1.dcd (opt.)  : trajectory; atoms are given by their serial numbers (1-based, as in CP2K)

Outputs (PREFIX = plot_smd by default)
  PREFIX.png       figure, shown on screen and saved at --dpi (600)
  PREFIX.dat       time, lambda, CV, restraint energy and work on the time grid of the metadynLog
  PREFIX_dist.dat  time and distances of every trajectory frame (with --dcd)

Usage examples
  python3 plot_smd.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out
  python3 plot_smd.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out \\
      --dcd QMMM_MD_SMD-pos-1.dcd --pairs "7699:7700 3607:7699 3026:7700 3607:7783 7785:7699" \\
      --pair-labels "O4-C1,HE2-O4,OD1-C1,HE2-Ow,H2w-O4"
  python3 plot_smd.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out --xunit ps --no-show
  python3 plot_smd.py -l QMMM_MD_SMD-COLVAR.metadynLog -o qmmm_md_smd.out --theme light

The figure uses a dark background by default (--theme dark, matching the tutorial pages);
--theme light gives a white background for documents and slides.
"""

import argparse
import os
import re
import struct
import sys

import numpy as np

HARTREE_TO_KCAL = 627.509474
HARTREE_TO_KJ = 2625.499639
BOHR_TO_ANG = 0.529177210903
AUT_PER_FS = 41.341373336  # atomic time units per fs

# background, text, axes/ticks, grid and categorical colours of each theme (series assigned in this
# order, at most 8 distances); dark: same surface and text colours as the tutorial pages
THEMES = {
    "dark": {"bg": "#0e1014", "ink": "#f3f5f8", "muted": "#9ba3b1", "grid": "#262a32",
             "series": ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]},
    "light": {"bg": "#ffffff", "ink": "#1f1f1e", "muted": "#6b6a64", "grid": "#e4e3de",
              "series": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]},
}
BG, INK, MUTED, GRID = (THEMES["dark"][k] for k in ("bg", "ink", "muted", "grid"))
SERIES = list(THEMES["dark"]["series"])


def set_theme(name, matplotlib):
    """Colours of the chosen theme, also applied to the figure, axes, legend and saved file."""
    global BG, INK, MUTED, GRID, SERIES
    th = THEMES[name]
    BG, INK, MUTED, GRID, SERIES = th["bg"], th["ink"], th["muted"], th["grid"], list(th["series"])
    matplotlib.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG, "savefig.edgecolor": BG,
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
        "xtick.color": MUTED, "ytick.color": MUTED, "grid.color": GRID,
    })


# ------------------------------------------------------------------------------------------ input
def parse_args():
    ap = argparse.ArgumentParser(
        description="Plot the CV, the work and selected distances of a CP2K steered MD run against time.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-l", "--log", required=True, help="CP2K metadynLog file with the CV")
    ap.add_argument("-o", "--out", default=None,
                    help="CP2K output file: restraint parameters, time step and trajectory stride from "
                         "the input echo")
    ap.add_argument("--target", type=float, default=None, help="Initial TARGET (internal units, e.g. bohr)")
    ap.add_argument("--growth", type=float, default=None,
                    help="TARGET_GROWTH (internal units per atomic time unit)")
    ap.add_argument("--limit", type=float, default=None, help="TARGET_LIMIT (internal units)")
    ap.add_argument("-k", "--kforce", type=float, default=None, help="Restraint K (Hartree / internal unit^2)")
    ap.add_argument("--kfac", type=float, default=2.0,
                    help="Prefactor of K in dU/dlambda: 2 for CP2K's U = K(s-lambda)^2, 1 for U = K/2(s-lambda)^2")
    ap.add_argument("--cv-col", type=int, default=1, help="0-based column of the CV in the metadynLog")
    ap.add_argument("--cv-unit", choices=["bohr", "angstrom", "none"], default="bohr",
                    help="Unit of the CV in the log (internal CP2K units are bohr for distances)")
    ap.add_argument("--cv-label", default="CV", help="Label of the CV axis")
    ap.add_argument("--energy-unit", choices=["kcal", "kj", "hartree"], default="kcal",
                    help="Unit of the work and of the restraint energy")
    ap.add_argument("--dcd", default=None, help="DCD trajectory of the run, to plot distances")
    ap.add_argument("--pairs", default=None,
                    help='Atom pairs for the distance panel, serial numbers (1-based) as "i:j i:j ..." (at most 8)')
    ap.add_argument("--pair-labels", default=None,
                    help='Comma-separated labels of the pairs (default: the "i:j" strings)')
    ap.add_argument("--dcd-stride", type=int, default=None,
                    help="MD steps between trajectory frames (default: &TRAJECTORY/&EACH/MD of the output, else 5)")
    ap.add_argument("--timestep", type=float, default=None,
                    help="MD time step in fs (default: &MD/TIMESTEP of the output, else 1.0)")
    ap.add_argument("--xunit", choices=["fs", "ps"], default="fs", help="Unit of the time axis")
    ap.add_argument("--window", type=float, default=150.0,
                    help="Half width (fs) of the window used to report the local maxima and minima of W")
    ap.add_argument("--prefix", default="plot_smd", help="Prefix of the output .png and .dat files")
    ap.add_argument("--dpi", type=int, default=600, help="Resolution of the saved figure")
    ap.add_argument("--no-show", action="store_true", help="Do not open a window with the plot (only save it)")
    ap.add_argument("--theme", choices=sorted(THEMES), default="dark",
                    help="Colour theme: dark background (tutorial pages) or white background")
    return ap.parse_args()


def section(text, name):
    """Body of the first &NAME ... &END NAME block of the input echo (None if absent)."""
    m = re.search(rf"&{name}\b(.*?)&END {name}\b", text, re.S)
    return m.group(1) if m else None


def keyword(block, key):
    if block is None:
        return None
    m = re.search(rf"^\s*{key}\s+(\S+)", block, re.M)
    return float(m.group(1)) if m else None


def parse_cp2k_echo(path):
    """Restraint parameters, MD time step and trajectory stride from the input echo of a CP2K output."""
    with open(path, "r", errors="replace") as fh:
        text = fh.read()
    p = {}
    coll = section(text, "COLLECTIVE")
    if coll is None:
        print(f"WARNING: no &COLLECTIVE section in {path}")
    for key, kw in (("target", "TARGET"), ("growth", "TARGET_GROWTH"), ("limit", "TARGET_LIMIT"), ("k", "K")):
        val = keyword(coll, kw)
        if val is not None:
            p[key] = val
    md = section(text, "MD")
    if md is not None:
        p["timestep"] = keyword(md, "TIMESTEP")
    traj = section(text, "TRAJECTORY")
    each = section(traj, "EACH") if traj else None
    stride = keyword(each, "MD")
    if stride is not None:
        p["stride"] = int(stride)
    return p


def read_metadynlog(path, cv_col):
    data = np.loadtxt(path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    return data[:, 0], data[:, cv_col]


# -------------------------------------------------------------------------------------------- DCD
def read_dcd_atoms(path, indices):
    """
    Coordinates (angstrom) of the atoms with the given 0-based indices in every frame of a CHARMM/CP2K
    DCD file, and the unit cell of every frame (a, b, c, alpha, beta, gamma) or None.
    The number of frames is taken from the file size: CP2K leaves it at 0 in the header.
    """
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        head = fh.read(4)
        if struct.unpack("<i", head)[0] == 84:
            e = "<"
        elif struct.unpack(">i", head)[0] == 84:
            e = ">"
        else:
            sys.exit(f"ERROR: {path} does not look like a DCD file")
        rec = fh.read(84)
        if rec[:4] != b"CORD":
            sys.exit(f"ERROR: {path} is not a coordinate DCD file (no CORD tag)")
        icntrl = struct.unpack(f"{e}9if10i", rec[4:])
        has_cell, has_4d = icntrl[10] != 0, icntrl[11] != 0
        if has_4d:
            sys.exit("ERROR: DCD files with a 4th dimension are not supported")
        fh.read(4)
        ntitle_bytes = struct.unpack(f"{e}i", fh.read(4))[0]
        fh.seek(ntitle_bytes + 4, 1)
        if struct.unpack(f"{e}i", fh.read(4))[0] != 4:
            sys.exit(f"ERROR: unexpected DCD layout in {path}")
        natom = struct.unpack(f"{e}i", fh.read(4))[0]
        fh.read(4)
        start = fh.tell()
        cell_bytes = 56 if has_cell else 0
        coord_bytes = 8 + 4 * natom
        frame_bytes = cell_bytes + 3 * coord_bytes
        nframes = (size - start) // frame_bytes
        if nframes == 0:
            sys.exit(f"ERROR: no complete frame in {path}")
        fh.seek(start)
        first = fh.read(frame_bytes)
        markers = [struct.unpack(f"{e}i", first[o:o + 4])[0]
                   for o in ([0] if has_cell else []) + [cell_bytes + n * coord_bytes for n in range(3)]]
        expected = ([48] if has_cell else []) + [4 * natom] * 3
        if markers != expected:
            sys.exit(f"ERROR: unexpected frame layout in {path} (record markers {markers}, expected {expected})")
    if (size - start) % frame_bytes:
        print(f"WARNING: {path} ends with an incomplete frame, which is ignored")
    if np.max(indices) >= natom:
        sys.exit(f"ERROR: atom serial {np.max(indices) + 1} is beyond the {natom} atoms of {path}")

    raw = np.memmap(path, dtype=np.uint8, mode="r", offset=start, shape=(nframes, frame_bytes))
    xyz = np.empty((nframes, len(indices), 3))
    for dim in range(3):
        base = cell_bytes + dim * coord_bytes + 4
        cols = (base + 4 * np.asarray(indices))[:, None] + np.arange(4)[None, :]
        xyz[:, :, dim] = raw[:, cols.ravel()].reshape(nframes, len(indices), 4).copy().view(f"{e}f4")[..., 0]
    cell = None
    if has_cell:
        c = raw[:, 4:52].copy().view(f"{e}f8")  # CHARMM order: a, gamma, b, beta, alpha, c
        cell = np.column_stack([c[:, 0], c[:, 2], c[:, 5], c[:, 4], c[:, 3], c[:, 1]])
    return xyz, cell, natom


def pair_distances(xyz, cell, pairs_idx, atom_pos):
    """Distances of the pairs in every frame, with the minimum image for orthorhombic cells."""
    d = np.empty((xyz.shape[0], len(pairs_idx)))
    ortho = False
    if cell is not None:
        ang = cell[:, 3:]
        # angles are written in degrees by CP2K; some programs write cosines instead
        ortho = bool(np.all(np.abs(ang - 90.0) < 1e-3) or np.all(np.abs(ang) < 1e-6))
    for n, (i, j) in enumerate(pairs_idx):
        v = xyz[:, atom_pos[j]] - xyz[:, atom_pos[i]]
        if ortho:
            box = cell[:, :3]
            v = v - box * np.round(v / box)
        d[:, n] = np.linalg.norm(v, axis=1)
    return d, ortho


# ------------------------------------------------------------------------------------------- work
def restraint_and_work(time_fs, s, lam0, v_aut, lam_lim, k, kfac):
    """lambda(t), restraint energy K(s-lambda)^2 and trapezoidal work (internal units, Hartree)."""
    lam = lam0 + v_aut * time_fs * AUT_PER_FS
    if lam_lim is not None:
        lam = np.minimum(lam, lam_lim) if v_aut > 0 else np.maximum(lam, lam_lim)
    force = -kfac * k * (s - lam)  # dU/dlambda
    w = np.concatenate(([0.0], np.cumsum(0.5 * (force[1:] + force[:-1]) * np.diff(lam))))
    u = (kfac / 2.0) * k * (s - lam) ** 2
    return lam, u, w


def turning_points(t, w, half_width):
    """Local maxima and minima of w: points that are the extreme of w within +/- half_width (fs)."""
    out = []
    for i in range(1, len(t) - 1):
        if t[i] - half_width < t[0] or t[i] + half_width > t[-1]:
            continue
        win = (t >= t[i] - half_width) & (t <= t[i] + half_width)
        if w[i] == w[win].max():
            out.append(("max", i))
        elif w[i] == w[win].min():
            out.append(("min", i))
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
    ax.yaxis.label.set_color(INK)


def main():
    a = parse_args()

    # ---- parameters
    p = parse_cp2k_echo(a.out) if a.out else {}
    for key, val in (("target", a.target), ("growth", a.growth), ("limit", a.limit), ("k", a.kforce),
                     ("timestep", a.timestep), ("stride", a.dcd_stride)):
        if val is not None:
            p[key] = val
    missing = [x for x in ("target", "growth", "k") if p.get(x) is None]
    if missing:
        sys.exit(f"ERROR: missing restraint parameters {missing}; give -o/--out or set them explicitly")
    dt = p.get("timestep") or 1.0
    stride = int(p.get("stride") or 5)
    if bool(a.dcd) != bool(a.pairs):
        sys.exit("ERROR: --dcd and --pairs go together")

    e_conv = {"kcal": HARTREE_TO_KCAL, "kj": HARTREE_TO_KJ, "hartree": 1.0}[a.energy_unit]
    e_lab = {"kcal": "kcal/mol", "kj": "kJ/mol", "hartree": "Ha"}[a.energy_unit]
    x_conv = BOHR_TO_ANG if a.cv_unit == "bohr" else 1.0
    x_lab = "Å" if a.cv_unit in ("bohr", "angstrom") else "CV units"
    s_to_internal = 1.0 / BOHR_TO_ANG if a.cv_unit == "angstrom" else 1.0
    t_conv, t_lab = (1.0, "fs") if a.xunit == "fs" else (1e-3, "ps")

    # ---- CV, restraint and work
    t, s = read_metadynlog(a.log, a.cv_col)
    s = s * s_to_internal
    lam, u, w = restraint_and_work(t, s, p["target"], p["growth"], p.get("limit"), p["k"], a.kfac)
    lag = (s - lam) * x_conv
    k_user = p["k"] * HARTREE_TO_KCAL / (x_conv ** 2 if a.cv_unit == "bohr" else 1.0)
    print(f"Restraint: TARGET {p['target'] * x_conv:.3f} {x_lab}, TARGET_LIMIT "
          f"{(p['limit'] * x_conv if p.get('limit') is not None else float('nan')):.3f} {x_lab}, "
          f"K {k_user:.2f} kcal/mol/{x_lab}^2 (U = K (s-lambda)^2), "
          f"speed {p['growth'] * AUT_PER_FS * x_conv * 1000:.3f} {x_lab}/ps; time step {dt:g} fs")
    print(f"{a.log}: {len(t)} points, t = {t[0]:.1f} to {t[-1]:.1f} fs")
    print(f"  CV     start {s[0] * x_conv:8.3f}  end {s[-1] * x_conv:8.3f}  min {s.min() * x_conv:8.3f}  "
          f"max {s.max() * x_conv:8.3f} {x_lab}")
    i_lag = int(np.argmax(np.abs(lag)))
    print(f"  |CV - lambda|  mean {np.abs(lag).mean():.3f} {x_lab}, max {np.abs(lag[i_lag]):.3f} {x_lab} "
          f"at t = {t[i_lag]:.1f} fs")
    print(f"  work   end {w[-1] * e_conv:8.2f} {e_lab}")
    for kind, i in turning_points(t, w, a.window):
        print(f"  work   local {kind} {w[i] * e_conv:8.2f} {e_lab} at t = {t[i]:7.1f} fs, "
              f"lambda = {lam[i] * x_conv:6.2f} {x_lab}")

    np.savetxt(f"{a.prefix}.dat", np.column_stack([t, lam * x_conv, s * x_conv, u * e_conv, w * e_conv]),
               fmt="%12.5f", header=f"time[fs] lambda[{x_lab}] CV[{x_lab}] U_restraint[{e_lab}] W[{e_lab}]")
    written = [f"{a.prefix}.dat"]

    # ---- distances
    td = dist = labels = None
    if a.dcd:
        try:
            pairs = [tuple(int(x) for x in tok.split(":")) for tok in a.pairs.split()]
            assert pairs and all(len(pr) == 2 for pr in pairs)
        except (ValueError, AssertionError):
            sys.exit('ERROR: --pairs must look like "7699:7700 3607:7699"')
        if len(pairs) > len(SERIES):
            sys.exit(f"ERROR: at most {len(SERIES)} pairs can be plotted together; run the script again for the others")
        labels = a.pair_labels.split(",") if a.pair_labels else [f"{i}:{j}" for i, j in pairs]
        if len(labels) != len(pairs):
            sys.exit("ERROR: --pair-labels must have one label per pair")
        serials = sorted({x for pr in pairs for x in pr})
        xyz, cell, natom = read_dcd_atoms(a.dcd, np.array(serials) - 1)
        pos = {sr - 1: n for n, sr in enumerate(serials)}
        dist, ortho = pair_distances(xyz, cell, [(i - 1, j - 1) for i, j in pairs], pos)
        nfr = xyz.shape[0]
        # the last frame is the last MD step (ADD_LAST); earlier frames are 'stride' steps apart
        td = t[-1] - (nfr - 1 - np.arange(nfr)) * stride * dt
        print(f"{a.dcd}: {natom} atoms, {nfr} frames every {stride} steps "
              f"(t = {td[0]:.1f} to {td[-1]:.1f} fs){', minimum image' if ortho else ''}")
        for n, lab in enumerate(labels):
            print(f"  {lab:>12s}  start {dist[0, n]:6.2f}  end {dist[-1, n]:6.2f}  "
                  f"min {dist[:, n].min():6.2f}  max {dist[:, n].max():6.2f} Å")
        np.savetxt(f"{a.prefix}_dist.dat", np.column_stack([td, dist]), fmt="%12.4f",
                   header="time[fs] " + " ".join(f"{lab}[Å]" for lab in labels))
        written.append(f"{a.prefix}_dist.dat")

    # ---- figure
    import matplotlib
    if a.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    set_theme(a.theme, matplotlib)

    panels = 3 if dist is not None else 2
    fig, axes = plt.subplots(panels, 1, figsize=(7.0, 2.4 * panels), sharex=True, constrained_layout=True)

    ax = axes[0]
    ax.plot(t * t_conv, lam * x_conv, color=MUTED, linewidth=1.5, linestyle=(0, (4, 3)),
            label="restraint centre λ")
    ax.plot(t * t_conv, s * x_conv, color=SERIES[0], linewidth=2.0, label=a.cv_label)
    ax.set_ylabel(f"{a.cv_label} ({x_lab})" if x_lab != "CV units" else a.cv_label)
    ax.legend(frameon=False, fontsize=9, loc="best")
    style(ax)

    ax = axes[1]
    ax.plot(t * t_conv, w * e_conv, color=SERIES[0], linewidth=2.0)
    ax.set_ylabel(f"Work ({e_lab})")
    style(ax)

    if dist is not None:
        ax = axes[2]
        for m, lab in enumerate(labels):
            ax.plot(td * t_conv, dist[:, m], color=SERIES[m], linewidth=1.5, label=lab)
        ax.set_ylabel("Distance (Å)")
        ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
        style(ax)

    axes[-1].set_xlabel(f"Time ({t_lab})", color=INK)
    fig.savefig(f"{a.prefix}.png", dpi=a.dpi)
    written.append(f"{a.prefix}.png ({a.dpi} dpi)")
    for f in written:
        print(f"wrote {f}")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
