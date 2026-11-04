#!/usr/bin/env python3
"""
smd_jarzynski.py

Free energy profile from CP2K steered MD (moving harmonic restraint on a COLVAR,
&CONSTRAINT/&COLLECTIVE with TARGET_GROWTH) using the Jarzynski equality.

For each trajectory the accumulated work is

    W(t) = int_0^t  dU/dlambda * (dlambda/dt') dt'
         = int_0^t  -c*K*(s(t') - lambda(t')) * v  dt'

where lambda(t) = TARGET + v*t (clipped at TARGET_LIMIT), v = TARGET_GROWTH and
CP2K's restraint is U = K*(s - lambda)^2 (no 1/2, so c = 2 by default).

With N trajectories pulled under the same protocol:
    Jarzynski (exponential average):  F = -kT ln < exp(-W/kT) >
    2nd-order cumulant expansion:     F = <W> - var(W)/(2kT)
With a single trajectory both reduce to W itself (an upper bound, not a PMF).

Input: CP2K *.metadynLog files (col 1 = time in fs, col 2 = CV in internal units,
i.e. bohr for distance CVs). Restraint parameters are read from the CP2K output
file (--out) or given explicitly on the command line (internal CP2K units).
"""

import argparse
import re
import sys

import numpy as np

HARTREE_TO_KCAL = 627.509474
HARTREE_TO_KJ = 2625.499639
BOHR_TO_ANG = 0.529177210903
AUT_PER_FS = 41.341373336  # atomic time units per fs
KB_HARTREE = 3.166811563e-6  # Boltzmann constant in Ha/K


def parse_cp2k_out(path):
    """Read restraint/MD parameters from the input echo of a CP2K output file."""
    params = {}
    with open(path, "r", errors="replace") as fh:
        text = fh.read()
    coll = re.search(r"&COLLECTIVE(.*?)&END COLLECTIVE", text, re.S)
    if coll is None:
        sys.exit(f"ERROR: no &COLLECTIVE block found in {path}")
    block = coll.group(1)
    pats = {
        "target": r"^\s*TARGET\s+(\S+)",
        "growth": r"^\s*TARGET_GROWTH\s+(\S+)",
        "limit": r"^\s*TARGET_LIMIT\s+(\S+)",
        "k": r"^\s*K\s+(\S+)",
    }
    for key, pat in pats.items():
        m = re.search(pat, block, re.M)
        if m:
            params[key] = float(m.group(1))
    md = re.search(r"&MD(.*?)&END MD", text, re.S)
    if md:
        m = re.search(r"^\s*TEMPERATURE\s+(\S+)", md.group(1), re.M)
        if m:
            params["temp"] = float(m.group(1))
    return params


def read_metadynlog(path, cv_col):
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data[None, :]
    return data[:, 0], data[:, cv_col]


def cumulative_work(time_fs, s, lam0, v_aut, lam_lim, k, kfac):
    """Trapezoidal integration of dW = -kfac*K*(s - lambda)*dlambda (Hartree)."""
    t_aut = time_fs * AUT_PER_FS
    lam = lam0 + v_aut * t_aut
    if lam_lim is not None:
        lam = np.minimum(lam, lam_lim) if v_aut > 0 else np.maximum(lam, lam_lim)
    force = -kfac * k * (s - lam)  # dU/dlambda
    dlam = np.diff(lam)
    w = np.concatenate(([0.0], np.cumsum(0.5 * (force[1:] + force[:-1]) * dlam)))
    return lam, w


def bootstrap(wmat, beta, nboot, rng):
    """Bootstrap std of the Jarzynski and cumulant estimators over replicas."""
    n = wmat.shape[0]
    jar = np.empty((nboot, wmat.shape[1]))
    cum = np.empty_like(jar)
    for b in range(nboot):
        idx = rng.integers(0, n, n)
        sample = wmat[idx]
        jar[b] = jarzynski(sample, beta)
        cum[b] = sample.mean(0) - 0.5 * beta * sample.var(0)
    return jar.std(0), cum.std(0)


def jarzynski(wmat, beta):
    x = -beta * wmat
    xmax = x.max(axis=0)
    return -(xmax + np.log(np.mean(np.exp(x - xmax), axis=0))) / beta


def main():
    ap = argparse.ArgumentParser(
        description="Jarzynski free energy profile from CP2K steered MD metadynLog files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("-l", "--logs", nargs="+", required=True,
                    help="One or more CP2K *.metadynLog files (replicas of the same pulling protocol)")
    ap.add_argument("-o", "--out", default=None,
                    help="CP2K output file to read TARGET, TARGET_GROWTH, TARGET_LIMIT, K and TEMPERATURE from")
    ap.add_argument("--target", type=float, default=None, help="Initial TARGET (internal units, e.g. bohr)")
    ap.add_argument("--growth", type=float, default=None,
                    help="TARGET_GROWTH (internal units per atomic time unit)")
    ap.add_argument("--limit", type=float, default=None, help="TARGET_LIMIT (internal units)")
    ap.add_argument("-k", "--kforce", type=float, default=None, help="Restraint K (Hartree / internal unit^2)")
    ap.add_argument("-T", "--temp", type=float, default=None, help="Temperature (K)")
    ap.add_argument("--kfac", type=float, default=2.0,
                    help="Prefactor of K in dU/dlambda: 2 for CP2K's U = K(s-lambda)^2, 1 for U = K/2(s-lambda)^2")
    ap.add_argument("--cv-col", type=int, default=1, help="0-based column of the CV in the metadynLog")
    ap.add_argument("--cv-unit", choices=["bohr", "angstrom", "none"], default="bohr",
                    help="Unit of the CV in the log (internal CP2K units are bohr for distances)")
    ap.add_argument("--energy-unit", choices=["kcal", "kj", "hartree"], default="kcal", help="Output energy unit")
    ap.add_argument("--xaxis", choices=["lambda", "cv"], default="lambda",
                    help="Plot vs restraint centre lambda (proper Jarzynski variable) or instantaneous CV")
    ap.add_argument("--nboot", type=int, default=200, help="Bootstrap resamples for error bars (N>1 only)")
    ap.add_argument("--seed", type=int, default=1, help="RNG seed for the bootstrap")
    ap.add_argument("--prefix", default="smd_jarzynski", help="Prefix for output .dat and .png files")
    ap.add_argument("--dpi", type=int, default=600, help="Figure resolution")
    ap.add_argument("--no-show", action="store_true", help="Do not open an interactive plot window")
    args = ap.parse_args()

    p = parse_cp2k_out(args.out) if args.out else {}
    for key, val in (("target", args.target), ("growth", args.growth), ("limit", args.limit),
                     ("k", args.kforce), ("temp", args.temp)):
        if val is not None:
            p[key] = val
    missing = [x for x in ("target", "growth", "k", "temp") if x not in p]
    if missing:
        sys.exit(f"ERROR: missing parameters {missing}; give --out or set them explicitly")
    lam_lim = p.get("limit")

    beta = 1.0 / (KB_HARTREE * p["temp"])
    e_conv = {"kcal": HARTREE_TO_KCAL, "kj": HARTREE_TO_KJ, "hartree": 1.0}[args.energy_unit]
    e_lab = {"kcal": "kcal/mol", "kj": "kJ/mol", "hartree": "Ha"}[args.energy_unit]
    x_conv = BOHR_TO_ANG if args.cv_unit == "bohr" else 1.0
    x_lab = "Å" if args.cv_unit in ("bohr", "angstrom") else "CV units"
    # restraint parameters are internal (bohr); convert CV to bohr if log is in angstrom
    s_to_internal = 1.0 / BOHR_TO_ANG if args.cv_unit == "angstrom" else 1.0

    print(f"# TARGET = {p['target']:.6f}  TARGET_GROWTH = {p['growth']:.6e}  "
          f"TARGET_LIMIT = {lam_lim}  K = {p['k']:.6e}  T = {p['temp']:.1f} K")
    k_kcal = p['k'] * HARTREE_TO_KCAL / (x_conv**2 if args.cv_unit == 'bohr' else 1)
    print(f"# K = {k_kcal:.2f} kcal/mol/{x_lab}^2 as used by CP2K (U = K (s-lambda)^2); "
          f"equivalent k = {args.kfac*k_kcal:.2f} in the U = k/2 (s-lambda)^2 convention")
    print(f"# pulling speed = {p['growth']*AUT_PER_FS*x_conv*1000:.3f} {x_lab}/ps")

    works, lams, cvs, ref_t = [], [], [], None
    for f in args.logs:
        t, s = read_metadynlog(f, args.cv_col)
        s = s * s_to_internal
        lam, w = cumulative_work(t, s, p["target"], p["growth"], lam_lim, p["k"], args.kfac)
        if ref_t is None:
            ref_t = t
        n = min(len(ref_t), len(t))
        if len(t) != len(ref_t) or not np.allclose(t[:n], ref_t[:n]):
            print(f"# WARNING: {f} has a different time grid; truncating to common length")
        ref_t = ref_t[:n]
        works.append(w[:n]); lams.append(lam[:n]); cvs.append(s[:n])
        print(f"# {f}: {len(t)} points, final W = {w[-1]*e_conv:.2f} {e_lab}")

    n = len(ref_t)
    wmat = np.array([w[:n] for w in works])
    lam = lams[0][:n]
    smean = np.mean([c[:n] for c in cvs], axis=0)
    nrep = wmat.shape[0]

    f_jar = jarzynski(wmat, beta)
    f_cum = wmat.mean(0) - 0.5 * beta * wmat.var(0)
    w_mean = wmat.mean(0)
    if nrep > 1:
        rng = np.random.default_rng(args.seed)
        e_jar, e_cum = bootstrap(wmat, beta, args.nboot, rng)
    else:
        e_jar = e_cum = np.zeros(n)
        print("# NOTE: single trajectory -> F(lambda) = W(lambda), an upper bound that includes dissipation")

    x = (lam if args.xaxis == "lambda" else smean) * x_conv
    xname = "lambda" if args.xaxis == "lambda" else "<s>"

    datfile = f"{args.prefix}.dat"
    with open(datfile, "w") as fh:
        fh.write(f"# N_replicas = {nrep}, T = {p['temp']} K, energies in {e_lab}, distances in {x_lab}\n")
        fh.write("# time_fs  lambda  <s>  <W>  std(W)  F_Jarzynski  err_Jar  F_cumulant2  err_cum\n")
        for i in range(n):
            fh.write(f"{ref_t[i]:10.2f} {lam[i]*x_conv:10.5f} {smean[i]*x_conv:10.5f} "
                     f"{w_mean[i]*e_conv:12.5f} {wmat[:, i].std()*e_conv:10.5f} "
                     f"{f_jar[i]*e_conv:12.5f} {e_jar[i]*e_conv:10.5f} "
                     f"{f_cum[i]*e_conv:12.5f} {e_cum[i]*e_conv:10.5f}\n")
    print(f"# wrote {datfile}")

    import matplotlib
    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    for w in wmat:
        ax1.plot(x, w * e_conv, color="0.75", lw=0.8, zorder=1)
    if nrep > 1:
        ax1.plot(x, f_jar * e_conv, color="#1f5aa6", lw=2, label="Jarzynski")
        ax1.fill_between(x, (f_jar - e_jar) * e_conv, (f_jar + e_jar) * e_conv, color="#1f5aa6", alpha=0.2)
        ax1.plot(x, f_cum * e_conv, color="#c0392b", lw=1.5, ls="--", label="2nd-order cumulant")
        ax1.legend(frameon=False)
    else:
        ax1.plot(x, wmat[0] * e_conv, color="#1f5aa6", lw=2, label="W (single trajectory)")
        ax1.legend(frameon=False)
    ax1.set_xlabel(f"{xname} ({x_lab})")
    ax1.set_ylabel(f"Free energy / work ({e_lab})")
    ax1.set_title(f"SMD work, N = {nrep}")

    ax2.plot(ref_t / 1000, smean * x_conv, color="#1f5aa6", lw=1.2, label="CV  s(t)")
    ax2.plot(ref_t / 1000, lam * x_conv, color="k", lw=1, ls="--", label="restraint centre  λ(t)")
    ax2.set_xlabel("time (ps)")
    ax2.set_ylabel(f"CV ({x_lab})")
    ax2.set_title("CV lag behind the restraint")
    ax2.legend(frameon=False)

    fig.tight_layout()
    png = f"{args.prefix}.png"
    fig.savefig(png, dpi=args.dpi)
    print(f"# wrote {png} ({args.dpi} dpi)")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
