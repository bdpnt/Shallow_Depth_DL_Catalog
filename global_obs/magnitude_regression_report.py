#!/usr/bin/env python3
"""
magnitude_regression_report.py
==============================
Report-quality version of the magnitude-conversion figure produced by
`global_obs/generate_magnitude_models.py`.

The original script writes one figure per magnitude pair.  This one builds the
three pairs into a single stacked three-panel figure (a, b, c) sharing the same
axis ranges, so the three regressions can be read against one another.

The matching and the piecewise ODR fit are reproduced verbatim from
`generate_magnitude_models.py` (KDTree candidate search, greedy closest-first
one-to-one matching, unconstrained ODR above M = 2, ODR constrained to be
continuous at M = 2 below it).  Only the plotting differs.

Input: per-catalogue event-header tables extracted from the PRE-CONVERSION
bulletins (columns: lat, lon, time, mag, magtype, magauthor).  The bulletins
in `obs/` have since been overwritten by `apply_magnitude_models.py`, whose
magnitudes are already ML LDG; the originals are the ones archived in
`obs/26-04-22_before.zip`.  Run with `--extract-from <dir>` once to build the
tables from those originals, then run without it to draw the figure.

Outputs: 300 dpi PNG + vector PDF (fonttype 42).
"""

import argparse
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# scipy / scikit-learn are imported lazily inside the fitting routines so that
# --extract-from works in a plain numpy + pandas environment.

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SRC = "mag_model/matched"   # override with --src
OUT = "mag_model/FIGURES"   # override with --out

DIST_THRESH = 10.0   # km   (default of generate_magnitude_models.py)
TIME_THRESH = 2.0    # s
BREAK = 2.0          # magnitude breakpoint

PANELS = [
    dict(tag="RESIF_MLv", label="ML$_\\mathrm{v}$ RESIF (BCSF-RéNaSS + OMP)",
         short="MLv RESIF"),
    dict(tag="IGN_mbLg",  label="m$_\\mathrm{b}$Lg IGN",        short="mb_Lg IGN"),
    dict(tag="ICGC_ML",   label="ML ICGC",              short="ML ICGC"),
]
TARGET = dict(tag="LDG_ML", label="ML LDG")

C_LOW  = "#DD8452"   # seaborn deep orange  — segment below the breakpoint
C_HIGH = "#4C72B0"   # seaborn deep blue    — segment at and above the breakpoint
C_GREY = "#8C8C8C"

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


# tag -> (bulletin basename, magnitude-type token searched in the header line)
SOURCES = {
    "RESIF_MLv": ("RESIF_20-25", "MLv"),
    "IGN_mbLg":  ("IGN_20-25",   "mb_Lg"),
    "ICGC_ML":   ("ICGC_20-25",  "ML"),
    "LDG_ML":    ("LDG_20-25",   "ML"),
}


def extract_headers(obs_dir, dest):
    """Write <tag>_headers.csv for each source bulletin in obs_dir."""
    import csv
    os.makedirs(dest, exist_ok=True)
    for tag, (base, mtype) in SOURCES.items():
        out = os.path.join(dest, f"{tag}_headers.csv")
        n = 0
        with open(os.path.join(obs_dir, base + ".obs"), "r", encoding="utf-8",
                  errors="ignore") as f, open(out, "w", newline="") as g:
            w = csv.writer(g)
            w.writerow(["lat", "lon", "time", "mag", "magtype", "magauthor"])
            for line in f:
                if (not line.startswith("#") or line.startswith("###")
                        or mtype not in line):
                    continue
                p = line.rstrip("\n").lstrip("# ").split()
                w.writerow([p[6], p[7],
                            f"{p[0]}-{p[1]}-{p[2]}T{p[3]}:{p[4]}:{p[5]}Z",
                            p[9], p[10], p[11]])
                n += 1
        print(f"{tag}: {n} events -> {out}")


def load(tag):
    df = pd.read_csv(os.path.join(SRC, f"{tag}_headers.csv"))
    df["time"] = pd.to_datetime(df["time"], format="ISO8601")
    return df.rename(columns={"lat": "latitude", "lon": "longitude",
                              "mag": "magnitude"})


def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dlat, dlon = np.radians(lat2 - lat1), np.radians(lon2 - lon1)
    a = (np.sin(dlat / 2) ** 2
         + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2))
         * np.sin(dlon / 2) ** 2)
    return 2 * R * np.arcsin(np.sqrt(a))


def match_events(cat1, cat2, dist_thresh=DIST_THRESH, time_thresh=TIME_THRESH):
    """Greedy one-to-one matching, identical in logic to the pipeline script."""
    from scipy.spatial import KDTree

    tree = KDTree(cat2[["latitude", "longitude"]].values)
    lat2 = cat2["latitude"].values
    lon2 = cat2["longitude"].values
    t2 = cat2["time"].values.astype("datetime64[ns]")
    m2 = cat2["magnitude"].values

    used2, pairs = set(), []
    for _, row in cat1.iterrows():
        _, cand = tree.query([row["latitude"], row["longitude"]], k=100)
        best, best_d, best_t = None, np.inf, np.inf
        for i in cand:
            if i in used2:
                continue
            d = haversine(row["latitude"], row["longitude"], lat2[i], lon2[i])
            if d > dist_thresh:
                continue
            dt = abs((np.datetime64(row["time"]) - t2[i]) / np.timedelta64(1, "s"))
            if dt > time_thresh:
                continue
            if d < best_d or (d == best_d and dt < best_t):
                best, best_d, best_t = i, d, dt
        if best is not None:
            used2.add(best)
            pairs.append((row["magnitude"], m2[best], best_d, best_t))
    return pd.DataFrame(pairs, columns=["m1", "m2", "dist_km", "dt_s"])


# ---------------------------------------------------------------------------
# Piecewise ODR fit
# ---------------------------------------------------------------------------


def _linear(p, x):
    return p[0] * x + p[1]


def fit_piecewise(m1, m2, xbreak=BREAK):
    from scipy.odr import Model, ODR, RealData
    from scipy.optimize import minimize
    from sklearn.metrics import r2_score

    hi = m1 >= xbreak
    lo = ~hi

    out = ODR(RealData(m1[hi], m2[hi]), Model(_linear), beta0=[1.0, 0.0]).run()
    s_hi, b_hi = out.beta
    y_break = s_hi * xbreak + b_hi

    X, Y = m1[lo], m2[lo]
    res = minimize(lambda p: np.sum((Y - (p[1] + p[0] * X)) ** 2),
                   x0=[1.0, 0.0], method="SLSQP",
                   constraints={"type": "eq",
                                "fun": lambda p: p[1] + p[0] * xbreak - y_break})
    s_lo, b_lo = res.x

    r2_hi = r2_score(m2[hi], s_hi * m1[hi] + b_hi)
    r2_lo = r2_score(m2[lo], s_lo * m1[lo] + b_lo)
    sd_hi = float(np.std(m2[hi] - (s_hi * m1[hi] + b_hi), ddof=1))
    sd_lo = float(np.std(m2[lo] - (s_lo * m1[lo] + b_lo), ddof=1))

    return dict(s_hi=s_hi, b_hi=b_hi, s_lo=s_lo, b_lo=b_lo,
                r2_hi=r2_hi, r2_lo=r2_lo, sd_hi=sd_hi, sd_lo=sd_lo,
                n_hi=int(hi.sum()), n_lo=int(lo.sum()), n=len(m1))


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------


def style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["TeX Gyre Heros", "DejaVu Sans"],
        "mathtext.fontset": "stixsans",
        "font.size": 8.5,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def main():
    global SRC, OUT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--src", default=SRC,
                    help="directory holding <TAG>_headers.csv tables")
    ap.add_argument("--out", default=OUT, help="output directory")
    ap.add_argument("--extract-from", default=None,
                    help="directory holding the pre-conversion .obs bulletins; "
                         "rebuild the header tables in --src and exit")
    args = ap.parse_args()
    SRC, OUT = args.src, args.out

    if args.extract_from:
        extract_headers(args.extract_from, SRC)
        return

    style()
    os.makedirs(OUT, exist_ok=True)

    target = load(TARGET["tag"])
    results = []
    for p in PANELS:
        src = load(p["tag"])
        mf = match_events(src, target)
        fit = fit_piecewise(mf["m1"].values, mf["m2"].values)
        results.append((p, mf, fit))
        print(f"{p['short']:>18s}  n={fit['n']:5d}  "
              f"hi: y={fit['s_hi']:.3f}x+{fit['b_hi']:+.3f} R2={fit['r2_hi']:.3f}  "
              f"lo: y={fit['s_lo']:.3f}x+{fit['b_lo']:+.3f} R2={fit['r2_lo']:.3f}")

    # Common axis limits across the three panels
    xs = np.concatenate([r[1]["m1"].values for r in results])
    ys = np.concatenate([r[1]["m2"].values for r in results])
    lo = np.floor(min(xs.min(), ys.min()) * 2) / 2 - 0.1
    hi = np.ceil(max(xs.max(), ys.max()) * 2) / 2 + 0.1

    fig, axes = plt.subplots(3, 1, figsize=(6.5, 7.1), sharex=True, sharey=True)

    for k, (ax, (p, mf, f)) in enumerate(zip(axes, results)):
        m1, m2 = mf["m1"].values, mf["m2"].values
        is_hi = m1 >= BREAK

        # 1:1 reference
        ax.plot([lo, hi], [lo, hi], color=C_GREY, lw=0.7, ls=(0, (4, 3)),
                zorder=1)
        # breakpoint
        ax.axvline(BREAK, color=C_GREY, lw=0.6, ls=(0, (1, 2.5)), zorder=1)

        # data
        ax.scatter(m1[~is_hi], m2[~is_hi], s=4, alpha=0.18, color=C_LOW,
                   linewidths=0, rasterized=True, zorder=2)
        ax.scatter(m1[is_hi], m2[is_hi], s=4, alpha=0.18, color=C_HIGH,
                   linewidths=0, rasterized=True, zorder=2)

        # fitted segments with ±1 sigma bands
        x_lo = np.linspace(m1.min(), BREAK, 50)
        x_hi = np.linspace(BREAK, m1.max(), 50)
        y_lo = f["s_lo"] * x_lo + f["b_lo"]
        y_hi = f["s_hi"] * x_hi + f["b_hi"]
        ax.fill_between(x_lo, y_lo - f["sd_lo"], y_lo + f["sd_lo"],
                        color=C_LOW, alpha=0.20, lw=0, zorder=3)
        ax.fill_between(x_hi, y_hi - f["sd_hi"], y_hi + f["sd_hi"],
                        color=C_HIGH, alpha=0.20, lw=0, zorder=3)
        ax.plot(x_lo, y_lo, color=C_LOW, lw=1.5, solid_capstyle="round", zorder=4)
        ax.plot(x_hi, y_hi, color=C_HIGH, lw=1.5, solid_capstyle="round", zorder=4)

        # equations, bottom right (empty corner)
        txt_hi = (f"$M \\geq 2$:  $y = {f['s_hi']:.3f}\\,x {f['b_hi']:+.3f}$\n"
                  f"$R^2 = {f['r2_hi']:.3f}$,  $\\sigma = {f['sd_hi']:.3f}$,  "
                  f"$n = {f['n_hi']}$")
        txt_lo = (f"$M < 2$:  $y = {f['s_lo']:.3f}\\,x {f['b_lo']:+.3f}$\n"
                  f"$R^2 = {f['r2_lo']:.3f}$,  $\\sigma = {f['sd_lo']:.3f}$,  "
                  f"$n = {f['n_lo']}$")
        ax.text(0.97, 0.30, txt_hi, transform=ax.transAxes, ha="right",
                va="bottom", color=C_HIGH, fontsize=9.2, linespacing=1.45)
        ax.text(0.97, 0.03, txt_lo, transform=ax.transAxes, ha="right",
                va="bottom", color=C_LOW, fontsize=9.2, linespacing=1.45)

        # panel letter only
        ax.text(0.012, 0.975, "abc"[k], transform=ax.transAxes,
                ha="left", va="top", fontsize=11, fontweight="bold")

        ax.set_xlabel(p["label"])
        ax.set_ylabel(TARGET["label"])
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xticks(np.arange(0, 5.5, 1.0))
        ax.set_yticks(np.arange(0, 5.5, 1.0))
        ax.tick_params(labelbottom=True, labelleft=True)

    fig.subplots_adjust(left=0.082, right=0.995, top=0.985, bottom=0.075,
                        hspace=0.24)

    png = os.path.join(OUT, "magnitude_conversion_regressions.png")
    pdf = os.path.join(OUT, "magnitude_conversion_regressions.pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)
    print("wrote", png, "and", pdf)


if __name__ == "__main__":
    main()
