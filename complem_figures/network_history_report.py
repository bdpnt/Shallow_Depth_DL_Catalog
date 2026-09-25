#!/usr/bin/env python3
"""Report figure: seismicity, station availability and temporary deployments through time.

Three panels on a shared time axis:
  (a) relocated seismicity, events per year
  (b) number of stations seen by the catalog each year, with the permanent
      backbone separated out, and the inventory-declared count as a lower bound
  (c) time span of the temporary / campaign deployments of the Pyrenees

New script (see pipeline/10_figures.md conventions); no existing figure script is modified.

CLI:
  python complem_figures/network_history_report.py \
      --catalog RESULT/SSST_result.csv \
      --obs obs/SSST_result.obs \
      --inventory stations/GLOBAL_inventory.xml \
      --output complem_figures/network_history/network_history_report.png
"""

import argparse
import collections
import datetime as dt
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors
import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.size"] = 11.5
matplotlib.rcParams["axes.labelsize"] = 12
matplotlib.rcParams["xtick.labelsize"] = 11
matplotlib.rcParams["ytick.labelsize"] = 11
matplotlib.rcParams["ps.fonttype"] = 42

BLUE = "#4C72B0"        # repo-wide seaborn deep blue
DARK = "#2F3B52"
GREY = "#9AA3B2"
ORANGE = "#DD8452"
GREEN = "#55A868"

PERMANENT = {"FR", "RA", "RD", "ES", "CA", "LC", "G", "GE", "IU", "MN", "WM"}
PLACEHOLDER_START = 1900      # startDate="1900-01-01" means "unknown" in this inventory
OPEN_END = 2400               # endDate 2500-12-31 / 2599-12-31 means "still running"

# ---------------------------------------------------------------------------
# Temporary / campaign deployments of the Pyrenees.
#   prec = "day"   -> window taken from station epochs or from a dated publication
#   prec = "year"  -> only the calendar year(s) are known (Tab. 2.1)
#   family: "largeN" | "classic" | "semiperm"
# Deployments outside the range (MACIV XF, AdriaArray XP, K3 XI/XR, TRACK YO,
# AlpArray Z3, OMIV 1N/3I) are deliberately excluded.
# ---------------------------------------------------------------------------
DEPLOYMENTS = [
    # label,                                   start,        end,          prec,   family,     net
    ("Bigorre (XX, 25 st.)",                   "2001-01-01", "2002-12-31", "year", "classic",  "XX"),
    ("Saint-Paul-de-Fenouillet (XX, 2 st.)",   "2004-01-01", "2004-12-31", "year", "classic",  "XX"),
    ("Argelès-Gazost (XX, 12 st.)",            "2006-01-01", "2006-12-31", "year", "classic",  "XX"),
    ("IberArray (IB, 16 st.)",                 "2010-11-08", "2014-04-24", "day",  "classic",  "IB"),
    ("PYROPE (X7, 113 st.)",                   "2010-11-09", "2014-10-22", "day",  "classic",  "X7"),
    ("Aucun (XX, 3 st.)",                      "2013-01-01", "2013-12-31", "year", "classic",  "XX"),
    ("MISTERIOS (2M, 27 st.)",                 "2015-01-26", "2018-11-20", "day",  "classic",  "2M"),
    ("OROGEN-X (ZU, 47 st.)",                  "2015-03-30", "2017-02-12", "day",  "classic",  "ZU"),
    ("MAUPASACQ (XD, 432 st.)",                "2017-04-01", "2017-10-31", "day",  "largeN",   "XD"),
    ("RaspberryShake (AM, 68 st.)",            "2019-07-01", None,         "day",  "semiperm", "AM"),
    ("POCRISC (YK, 14 st.)",                   "2019-04-03", "2020-12-16", "day",  "classic",  "YK"),
    ("SANIMS / IMAGYN (YS, 211 st.)",          "2019-09-04", "2022-10-19", "day",  "classic",  "YS"),
    ("Campan–Gripp (XX, 9 st.)",               "2020-01-01", "2024-12-31", "year", "classic",  "XX"),
    ("SISLACQ (5M, 12 st.)",                   "2021-01-11", "2023-12-31", "day",  "classic",  "5M"),
    ("RASPY (4I, 14 st.)",                     "2021-03-23", None,         "day",  "semiperm", "4I"),
    ("Arette cluster / Stryde (8M, 217 st.)",  "2022-05-21", "2022-06-18", "day",  "largeN",   "8M"),
    ("Vielha–Maladeta / Aneto (20, 12 st.)",   "2023-08-01", "2025-01-31", "day",  "classic",  "20"),
    ("TEMPYR (5I, 7 st.)",                     "2023-09-01", None,         "day",  "classic",  "5I"),
    ("STGO, pre-PISCO (23, 41 st.)",           "2023-10-01", "2023-12-31", "day",  "classic",  "23"),
    ("PISCO (24, 323 st.)",                    "2024-05-21", "2024-07-13", "day",  "largeN",   "24"),
    ("NAFOSSA / Arette shallow (6E, 58 st.)",  "2025-02-12", "2025-03-19", "day",  "largeN",   "6E"),
    # Permanent networks, drawn as a separate block at the bottom of panel c.
    # Start years are the network start years of Tab. 2.1 (all still running);
    # FR and RD begin before the axis and are marked with a left arrow.
    ("RéNaSS–RLBP (FR, 71 st.)",               "1962-01-01", None,         "year", "perm",     "FR"),
    ("CEA-DASE RLBP (RD, 39 st.)",             "1962-01-01", None,         "year", "perm",     "RD"),
    ("Catalan network (CA, 73 st.)",           "1984-01-01", None,         "year", "perm",     "CA"),
    ("RAP accelerometric (RA, 30 st.)",        "1995-01-01", None,         "year", "perm",     "RA"),
    ("Spanish network (ES, 95 st.)",           "1999-01-01", None,         "year", "perm",     "ES"),
    ("LSC Canfranc (LC, 1 st.)",               "2011-01-01", None,         "year", "perm",     "LC"),
]
N_PERM = 6      # trailing rows of DEPLOYMENTS that are permanent networks

FAMILY_COLOR = {"largeN": ORANGE, "classic": BLUE, "semiperm": GREEN, "perm": DARK}
FAMILY_LABEL = {
    "largeN": "large-N nodal campaign",
    "classic": "temporary campaign",
    "semiperm": "semi-permanent / low-cost",
    "perm": "permanent network",
}
# deployments highlighted by a guide line across the three panels
GUIDES = {"XD", "8M", "24"}


def dec_year(date):
    """Decimal year of a datetime.date."""
    y = date.year
    start = dt.date(y, 1, 1).toordinal()
    end = dt.date(y + 1, 1, 1).toordinal()
    return y + (date.toordinal() - start) / (end - start)


def events_per_year(path):
    df = pd.read_csv(path, usecols=["date-time"])
    yr = pd.to_datetime(df["date-time"], format="mixed").dt.year
    return yr.value_counts().sort_index(), len(df)


def stations_per_year_from_obs(path):
    """Distinct station codes carrying at least one arrival, per event year."""
    year = None
    allsta = collections.defaultdict(set)
    perm = collections.defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith("# 1") or line.startswith("# 2"):
                year = int(line.split()[1])
                continue
            if line.startswith("#") or line.startswith("PUBLIC_ID") or not line.strip():
                continue
            code = line.split(None, 1)[0]
            if "." not in code:
                continue
            allsta[year].add(code)
            if code.split(".")[0] in PERMANENT:
                perm[year].add(code)
    return ({y: len(s) for y, s in allsta.items()},
            {y: len(perm.get(y, ())) for y in allsta})


def declared_per_year(path, last_year):
    """Stations declared active per year, counted on unified codes.

    Only codes carrying a real start date are counted: 1900-01-01 is the
    inventory's placeholder for an unknown date, and several deployments carry
    no date at all, so this curve is a lower bound (see figure caption).
    """
    txt = Path(path).read_text()
    codes = {}
    for m in re.finditer(r"<Station ([^>]+)>", txt):
        attrs = m.group(1)
        ac = re.search(r'alternateCode="([^"]+)"', attrs)
        if not ac:
            continue
        ac = ac.group(1)
        s = re.search(r'startDate="(\d{4})', attrs)
        e = re.search(r'endDate="(\d{4})', attrs)
        s = int(s.group(1)) if s else None
        if s is not None and s <= PLACEHOLDER_START:
            s = None
        e = int(e.group(1)) if e else None
        if e is not None and e >= OPEN_END:
            e = None
        cur = codes.get(ac)
        if cur is None:
            codes[ac] = [s, e, e is None]
        else:  # widen the epoch, as the merge itself does
            if s is not None:
                cur[0] = s if cur[0] is None else min(cur[0], s)
            if e is None:
                cur[2] = True
            elif not cur[2]:
                cur[1] = e if cur[1] is None else max(cur[1], e)
    counts = collections.Counter()
    n_undated = 0
    for s, e, is_open in codes.values():
        if s is None:
            n_undated += 1
            continue
        end = last_year if is_open or e is None else min(e, last_year)
        for y in range(s, end + 1):
            counts[y] += 1
    return counts, len(codes), n_undated


def step_xy(counts, y0, y1):
    """Year-binned counts as a staircase (value held over the calendar year)."""
    xs, ys = [], []
    for y in range(y0, y1 + 1):
        v = counts.get(y, 0)
        xs += [y, y + 1]
        ys += [v, v]
    return np.array(xs), np.array(ys)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalog", default="RESULT/SSST_result.csv")
    ap.add_argument("--obs", default="obs/SSST_result.obs")
    ap.add_argument("--inventory", default="stations/GLOBAL_inventory.xml")
    ap.add_argument("--output", default="complem_figures/network_history/network_history_report.png")
    ap.add_argument("--year-min", type=int, default=1978)
    ap.add_argument("--year-max", type=int, default=2025)
    ap.add_argument("--no-declared", action="store_true",
                    help="drop the inventory-declared curve from panel (b)")
    ap.add_argument("--no-guides", action="store_true",
                    help="drop the vertical guide lines of the large-N campaigns")
    args = ap.parse_args()

    y0, y1 = args.year_min, args.year_max

    ev, n_ev = events_per_year(args.catalog)
    contrib, contrib_perm = stations_per_year_from_obs(args.obs)
    declared, n_codes, n_undated = declared_per_year(args.inventory, y1)

    nets_with_picks = set()
    with open(args.obs) as fh:
        for line in fh:
            if line[:1] in "#P" or not line.strip():
                continue
            code = line.split(None, 1)[0]
            if "." in code:
                nets_with_picks.add(code.split(".")[0])

    print(f"events {n_ev}  years {ev.index.min()}-{ev.index.max()}")
    print(f"stations contributing arrivals: max {max(contrib.values())} "
          f"in {max(contrib, key=contrib.get)}")
    print(f"unified codes {n_codes}, of which undated (or placeholder-dated) {n_undated}")

    # ---------------------------------------------------------------- layout
    nrow = len(DEPLOYMENTS)
    fig = plt.figure(figsize=(9.2, 10.7))
    gs = fig.add_gridspec(3, 1, height_ratios=[2.0, 2.0, 0.25 * nrow],
                          hspace=0.055, left=0.325, right=0.975, top=0.975, bottom=0.105)
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1], sharex=axA)
    axC = fig.add_subplot(gs[2], sharex=axA)

    for ax in (axA, axB, axC):
        ax.set_xlim(y0, y1 + 1)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.set_facecolor("white")
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
        ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
        ax.xaxis.grid(True, color="#EEF0F3", lw=0.7)
    fig.patch.set_facecolor("white")

    # guide lines for the large-N campaigns
    if not args.no_guides:
        for label, s, e, prec, fam, net in DEPLOYMENTS:
            if net not in GUIDES:
                continue
            xs = dec_year(dt.date.fromisoformat(s))
            xe = dec_year(dt.date.fromisoformat(e)) if e else y1 + 1
            xe = max(xe, xs + 0.32)      # keep a one-month campaign visible
            for ax in (axA, axB, axC):
                ax.axvspan(xs, xe, color=ORANGE, alpha=0.16, lw=0, zorder=0)

    # ------------------------------------------------------------- panel (a)
    years = np.arange(y0, y1 + 1)
    counts = np.array([ev.get(y, 0) for y in years])
    axA.bar(years + 0.5, counts, width=0.86, color=BLUE, edgecolor="white", lw=0.6, zorder=3)
    axA.set_ylabel("Events per year")
    axA.yaxis.grid(True, color="#E3E6EB", lw=0.8)
    axA.set_axisbelow(True)
    axA.tick_params(labelbottom=False)

    # ------------------------------------------------------------- panel (b)
    x, y = step_xy(contrib, y0, y1)
    axB.fill_between(x, y, step=None, color=DARK, alpha=0.10, lw=0, zorder=2)
    axB.plot(x, y, color=DARK, lw=1.6, zorder=4,
             label="stations carrying ≥ 1 arrival in the catalog")
    x, y = step_xy(contrib_perm, y0, y1)
    axB.plot(x, y, color=BLUE, lw=1.4, zorder=4,
             label="of which permanent networks (FR, RA, RD, ES, CA, LC)")
    if not args.no_declared:
        x, y = step_xy(declared, y0, y1)
        axB.plot(x, y, color=GREY, lw=1.2, ls="--", zorder=3,
                 label="declared active in the unified inventory")
    axB.set_ylabel("Number of stations")
    axB.yaxis.grid(True, color="#E3E6EB", lw=0.8)
    axB.set_axisbelow(True)
    axB.tick_params(labelbottom=False)
    axB.legend(loc="upper left", frameon=False, fontsize=9.6, handlelength=2.4)

    # ------------------------------------------------------------- panel (c)
    for i, (label, s, e, prec, fam, net) in enumerate(DEPLOYMENTS):
        yi = nrow - 1 - i
        xs = dec_year(dt.date.fromisoformat(s))
        starts_before = xs < y0
        xs = max(xs, y0)                     # never extend the axis leftwards
        open_ended = e is None
        xe = dec_year(dt.date.fromisoformat(e)) if e else y1 + 0.55
        width = max(xe - xs, 0.22)          # keep a one-month campaign visible
        contributes = net in nets_with_picks
        col = FAMILY_COLOR[fam]
        # `prec` ("day" vs "year") stays in the table and is documented in the
        # caption, but is no longer drawn: the five year-resolution deployments
        # (Bigorre, Saint-Paul-de-Fenouillet, Argeles-Gazost, Aucun,
        # Campan-Gripp) are drawn like the others.
        face = col if contributes else "white"
        axC.barh(yi, width, left=xs, height=0.58,
                 color=face, edgecolor=col, lw=1.1, zorder=3)
        if open_ended:      # still running: arrow past the right-hand edge
            axC.annotate("", xy=(y1 + 1.5, yi), xytext=(y1 + 0.6, yi),
                         arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2),
                         annotation_clip=False, zorder=4)
        if starts_before:   # started before 1978: arrow past the left-hand edge
            axC.annotate("", xy=(y0 - 0.55, yi), xytext=(y0 + 0.15, yi),
                         arrowprops=dict(arrowstyle="-|>", color=col, lw=1.2),
                         annotation_clip=False, zorder=4)
    axC.axhline(N_PERM - 0.5, color="#C9CFD8", lw=0.9, zorder=2)
    axC.set_yticks(range(nrow))
    axC.set_yticklabels([d[0] for d in DEPLOYMENTS][::-1], fontsize=9.6)
    axC.set_ylim(-0.7, nrow - 0.3)
    axC.set_xlabel("Year")
    axC.xaxis.grid(True, color="#E3E6EB", lw=0.8)
    axC.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
    axC.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1))
    axC.set_axisbelow(True)
    axC.spines["left"].set_visible(False)
    axC.tick_params(axis="y", length=0)

    handles = [Patch(facecolor=FAMILY_COLOR[f], edgecolor=FAMILY_COLOR[f], label=FAMILY_LABEL[f])
               for f in ("largeN", "classic", "semiperm", "perm")]
    handles += [
        Patch(facecolor="white", edgecolor=DARK, label="no arrival under this network's codes"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.004),
               ncol=3, frameon=False, fontsize=9.6, handlelength=1.9, columnspacing=2.2)

    # panel letters (journal style): bold, no parentheses, no panel titles,
    # placed at the far left of the figure, level with each panel's top edge
    for ax, letter in ((axA, "a"), (axB, "b"), (axC, "c")):
        fig.text(0.018, ax.get_position().y1, letter, fontsize=14,
                 fontweight="bold", ha="left", va="top")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
    print("wrote", out, "and", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
