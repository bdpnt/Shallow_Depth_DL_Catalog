"""
depth_histogram_compare_report.py
=================================
Report figure: depth distribution of the same events before and after
relocation, on two stacked panels sharing both axes.

  (a) depths published by the source agencies (merged bulletin, .obs), blue;
  (b) relocated depths from a RESULT/*.csv (default column: expect_z, the
      PDF expectation published as the preferred origin), orange;
  (c) same relocation, maximum-likelihood hypocentres (default column:
      depth), green.
No titles; panel letters only.

Only events present in both files (matched on the permanent identifier)
are drawn, so that the two panels describe the same earthquakes.
Sibling of depth_histogram_report.py; that module is left untouched.

Usage
-----
    python complem_figures/depth_histogram_compare_report.py \\
        --before obs/GLOBAL.obs --after RESULT/SSST_result.csv \\
        --output complem_figures/depth_histogram/GLOBAL_vs_SSST_report.png \\
        --bin-size 0.2 --min-depth -3 --max-depth 25
"""

import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

_BAR_COLOR  = '#4C72B0'   # seaborn deep blue, as in the other report figures
_ROUND_COL  = '#DD8452'   # seaborn deep orange: relocated depths, panel (b)
_ML_COL     = '#55A868'   # seaborn deep green: maximum likelihood, panel (c)
_INK        = '#1A1A1A'
_INK_MUTED  = '#5A5A5A'
_GRID       = '#D9D9D9'


def read_obs_depths(path):
    """(publicId, depth) for every event of a .obs bulletin."""
    rows, dep = [], None
    with open(path) as f:
        for line in f:
            if line.startswith('# '):
                dep = float(line.split()[9])
            elif line.startswith('PUBLIC_ID') and dep is not None:
                rows.append((line.split()[1], dep))
                dep = None
    return pd.DataFrame(rows, columns=['publicId', 'dep_before'])


def fmt_int(n):
    return f'{n:,}'.replace(',', ' ')   # thin space as thousands separator


def generate_figure(before, after, output, column='expect_z', column_ml='depth',
                    max_err=None,
                    bin_size=0.2, min_depth=-3.0, max_depth=25.0,
                    also_pdf=True):
    b = read_obs_depths(before)
    a = pd.read_csv(after, usecols=['publicId', column, column_ml,
                                    'true_erh', 'true_erz']).rename(
        columns={column: 'dep_after', column_ml: 'dep_ml'})
    m = b.merge(a, on='publicId', how='inner').dropna()
    n = len(m)

    z0, z1 = m['dep_before'].to_numpy(), m['dep_after'].to_numpy()
    z2 = m['dep_ml'].to_numpy()
    if max_err is not None:
        # panels (b) and (c) only: events whose published horizontal AND
        # vertical uncertainties are both <= max_err km; (a) keeps them all
        keep = ((m['true_erh'] <= max_err) & (m['true_erz'] <= max_err)).to_numpy()
        z1, z2 = z1[keep], z2[keep]
        print(f'max_err {max_err:g} km: {keep.sum()} of {n} events kept in (b) and (c)')
    whole  = np.isclose(z0, np.round(z0))
    # rounded so that whole-km depths fall at the left of their own bin,
    # not on a floating-point edge that sends some of them to the bin below
    bins   = np.round(np.arange(min_depth, max_depth + bin_size / 2, bin_size), 6)

    mpl.rcParams.update({
        'font.size':        10,
        'axes.edgecolor':   _INK_MUTED,
        'axes.labelcolor':  _INK,
        'text.color':       _INK,
        'xtick.color':      _INK_MUTED,
        'ytick.color':      _INK_MUTED,
        'figure.facecolor': 'white',
        'axes.facecolor':   'white',
        'pdf.fonttype':     42,
        'legend.frameon':   False,
    })

    fig, axes = plt.subplots(3, 1, figsize=(6.5, 7.2), sharex=True,
                             sharey=True, layout='constrained')

    ax = axes[0]
    ax.hist(z0, bins=bins, color=_BAR_COLOR, edgecolor='white', linewidth=0.3)

    ax = axes[1]
    ax.hist(z1, bins=bins, color=_ROUND_COL, edgecolor='white', linewidth=0.3)

    ax = axes[2]
    ax.hist(z2, bins=bins, color=_ML_COL, edgecolor='white', linewidth=0.3)
    ax.set_xlabel('Depth (km)')

    for ax, letter in zip(axes, 'abc'):
        ax.set_ylabel('Number of events')
        ax.set_xlim(min_depth, max_depth)
        ax.set_xticks(np.arange(0, max_depth + 0.1, 5))
        ax.axvline(0, color=_INK_MUTED, linewidth=0.8,
                   linestyle=(0, (4, 3)), zorder=3)
        ax.yaxis.grid(True, color=_GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.text(-0.09, 1.0, letter, transform=ax.transAxes,
                fontsize=12, fontweight='bold', va='bottom', ha='left')
    axes[0].margins(y=0.02)

    os.makedirs(os.path.dirname(output) or '.', exist_ok=True)
    fig.savefig(output, dpi=300)
    outs = [output]
    if also_pdf:
        pdf = os.path.splitext(output)[0] + '.pdf'
        fig.savefig(pdf)
        outs.append(pdf)
    plt.close(fig)

    out_b = int(((z0 < min_depth) | (z0 > max_depth)).sum())
    out_a = int(((z1 < min_depth) | (z1 > max_depth)).sum())
    out_m = int(((z2 < min_depth) | (z2 > max_depth)).sum())
    print(f'N matched = {n}; medians {np.median(z0):.2f} / {np.median(z1):.2f} km; '
          f'whole-km depths before = {whole.sum()} ({100*whole.mean():.1f} %); '
          f'not drawn (outside {min_depth:g}..{max_depth:g} km): '
          f'before {out_b}, after {out_a}, max-likelihood {out_m}; '
          f'median max-likelihood {np.median(z2):.2f} km')
    for o in outs:
        print('Figure saved @', o)
    return m


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    p.add_argument('--before', default=os.path.join(_PROJECT_ROOT, 'obs', 'GLOBAL.obs'))
    p.add_argument('--after',  default=os.path.join(_PROJECT_ROOT, 'RESULT', 'SSST_result.csv'))
    p.add_argument('--column', default='expect_z')
    p.add_argument('--column-ml', default='depth',
                   help='maximum-likelihood depth column for panel (c)')
    p.add_argument('--output', default=os.path.join(
        _MODULE_DIR, 'depth_histogram', 'GLOBAL_vs_SSST_report.png'))
    p.add_argument('--bin-size',  type=float, default=0.2)
    p.add_argument('--min-depth', type=float, default=-3.0)
    p.add_argument('--max-depth', type=float, default=25.0)
    p.add_argument('--max-err', type=float, default=None,
                   help='keep in panels (b) and (c) only events with ERH and '
                        'ERZ both <= this value (km); panel (a) is unfiltered')
    p.add_argument('--no-pdf', action='store_true')
    a = p.parse_args()
    generate_figure(a.before, a.after, a.output, a.column, a.column_ml, a.max_err, a.bin_size,
                    a.min_depth, a.max_depth, not a.no_pdf)


if __name__ == '__main__':
    main()
