"""
depth_histogram_report.py
=========================
Report-ready variant of ``depth_histogram.py``.

Differences with the original module (which is left untouched):
  * depth on the x-axis, count on the y-axis (no inverted axis);
  * N and median are carried by the title, so nothing is drawn on top
    of the bars;
  * light figure style (no seaborn grey panel), thin bars, horizontal
    grid only, no top/right spines;
  * writes both a 300 dpi PNG and a vector PDF for inclusion in the report.

Usage
-----
    python complem_figures/depth_histogram_report.py \\
        --bulletin obs/GLOBAL.obs \\
        --output   complem_figures/depth_histogram/GLOBAL_report.png
"""

import argparse
import os
from dataclasses import dataclass

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

_BAR_COLOR  = '#4C72B0'   # same blue as the seaborn deep palette used elsewhere
_INK        = '#1A1A1A'
_INK_MUTED  = '#5A5A5A'
_GRID       = '#D9D9D9'


@dataclass
class DepthHistogramParams:
    file_bulletin: str
    fig_save:      str
    bin_size:      float = 1.0
    max_depth:     float = 40.0
    min_depth:     float = 0.0
    column:        str   = None   # csv input only, e.g. 'expect_z'
    also_pdf:      bool  = True


def _read_depths(file_bulletin, column=None):
    """Depths from a .obs bulletin, or from a RESULT/*.csv column.

    ``column`` (e.g. 'expect_z') only applies to a .csv input; the .obs
    files carry whichever solution was published when they were written.
    """
    if file_bulletin.lower().endswith('.csv'):
        col = column or 'depth'
        return pd.read_csv(file_bulletin, usecols=[col])[col].dropna()

    with open(file_bulletin, 'r') as f:
        bull = [line.lstrip('# ').rstrip('\n').split()
                for line in f if line.startswith('# ')]
    bull_df = pd.DataFrame(bull, columns=[
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'Mag', 'MagType', 'MagAuthor',
        'PhaseCount', 'HorUncer', 'VerUncer', 'AzGap', 'RMS',
    ])
    return pd.to_numeric(bull_df['Dep'], errors='coerce').dropna()


def generate_figure(parameters):
    depths = _read_depths(parameters.file_bulletin, parameters.column)

    n_total  = len(depths)
    median   = depths.median()
    shown    = depths[(depths >= parameters.min_depth)
                      & (depths <= parameters.max_depth)]
    n_beyond = n_total - len(shown)

    bins = np.arange(parameters.min_depth,
                     parameters.max_depth + parameters.bin_size,
                     parameters.bin_size)

    os.makedirs(os.path.dirname(parameters.fig_save) or '.', exist_ok=True)

    mpl.rcParams.update({
        'font.size':        11,
        'axes.edgecolor':   _INK_MUTED,
        'axes.labelcolor':  _INK,
        'text.color':       _INK,
        'xtick.color':      _INK_MUTED,
        'ytick.color':      _INK_MUTED,
        'figure.facecolor': 'white',
        'axes.facecolor':   'white',
        'pdf.fonttype':     42,
    })

    fig, ax = plt.subplots(figsize=(7.0, 3.6), layout='constrained')

    ax.hist(shown, bins=bins, color=_BAR_COLOR,
            edgecolor='white', linewidth=0.6)

    ax.set_xlabel('Depth (km)')
    ax.set_ylabel('Number of events')
    ax.set_xlim(parameters.min_depth, parameters.max_depth)
    ax.set_xticks(np.arange(0, parameters.max_depth + 1, 5))
    if parameters.min_depth < 0:
        # sea level, the reference the search grid is hung from
        ax.axvline(0, color=_INK_MUTED, linewidth=0.8, linestyle=(0, (4, 3)),
                   zorder=3)
    ax.margins(y=0.02)

    ax.yaxis.grid(True, color=_GRID, linewidth=0.6)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_linewidth(0.8)

    ax.set_title(
        f'Depth distribution  —  N = {n_total:,}'.replace(',', ' ')
        + f', median = {median:.1f} km',
        pad=10)

    fig.savefig(parameters.fig_save, dpi=300)
    outputs = [parameters.fig_save]
    if parameters.also_pdf:
        pdf_path = os.path.splitext(parameters.fig_save)[0] + '.pdf'
        fig.savefig(pdf_path)
        outputs.append(pdf_path)
    plt.close(fig)

    print(f'N = {n_total}, median = {median:.2f} km, '
          f'events outside {parameters.min_depth:g}-{parameters.max_depth:g} km '
          f'(not drawn) = {n_beyond}')
    for o in outputs:
        print(f'Figure saved @ {o}')
    return {'output': outputs, 'n': n_total, 'median': median,
            'n_outside': n_beyond}


def main():
    parser = argparse.ArgumentParser(
        description='Report figure: depth histogram, depth on x-axis.')
    parser.add_argument('--bulletin',
                        default=os.path.join(_PROJECT_ROOT, 'obs', 'GLOBAL.obs'))
    parser.add_argument('--output',
                        default=os.path.join(_MODULE_DIR, 'depth_histogram',
                                             'GLOBAL_report.png'))
    parser.add_argument('--bin-size',  type=float, default=1.0)
    parser.add_argument('--max-depth', type=float, default=40.0)
    parser.add_argument('--min-depth', type=float, default=0.0,
                        help='Left axis limit in km; use -3 for relocated '
                             'bulletins, whose search grid tops at -3 km')
    parser.add_argument('--column', default=None,
                        help="Column to read when --bulletin is a RESULT/*.csv "
                             "(e.g. expect_z for the PDF expectation)")
    parser.add_argument('--no-pdf', action='store_true')
    args = parser.parse_args()

    generate_figure(DepthHistogramParams(
        file_bulletin=args.bulletin,
        fig_save=args.output,
        bin_size=args.bin_size,
        max_depth=args.max_depth,
        min_depth=args.min_depth,
        column=args.column,
        also_pdf=not args.no_pdf,
    ))


if __name__ == '__main__':
    main()
