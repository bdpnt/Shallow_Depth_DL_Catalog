"""
depth_histogram.py
============================
Plot a histogram of earthquake depths from a .obs bulletin or a RESULT/*.csv.

Depth on the x-axis, count on the y-axis; N and median are carried by the
title, so nothing is drawn on top of the bars. Light figure style (thin bars,
horizontal grid only, no top/right spines). Writes both a 300 dpi PNG and a
vector PDF. Depths recorded as 'None' in a .obs bulletin are ignored.

For a RESULT/*.csv, `--solution expect` (default) reads the PDF expectation
depth `expect_z`; `--solution maxlike` the maximum-likelihood depth (see
`event_maps.hypocentre_columns`).

Usage
-----
    python complem_figures/depth_histogram.py \\
        --bulletin obs/GLOBAL.obs \\
        --output   complem_figures/depth_histogram/GLOBAL.png

    python complem_figures/depth_histogram.py \\
        --bulletin RESULT/SSST_result.csv --min-depth -3 \\
        --output   complem_figures/depth_histogram/SSST_result.png

    python complem_figures/depth_histogram.py \\
        --bulletin RESULT/SSST_result.csv --min-depth -3 --before-year 2020 \\
        --output   complem_figures/depth_histogram/SSST_result_before_2020.png
"""

import argparse
import os
import sys
from dataclasses import dataclass

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
from complem_figures.event_maps import hypocentre_columns  # noqa: E402

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
    solution:      str   = 'expect'   # csv input only: 'expect' | 'maxlike'
    before_year:   int | None = None  # keep only events with year < before_year
    also_pdf:      bool  = True


def _read_depths(file_bulletin, solution='expect', before_year=None):
    """Depths from a .obs bulletin, or from a RESULT/*.csv for the chosen solution.

    ``solution`` only applies to a .csv input; the .obs files carry whichever
    solution was published when they were written. ``before_year`` keeps only
    events strictly earlier than that year.
    """
    if file_bulletin.lower().endswith('.csv'):
        df = pd.read_csv(file_bulletin)
        if before_year is not None:
            df = df[pd.to_datetime(df['date-time']).dt.year < before_year]
        _, _, col = hypocentre_columns(df, solution)
        return df[col].dropna()

    with open(file_bulletin, 'r') as f:
        bull = [line.lstrip('# ').rstrip('\n').split()
                for line in f if line.startswith('# ')]
    bull_df = pd.DataFrame(bull, columns=[
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'Mag', 'MagType', 'MagAuthor',
        'PhaseCount', 'HorUncer', 'VerUncer', 'AzGap', 'RMS',
    ])
    if before_year is not None:
        bull_df = bull_df[pd.to_numeric(bull_df['Year']) < before_year]
    return pd.to_numeric(bull_df['Dep'], errors='coerce').dropna()


def generate_figure(parameters):
    depths = _read_depths(parameters.file_bulletin, parameters.solution,
                          parameters.before_year)

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
        + f', median = {median:.1f} km'
        + (f', before {parameters.before_year}'
           if parameters.before_year is not None else ''),
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
        description='Plot a histogram of earthquake depths.')
    parser.add_argument('--bulletin',
                        default=os.path.join(_PROJECT_ROOT, 'obs', 'GLOBAL.obs'))
    parser.add_argument('--output',
                        default=os.path.join(_MODULE_DIR, 'depth_histogram',
                                             'GLOBAL.png'))
    parser.add_argument('--bin-size',  type=float, default=1.0)
    parser.add_argument('--max-depth', type=float, default=40.0)
    parser.add_argument('--min-depth', type=float, default=0.0,
                        help='Left axis limit in km; use -3 for relocated '
                             'bulletins, whose search grid tops at -3 km')
    parser.add_argument('--solution', default='expect', choices=['expect', 'maxlike'],
                        help='Depth of a RESULT/*.csv: PDF expectation (default) '
                             'or maximum likelihood')
    parser.add_argument('--before-year', type=int, default=None,
                        help='Keep only events strictly before this year '
                             '(2020 keeps 1978-2019)')
    parser.add_argument('--no-pdf', action='store_true')
    args = parser.parse_args()

    generate_figure(DepthHistogramParams(
        file_bulletin=args.bulletin,
        fig_save=args.output,
        bin_size=args.bin_size,
        max_depth=args.max_depth,
        min_depth=args.min_depth,
        solution=args.solution,
        before_year=args.before_year,
        also_pdf=not args.no_pdf,
    ))


if __name__ == '__main__':
    main()
