"""
error_vs_quality_report.py
==========================
Report-ready version of the `<run>_error_vs_quality.pdf` figure written by
`NLL_run/pdf_metrics.py`, restricted to the four classical quality
indicators: published ERH (top row) and ERZ (bottom row) against RMS, azimuthal
gap, phase count and nearest-station distance.

Differences from `pdf_metrics._generate_error_figure`
-----------------------------------------------------
* Four columns (RMS, Gap, Nphs, Dist) instead of six plus the dip-test boxes.
* Drawn at its printed size (`--width-cm`, default 16 cm = `\\linewidth` of an
  A4 page with 2.5 cm margins) with 6.5-8 pt fonts, so LaTeX includes it at
  scale 1 — the original is ~29 x 7 in.
* Light style, no suptitle (the legend belongs in the caption), thinner
  median line, smaller statistics box; y labels on the first column only,
  x labels on the bottom row only; the log error axis is shared along each row
  as in the original.
* Everything computed is the original's: the 0.5-99.5 % x clip, the shared
  row y-range (`_error_ylim`), the log-count hexbin, the median + IQR over 20
  equal-count bins (`_quantile_bins`), Spearman rho and the Nphs-controlled
  partial rho (`_partial_spearman`) — imported, not copied.
* Writes both a 300 dpi PNG and a vector PDF (`pdf.fonttype = 42`).

Usage
-----
    python complem_figures/error_vs_quality_report.py \\
        --csv RESULT/SSST_result.csv \\
        --output complem_figures/pdf_metrics/ssst_run1_error_vs_quality_report.png
"""

import argparse
import importlib.util
import os
from dataclasses import dataclass
from typing import Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

# Loaded by path: importing it through the NLL_run package would run
# NLL_run/__init__.py and pull in the whole relocation stack.
_spec = importlib.util.spec_from_file_location(
    '_pdf_metrics', os.path.join(_PROJECT_ROOT, 'NLL_run', 'pdf_metrics.py'))
pm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pm)

_ROWS = pm._ERROR_ROWS                                   # ERH, ERZ
_DEFAULT_PREDICTORS = ('RMS', 'Gap', 'Nphs', 'Dist')

# the original's labels, 'Dist' shortened so it fits under a 3.5 cm panel
_LABELS = {p: ('Nearest station (km)' if p == 'Dist' else label, scale)
           for p, label, scale in pm._ERROR_PREDICTORS}


def _columns(predictors):
    unknown = [p for p in predictors if p not in _LABELS]
    if unknown:
        raise SystemExit(f'unknown predictor(s) {unknown}; available: {sorted(_LABELS)}')
    return [(p, _LABELS[p][0], _LABELS[p][1]) for p in predictors]

_INK       = '#1A1A1A'
_INK_MUTED = '#5A5A5A'
_BLUE      = '#1f77b4'   # tab:blue, as the original
_CM        = 1 / 2.54


@dataclass
class ErrorVsQualityReportParams:
    csv:       str = 'RESULT/SSST_result.csv'
    figSave:   Optional[str] = None
    width_cm:  float = 16.0
    height_cm: float = 8.0
    predictors: tuple = _DEFAULT_PREDICTORS


def _style():
    mpl.rcParams.update({
        'font.family':       'sans-serif',
        'font.size':         7,
        'axes.labelsize':    8,
        'xtick.labelsize':   7,
        'ytick.labelsize':   7,
        'axes.edgecolor':    _INK_MUTED,
        'axes.linewidth':    0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.minor.width': 0.4,
        'ytick.minor.width': 0.4,
        'xtick.major.size':  2.5,
        'ytick.major.size':  2.5,
        'xtick.minor.size':  1.5,
        'ytick.minor.size':  1.5,
        'xtick.major.pad':   2,
        'ytick.major.pad':   2,
        'text.color':        _INK,
        'axes.labelcolor':   _INK,
        'xtick.color':       _INK,
        'ytick.color':       _INK,
        'pdf.fonttype':      42,
        'ps.fonttype':       42,
        'savefig.dpi':       300,
    })


def _emptiest_corner(ax, x, y, curve, box=(0.45, 0.30)):
    """
    The axes corner ('TR', 'TL', 'BR', 'BL') where the statistics box hides the
    least: fewest events inside a box-sized rectangle, with any point of the
    median curve or of its IQR band counted as heavily as the whole catalogue,
    so the box never lands on the curve if a curve-free corner exists.
    """
    to_axes = ax.transData + ax.transAxes.inverted()
    events = to_axes.transform(np.column_stack((x, y)))
    line = (to_axes.transform(np.column_stack(curve)) if len(curve[0])
            else np.empty((0, 2)))
    bw, bh = box
    best, best_score = 'TR', np.inf
    for corner in ('TR', 'TL', 'BR', 'BL'):          # ties keep the original top-right
        x0, x1 = (1 - bw, 1) if corner[1] == 'R' else (0, bw)
        y0, y1 = (1 - bh, 1) if corner[0] == 'T' else (0, bh)
        def inside(pts):
            return int(np.sum((pts[:, 0] >= x0) & (pts[:, 0] <= x1)
                              & (pts[:, 1] >= y0) & (pts[:, 1] <= y1)))
        score = inside(events) + inside(line) * len(x)
        if score < best_score:
            best, best_score = corner, score
    return best


def _panel(ax, data, error_col, predictor, log_x, y_lim):
    """pdf_metrics._plot_error_panel at report size: same selection, statistics and marks."""
    control_col = pm._ERROR_CONTROL
    columns = list(dict.fromkeys([predictor, error_col, control_col]))
    panel = data[columns].dropna(subset=[predictor, error_col])
    panel = panel[panel[error_col] > 0]
    if log_x:
        panel = panel[panel[predictor] > 0]

    x_lo, x_hi = np.nanpercentile(panel[predictor], pm._CLIP_PCT)
    panel = panel[(panel[predictor] >= x_lo) & (panel[predictor] <= x_hi)]
    x = panel[predictor].to_numpy(dtype=float)
    y = panel[error_col].to_numpy(dtype=float)
    y_lo, y_hi = y_lim
    extent = (((np.log10(x_lo), np.log10(x_hi)) if log_x else (x_lo, x_hi))
              + (np.log10(y_lo), np.log10(y_hi)))

    ax.hexbin(x, y, gridsize=45, bins='log', cmap='Greys', mincnt=1,
              xscale='log' if log_x else 'linear', yscale='log',
              extent=extent, zorder=1, linewidths=0, rasterized=True)

    curve = (np.array([]), np.array([]))
    bin_index, centres = pm._quantile_bins(x)
    if bin_index is not None:
        binned = pd.DataFrame({'bin': bin_index, 'error': y}).groupby('bin')['error']
        centre = centres[binned.median().index.to_numpy()]
        ax.fill_between(centre, binned.quantile(0.25), binned.quantile(0.75),
                        color=_BLUE, alpha=0.25, lw=0, zorder=3)
        ax.plot(centre, binned.median(), color=_BLUE, lw=1.1, zorder=4)
        # the median curve and both IQR edges, which the statistics box must not cover
        curve = (np.concatenate([centre] * 3),
                 np.concatenate([binned.median(), binned.quantile(0.25), binned.quantile(0.75)]))

    rho = spearmanr(x, y).statistic
    lines = [f'$\\rho$ = {rho:+.2f}']
    stats = {'rho': rho, 'n': len(x)}
    if predictor != control_col:
        control = panel[control_col].to_numpy(dtype=float)
        finite = np.isfinite(control)
        if finite.sum() > 2:
            stats['rho_partial'] = pm._partial_spearman(x[finite], y[finite], control[finite])
            lines.append(f'$\\rho\\,|\\,${control_col} = {stats["rho_partial"]:+.2f}')
    lines.append(f'N = {len(x):,}'.replace(',', ' '))
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    corner = _emptiest_corner(ax, x, y, curve)
    ha = 'right' if corner[1] == 'R' else 'left'
    va = 'top' if corner[0] == 'T' else 'bottom'
    ax.text(0.97 if ha == 'right' else 0.03, 0.97 if va == 'top' else 0.03,
            '\n'.join(lines), transform=ax.transAxes, ha=ha, va=va,
            fontsize=5.5, linespacing=1.15,
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.8, lw=0))

    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    return stats


def generate_figure(parameters):
    csv = parameters.csv if os.path.isabs(parameters.csv) else os.path.join(_PROJECT_ROOT, parameters.csv)
    data = pd.read_csv(csv)
    columns = _columns(parameters.predictors)
    missing = [c for c in [r[0] for r in _ROWS] + [c[0] for c in columns] if c not in data.columns]
    if missing:
        raise SystemExit(f'{csv} lacks {missing}')
    _style()

    fig, axes = plt.subplots(len(_ROWS), len(columns), squeeze=False, layout='constrained',
                             figsize=(parameters.width_cm * _CM, parameters.height_cm * _CM))
    fig.get_layout_engine().set(w_pad=0.02, h_pad=0.03, wspace=0.03, hspace=0.04)

    for r, (error_col, ylabel) in enumerate(_ROWS):
        y_lim = pm._error_ylim(data[error_col].to_numpy(dtype=float))
        for c, (predictor, xlabel, scale) in enumerate(columns):
            ax = axes[r][c]
            stats = _panel(ax, data, error_col, predictor, scale == 'log', y_lim)
            if r == len(_ROWS) - 1:
                ax.set_xlabel(xlabel)
            if c == 0:
                ax.set_ylabel(ylabel)
            else:
                ax.tick_params(labelleft=False)   # the y-range is shared along the row
            print(f'{ylabel:9s} vs {predictor:5s}: rho {stats["rho"]:+.2f}'
                  + (f', rho|Nphs {stats["rho_partial"]:+.2f}' if 'rho_partial' in stats else '')
                  + f', N {stats["n"]}')

    if not parameters.figSave:
        parameters.figSave = os.path.join(_MODULE_DIR, 'pdf_metrics',
                                          'ssst_run1_error_vs_quality_report.png')
    base, _ = os.path.splitext(parameters.figSave)
    outputs = []
    for path in (f'{base}.png', f'{base}.pdf'):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path)
        outputs.append(path)
        print(f'Figure saved @ {path}')
    plt.close(fig)
    return {'output': outputs[0], 'outputs': outputs}


def main():
    parser = argparse.ArgumentParser(
        description='Report-ready ERH/ERZ vs RMS, gap, phase count and nearest-station distance.')
    parser.add_argument('--csv', default='RESULT/SSST_result.csv',
                        help='result CSV annotated by NLL_run/pdf_metrics.py (default: %(default)s)')
    parser.add_argument('--output', default=None,
                        help='default: complem_figures/pdf_metrics/ssst_run1_error_vs_quality_report.png')
    parser.add_argument('--predictors', default=','.join(_DEFAULT_PREDICTORS),
                        help='comma-separated columns, one panel each '
                             f'(default: %(default)s; also Psi, C68)')
    parser.add_argument('--width-cm',  type=float, default=16.0,
                        help='printed width, i.e. your \\linewidth in cm (default: 16)')
    parser.add_argument('--height-cm', type=float, default=8.0,
                        help='printed height in cm (default: 8)')
    args = parser.parse_args()
    generate_figure(ErrorVsQualityReportParams(
        csv=args.csv, figSave=args.output, width_cm=args.width_cm, height_cm=args.height_cm,
        predictors=tuple(p.strip() for p in args.predictors.split(',') if p.strip())))


if __name__ == '__main__':
    main()
