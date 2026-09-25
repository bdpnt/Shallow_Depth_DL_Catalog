"""
temporary_network_impact_report.py
==================================
Report-ready version of the `phase_count_timeline` figure of
`temporary_network_impact.py`: two chronological strips, 1978 -> 2025,
coloured by the monthly median number of phases per event, before (a) and
after (b) the temporary-network picks, on a shared colour scale.

Differences from `temporary_network_impact.py`
----------------------------------------------
* Drawn at its printed size (`--width-cm`, default 16 cm = `\\linewidth` of an
  A4 page with 2.5 cm margins) with 8-9 pt fonts, so LaTeX includes it at
  scale 1 — the original is 13 x 4.2 in and loses half its text size.
* Light style (no seaborn grey panel), no suptitle, no row labels: a bold
  "a" / "b" in front of each strip instead.
* Timeline only; the parsing (`_load_station_coords`, `_scan_bulletin`,
  `_monthly_medians`) and the colour-scale rule (shared, clipped at the pooled
  p99, `extend='max'`) are imported / reproduced unchanged.
* Writes both a 300 dpi PNG and a vector PDF (`pdf.fonttype = 42`).

Usage
-----
    python complem_figures/temporary_network_impact_report.py

    python complem_figures/temporary_network_impact_report.py \\
        --output complem_figures/temporary_network_impact/phase_count_timeline_report.png \\
        --width-cm 16 --bin-months 1
"""

import argparse
import os
import sys
from dataclasses import dataclass
from typing import Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
from complem_figures.temporary_network_impact import (  # noqa: E402
    _COLOR_EMPTY, _load_station_coords, _monthly_medians, _scan_bulletin,
)

_INK       = '#1A1A1A'
_INK_MUTED = '#5A5A5A'
_CM        = 1 / 2.54


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class TimelineReportParams:
    file_before:   Optional[str]   = None
    file_after:    Optional[str]   = None
    stations_glob: Optional[str]   = None
    figSave:       Optional[str]   = None
    bin_months:    int             = 1
    vmin:          Optional[float] = None
    vmax:          Optional[float] = None
    width_cm:      float           = 16.0
    height_cm:     float           = 5.0
    tags:          tuple           = ('a', 'b')


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _style():
    mpl.rcParams.update({
        'font.family':       'sans-serif',
        'font.size':         8,
        'axes.labelsize':    9,
        'xtick.labelsize':   8,
        'ytick.labelsize':   8,
        'axes.edgecolor':    _INK_MUTED,
        'axes.linewidth':    0.6,
        'xtick.color':       _INK,
        'ytick.color':       _INK,
        'xtick.major.width': 0.6,
        'xtick.major.size':  3,
        'ytick.major.width': 0.6,
        'ytick.major.size':  3,
        'text.color':        _INK,
        'axes.labelcolor':   _INK,
        'pdf.fonttype':      42,
        'ps.fonttype':       42,
        'savefig.dpi':       300,
    })


def plot_timeline(before, after, parameters):
    """Two strips on a shared colour scale, sized for a \\linewidth figure."""
    _style()

    edges, med_before = _monthly_medians(before, parameters.bin_months)
    _,     med_after  = _monthly_medians(after,  parameters.bin_months)

    # Same colour-scale rule as the original: shared, clipped at the pooled p99.
    finite = np.concatenate([med_before[np.isfinite(med_before)],
                             med_after[np.isfinite(med_after)]])
    vmin = parameters.vmin if parameters.vmin is not None else float(finite.min())
    vmax = parameters.vmax if parameters.vmax is not None else float(np.percentile(finite, 99))

    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad(_COLOR_EMPTY)

    fig, axes = plt.subplots(2, 1, sharex=True, layout='constrained',
                             figsize=(parameters.width_cm * _CM, parameters.height_cm * _CM))
    fig.get_layout_engine().set(h_pad=0.02, hspace=0.04, w_pad=0.02)

    for axis, medians, tag in zip(axes, (med_before, med_after), parameters.tags):
        mesh = axis.pcolormesh(edges, [0, 1], np.ma.masked_invalid(medians)[None, :],
                               cmap=cmap, shading='flat', vmin=vmin, vmax=vmax,
                               rasterized=True)
        axis.set_yticks([])
        for side in ('top', 'right', 'left', 'bottom'):
            axis.spines[side].set_visible(False)
        if tag:
            # in the left margin, level with the top of the strip — nothing over the data
            axis.text(-0.012, 1.0, tag, transform=axis.transAxes, ha='right', va='top',
                      fontsize=10, fontweight='bold')

    axes[0].tick_params(axis='x', length=0)
    axes[1].set_xlabel('Year')
    axes[1].set_xlim(edges[0], edges[-1])
    axes[1].set_xticks(np.arange(np.ceil(edges[0] / 5) * 5, edges[-1], 5))
    axes[1].xaxis.set_major_formatter(lambda value, _: f'{value:.0f}')

    bin_label = 'monthly' if parameters.bin_months == 1 else f'{parameters.bin_months}-month'
    cbar = fig.colorbar(mesh, ax=axes, pad=0.01, fraction=0.03, aspect=18,
                        extend='max' if finite.max() > vmax else 'neither')
    cbar.set_label(f'Median phases\nper event ({bin_label})')
    cbar.outline.set_linewidth(0.6)
    cbar.ax.tick_params(width=0.6, length=3)
    cbar.ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(5))

    base, _ = os.path.splitext(parameters.figSave)
    outputs = []
    for path in (f'{base}.png', f'{base}.pdf'):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path)
        outputs.append(path)
        print(f'Figure saved @ {path}')
    plt.close(fig)
    return outputs, vmin, vmax


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    file_before   = parameters.file_before   or os.path.join(_PROJECT_ROOT, 'obs', 'NLL_result.obs')
    file_after    = parameters.file_after    or os.path.join(_PROJECT_ROOT, 'obs', 'NLL_result_augmented.obs')
    stations_glob = parameters.stations_glob or os.path.join(_PROJECT_ROOT, 'stations', 'GTSRCE_*.txt')
    if not parameters.figSave:
        parameters.figSave = os.path.join(_MODULE_DIR, 'temporary_network_impact',
                                          'phase_count_timeline_report.png')

    station_coords = _load_station_coords(stations_glob)
    before, s_before = _scan_bulletin(file_before, station_coords)
    after,  s_after  = _scan_bulletin(file_after,  station_coords)
    print(f'Before : {s_before["n_events"]} events / {s_before["n_picks"]} picks')
    print(f'After  : {s_after["n_events"]} events / {s_after["n_picks"]} picks')

    outputs, vmin, vmax = plot_timeline(before, after, parameters)
    print(f'Colour scale : {vmin:g} - {vmax:g} phases')
    return {'output': outputs[0], 'outputs': outputs, 'vmin': vmin, 'vmax': vmax}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Report-ready phase-count timeline (before / after temporary picks).'
    )
    parser.add_argument('--before',   default=None,
                        help='Bulletin before augmentation (default: obs/NLL_result.obs)')
    parser.add_argument('--after',    default=None,
                        help='Bulletin after augmentation (default: obs/NLL_result_augmented.obs)')
    parser.add_argument('--stations', default=None,
                        help='Glob of the GTSRCE station files (default: stations/GTSRCE_*.txt)')
    parser.add_argument('--output',   default=None,
                        help='Output path; both .png and .pdf are written '
                             '(default: complem_figures/temporary_network_impact/phase_count_timeline_report.png)')
    parser.add_argument('--bin-months', type=int, default=1,
                        help='Time bin width in months (default: 1)')
    parser.add_argument('--vmin', type=float, default=None, help='Colour scale minimum (default: data)')
    parser.add_argument('--vmax', type=float, default=None, help='Colour scale maximum (default: pooled p99)')
    parser.add_argument('--width-cm',  type=float, default=16.0,
                        help='Printed width, i.e. your \\linewidth in cm (default: 16)')
    parser.add_argument('--height-cm', type=float, default=5.0,
                        help='Printed height in cm (default: 5)')
    parser.add_argument('--no-tags', action='store_true', help='Omit the "a" / "b" tags')
    args = parser.parse_args()

    generate_figure(TimelineReportParams(
        file_before=args.before, file_after=args.after, stations_glob=args.stations,
        figSave=args.output, bin_months=args.bin_months, vmin=args.vmin, vmax=args.vmax,
        width_cm=args.width_cm, height_cm=args.height_cm,
        tags=(None, None) if args.no_tags else ('a', 'b'),
    ))


if __name__ == '__main__':
    main()
