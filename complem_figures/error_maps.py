"""
error_maps.py
============================
Windowed-median maps of the horizontal (ERH) and vertical (ERV) location error,
with the events overlaid — one page per time period.

Two modes:
* batch (default): consecutive `--time-range`-year periods from 1976 to 2025,
  one page each in `--map-folder` (this is what `generate_complem_figures.py`
  runs);
* single period: `--start` / `--end` (whole calendar years, both inclusive),
  e.g. 2020-2025 as one page.

Grid: 400 x 860 cells over 42-44 N, 2.25 W-3.5 E; a cell's value is the median
over its 9 x 9-cell window (the original's exact inclusive edge test), masked
where the window holds fewer than 10 events; 0-5 km `rocket_r` scale. The
median is computed from the events binned once rather than by 344 000
DataFrame masks — identical numbers, seconds instead of hours. Events are
placed at the PDF expectation (`expect_lat` / `expect_lon`), the point where
the published errors are evaluated.

Pages are drawn at their printed size (`--width-cm`, default 16 cm =
`\\linewidth` of an A4 page with 2.5 cm margins) with 7-9 pt fonts, at the
true geographic aspect, one shared colour bar to the right. National borders
and shorelines are drawn *under* the data, from the same GMT database as the
PyGMT maps, through the `.xy` caches of `basemap_lines.py`; `--no-borders` /
`--no-coastlines` turn them off. The caches must be dumped once, in an
environment with GMT: `conda run -n pygmt_env python complem_figures/basemap_lines.py --dump`.
Without them the figure is drawn as before, with a warning. Writes both a
300 dpi PNG and a vector PDF (`pdf.fonttype = 42`).

Usage
-----
    python complem_figures/error_maps.py \\
        --file RESULT/NLL_result.csv --map-folder complem_figures/error_maps/ --time-range 5

    python complem_figures/error_maps.py \\
        --file RESULT/SSST_result.csv --start 2020 --end 2025 \\
        --output complem_figures/error_maps_ssst/2020-2025.png
"""

import argparse
import os
from dataclasses import dataclass
from typing import Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402  (rocket_r colormap)

try:                                                    # noqa: E402
    from complem_figures.basemap_lines import draw_basemap_lines
except ImportError:                                     # run as a script
    from basemap_lines import draw_basemap_lines

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

# Map grid
_LAT_MIN, _LAT_MAX = 42.0, 44.0
_LON_MIN, _LON_MAX = -2.25, 3.5
_BINS_LAT, _BINS_LON = 400, 860
_WINDOW   = 4      # cells on each side -> 9 x 9
_MIN_N    = 10
_VMIN, _VMAX = 0.0, 5.0

_INK       = '#1A1A1A'
_INK_MUTED = '#5A5A5A'
_GRID      = '#D9D9D9'
_CM        = 1 / 2.54


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ErrorMapsParams:
    file:       str
    mapFolder:  Optional[str] = None   # batch mode: output folder
    time_range: int           = 5      # batch mode: period length in years
    start:      Optional[int] = None   # single-period mode (with end)
    end:        Optional[int] = None
    figSave:    Optional[str] = None   # single-period mode
    width_cm:  float = 16.0
    borders:    bool = True
    coastlines: bool = True
    basemap_resolution: str = 'i'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_events(file):
    """RESULT/*.csv -> time, latitude, longitude (PDF expectation), erh, erv."""
    df = pd.read_csv(file)
    return pd.DataFrame({'time':      pd.to_datetime(df['date-time']),
                         'latitude':  df['expect_lat'],
                         'longitude': df['expect_lon'],
                         'erh':       df['true_erh'],
                         'erv':       df['true_erz']})


def _select_period(events, start, end):
    """Events of [start-01-01, end-12-31]."""
    mask = ((events['time'] >= pd.Timestamp(f'{start}-01-01'))
            & (events['time'] <= pd.Timestamp(f'{end}-12-31 23:59:59.999999')))
    return events.loc[mask].reset_index(drop=True)


def windowed_median(events, column):
    """
    Median of `column` over the 9 x 9-cell window of every cell; NaN where the
    window holds fewer than 10 events.

    Events are binned once; a cell's window is the union of the bins within
    +/-4 cells, clipped to the grid — the index form of the original's
    edge-based window.
    """
    lat_step = (_LAT_MAX - _LAT_MIN) / _BINS_LAT
    lon_step = (_LON_MAX - _LON_MIN) / _BINS_LON
    lat, lon = events['latitude'].to_numpy(), events['longitude'].to_numpy()
    val      = events[column].to_numpy()
    inside   = ((lat >= _LAT_MIN) & (lat <= _LAT_MAX) & (lon >= _LON_MIN) & (lon <= _LON_MAX)
                & np.isfinite(val))
    lat, lon, val = lat[inside], lon[inside], val[inside]
    bi = np.clip(((lat - _LAT_MIN) / lat_step).astype(int), 0, _BINS_LAT - 1)
    bj = np.clip(((lon - _LON_MIN) / lon_step).astype(int), 0, _BINS_LON - 1)

    counts = np.zeros((_BINS_LAT, _BINS_LON), dtype=int)
    np.add.at(counts, (bi, bj), 1)
    # window counts by a 2-D prefix sum, to skip cells with no event nearby cheaply
    csum = np.pad(counts.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
    i0 = np.clip(np.arange(_BINS_LAT) - _WINDOW, 0, _BINS_LAT)
    i1 = np.clip(np.arange(_BINS_LAT) + _WINDOW + 1, 0, _BINS_LAT)
    j0 = np.clip(np.arange(_BINS_LON) - _WINDOW, 0, _BINS_LON)
    j1 = np.clip(np.arange(_BINS_LON) + _WINDOW + 1, 0, _BINS_LON)
    win = (csum[i1][:, j1] - csum[i0][:, j1] - csum[i1][:, j0] + csum[i0][:, j0])

    order = np.lexsort((bj, bi))
    bi_s, bj_s = bi[order], bj[order]
    lat_s, lon_s, val_s = lat[order], lon[order], val[order]
    key   = bi_s * _BINS_LON + bj_s
    cells = np.arange(_BINS_LAT * _BINS_LON)
    first = np.searchsorted(key, cells, side='left')
    last  = np.searchsorted(key, cells, side='right')

    lat_edges = np.linspace(_LAT_MIN, _LAT_MAX, _BINS_LAT + 1)
    lon_edges = np.linspace(_LON_MIN, _LON_MAX, _BINS_LON + 1)
    dlat, dlon = lat_edges[1] - lat_edges[0], lon_edges[1] - lon_edges[0]

    median = np.full((_BINS_LAT, _BINS_LON), np.nan)
    # candidates: one extra ring of bins around the 9 x 9 block, then the
    # original's exact inclusive float bounds — catalogue coordinates are
    # rounded to 1e-4 deg, so many events sit exactly on a 0.005 deg cell edge
    # and the edge convention changes the median in ~4 % of cells.
    reach = (win >= 1)
    for i, j in zip(*np.nonzero(reach)):
        r0, r1 = max(i - _WINDOW - 1, 0), min(i + _WINDOW + 2, _BINS_LAT)
        c0, c1 = max(j - _WINDOW - 1, 0), min(j + _WINDOW + 2, _BINS_LON)
        idx = [np.arange(first[r * _BINS_LON + c], last[r * _BINS_LON + c])
               for r in range(r0, r1) for c in range(c0, c1) if counts[r, c]]
        if not idx:
            continue
        idx = np.concatenate(idx)
        lat_low  = max(lat_edges[i]     - _WINDOW * dlat, _LAT_MIN)
        lat_high = min(lat_edges[i + 1] + _WINDOW * dlat, _LAT_MAX)
        lon_low  = max(lon_edges[j]     - _WINDOW * dlon, _LON_MIN)
        lon_high = min(lon_edges[j + 1] + _WINDOW * dlon, _LON_MAX)
        keep = ((lat_s[idx] >= lat_low) & (lat_s[idx] <= lat_high)
                & (lon_s[idx] >= lon_low) & (lon_s[idx] <= lon_high))
        if keep.sum() >= _MIN_N:
            median[i, j] = np.median(val_s[idx[keep]])
    return median


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
        'xtick.major.size':  2.5,
        'ytick.major.size':  2.5,
        'text.color':        _INK,
        'axes.labelcolor':   _INK,
        'xtick.color':       _INK,
        'ytick.color':       _INK,
        'pdf.fonttype':      42,
        'ps.fonttype':       42,
        'savefig.dpi':       300,
    })


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def plot_period(events, fig_save, parameters):
    """Draw the ERH / ERV page of one period; writes <fig_save>.png and .pdf."""
    _style()

    lat_edges = np.linspace(_LAT_MIN, _LAT_MAX, _BINS_LAT + 1)
    lon_edges = np.linspace(_LON_MIN, _LON_MAX, _BINS_LON + 1)
    aspect    = 1.0 / np.cos(np.deg2rad(0.5 * (_LAT_MIN + _LAT_MAX)))

    width = parameters.width_cm
    left, right_block = 1.45, 1.6                       # cm: lat labels | colour bar
    map_w = width - left - right_block
    map_h = map_w * (_LAT_MAX - _LAT_MIN) * aspect / (_LON_MAX - _LON_MIN)
    bottom, gap, top = 0.9, 0.35, 0.15                 # cm
    height = bottom + 2 * map_h + gap + top

    fig  = plt.figure(figsize=(width * _CM, height * _CM))
    cmap = sns.color_palette('rocket_r', as_cmap=True)
    mesh = None
    for k, column in enumerate(('erh', 'erv')):
        y  = bottom + (1 - k) * (map_h + gap)
        ax = fig.add_axes([left / width, y / height, map_w / width, map_h / height])
        ax.set_axisbelow(True)
        ax.grid(color=_GRID, linewidth=0.4)
        median = windowed_median(events, column)
        mesh = ax.pcolormesh(lon_edges, lat_edges, np.ma.masked_invalid(median),
                             vmin=_VMIN, vmax=_VMAX, cmap=cmap, shading='auto',
                             alpha=0.9, rasterized=True)
        ax.scatter(events['longitude'], events['latitude'], s=0.08, c='black',
                   linewidths=0, rasterized=True)
        draw_basemap_lines(ax, borders=parameters.borders,
                           coastlines=parameters.coastlines,
                           resolution=parameters.basemap_resolution,
                           verbose=(k == 0))
        ax.text(0.01, 0.97, column.upper(), transform=ax.transAxes,
                ha='left', va='top', fontsize=9, fontweight='bold')
        ax.set_xlim(_LON_MIN, _LON_MAX)
        ax.set_ylim(_LAT_MIN, _LAT_MAX)
        ax.set_aspect(aspect)
        ax.set_xticks(np.arange(-2, 4, 1))
        ax.set_yticks(np.arange(42.0, 44.01, 0.5))
        ax.xaxis.set_major_formatter(lambda v, _: f'{v:.0f}°'.replace('-', '\u2212'))
        ax.yaxis.set_major_formatter(lambda v, _: f'{v:.1f}°')
        ax.set_ylabel('Latitude')
        if k == 0:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel('Longitude')
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        valid = median[np.isfinite(median)]
        print(f'  {column.upper()}: {valid.size} cells, mean of medians {np.mean(valid):.2f} km, '
              f'max {np.max(valid):.2f} km, median over events {np.nanmedian(events[column]):.2f} km')

    cb_h = 0.7 * (2 * map_h + gap)
    cax  = fig.add_axes([(left + map_w + 0.35) / width,
                         (bottom + (2 * map_h + gap - cb_h) / 2) / height,
                         0.25 / width, cb_h / height])
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_label('Median error (km), 9 × 9-cell window', labelpad=4)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2.5)
    cb.solids.set(alpha=1)

    base, _ = os.path.splitext(fig_save)
    outputs = []
    for path in (f'{base}.png', f'{base}.pdf'):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path)
        outputs.append(path)
        print(f'Figure saved @ {path}')
    plt.close(fig)
    return outputs


def generate_figure(parameters):
    """
    Draw the error maps: one page for [start, end] when both are given,
    otherwise one page per `time_range`-year period from 1976 in `mapFolder`.

    Returns
    -------
    dict with keys: output, outputs, n_figures
    """
    events = read_events(parameters.file)

    if parameters.start is not None and parameters.end is not None:
        periods = [(parameters.start, parameters.end,
                    parameters.figSave or os.path.join(
                        _MODULE_DIR, 'error_maps_ssst',
                        f'{parameters.start}-{parameters.end}.png'))]
    else:
        folder  = parameters.mapFolder or os.path.join(_MODULE_DIR, 'error_maps')
        periods = [(y, y + parameters.time_range - 1,
                    os.path.join(folder, f'{y}-{y + parameters.time_range - 1}.png'))
                   for y in range(1976, 2026, parameters.time_range)]

    outputs = []
    for start, end, fig_save in periods:
        period = _select_period(events, start, end)
        print(f'{parameters.file}: {len(period)} events in {start}-{end}')
        if period.empty:
            continue
        outputs += plot_period(period, fig_save, parameters)
    return {'output': outputs[0] if outputs else None, 'outputs': outputs,
            'n_figures': len(outputs) // 2}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='ERH/ERV windowed-median error maps.')
    parser.add_argument('--file', default='RESULT/SSST_result.csv',
                        help='result CSV (default: RESULT/SSST_result.csv)')
    parser.add_argument('--map-folder', default=None,
                        help='batch mode output folder (default: complem_figures/error_maps/)')
    parser.add_argument('--time-range', type=int, default=5,
                        help='batch mode period length in years (default: 5)')
    parser.add_argument('--start', type=int, default=None,
                        help='single period: first year, inclusive (with --end)')
    parser.add_argument('--end',   type=int, default=None,
                        help='single period: last year, inclusive (with --start)')
    parser.add_argument('--output', default=None,
                        help='single period output (default: complem_figures/error_maps_ssst/<start>-<end>.png)')
    parser.add_argument('--width-cm', type=float, default=16.0,
                        help='printed width, i.e. your \\linewidth in cm (default: 16)')
    parser.add_argument('--no-borders', action='store_true',
                        help='do not draw the national borders')
    parser.add_argument('--no-coastlines', action='store_true',
                        help='do not draw the shorelines')
    parser.add_argument('--basemap-resolution', default='i', choices=['c', 'l', 'i', 'h', 'f'],
                        help='which basemap_lines.py cache to read (default: i)')
    args = parser.parse_args()
    if (args.start is None) != (args.end is None):
        parser.error('--start and --end go together')
    generate_figure(ErrorMapsParams(file=args.file, mapFolder=args.map_folder,
                                    time_range=args.time_range,
                                    start=args.start, end=args.end,
                                    figSave=args.output, width_cm=args.width_cm,
                                    borders=not args.no_borders,
                                    coastlines=not args.no_coastlines,
                                    basemap_resolution=args.basemap_resolution))


if __name__ == '__main__':
    main()
