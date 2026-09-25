"""
error_maps_report.py
====================
Report-ready version of one page of `error_maps.py`: windowed-median maps of
the horizontal (ERH) and vertical (ERV) location error over a chosen period,
with the events overlaid.

Differences from `error_maps.py`
--------------------------------
* Any period (`--start` / `--end`, whole calendar years, both inclusive)
  instead of the fixed 5-year windows starting in 1976, so 2020-2025 can be
  drawn as one page.
* Drawn at its printed size (`--width-cm`, default 16 cm = `\\linewidth` of an
  A4 page with 2.5 cm margins) with 7-9 pt fonts, so LaTeX includes it at
  scale 1 — the original page is 12 x 10 in.
* Maps at the true geographic aspect (1 deg lon = cos(lat) deg lat), light
  style, no suptitle, no statistics block over the data; "ERH" / "ERV" tag in
  each map's top-left corner; one shared colour bar to the right.
* Same grid (400 x 860 cells over 42-44 N, 2.25 W-3.5 E), same 9 x 9-cell
  window, same `count >= 10` mask, same 0-5 km `rocket_r` scale as the
  original. The median is computed per cell from the events binned into the
  9 x 9 neighbourhood (then the original's exact inclusive edge test) rather
  than by 344 000 DataFrame masks — identical numbers, seconds instead of hours.
* National borders and shorelines drawn *under* the data (they show through
  where the map is empty, the median field and the epicentres cover them
  elsewhere), from the same GMT database as the PyGMT maps (`aoi_map_report.py`'s `fig.coast(..., borders='1/...')`), through
  the `.xy` caches of `basemap_lines.py`; `--no-borders` / `--no-coastlines`
  turn them off.  The caches must be dumped once, in an environment with GMT:
  `conda run -n pygmt_env python complem_figures/basemap_lines.py --dump`.
  Without them the figure is drawn as before, with a warning.
* Writes both a 300 dpi PNG and a vector PDF (`pdf.fonttype = 42`).

Usage
-----
    python complem_figures/error_maps_report.py \\
        --file RESULT/SSST_result.csv --start 2020 --end 2025 \\
        --output complem_figures/error_maps_ssst/2020-2025_report.png
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

# Grid of error_maps._add_subplot
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
class ErrorMapsReportParams:
    file:      str
    start:     int
    end:       int
    figSave:   Optional[str] = None
    width_cm:  float = 16.0
    borders:    bool = True
    coastlines: bool = True
    basemap_resolution: str = 'i'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_events(file, start, end):
    """RESULT/*.csv -> events of [start-01-01, end-12-31], as error_maps._read_file."""
    df = pd.read_csv(file)
    df['time'] = pd.to_datetime(df['date-time'])
    df = df.rename(columns={'true_erh': 'erh', 'true_erz': 'erv'})
    mask = ((df['time'] >= pd.Timestamp(f'{start}-01-01'))
            & (df['time'] <= pd.Timestamp(f'{end}-12-31 23:59:59.999999')))
    return df.loc[mask, ['time', 'latitude', 'longitude', 'erh', 'erv']].reset_index(drop=True)


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

def generate_figure(parameters):
    events = read_events(parameters.file, parameters.start, parameters.end)
    print(f'{parameters.file}: {len(events)} events in {parameters.start}-{parameters.end}')
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

    if not parameters.figSave:
        parameters.figSave = os.path.join(_MODULE_DIR, 'error_maps_ssst',
                                          f'{parameters.start}-{parameters.end}_report.png')
    base, _ = os.path.splitext(parameters.figSave)
    outputs = []
    for path in (f'{base}.png', f'{base}.pdf'):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path)
        outputs.append(path)
        print(f'Figure saved @ {path}')
    plt.close(fig)
    return {'output': outputs[0], 'outputs': outputs, 'n_events': len(events)}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Report-ready ERH/ERV error maps for one period.')
    parser.add_argument('--file', default='RESULT/SSST_result.csv',
                        help='result CSV (default: RESULT/SSST_result.csv)')
    parser.add_argument('--start', type=int, default=2020, help='first year, inclusive (default: 2020)')
    parser.add_argument('--end',   type=int, default=2025, help='last year, inclusive (default: 2025)')
    parser.add_argument('--output', default=None,
                        help='default: complem_figures/error_maps_ssst/<start>-<end>_report.png')
    parser.add_argument('--width-cm', type=float, default=16.0,
                        help='printed width, i.e. your \\linewidth in cm (default: 16)')
    parser.add_argument('--no-borders', action='store_true',
                        help='do not draw the national borders')
    parser.add_argument('--no-coastlines', action='store_true',
                        help='do not draw the shorelines')
    parser.add_argument('--basemap-resolution', default='i', choices=['c', 'l', 'i', 'h', 'f'],
                        help='which basemap_lines.py cache to read (default: i)')
    args = parser.parse_args()
    generate_figure(ErrorMapsReportParams(file=args.file, start=args.start, end=args.end,
                                          figSave=args.output, width_cm=args.width_cm,
                                          borders=not args.no_borders,
                                          coastlines=not args.no_coastlines,
                                          basemap_resolution=args.basemap_resolution))


if __name__ == '__main__':
    main()
