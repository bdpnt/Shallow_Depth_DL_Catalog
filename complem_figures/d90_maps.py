r"""
d90_maps.py
============================
D90 maps — the depth above which 90 % of the seismicity lies — of a result CSV,
for the two hypocentre solutions and with / without the `pyr:usable` filter.

D90 is the usual proxy for the base of the seismogenic layer. It is a tail
statistic, so it is the first thing a depth error or a mode-vs-mean shift
moves; the four panels put those two effects side by side:

    rows    : PDF expectation (published) | maximum-likelihood point
    columns : all events                  | pyr:usable events only

`pyr:usable` is `NLL_run.export_quakeml.classify` itself, so it needs the
location-PDF columns of `RESULT/SSST_result.csv`. Each solution supplies its
own epicentre and depth (`event_maps.hypocentre_columns`).

Grid: nodes every `--step` degrees over 42-44 N, 2.25 W-3.5 E (the
`error_maps.py` extent); a node's D90 is the 90th percentile of the depths of
the events within `--radius` km of it (epicentral distance, local
equirectangular km), masked where fewer than `--min-n` events fall inside.
Depths are km below sea level, as in the CSV. Same page conventions as
`error_maps.py`: printed size `--width-cm`, shaded relief under opaque cells
(`--no-dem`), borders / shorelines from the `basemap_lines.py` caches,
300 dpi PNG and vector PDF.

Usage
-----
    python complem_figures/d90_maps.py

    python complem_figures/d90_maps.py --file RESULT/SSST_result.csv \
        --radius 10 --min-n 50 --step 0.05 --output complem_figures/d90_maps/d90_maps.png
"""

import argparse
import os
import sys
from dataclasses import dataclass

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402  (mako_r colormap)
from scipy.spatial import cKDTree  # noqa: E402

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
from complem_figures.basemap_lines import draw_basemap_lines, draw_relief  # noqa: E402
from complem_figures.event_maps import hypocentre_columns                  # noqa: E402
from NLL_run.export_quakeml import classify                                # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_LAT_MIN, _LAT_MAX = 42.0, 44.0
_LON_MIN, _LON_MAX = -2.25, 3.5
_KM_PER_DEG = 111.19
_PERCENTILE = 90
_VMIN, _VMAX = 0.0, 20.0

_SOLUTIONS = (('expect',  'PDF expectation'),
              ('maxlike', 'Maximum likelihood'))
_SUBSETS   = (('all',    'All events'),
              ('usable', 'pyr:usable only'))

_INK       = '#1A1A1A'
_INK_MUTED = '#5A5A5A'
_GRID      = '#D9D9D9'
_CM        = 1 / 2.54


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class D90MapsParams:
    file:       str   = os.path.join(_PROJECT_ROOT, 'RESULT', 'SSST_result.csv')
    figSave:    str   = os.path.join(_MODULE_DIR, 'd90_maps', 'd90_maps.png')
    radius_km:  float = 10.0
    min_n:      int   = 50
    step_deg:   float = 0.05
    width_cm:   float = 16.0
    borders:    bool  = True
    coastlines: bool  = True
    basemap_resolution: str = 'i'
    dem:        bool  = True


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def read_events(file):
    """
    RESULT/*.csv -> one (lat, lon, depth) table per solution, plus the
    `pyr:usable` mask shared by both.
    """
    df = pd.read_csv(file, skipinitialspace=True)
    if 'C68' not in df.columns:
        raise ValueError(f'{file} has no location-PDF columns, so pyr:usable cannot be '
                         f'evaluated; use RESULT/SSST_result.csv after NLL_run/pdf_metrics.py.')
    usable = np.array([classify(row)[0] for _, row in df.iterrows()])
    events = {}
    for solution, _ in _SOLUTIONS:
        lat_col, lon_col, z_col = hypocentre_columns(df, solution)
        events[solution] = pd.DataFrame({'latitude':  df[lat_col].to_numpy(),
                                         'longitude': df[lon_col].to_numpy(),
                                         'depth':     df[z_col].to_numpy()})
    return events, usable


def grid_nodes(step_deg):
    """Node longitudes and latitudes (1-D) of the map grid."""
    lon = np.arange(_LON_MIN, _LON_MAX + 1e-9, step_deg)
    lat = np.arange(_LAT_MIN, _LAT_MAX + 1e-9, step_deg)
    return lon, lat


def d90_grid(events, lon_nodes, lat_nodes, radius_km, min_n):
    """
    D90 and event count at every node: 90th percentile of the depths of the
    events within `radius_km` of the node; NaN where fewer than `min_n`.
    """
    lat0 = np.deg2rad(0.5 * (_LAT_MIN + _LAT_MAX))
    ok   = np.isfinite(events[['latitude', 'longitude', 'depth']].to_numpy()).all(axis=1)
    ev   = events.loc[ok]
    xy   = np.column_stack((ev['longitude'].to_numpy() * _KM_PER_DEG * np.cos(lat0),
                            ev['latitude'].to_numpy() * _KM_PER_DEG))
    depth = ev['depth'].to_numpy()

    glon, glat = np.meshgrid(lon_nodes, lat_nodes)
    nodes = np.column_stack((glon.ravel() * _KM_PER_DEG * np.cos(lat0),
                             glat.ravel() * _KM_PER_DEG))
    neighbours = cKDTree(xy).query_ball_point(nodes, r=radius_km)

    d90   = np.full(len(nodes), np.nan)
    count = np.array([len(idx) for idx in neighbours])
    for k in np.nonzero(count >= min_n)[0]:
        d90[k] = np.percentile(depth[neighbours[k]], _PERCENTILE)
    return d90.reshape(glon.shape), count.reshape(glon.shape)


def _cell_edges(nodes, step_deg):
    """pcolormesh edges centred on the nodes."""
    return np.append(nodes - step_deg / 2, nodes[-1] + step_deg / 2)


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
    """
    Compute the four D90 grids and draw them on one 2 x 2 page; writes
    <figSave>.png and .pdf.

    Returns
    -------
    dict with keys: output, outputs, lon, lat, d90 and count (both dicts keyed
    by (solution, subset)), n_events, n_usable
    """
    events, usable = read_events(parameters.file)
    print(f'{parameters.file}: {len(usable)} events, {int(usable.sum())} pyr:usable')

    lon_nodes, lat_nodes = grid_nodes(parameters.step_deg)
    d90, count = {}, {}
    for solution, _ in _SOLUTIONS:
        for subset, _ in _SUBSETS:
            sel = events[solution] if subset == 'all' else events[solution].loc[usable]
            key = (solution, subset)
            d90[key], count[key] = d90_grid(sel, lon_nodes, lat_nodes,
                                            parameters.radius_km, parameters.min_n)
            valid = d90[key][np.isfinite(d90[key])]
            print(f'  {solution:7s} / {subset:6s}: {len(sel)} events, {valid.size} nodes, '
                  f'median D90 {np.median(valid):.2f} km, range {valid.min():.2f}-{valid.max():.2f} km'
                  if valid.size else
                  f'  {solution:7s} / {subset:6s}: no node reaches {parameters.min_n} events')
    reference = d90[('expect', 'all')]
    for key, grid in d90.items():
        if key == ('expect', 'all'):
            continue
        diff = (grid - reference)[np.isfinite(grid) & np.isfinite(reference)]
        if diff.size:
            print(f'  {key[0]}/{key[1]} - expect/all: median {np.median(diff):+.2f} km, '
                  f'5-95 % {np.percentile(diff, 5):+.2f} / {np.percentile(diff, 95):+.2f} km '
                  f'over {diff.size} shared nodes')

    _style()
    lon_edges = _cell_edges(lon_nodes, parameters.step_deg)
    lat_edges = _cell_edges(lat_nodes, parameters.step_deg)
    aspect    = 1.0 / np.cos(np.deg2rad(0.5 * (_LAT_MIN + _LAT_MAX)))

    width = parameters.width_cm
    left, right, hgap = 1.3, 0.15, 0.3                  # cm
    map_w = (width - left - right - hgap) / 2
    map_h = map_w * (_LAT_MAX - _LAT_MIN) * aspect / (_LON_MAX - _LON_MIN)
    cb_block, bottom, vgap, top = 1.2, 0.8, 0.3, 0.45   # cm
    height = cb_block + bottom + 2 * map_h + vgap + top

    fig  = plt.figure(figsize=(width * _CM, height * _CM))
    cmap = sns.color_palette('mako_r', as_cmap=True)
    mesh = None
    first = True
    for r, (solution, solution_label) in enumerate(_SOLUTIONS):
        for c, (subset, subset_label) in enumerate(_SUBSETS):
            x  = left + c * (map_w + hgap)
            y  = cb_block + bottom + (1 - r) * (map_h + vgap)
            ax = fig.add_axes([x / width, y / height, map_w / width, map_h / height])
            ax.set_axisbelow(True)
            ax.grid(color=_GRID, linewidth=0.4)
            # opaque cells over the relief: mako_r and the hillshade both encode lightness
            relief = parameters.dem and draw_relief(ax, verbose=first)
            mesh = ax.pcolormesh(lon_edges, lat_edges,
                                 np.ma.masked_invalid(d90[(solution, subset)]),
                                 vmin=_VMIN, vmax=_VMAX, cmap=cmap, shading='flat',
                                 alpha=1.0 if relief else 0.9, rasterized=True)
            draw_basemap_lines(ax, borders=parameters.borders,
                               coastlines=parameters.coastlines,
                               resolution=parameters.basemap_resolution,
                               verbose=first)
            first = False
            ax.text(0.015, 0.96, f'{"ab"[r]}{c + 1}  {solution_label}', transform=ax.transAxes,
                    ha='left', va='top', fontsize=8, fontweight='bold')
            if r == 0:
                ax.set_title(subset_label, fontsize=8, pad=3)
            ax.set_xlim(_LON_MIN, _LON_MAX)
            ax.set_ylim(_LAT_MIN, _LAT_MAX)
            ax.set_aspect(aspect)
            ax.set_xticks(np.arange(-2, 4, 1))
            ax.set_yticks(np.arange(42.0, 44.01, 0.5))
            ax.xaxis.set_major_formatter(lambda v, _: f'{v:.0f}°'.replace('-', '−'))
            ax.yaxis.set_major_formatter(lambda v, _: f'{v:.1f}°')
            if c == 0:
                ax.set_ylabel('Latitude')
            else:
                ax.tick_params(labelleft=False)
            if r == 0:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel('Longitude')
            for side in ('top', 'right'):
                ax.spines[side].set_visible(False)

    cb_w = 0.5 * (2 * map_w + hgap)
    cax  = fig.add_axes([(left + (2 * map_w + hgap - cb_w) / 2) / width,
                         (cb_block - 0.4) / height, cb_w / width, 0.22 / height])
    cb = fig.colorbar(mesh, cax=cax, orientation='horizontal', extend='max')
    cb.set_label(f'D90 (km b.s.l.), events within {parameters.radius_km:g} km, '
                 f'n ≥ {parameters.min_n}', labelpad=3)
    cb.outline.set_linewidth(0.5)
    cb.ax.tick_params(width=0.5, length=2.5)
    cb.solids.set(alpha=1)

    base, _ = os.path.splitext(parameters.figSave)
    outputs = []
    for path in (f'{base}.png', f'{base}.pdf'):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path)
        outputs.append(path)
        print(f'Figure saved @ {path}')
    plt.close(fig)
    return {'output': outputs[0], 'outputs': outputs, 'lon': lon_nodes, 'lat': lat_nodes,
            'd90': d90, 'count': count, 'n_events': len(usable),
            'n_usable': int(usable.sum())}


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    defaults = D90MapsParams()
    parser = argparse.ArgumentParser(description='D90 maps: expectation / maximum likelihood, '
                                                 'all / pyr:usable events.')
    parser.add_argument('--file', default=defaults.file,
                        help='result CSV with the location-PDF columns (default: RESULT/SSST_result.csv)')
    parser.add_argument('--output', default=defaults.figSave,
                        help='output figure, .png and .pdf written (default: complem_figures/d90_maps/d90_maps.png)')
    parser.add_argument('--radius', type=float, default=defaults.radius_km,
                        help='search radius around each node in km (default: 10)')
    parser.add_argument('--min-n', type=int, default=defaults.min_n,
                        help='minimum events within the radius for a D90 (default: 50)')
    parser.add_argument('--step', type=float, default=defaults.step_deg,
                        help='node spacing in degrees (default: 0.05)')
    parser.add_argument('--width-cm', type=float, default=defaults.width_cm,
                        help='printed width, i.e. your \\linewidth in cm (default: 16)')
    parser.add_argument('--no-borders', action='store_true',
                        help='do not draw the national borders')
    parser.add_argument('--no-coastlines', action='store_true',
                        help='do not draw the shorelines')
    parser.add_argument('--basemap-resolution', default='i', choices=['c', 'l', 'i', 'h', 'f'],
                        help='which basemap_lines.py cache to read (default: i)')
    parser.add_argument('--no-dem', action='store_true',
                        help='do not draw the shaded relief under the D90 cells')
    args = parser.parse_args()
    generate_figure(D90MapsParams(file=args.file, figSave=args.output,
                                  radius_km=args.radius, min_n=args.min_n,
                                  step_deg=args.step, width_cm=args.width_cm,
                                  borders=not args.no_borders,
                                  coastlines=not args.no_coastlines,
                                  basemap_resolution=args.basemap_resolution,
                                  dem=not args.no_dem))


if __name__ == '__main__':
    main()
