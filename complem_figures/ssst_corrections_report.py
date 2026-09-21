"""
ssst_corrections_report.py
==========================
Report-ready version of one page of the station atlas drawn by
`ssst_corrections.py`: the five SSST travel-time correction increments
(L = 9999 -> 50 -> 15 -> 5 -> 1 km) and their cumulative sum, for one
station and phase, as depth-slice maps.

Differences from the atlas page
-------------------------------
* Drawn at its printed size (`--width-cm`, default 16 cm = `\\linewidth` of an
  A4 page with 2.5 cm margins) with 7-9 pt fonts, so LaTeX includes it at
  scale 1 — the atlas page is 15 x 9 in and loses ~60 % of its text size.
* Three rows of two maps instead of two rows of three, so each map is twice
  as wide on the page.
* No suptitle and no footnote — that information belongs in the caption.
  Each panel carries a bold tag (a-f), its smoothing length, and on the right
  its arrival count and rms. Longitude labels on the bottom row only,
  latitude labels on the left column only.
* The reconstruction (cache, station_fields, LSGRID clipping, the 98th
  percentile colour scales, the grey support mask) is imported unchanged.
* Writes both a 300 dpi PNG and a vector PDF (`pdf.fonttype = 42`).

Usage
-----
    python complem_figures/ssst_corrections_report.py --station CA.0051:S

    # finer node spacing, which resolves the L = 1 km panel (~4x slower per 1/2 spacing)
    python complem_figures/ssst_corrections_report.py --station CA.0051:S --map-spacing 0.01
"""

import argparse
import os
import sys
from dataclasses import dataclass
from typing import Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
from complem_figures import ssst_corrections as sc  # noqa: E402

_INK       = '#1A1A1A'
_INK_MUTED = '#5A5A5A'
_CM        = 1 / 2.54

# Shared frame of the atlas page supplied for the report
DEFAULT_EXTENT = (-2.25, 3.5, 41.75, 43.75)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SsstCorrReportParams:
    station:     str                   # 'STA:PHASE'
    figSave:     Optional[str]  = None
    run_name:    str            = 'ssst_run1'
    zone:        Optional[int]  = None   # None = the station's best-sampled zone, as the atlas
    depth:       Optional[float] = None  # None = the station's median event depth, as the atlas
    map_spacing: float          = sc._MAP_SPACING
    extent:      Optional[tuple] = DEFAULT_EXTENT
    width_cm:    float          = 16.0


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _style():
    mpl.rcParams.update({
        'font.family':       'sans-serif',
        'font.size':         7,
        'axes.titlesize':    8,
        'axes.labelsize':    8,
        'xtick.labelsize':   7,
        'ytick.labelsize':   7,
        'axes.edgecolor':    _INK_MUTED,
        'axes.linewidth':    0.5,
        'xtick.major.width': 0.5,
        'ytick.major.width': 0.5,
        'xtick.major.size':  2.5,
        'ytick.major.size':  2.5,
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


def _draw(ax, lats, lons, values, mask, vmax, sta_lat, sta_lon):
    mesh = ax.pcolormesh(lons, lats, np.ma.masked_where(mask, values), cmap='RdBu_r',
                         vmin=-vmax, vmax=vmax, shading='auto', rasterized=True)
    ax.set_facecolor('0.88')
    ax.plot(sta_lon, sta_lat, marker='v', color='k', markersize=5,
            markeredgecolor='w', markeredgewidth=0.6, zorder=5)
    ax.set_aspect(1.0 / np.cos(np.deg2rad(float(np.mean(lats)))))
    ax.set_xlim(lons[0] - 0.5 * (lons[1] - lons[0]), lons[-1] + 0.5 * (lons[1] - lons[0]))
    ax.set_ylim(lats[0] - 0.5 * (lats[1] - lats[0]), lats[-1] + 0.5 * (lats[1] - lats[0]))
    return mesh


def plot_page(fields, sta_lat, sta_lon, parameters):
    """Five increments + cumulative, 3 rows x 2 columns, at printed size."""
    _style()

    increments, weights = fields['increments'], fields['weights']
    lats, lons = fields['lats'], fields['lons']

    # colour scales and support mask exactly as atlas_page
    masks = [w < sc._MIN_EFF_N for w in weights]
    stack = np.concatenate([np.asarray(v)[~m].ravel()
                            for v, m in zip(increments, masks) if (~m).any()])
    vmax_inc = float(np.percentile(np.abs(stack), 98)) if stack.size else 1e-3
    inside_total = fields['total'][~fields['outside']]
    vmax_tot = (float(np.percentile(np.abs(inside_total), 98))
                if inside_total.size else 1e-3) or 1e-3

    width = parameters.width_cm
    lon_span = (lons[-1] - lons[0]) * np.cos(np.deg2rad(float(np.mean(lats))))
    map_ratio = (lats[-1] - lats[0]) / lon_span
    # margins in cm: left for lat labels, right gap between columns, title strip per row
    left, gap, right = 1.0, 0.35, 0.1
    map_w = (width - left - gap - right) / 2
    map_h = map_w * map_ratio
    title_h, bottom_axis, cbar_block = 0.45, 0.55, 1.35
    row_h = title_h + map_h
    height = 3 * row_h + 2 * 0.12 + bottom_axis + cbar_block

    fig = plt.figure(figsize=(width * _CM, height * _CM))

    def rect(row, col):
        x = left + col * (map_w + gap)
        y = height - (row + 1) * row_h - row * 0.12
        return [x / width, y / height, map_w / width, map_h / height]

    tags = 'abcdef'
    axes = []
    for k in range(6):
        row, col = divmod(k, 2)
        ax = fig.add_axes(rect(row, col))
        axes.append(ax)
        if k < 5:
            mesh_inc = _draw(ax, lats, lons, increments[k], masks[k], vmax_inc, sta_lat, sta_lon)
            char_dist = sc._CHAR_DISTS[k]
            label = 'L = 9999 km (static)' if char_dist > 1000 else f'L = {char_dist:g} km'
            rms = (float(np.sqrt((np.asarray(increments[k])[~masks[k]] ** 2).mean()))
                   if (~masks[k]).any() else 0.0)
            n_arr = f'{fields["n_arrivals"][k]:,}'.replace(',', '\u2009')
            right_txt = f'{n_arr} arrivals, rms {rms * 1000:.0f} ms'
        else:
            mesh_tot = _draw(ax, lats, lons, fields['total'], fields['outside'],
                             vmax_tot, sta_lat, sta_lon)
            label, right_txt = 'Cumulative (final grids)', ''
        ax.set_title(f'$\\bf{{{tags[k]}}}$   {label}', loc='left', pad=3)
        if right_txt:
            ax.set_title(right_txt, loc='right', pad=3, color=_INK_MUTED, fontsize=7)

        ax.set_xticks(np.arange(np.ceil(lons[0]), lons[-1] + 1e-9, 1.0))
        ax.set_yticks(np.arange(np.ceil(lats[0] * 2) / 2, lats[-1] + 1e-9, 0.5))
        ax.xaxis.set_major_formatter(lambda v, _: f'{v:.0f}°')
        ax.yaxis.set_major_formatter(lambda v, _: f'{v:.1f}°')
        if row < 2:
            ax.tick_params(labelbottom=False)
        if col == 1:
            ax.tick_params(labelleft=False)

    # two colour bars under the bottom row: increments under column a, total under column b
    cb_y = (cbar_block - 0.45) / height
    for col, mesh, label in ((0, mesh_inc, 'Increment (s), shared by a–e'),
                             (1, mesh_tot, 'Total correction (s)')):
        x = left + col * (map_w + gap) + 0.1 * map_w
        cax = fig.add_axes([x / width, cb_y, 0.8 * map_w / width, 0.2 / height])
        cb = fig.colorbar(mesh, cax=cax, orientation='horizontal')
        cb.set_label(label, fontsize=8, labelpad=2)
        cb.outline.set_linewidth(0.5)
        cb.ax.tick_params(width=0.5, length=2.5, labelsize=7)
        cb.ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(5, symmetric=True))

    base, _ = os.path.splitext(parameters.figSave)
    outputs = []
    for path in (f'{base}.png', f'{base}.pdf'):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path)
        outputs.append(path)
        print(f'Figure saved @ {path}')
    plt.close(fig)
    return outputs


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    station, _, phase = parameters.station.partition(':')
    if phase not in ('P', 'S'):
        raise SystemExit('--station must be STA:PHASE, e.g. CA.0051:S')

    params = sc.SsstCorrParams(run_name=parameters.run_name, depth=parameters.depth,
                               map_spacing=parameters.map_spacing, extent=parameters.extent)
    zone = parameters.zone
    if zone is None:
        # same choice as the atlas: the zone where this field has the most usable arrivals
        census = sc.station_census(sc.load_all(params, iterations=[0]), params.zones)
        counts = {z: n for (s, p, z), n in census.items() if s == station and p == phase}
        if not counts:
            raise SystemExit(f'{station}:{phase} has no usable arrival in {parameters.run_name}')
        zone = max(counts, key=counts.get)
    params.zones = [zone]

    cache = sc.load_all(params, iterations=range(len(sc._CHAR_DISTS)))
    fields = sc.station_fields(cache, zone, station, phase, params.depth,
                               params.map_spacing, params.extent)
    if fields is None:
        raise SystemExit(f'{station}:{phase} has no field in zone {zone}')
    sta_lat, sta_lon = sc.read_station_latlon(zone).get(station, (np.nan, np.nan))

    if not parameters.figSave:
        parameters.figSave = os.path.join(
            _MODULE_DIR, 'ssst_corrections',
            f'{parameters.run_name}_{station}_{phase}_report.png')
    outputs = plot_page(fields, sta_lat, sta_lon, parameters)
    print(f'{station}:{phase}  zone {zone}  depth slice {fields["depth"]:.1f} km  '
          f'arrivals {fields["n_arrivals"]}')
    return {'output': outputs[0], 'outputs': outputs, 'zone': zone, 'depth': fields['depth']}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Report-ready SSST correction page for one station and phase.')
    parser.add_argument('--station', required=True, help="'STA:PHASE', e.g. CA.0051:S")
    parser.add_argument('--run-name', default='ssst_run1')
    parser.add_argument('--zone', type=int, default=None,
                        help="default: the station's best-sampled zone, as in the atlas")
    parser.add_argument('--depth', type=float, default=None,
                        help="slice depth in km; default: the station's median event depth")
    parser.add_argument('--map-spacing', type=float, default=sc._MAP_SPACING,
                        help='node spacing in degrees (default %(default)s, as the atlas)')
    parser.add_argument('--extent', default=','.join(str(v) for v in DEFAULT_EXTENT),
                        help='lon0,lon1,lat0,lat1 — needs the = form when lon0 is negative '
                             '(default: %(default)s); "auto" frames on the station\'s events')
    parser.add_argument('--output', default=None,
                        help='default: complem_figures/ssst_corrections/<run>_<STA>_<PHASE>_report.png')
    parser.add_argument('--width-cm', type=float, default=16.0,
                        help='Printed width, i.e. your \\linewidth in cm (default: 16)')
    args = parser.parse_args()

    extent = None if args.extent == 'auto' else tuple(float(v) for v in args.extent.split(','))
    generate_figure(SsstCorrReportParams(
        station=args.station, figSave=args.output, run_name=args.run_name,
        zone=args.zone, depth=args.depth, map_spacing=args.map_spacing,
        extent=extent, width_cm=args.width_cm))


if __name__ == '__main__':
    main()
