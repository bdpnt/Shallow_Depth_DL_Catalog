"""
aoi_map_report.py
============================
Report-ready map of the Area of Interest applied to the source bulletins
(report 03, D-11) -- the figure for the AOI paragraph of the deliverable.

This is a *new sibling* of `global_obs/filter_events_by_aoi.py`, not a change
to it: that module filters the bulletins in place as well as drawing its map,
so it must not be re-run merely to regenerate a figure.  Nothing here writes
to any bulletin; the geometry below is copied verbatim from it.

A single panel carries the whole boundary.  The two lines applied to every
bulletin are drawn in grey; the two source-specific lines -- RESIF admitted
only north of its own, IGN and ICGC only south of theirs -- in colour.  The
shaded polygon is their intersection, the overlap band in which all five
agencies contribute.  No legend is drawn: the colours are explained in the
caption (`deliverable_reports/drafts/caption_aoi_map.txt`).  National borders
and shorelines are drawn by the PyGMT backend only.  The RESIF and the IGN + ICGC epicentres are plotted as
they stand after filtering, so that each catalog is seen to terminate on its
own line; LDG and OMP have no source-specific line and are not drawn.

Backends
--------
`--backend pygmt` (default) draws the repo's usual basemap (`fig.coast`,
land `#777777`, water `skyblue`, national borders) and writes a 300 dpi PNG
and a vector PDF.  Needs `pygmt_env`.

`--backend mpl` draws the same geometry with matplotlib and **no coastline**
-- a layout preview for environments without GMT, not the report figure.

Usage
-----
    conda run -n pygmt_env python complem_figures/aoi_map_report.py \
        --output complem_figures/aoi_map/aoi_map_report.png
"""

import argparse
import os
from dataclasses import dataclass, field
from typing import List

import numpy as np

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

# ---------------------------------------------------------------------------
# The AOI geometry -- copied verbatim from global_obs/filter_events_by_aoi.py
#
# Each entry is ((lat1, lon1), (lat2, lon2), aoi_above), where `aoi_above`
# has the meaning of `_is_in_aoi`: True keeps the cross < 0 side, False the
# cross > 0 side.  Verified against the filtered bulletins on 2026-09-21:
# all six satisfy both universal lines, RESIF is 100 % on the kept side of
# its own line and IGN/ICGC 100 % on theirs.
# ---------------------------------------------------------------------------

LINE_NORTH_ALL = ((44.00, -0.25), (43.25,  3.50), False)   # every source: south of it
LINE_SOUTH_ALL = ((42.50, -2.25), (42.00,  0.25), True)    # every source: north of it
LINE_RESIF     = ((43.00, -2.25), (42.00,  2.25), True)    # RESIF: north of it
LINE_IGN_ICGC  = ((43.75, -2.25), (42.00,  6.25), False)   # IGN + ICGC: south of it

UNIVERSAL = [LINE_NORTH_ALL, LINE_SOUTH_ALL]

# Map extent of the report figure: lon_min lon_max lat_min lat_max
DEFAULT_REGION = [-3.25, 4.35, 41.25, 44.45]

# Source bulletins drawn
SRC_RESIF   = ['obs/RESIF_20-25.obs']
SRC_IGNICGC = ['obs/IGN_20-25.obs', 'obs/ICGC_20-25.obs']

# Colours
C_ALL     = '#4D4D4D'
C_RESIF   = '#0072B2'
C_IGNICGC = '#D55E00'
C_BAND    = '#EADFC8'   # warm, so it stays distinct from land and water

# Basemap defaults.  Land is lighter than the repo's usual '#777777' so that
# the two transparent event clouds stay readable over it; pass --land 777777
# to go back to the pipeline maps' shade.
LAND_DEFAULT       = '#F5F5F5'
WATER_DEFAULT      = '#CFE4F0'
SHORELINE_PEN      = '0.4p,gray35'
BORDER_PEN         = '0.9p,black'


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class AOIMapParams:
    figSave:     str         = 'complem_figures/aoi_map/aoi_map_report.png'
    map_region:  List[float] = field(default_factory=lambda: list(DEFAULT_REGION))
    width:       str         = '6i'
    backend:     str         = 'pygmt'
    event_size:  float       = 0.016
    max_events:  int         = 0          # 0 = draw them all
    no_events:   bool        = False
    land:        str         = LAND_DEFAULT
    water:       str         = WATER_DEFAULT


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _cross(lon, lat, line):
    """
    Signed side of `line` for the point (lon, lat).

    Identical to the cross product of `filter_events_by_aoi._is_in_aoi`,
    written for arrays: cross < 0 is the side kept when aoi_above is True.
    """
    (y1, x1), (y2, x2) = line[0], line[1]
    return (np.asarray(lon) - x1) * (y2 - y1) - (np.asarray(lat) - y1) * (x2 - x1)


def _kept(lon, lat, line):
    """Boolean mask: points on the side of `line` that the filter keeps."""
    c = _cross(lon, lat, line)
    return c < 0 if line[2] else c > 0


def _line_across(line, region):
    """
    The boundary line extended to both edges of the map frame.

    Returns (lons, lats), two-point arrays.  Every AOI line is oblique in
    longitude, so no vertical case has to be handled.
    """
    (y1, x1), (y2, x2) = line[0], line[1]
    slope = (y2 - y1) / (x2 - x1)
    lons  = np.array([region[0], region[1]])
    lats  = y1 + slope * (lons - x1)
    return lons, lats


def _clip_halfplane(poly, line):
    """Sutherland-Hodgman clip of a closed (lon, lat) polygon to the kept side."""
    if not poly:
        return []
    out = []
    n   = len(poly)
    for i in range(n):
        cur, nxt = poly[i], poly[(i + 1) % n]
        c_cur = float(_cross(cur[0], cur[1], line))
        c_nxt = float(_cross(nxt[0], nxt[1], line))
        in_cur = c_cur < 0 if line[2] else c_cur > 0
        in_nxt = c_nxt < 0 if line[2] else c_nxt > 0
        if in_cur:
            out.append(cur)
        if in_cur != in_nxt:
            t = c_cur / (c_cur - c_nxt)
            out.append((cur[0] + t * (nxt[0] - cur[0]),
                        cur[1] + t * (nxt[1] - cur[1])))
    return out


def _admitted_polygon(region, lines):
    """The map frame clipped successively to the kept side of every line."""
    poly = [(region[0], region[2]), (region[1], region[2]),
            (region[1], region[3]), (region[0], region[3])]
    for line in lines:
        poly = _clip_halfplane(poly, line)
    return poly


# ---------------------------------------------------------------------------
# Bulletin reading
# ---------------------------------------------------------------------------

def _read_obs(path):
    """
    Epicentres of an NLL `.obs` bulletin.

    Event headers start with '# '; fields 6 and 7 are latitude and longitude,
    the same columns `filter_events_by_aoi` reads.
    """
    lat, lon = [], []
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('# ') and not line.startswith('###'):
                fields = line[2:].split()
                try:
                    lat.append(float(fields[6]))
                    lon.append(float(fields[7]))
                except (IndexError, ValueError):
                    continue
    return np.array(lon), np.array(lat)


def _resolve(path):
    return path if os.path.isabs(path) else os.path.join(_PROJECT_ROOT, path)


def _read_group(paths, max_events):
    """Concatenated epicentres of several bulletins; missing files are skipped."""
    lons, lats = [], []
    for rel in paths:
        full = _resolve(rel)
        if not os.path.exists(full):
            print(f'  (skipped, not found: {rel})')
            continue
        lo, la = _read_obs(full)
        print(f'  {len(lo):6d} epicentres @ {rel}')
        lons.append(lo)
        lats.append(la)
    if not lons:
        return np.array([]), np.array([])
    lon, lat = np.concatenate(lons), np.concatenate(lats)
    if max_events and len(lon) > max_events:
        idx = np.linspace(0, len(lon) - 1, max_events).astype(int)
        lon, lat = lon[idx], lat[idx]
    return lon, lat


# ---------------------------------------------------------------------------
# PyGMT backend
# ---------------------------------------------------------------------------

_PEN_ALL     = f'0.9p,{C_ALL}'
_PEN_RESIF   = f'1.6p,{C_RESIF},8p_4p'
_PEN_IGNICGC = f'1.6p,{C_IGNICGC},8p_4p'


def _render_pygmt(parameters, ev_resif, ev_ignicgc):
    import pygmt as pg

    region = parameters.map_region
    fig    = pg.Figure()

    with pg.config(MAP_FRAME_TYPE='fancy+', FONT_ANNOT_PRIMARY='9p',
                   FONT_LABEL='10p'):
        fig.basemap(region=region, projection=f'M{parameters.width}',
                    frame=['af', 'xa1f0.5+lLongitude', 'ya1f0.5+lLatitude'])
    fig.coast(water=parameters.water, land=parameters.land, resolution='i',
              area_thresh='0/0/1', shorelines=SHORELINE_PEN,
              borders=[f'1/{BORDER_PEN}'])

    poly = _admitted_polygon(region, UNIVERSAL + [LINE_RESIF, LINE_IGN_ICGC])
    fig.plot(x=[p[0] for p in poly], y=[p[1] for p in poly],
             fill=C_BAND, transparency=35, close=True)

    for ev, colour in ((ev_ignicgc, C_IGNICGC), (ev_resif, C_RESIF)):
        if len(ev[0]):
            fig.plot(x=ev[0], y=ev[1], style=f'c{parameters.event_size}c',
                     fill=colour, transparency=72)

    for line in UNIVERSAL:
        xs, ys = _line_across(line, region)
        fig.plot(x=xs, y=ys, pen=_PEN_ALL)
    for line, pen in ((LINE_IGN_ICGC, _PEN_IGNICGC), (LINE_RESIF, _PEN_RESIF)):
        xs, ys = _line_across(line, region)
        fig.plot(x=xs, y=ys, pen=pen)

    out  = _resolve(parameters.figSave)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    root = os.path.splitext(out)[0]
    fig.savefig(root + '.png', dpi=300)
    fig.savefig(root + '.pdf')
    print(f'Figure successfully saved @ {root}.png and {root}.pdf')
    return root + '.png'


# ---------------------------------------------------------------------------
# Matplotlib backend -- layout preview, no coastline
# ---------------------------------------------------------------------------

def _render_mpl(parameters, ev_resif, ev_ignicgc):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    matplotlib.rcParams.update({'pdf.fonttype': 42, 'font.size': 9})

    region = parameters.map_region
    aspect = 1.0 / np.cos(np.radians(np.mean(region[2:])))

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.set_xlim(region[0], region[1])
    ax.set_ylim(region[2], region[3])
    ax.set_aspect(aspect)
    ax.grid(True, lw=0.4, color='0.90', zorder=0)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)

    poly = _admitted_polygon(region, UNIVERSAL + [LINE_RESIF, LINE_IGN_ICGC])
    ax.add_patch(Polygon(poly, closed=True, facecolor=C_BAND, alpha=0.55,
                         edgecolor='none', zorder=1))

    for ev, colour in ((ev_ignicgc, C_IGNICGC), (ev_resif, C_RESIF)):
        if len(ev[0]):
            ax.scatter(ev[0], ev[1], s=1.0, c=colour, alpha=0.25,
                       linewidths=0, zorder=2)

    for line in UNIVERSAL:
        xs, ys = _line_across(line, region)
        ax.plot(xs, ys, color=C_ALL, lw=1.0, zorder=4)
    for line, colour in ((LINE_IGN_ICGC, C_IGNICGC), (LINE_RESIF, C_RESIF)):
        xs, ys = _line_across(line, region)
        ax.plot(xs, ys, color=colour, lw=1.7, ls=(0, (5, 2.5)), zorder=5)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    fig.tight_layout()
    out  = _resolve(parameters.figSave)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    root = os.path.splitext(out)[0]
    fig.savefig(root + '.png', dpi=300)
    fig.savefig(root + '.pdf')
    plt.close(fig)
    print(f'Preview saved @ {root}.png and {root}.pdf')
    return root + '.png'


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def aoi_map_report(parameters):
    ev_resif = ev_ignicgc = (np.array([]), np.array([]))
    if not parameters.no_events:
        print('RESIF')
        ev_resif = _read_group(SRC_RESIF, parameters.max_events)
        print('IGN + ICGC')
        ev_ignicgc = _read_group(SRC_IGNICGC, parameters.max_events)

    if parameters.backend == 'pygmt':
        path = _render_pygmt(parameters, ev_resif, ev_ignicgc)
    else:
        path = _render_mpl(parameters, ev_resif, ev_ignicgc)
    return {'fig_save': path,
            'n_resif': int(len(ev_resif[0])),
            'n_ign_icgc': int(len(ev_ignicgc[0]))}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Map of the per-source Area of Interest (report 03, D-11).')
    parser.add_argument('--output', default='complem_figures/aoi_map/aoi_map_report.png')
    parser.add_argument('--backend', choices=['pygmt', 'mpl'], default='pygmt')
    parser.add_argument('--map-region', nargs=4, type=float, default=DEFAULT_REGION,
                        metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'))
    parser.add_argument('--width', default='6i', help='PyGMT map width')
    parser.add_argument('--event-size', type=float, default=0.016,
                        help='PyGMT circle diameter [cm]')
    parser.add_argument('--max-events', type=int, default=0,
                        help='decimate each group to at most N epicentres (0 = all)')
    parser.add_argument('--no-events', action='store_true',
                        help='draw the lines alone, without the source catalogs')
    parser.add_argument('--land', default=LAND_DEFAULT,
                        help="land fill; '#777777' restores the pipeline maps' shade")
    parser.add_argument('--water', default=WATER_DEFAULT, help='water fill')
    args = parser.parse_args()

    params = AOIMapParams(
        figSave=args.output,
        map_region=list(args.map_region),
        width=args.width,
        backend=args.backend,
        event_size=args.event_size,
        max_events=args.max_events,
        no_events=args.no_events,
        land=args.land,
        water=args.water,
    )
    aoi_map_report(params)


if __name__ == '__main__':
    main()
