"""
zone_map.py
============================
Generate a PyGMT overview map of all 6 NLL relocation zones.

Each zone is drawn with two nested rectangles:
- Outer box (dashed outline): the NLL velocity-model grid, extended 100 km
  beyond the event selection box.
- Inner box (solid outline + semi-transparent fill): the event selection box
  used to filter obs/GLOBAL.obs for each zone run.

Optionally overlays seismicity from a .obs bulletin, coloured by depth.
A quality filter (erh ≤ 3 km, erv ≤ 3 km, gap ≤ 300°, rms ≤ 0.5 s) is
applied by default; use --no-filter to skip it (e.g. for pre-relocation .obs
files where erh/erv are not available).

Usage
-----
    # Zones only (no seismicity)
    python complem_figures/zone_map.py \\
        --output complem_figures/zone_map/zones_only.pdf

    # With seismicity overlay, no quality filter (pre-relocation catalog)
    python complem_figures/zone_map.py \\
        --bulletin obs/GLOBAL.obs \\
        --output   complem_figures/zone_map/zones_GLOBAL.pdf \\
        --no-filter

    # With seismicity overlay, quality filter applied (post-relocation)
    python complem_figures/zone_map.py \\
        --bulletin obs/FINAL.obs \\
        --output   complem_figures/zone_map/zones_FINAL.pdf
"""

import argparse
import os
from dataclasses import dataclass
from typing import Optional

import pandas as pd
import pygmt as pg

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Inner: event selection box (from prepare_nll_inputs.py)
# Outer: NLL velocity-model grid (inner extended ~100 km, from event_maps.py)
_ZONES = {
    "1": {"inner": ((42.50, -2.00), (43.50, -0.75)), "outer": ((41.60, -3.22), (44.40,  0.46))},
    "2": {"inner": ((42.50, -1.00), (43.25,  0.50)), "outer": ((41.60, -2.22), (44.15,  1.71))},
    "3": {"inner": ((42.00,  0.25), (43.25,  1.00)), "outer": ((41.10, -0.96), (44.15,  2.20))},
    "4": {"inner": ((42.00,  0.75), (43.00,  2.25)), "outer": ((41.10, -0.46), (43.90,  3.45))},
    "5": {"inner": ((42.00,  2.00), (43.00,  3.50)), "outer": ((41.10,  0.79), (43.90,  4.70))},
    "6": {"inner": ((42.75,  2.25), (43.75,  3.50)), "outer": ((41.85,  1.03), (44.65,  4.75))},
}

_ZONE_COLORS = [
    "red", "dodgerblue", "forestgreen", "darkorange", "mediumpurple", "saddlebrown"
]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ZoneMapParams:
    figSave:      str
    fileBulletin: Optional[str]  = None  # .obs file; if None, no seismicity overlay
    no_filter:    bool           = False
    map_region:   Optional[list] = None  # [lon_min, lon_max, lat_min, lat_max]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _remove_high_err(df):
    """Filter out events with location errors or quality metrics above thresholds."""
    df = df[df.erh <= 3.0]
    df = df[df.erv <= 3.0]
    df = df[df.gap <= 300]
    df = df[df.rms <= 0.5]
    return df


def _box_coords(sw, ne):
    """Return (x, y) lists for a closed rectangle from SW and NE corners."""
    lat_sw, lon_sw = sw
    lat_ne, lon_ne = ne
    x = [lon_sw, lon_ne, lon_ne, lon_sw, lon_sw]
    y = [lat_sw, lat_sw, lat_ne, lat_ne, lat_sw]
    return x, y


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    """
    Draw all 6 NLL zones on a Pyrenees basemap and optionally overlay
    seismicity from a .obs bulletin.

    Parameters
    ----------
    parameters : ZoneMapParams

    Returns
    -------
    dict with keys: output
    """
    region = parameters.map_region if parameters.map_region else [-5.0, 6.0, 40.5, 45.5]

    fig = pg.Figure()
    with pg.config(MAP_FRAME_TYPE='fancy+'):
        fig.basemap(region=region, projection='M6i', frame='af')
    fig.coast(water='skyblue', land='#777777', resolution='i',
              area_thresh='0/0/1', borders='1/0.75p,black')

    for (zone_id, boxes), color in zip(_ZONES.items(), _ZONE_COLORS):
        sw_out, ne_out = boxes["outer"]
        sw_in,  ne_in  = boxes["inner"]

        x_out, y_out = _box_coords(sw_out, ne_out)
        x_in,  y_in  = _box_coords(sw_in,  ne_in)

        fig.plot(x=x_out, y=y_out, close=True, pen=f'1p,{color},-')
        fig.plot(x=x_in,  y=y_in,  close=True, pen=f'1.5p,{color}',
                 fill=color, transparency=75)

        label_lat = (sw_in[0] + ne_in[0]) / 2.0
        label_lon = (sw_in[1] + ne_in[1]) / 2.0
        fig.text(x=label_lon, y=label_lat, text=zone_id,
                 font=f'10p,Helvetica-Bold,{color}', justify='CM')

    if parameters.fileBulletin:
        with open(parameters.fileBulletin, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        print(f"Catalog read @ {parameters.fileBulletin}")
        events = [line.lstrip('# ').rstrip('\n').split()
                  for line in lines if line.startswith('# ')]
        events_df = (
            pd.DataFrame(events)
            .drop(columns=[0, 1, 2, 3, 4, 5, 10, 11, 12])
            .rename(columns={6: 'Latitude', 7: 'Longitude', 8: 'Depth',
                              9: 'Magnitude', 13: 'erh', 14: 'erv',
                              15: 'gap', 16: 'rms'})
            .replace('None', float('nan'))
            .astype(float)
        )
        if not parameters.no_filter:
            events_df = _remove_high_err(events_df)

        pg.makecpt(cmap='viridis', series=[0, 15, 1], reverse=True)
        fig.plot(
            x=events_df.Longitude,
            y=events_df.Latitude,
            style='c0.02c',
            fill=events_df.Depth,
            cmap=True,
            transparency=15,
        )
        fig.colorbar(frame=['a5f5+lDepth [km] (events above 15 are in black)'])

    os.makedirs(os.path.dirname(parameters.figSave), exist_ok=True)
    fig.savefig(parameters.figSave, dpi=300)

    print(f"Figure saved @ {parameters.figSave}")
    return {'output': parameters.figSave}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate a PyGMT overview map of all 6 NLL relocation zones.'
    )
    parser.add_argument('--output',    required=True,
                        help='Output figure path (PDF or PNG)')
    parser.add_argument('--bulletin',  default=None,
                        help='Optional .obs bulletin to overlay seismicity')
    parser.add_argument('--no-filter', action='store_true',
                        help='Skip the erh/erv/gap/rms quality filter on events')
    parser.add_argument('--map-region', nargs=4, type=float,
                        metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'),
                        default=None,
                        help='Map extent (default: -5 6 40.5 45.5)')
    args = parser.parse_args()

    generate_figure(ZoneMapParams(
        figSave      = args.output,
        fileBulletin = args.bulletin,
        no_filter    = args.no_filter,
        map_region   = args.map_region,
    ))


if __name__ == '__main__':
    main()
