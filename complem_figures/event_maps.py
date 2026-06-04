"""
event_maps.py
============================
Generate a PyGMT map of seismic events coloured by depth.

Reads a .txt (NLL result) or .obs bulletin, optionally filters high-error
events, and plots each event on a Pyrenees basemap coloured by depth.
Optionally overlays station positions and zone-boundary rectangles.

Quality filter (erh ≤ 3 km, erv ≤ 3 km, gap ≤ 300°, rms ≤ 0.5 s) is applied
by default. Use --no-filter to skip it (e.g. for pre-relocation .obs files
where erh/erv are not available).

Usage
-----
    # .obs bulletin — pre-relocation, no quality filter
    python complem_figures/event_maps.py \\
        --bulletin  obs/GLOBAL.obs \\
        --output    complem_figures/event_maps/GLOBAL.pdf \\
        --no-filter

    # .obs bulletin — post-relocation, with filter
    python complem_figures/event_maps.py \\
        --bulletin  obs/FINAL.obs \\
        --output    complem_figures/event_maps/FINAL.pdf

    # .txt result with stations and zone boxes
    python complem_figures/event_maps.py \\
        --bulletin   RESULT/GLOBAL_1_PR.txt \\
        --output     complem_figures/event_maps/GLOBAL_1_PR.pdf \\
        --stations   loc/GLOBAL_1/last.stations \\
        --region-in  42.50 -2.00 43.50 -0.75 \\
        --region-out 41.60 -3.22 44.40  0.46

    # custom map extent (zoom in on eastern Pyrenees)
    python complem_figures/event_maps.py \\
        --bulletin  obs/FINAL.obs \\
        --output    complem_figures/event_maps/FINAL_east.pdf \\
        --map-region 0.0 3.5 42.0 44.0
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
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class EventMapsParams:
    fileBulletin: str
    figSave:      str
    fileStations: Optional[str]   = None
    region_in:    Optional[tuple] = None  # ((lat_min, lon_min), (lat_max, lon_max))
    region_out:   Optional[tuple] = None  # ((lat_min, lon_min), (lat_max, lon_max))
    map_region:   Optional[list]  = None  # [lon_min, lon_max, lat_min, lat_max]
    no_filter:    bool            = False


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


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    """
    Read a .txt or .obs bulletin, filter high-error events, and save a
    PyGMT map coloured by depth.

    Parameters
    ----------
    parameters : EventMapsParams

    Returns
    -------
    dict with keys: output
    """
    ext = parameters.fileBulletin.split('.')[-1]

    if ext == 'txt':
        with open(parameters.fileBulletin, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        print(f"Catalog read @ {parameters.fileBulletin}")
        events = [event.split() for event in lines]
        events_df = (
            pd.DataFrame(events)
            .drop(columns=[0, 1, 2, 3, 4, 5, 11])
            .rename(columns={6: 'Latitude', 7: 'Longitude', 8: 'Depth',
                              9: 'Magnitude', 10: 'rms', 12: 'erh',
                              13: 'erv', 14: 'gap'})
            .astype(float)
        )
        if not parameters.no_filter:
            events_df = _remove_high_err(events_df)

    elif ext == 'obs':
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

    else:
        print(f'Unsupported format (expected "txt" or "obs"): {parameters.fileBulletin}')
        return {'output': None}

    region = parameters.map_region if parameters.map_region else [-4.0, 4, 41, 45]

    fig = pg.Figure()
    with pg.config(MAP_FRAME_TYPE='fancy+'):
        fig.basemap(region=region, projection='M6i', frame='af')
    fig.coast(water='skyblue', land='#777777', resolution='i',
              area_thresh='0/0/1', borders='1/0.75p,black')

    if parameters.region_out:
        ro = parameters.region_out
        fig.plot(
            x=[ro[0][1], ro[1][1], ro[1][1], ro[0][1], ro[0][1]],
            y=[ro[0][0], ro[0][0], ro[1][0], ro[1][0], ro[0][0]],
            close=True, pen='2p,red', transparency=50,
        )

    if parameters.region_in:
        ri = parameters.region_in
        fig.plot(
            x=[ri[0][1], ri[1][1], ri[1][1], ri[0][1], ri[0][1]],
            y=[ri[0][0], ri[0][0], ri[1][0], ri[1][0], ri[0][0]],
            close=True, pen='0.5p,blue', fill='blue', transparency=85,
        )

    if parameters.fileStations:
        stations = pd.read_csv(
            parameters.fileStations, header=0, delimiter=' ',
            names=['Code', 'x', 'y', 'z', 'Latitude', 'Longitude', 'Depth'],
        ).drop(columns=['x', 'y', 'z'])
        fig.plot(x=stations.Longitude, y=stations.Latitude,
                 style='i0.1c', fill='black', transparency=40)

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
    fig.savefig(parameters.figSave, dpi=300)

    print(f"Figure saved @ {parameters.figSave}")
    return {'output': parameters.figSave}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate a PyGMT depth-coloured event map.'
    )
    parser.add_argument('--bulletin',   required=True,
                        help='Input bulletin file (.txt NLL result or .obs)')
    parser.add_argument('--output',     required=True,
                        help='Output figure path (PDF or PNG)')
    parser.add_argument('--stations',   default=None,
                        help='Optional last.stations file to overlay station positions')
    parser.add_argument('--region-in',  nargs=4, type=float,
                        metavar=('LAT_MIN', 'LON_MIN', 'LAT_MAX', 'LON_MAX'),
                        default=None,
                        help='Inner zone box corners: lat_min lon_min lat_max lon_max')
    parser.add_argument('--region-out', nargs=4, type=float,
                        metavar=('LAT_MIN', 'LON_MIN', 'LAT_MAX', 'LON_MAX'),
                        default=None,
                        help='Outer zone box corners: lat_min lon_min lat_max lon_max')
    parser.add_argument('--map-region', nargs=4, type=float,
                        metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'),
                        default=None,
                        help='Map extent: lon_min lon_max lat_min lat_max (default: -4 4 41 45)')
    parser.add_argument('--no-filter', action='store_true',
                        help='Skip the erh/erv/gap/rms quality filter (useful for pre-relocation .obs)')
    args = parser.parse_args()

    ri = args.region_in
    ro = args.region_out

    generate_figure(EventMapsParams(
        fileBulletin = args.bulletin,
        figSave      = args.output,
        fileStations = args.stations,
        region_in    = ((ri[0], ri[1]), (ri[2], ri[3])) if ri else None,
        region_out   = ((ro[0], ro[1]), (ro[2], ro[3])) if ro else None,
        map_region   = args.map_region,
        no_filter    = args.no_filter,
    ))


if __name__ == '__main__':
    main()

    ''' ALL OPTIONS TO USE
    all_runs = {
        "1": ("RESULT/GLOBAL_1.txt", "loc/GLOBAL_1/last.stations", "complem_figures/event_maps/GLOBAL_1.pdf",
              ((42.50, -2.00), (43.50, -0.75)), ((41.60, -3.22), (44.40, 0.46))),
        "2": ("RESULT/GLOBAL_2.txt", "loc/GLOBAL_2/last.stations", "complem_figures/event_maps/GLOBAL_2.pdf",
              ((42.50, -1.00), (43.25, 0.50)), ((41.60, -2.22), (44.15, 1.71))),
        "3": ("RESULT/GLOBAL_3.txt", "loc/GLOBAL_3/last.stations", "complem_figures/event_maps/GLOBAL_3.pdf",
              ((42.00, 0.25), (43.25, 1.00)), ((41.10, -0.96), (44.15, 2.20))),
        "4": ("RESULT/GLOBAL_4.txt", "loc/GLOBAL_4/last.stations", "complem_figures/event_maps/GLOBAL_4.pdf",
              ((42.00, 0.75), (43.00, 2.25)), ((41.10, -0.46), (43.90, 3.45))),
        "5": ("RESULT/GLOBAL_5.txt", "loc/GLOBAL_5/last.stations", "complem_figures/event_maps/GLOBAL_5.pdf",
              ((42.00, 2.00), (43.00, 3.50)), ((41.10, 0.79), (43.90, 4.70))),
        "6": ("RESULT/GLOBAL_6.txt", "loc/GLOBAL_6/last.stations", "complem_figures/event_maps/GLOBAL_6.pdf",
              ((42.75, 2.25), (43.75, 3.50)), ((41.85, 1.03), (44.65, 4.75))),
        "Final": ("obs/FINAL.obs", None, "complem_figures/event_maps/FINAL.pdf",
                  None, None),
    }
    '''