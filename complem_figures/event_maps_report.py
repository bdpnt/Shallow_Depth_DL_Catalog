"""
event_maps_report.py
============================
Report-ready version of `event_maps.py`: a PyGMT map of the seismic catalogue
coloured by depth, with the cross-section profiles drawn on top.

Differences from `event_maps.py`
--------------------------------
* `--column` chooses the depth column of a result CSV (`depth`, the
  maximum-likelihood solution, or `expect_z`, the PDF expectation — see
  `NLL_run/pdf_metrics.py`). `--position` likewise chooses between the
  maximum-likelihood epicentre (`latitude`/`longitude`) and the expectation
  epicentre (`expect_lat`/`expect_lon`).
* `--profile` overlays a cross-section: a solid red axis line from the profile
  origin along its azimuth, and a dashed red rectangle showing the swath
  actually projected (origin ± half-width, perpendicular to the axis).
  Repeatable; `--profiles-default` draws the three profiles of the report
  (Aneto, Arette, Andorra). Both pens are tunable — `--box-pen` takes an
  absolute dash pattern (`3p_2p`) rather than GMT's `--` shorthand, whose dash
  and gap lengths are multiples of the pen width and so grow with it.
* Writes both a 300 dpi PNG and a vector PDF.

The swath geometry is computed with the same flat-Earth helper
(`_dest_point`, 1 deg = 111 km) that `complem_figures/cross_section.py` uses to
select the events, so the rectangle drawn here is exactly the strip projected
onto the section — not an approximation of it.

Usage
-----
    # the report map: SSST catalogue, expectation depths, no quality filter
    conda run -n pygmt_env python complem_figures/event_maps_report.py \\
        --bulletin RESULT/SSST_result.csv \\
        --output   complem_figures/event_maps/SSST_result_report.png \\
        --column   expect_z --no-filter --profiles-default

    # one custom profile: lon0 lat0 azimuth length half-width
    conda run -n pygmt_env python complem_figures/event_maps_report.py \\
        --bulletin RESULT/SSST_result.csv \\
        --output   complem_figures/event_maps/SSST_result_report.png \\
        --column   expect_z --no-filter \\
        --profile  0.6852 42.6 44 12 3
"""

import argparse
import os
from dataclasses import dataclass, field
from math import cos, radians, sin
from typing import List, Optional, Tuple

import pandas as pd
import pygmt as pg

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

# ---------------------------------------------------------------------------
# Cross-sections of the report
#   (lon0, lat0, azimuth [deg from N], length [km], half-width [km])
# ---------------------------------------------------------------------------

DEFAULT_PROFILES = [
    ( 0.6852, 42.6, 44.0, 12.0, 3.0),   # Aneto
    (-0.6275, 43.0,  0.0, 16.0, 8.0),   # Arette
    ( 1.4000, 42.48, 32.0, 10.0, 3.0),   # Andorra
]

# Map extent of the report figure: lon_min lon_max lat_min lat_max
DEFAULT_REGION = [-2.5, 3.4, 41.75, 43.75]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class EventMapsReportParams:
    fileBulletin:  str
    figSave:       str
    fileStations:  Optional[str]   = None
    map_region:    Optional[list]  = None
    no_filter:     bool            = False
    source_filter: Optional[str]   = None
    depth_column:  str             = 'depth'     # 'depth' | 'expect_z'
    position:      str             = 'maxlike'   # 'maxlike' | 'expect'
    max_depth:     float           = 15.0        # colour-bar saturation
    profiles:      List[Tuple[float, float, float, float, float]] = field(default_factory=list)
    axis_pen:      str             = '0.9p,red'
    box_pen:       str             = '0.4p,red,3p_2p'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _remove_high_err(df):
    """Quality filter of event_maps.py: erh/erv <= 3 km, gap <= 300 deg, rms <= 0.5 s."""
    df = df[df.erh <= 3.0]
    df = df[df.erv <= 3.0]
    df = df[df.gap <= 300]
    df = df[df.rms <= 0.5]
    return df


def _dest_point(lon, lat, azimut, dist_km):
    """
    Destination point from an origin, an azimuth and a distance.

    Flat-Earth approximation, 1 deg = 111 km — identical to the helper in
    `cross_section.py`, so the swath drawn here matches the swath projected.
    """
    R    = 111.0
    az   = radians(azimut)
    dlat = (dist_km * cos(az)) / R
    dlon = (dist_km * sin(az)) / (R * cos(radians(lat)))
    return lon + dlon, lat + dlat


def profile_geometry(lon0, lat0, azimut, length_km, half_width_km):
    """
    Axis and swath outline of a cross-section.

    Returns
    -------
    (axis, box) — two (xs, ys) tuples, the second closed on itself.
    """
    lon_end, lat_end = _dest_point(lon0, lat0, azimut, length_km)

    # Swath edges sit at +/- half-width perpendicular to the axis, i.e. along
    # azimut -/+ 90 -- not due east-west (cf. cross_section.py).
    lon1a, lat1a = _dest_point(lon0,  lat0,  azimut - 90.0, half_width_km)
    lon2a, lat2a = _dest_point(lon1a, lat1a, azimut,        length_km)
    lon1b, lat1b = _dest_point(lon0,  lat0,  azimut + 90.0, half_width_km)
    lon2b, lat2b = _dest_point(lon1b, lat1b, azimut,        length_km)

    axis = ([lon0, lon_end], [lat0, lat_end])
    box  = ([lon1a, lon2a, lon2b, lon1b, lon1a],
            [lat1a, lat2a, lat2b, lat1b, lat1a])
    return axis, box


def _read_catalog(parameters):
    """Read a .obs bulletin or a RESULT/*.csv into Longitude / Latitude / Depth."""
    ext = parameters.fileBulletin.split('.')[-1]

    if ext == 'obs':
        with open(parameters.fileBulletin, 'r', encoding='utf-8') as f:
            lines = f.readlines()
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
        if parameters.depth_column != 'depth':
            raise SystemExit('--column is meaningful only for a RESULT/*.csv; '
                             'a .obs bulletin carries the maximum-likelihood depth only.')

    elif ext == 'csv':
        df = pd.read_csv(parameters.fileBulletin, skipinitialspace=True)

        lat_col = 'expect_lat' if parameters.position == 'expect' else 'latitude'
        lon_col = 'expect_lon' if parameters.position == 'expect' else 'longitude'
        for col in (lat_col, lon_col, parameters.depth_column):
            if col not in df.columns:
                raise SystemExit(f'Column "{col}" absent from {parameters.fileBulletin}')

        events_df = df.rename(columns={
            lat_col:                  'Latitude',
            lon_col:                  'Longitude',
            parameters.depth_column:  'Depth',
            'RMS':                    'rms',
            'true_erh':               'erh',
            'true_erz':               'erv',
            'Gap':                    'gap',
            'source':                 'source',
        })
        if parameters.source_filter and 'source' in events_df.columns:
            events_df = events_df[events_df['source'] == parameters.source_filter]

    else:
        raise SystemExit(f'Unsupported format (expected "csv" or "obs"): {parameters.fileBulletin}')

    print(f"Catalog read @ {parameters.fileBulletin} — {len(events_df)} events")
    if not parameters.no_filter:
        events_df = _remove_high_err(events_df)
        print(f"Quality filter applied — {len(events_df)} events kept")

    return events_df.dropna(subset=['Latitude', 'Longitude', 'Depth'])


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    """
    Draw the depth-coloured event map and overlay the cross-section profiles.

    Returns
    -------
    dict with keys: output, outputs, n_events
    """
    events_df = _read_catalog(parameters)
    region    = parameters.map_region if parameters.map_region else list(DEFAULT_REGION)

    fig = pg.Figure()
    with pg.config(MAP_FRAME_TYPE='fancy+'):
        fig.basemap(region=region, projection='M6i', frame='af')
    fig.coast(water='skyblue', land='#777777', resolution='i',
              area_thresh='0/0/1', borders='1/0.75p,black')

    if parameters.fileStations:
        stations = pd.read_csv(
            parameters.fileStations, header=0, delimiter=' ',
            names=['Code', 'x', 'y', 'z', 'Latitude', 'Longitude', 'Depth'],
        ).drop(columns=['x', 'y', 'z'])
        fig.plot(x=stations.Longitude, y=stations.Latitude,
                 style='i0.1c', fill='black', transparency=40)

    pg.makecpt(cmap='viridis', series=[0, parameters.max_depth, 1], reverse=True)
    fig.plot(
        x=events_df.Longitude,
        y=events_df.Latitude,
        style='c0.02c',
        fill=events_df.Depth,
        cmap=True,
        transparency=15,
    )

    # -- Cross-sections, drawn last so they stay legible over the seismicity --
    for (lon0, lat0, azimut, length_km, half_width_km) in parameters.profiles:
        axis, box = profile_geometry(lon0, lat0, azimut, length_km, half_width_km)
        fig.plot(x=box[0],  y=box[1],  pen=parameters.box_pen)
        fig.plot(x=axis[0], y=axis[1], pen=parameters.axis_pen)

    label = 'Depth' if parameters.depth_column == 'depth' else 'Expectation depth'
    fig.colorbar(frame=[f'a5f5+l{label} [km] '
                        f'(events above {parameters.max_depth:g} are in black)'])

    base, _  = os.path.splitext(parameters.figSave)
    outputs  = []
    for path in (f'{base}.png', f'{base}.pdf'):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path, dpi=300)
        outputs.append(path)
        print(f"Figure saved @ {path}")

    return {'output': outputs[0], 'outputs': outputs, 'n_events': len(events_df)}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Report-ready PyGMT depth map with cross-section profiles overlaid.'
    )
    parser.add_argument('--bulletin', required=True,
                        help='Input bulletin (.csv result or .obs bulletin)')
    parser.add_argument('--output',   required=True,
                        help='Output path; the extension is replaced, both .png and .pdf are written')
    parser.add_argument('--stations', default=None,
                        help='Optional last.stations file to overlay station positions')
    parser.add_argument('--column',   default='depth', choices=['depth', 'expect_z'],
                        help='Depth column of a result CSV (default: depth, the maximum-likelihood solution)')
    parser.add_argument('--position', default='maxlike', choices=['maxlike', 'expect'],
                        help='Epicentre of a result CSV: maximum-likelihood (default) or PDF expectation')
    parser.add_argument('--map-region', nargs=4, type=float,
                        metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'),
                        default=None,
                        help=f'Map extent (default: {" ".join(str(v) for v in DEFAULT_REGION)})')
    parser.add_argument('--max-depth', type=float, default=15.0,
                        help='Colour-bar saturation depth in km (default: 15)')
    parser.add_argument('--no-filter', action='store_true',
                        help='Skip the erh/erv/gap/rms quality filter')
    parser.add_argument('--source-filter', default=None,
                        help='Keep only events from this source zone (e.g. "GLOBAL_1"); CSV only')
    parser.add_argument('--profile', action='append', nargs=5, type=float,
                        metavar=('LON0', 'LAT0', 'AZIMUTH', 'LENGTH_KM', 'HALF_WIDTH_KM'),
                        default=None,
                        help='Overlay a cross-section; repeatable')
    parser.add_argument('--profiles-default', action='store_true',
                        help='Overlay the three cross-sections of the report (Aneto, Arette, Andorra)')
    parser.add_argument('--axis-pen', default='0.9p,red',
                        help='GMT pen for the profile axis line (default: 0.9p,red)')
    parser.add_argument('--box-pen',  default='0.4p,red,3p_2p',
                        help='GMT pen for the swath outline; the dash pattern is given in '
                             'absolute points (dash_gap), not as "--", so it does not scale '
                             'with the pen width (default: 0.4p,red,3p_2p)')
    args = parser.parse_args()

    profiles = [tuple(p) for p in (args.profile or [])]
    if args.profiles_default:
        profiles = list(DEFAULT_PROFILES) + profiles

    generate_figure(EventMapsReportParams(
        fileBulletin  = args.bulletin,
        figSave       = args.output,
        fileStations  = args.stations,
        map_region    = args.map_region,
        no_filter     = args.no_filter,
        source_filter = args.source_filter,
        depth_column  = args.column,
        position      = args.position,
        max_depth     = args.max_depth,
        profiles      = profiles,
        axis_pen      = args.axis_pen,
        box_pen       = args.box_pen,
    ))


if __name__ == '__main__':
    main()
