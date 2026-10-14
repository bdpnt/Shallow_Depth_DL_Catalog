"""
event_maps.py
============================
Generate a PyGMT map of seismic events coloured by depth.

Reads a .csv (RESULT/*.csv) or .obs bulletin, optionally filters high-error
events, and plots each event on a Pyrenees basemap coloured by depth.
Optionally overlays station positions, zone-boundary rectangles and
cross-section profiles. Writes both a 300 dpi PNG and a vector PDF.

Quality filter (erh ≤ 3 km, erv ≤ 3 km, gap ≤ 300°, rms ≤ 0.5 s) is applied
by default. Use --no-filter to skip it (e.g. for pre-relocation .obs files
where erh/erv are not available).

Hypocentre of a result CSV
--------------------------
`--solution expect` (default) reads the PDF expectation (`expect_lat`,
`expect_lon`, `expect_z`), the published hypocentre (CLAUDE.md §4).
`--solution maxlike` reads the maximum-likelihood point from `maxlike_*`; a CSV
written before `SAVE_NLLOC_EXPECTATION` has no such columns, and there
`latitude`/`longitude`/`depth` are the maximum-likelihood point, so they are
read instead (a note is printed). A .obs bulletin carries one hypocentre only.

Cross-section profiles
----------------------
`--profile` overlays a cross-section: a solid red axis line from the profile
origin along its azimuth, and a dashed red rectangle showing the swath
actually projected (origin ± half-width, perpendicular to the axis).
Repeatable; `--profiles-default` draws the three profiles of the report
(Aneto, Arette, Andorra). `--box-pen` takes an absolute dash pattern (`3p_2p`)
rather than GMT's `--` shorthand, whose dash and gap lengths are multiples of
the pen width and so grow with it. The swath geometry is computed with the
same flat-Earth helper (`_dest_point`, 1 deg = 111 km) that `cross_section.py`
uses to select the events, so the rectangle drawn is exactly the strip
projected onto the section.

Usage
-----
    # .obs bulletin — pre-relocation, no quality filter
    python complem_figures/event_maps.py \\
        --bulletin  obs/GLOBAL.obs \\
        --output    complem_figures/event_maps/GLOBAL.png \\
        --no-filter --map-region -2.25 3.5 41.75 43.75

    # .csv NLL result — single zone from NLL_result.csv
    python complem_figures/event_maps.py \\
        --bulletin      RESULT/NLL_result.csv \\
        --output        complem_figures/event_maps/GLOBAL_1.png \\
        --stations      run/nll_loc/GLOBAL_1/last.stations \\
        --source-filter GLOBAL_1 \\
        --region-in     42.50 -2.00 43.50 -0.75 \\
        --region-out    41.60 -3.22 44.40  0.46

    # the report map: SSST catalogue, no quality filter, the three profiles
    conda run -n pygmt_env python complem_figures/event_maps.py \\
        --bulletin   RESULT/SSST_result.csv \\
        --output     complem_figures/event_maps/SSST_result.png \\
        --no-filter --profiles-default --map-region -2.5 3.4 41.75 43.75

    # one custom profile: lon0 lat0 azimuth length half-width
    conda run -n pygmt_env python complem_figures/event_maps.py \\
        --bulletin   RESULT/SSST_result.csv \\
        --output     complem_figures/event_maps/SSST_result.png \\
        --no-filter --profile 0.6852 42.6 44 12 3
"""

import argparse
import os
from dataclasses import dataclass, field
from math import cos, radians, sin
from typing import List, Optional, Tuple

import pandas as pd

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

# Default map extent (wide enough for the zone outer boxes), and the tighter
# extent of the report map: lon_min lon_max lat_min lat_max
DEFAULT_REGION = [-4.0, 4.0, 41.0, 45.0]
REPORT_REGION  = [-2.5, 3.4, 41.75, 43.75]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class EventMapsParams:
    fileBulletin:  str
    figSave:       str
    fileStations:  Optional[str]   = None
    region_in:     Optional[tuple] = None  # ((lat_min, lon_min), (lat_max, lon_max))
    region_out:    Optional[tuple] = None  # ((lat_min, lon_min), (lat_max, lon_max))
    map_region:    Optional[list]  = None  # [lon_min, lon_max, lat_min, lat_max]
    no_filter:     bool            = False
    source_filter: Optional[str]   = None  # e.g. "GLOBAL_1"; None → all sources
    solution:      str             = 'expect'   # 'expect' | 'maxlike' (result CSV only)
    max_depth:     float           = 15.0       # colour-bar saturation
    profiles:      List[Tuple[float, float, float, float, float]] = field(default_factory=list)
    axis_pen:      str             = '0.9p,red'
    box_pen:       str             = '0.4p,red,3p_2p'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def hypocentre_columns(df, solution):
    """
    (lat, lon, depth) column names of a result CSV for the chosen solution.

    'expect'  → expect_lat / expect_lon / expect_z (the published hypocentre).
    'maxlike' → maxlike_latitude / maxlike_longitude / maxlike_depth; a CSV
                written before SAVE_NLLOC_EXPECTATION has none, and there
                latitude / longitude / depth are the maximum-likelihood point.
    """
    if solution == 'expect':
        return 'expect_lat', 'expect_lon', 'expect_z'
    if solution != 'maxlike':
        raise ValueError(f'unknown solution {solution!r} (expect or maxlike)')
    if {'maxlike_latitude', 'maxlike_longitude', 'maxlike_depth'} <= set(df.columns):
        return 'maxlike_latitude', 'maxlike_longitude', 'maxlike_depth'
    print('Note: no maxlike_* columns (CSV predates SAVE_NLLOC_EXPECTATION) — '
          'latitude/longitude/depth are the maximum-likelihood point there, reading those')
    return 'latitude', 'longitude', 'depth'


def _remove_high_err(df):
    """Filter out events with location errors or quality metrics above thresholds."""
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
        if parameters.solution != 'expect':
            raise SystemExit('--solution is meaningful only for a RESULT/*.csv; '
                             'a .obs bulletin carries a single hypocentre.')
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

    elif ext == 'csv':
        df = pd.read_csv(parameters.fileBulletin, skipinitialspace=True)
        lat_col, lon_col, depth_col = hypocentre_columns(df, parameters.solution)
        for col in (lat_col, lon_col, depth_col):
            if col not in df.columns:
                raise SystemExit(f'Column "{col}" absent from {parameters.fileBulletin}')

        events_df = df.rename(columns={
            lat_col:    'Latitude',
            lon_col:    'Longitude',
            depth_col:  'Depth',
            'RMS':      'rms',
            'true_erh': 'erh',
            'true_erz': 'erv',
            'Gap':      'gap',
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
    Draw the depth-coloured event map, with the optional zone boxes, stations
    and cross-section profiles.

    Parameters
    ----------
    parameters : EventMapsParams

    Returns
    -------
    dict with keys: output, outputs, n_events
    """
    import pygmt as pg   # here, so seisbench_env modules can import hypocentre_columns

    events_df = _read_catalog(parameters)
    region    = parameters.map_region if parameters.map_region else list(DEFAULT_REGION)

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

    is_csv = parameters.fileBulletin.endswith('.csv')
    label  = {'expect': 'Expectation depth', 'maxlike': 'Max-likelihood depth'}[
        parameters.solution] if is_csv else 'Depth'
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
        description='Generate a PyGMT depth-coloured event map.'
    )
    parser.add_argument('--bulletin',   required=True,
                        help='Input bulletin file (.csv result or .obs bulletin)')
    parser.add_argument('--output',     required=True,
                        help='Output path; the extension is replaced, both .png and .pdf are written')
    parser.add_argument('--stations',   default=None,
                        help='Optional last.stations file to overlay station positions')
    parser.add_argument('--solution',   default='expect', choices=['expect', 'maxlike'],
                        help='Hypocentre of a result CSV: PDF expectation (default) '
                             'or maximum likelihood')
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
                        help=f'Map extent (default: {" ".join(str(v) for v in DEFAULT_REGION)}; '
                             f'report map: {" ".join(str(v) for v in REPORT_REGION)})')
    parser.add_argument('--max-depth', type=float, default=15.0,
                        help='Colour-bar saturation depth in km (default: 15)')
    parser.add_argument('--no-filter', action='store_true',
                        help='Skip the erh/erv/gap/rms quality filter (useful for pre-relocation .obs)')
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

    ri = args.region_in
    ro = args.region_out
    profiles = [tuple(p) for p in (args.profile or [])]
    if args.profiles_default:
        profiles = list(DEFAULT_PROFILES) + profiles

    generate_figure(EventMapsParams(
        fileBulletin  = args.bulletin,
        figSave       = args.output,
        fileStations  = args.stations,
        region_in     = ((ri[0], ri[1]), (ri[2], ri[3])) if ri else None,
        region_out    = ((ro[0], ro[1]), (ro[2], ro[3])) if ro else None,
        map_region    = args.map_region,
        no_filter     = args.no_filter,
        source_filter = args.source_filter,
        solution      = args.solution,
        max_depth     = args.max_depth,
        profiles      = profiles,
        axis_pen      = args.axis_pen,
        box_pen       = args.box_pen,
    ))


if __name__ == '__main__':
    main()
