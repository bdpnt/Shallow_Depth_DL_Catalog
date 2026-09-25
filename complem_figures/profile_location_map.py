"""
profile_location_map.py
=======================
Location map of cross-section profiles: the SSST-relocated seismicity of the
central Pyrenees coloured by depth, the profiles drawn on top and labelled
X-X', a few towns for orientation, and an inset of the whole range with the
frame of the main map.

Profile convention (same as `cross_section.py` and `event_maps_report.py`)
------------------------------------------------------------------------
    NAME LON0 LAT0 AZIMUTH LENGTH_KM HALF_WIDTH_KM
* LON0/LAT0 is the START of the profile (label X), not its centre;
* the swath is +/- HALF_WIDTH_KM perpendicular to the axis (`--width` of
  `cross_section.py`);
* the geometry is `event_maps_report.profile_geometry` (flat Earth,
  1 deg = 111 km), so the rectangle drawn is exactly the strip that
  `cross_section.py` projects.

Catalogue reading is `event_maps_report._read_catalog`, so the epicentre /
depth choices are the same flags. Defaults: PDF expectation for both
(`--position expect --column expect_z`) and no quality filter — the published
hypocentre, see CLAUDE.md §4.

Usage (pygmt_env, from the repository root)
-------------------------------------------
    conda run -n pygmt_env python complem_figures/profile_location_map.py \\
        --bulletin RESULT/SSST_result.csv \\
        --output   complem_figures/cross_section/profiles_location.png

    # other profiles / another frame
    conda run -n pygmt_env python complem_figures/profile_location_map.py \\
        --bulletin RESULT/SSST_result.csv \\
        --output   complem_figures/cross_section/profiles_location.png \\
        --profile A 0.3 42.7 45 30 5 --profile B 0.1 42.8 45 30 5 \\
        --map-region -1.2 1.8 42.2 43.45
"""

import argparse
import os
import tempfile

import pygmt as pg

from event_maps_report import (
    EventMapsReportParams,
    DEFAULT_REGION as PYRENEES_REGION,
    _dest_point,
    _read_catalog,
    profile_geometry,
)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

# name, lon0, lat0, azimuth [deg from N], length [km], half-width [km]
DEFAULT_PROFILES = [
    ('A', 0.3, 42.7, 45.0, 30.0, 5.0),
    ('B', 0.1, 42.8, 45.0, 30.0, 5.0),
]

# Central Pyrenees, ~250 x 140 km: Pau -> Andorra, Jaca -> Toulouse plain
DEFAULT_REGION = [-1.2, 1.8, 42.2, 43.45]

# Orientation towns: name, lon, lat, GMT justify of the label
TOWNS = [
    ('Pau',               -0.370, 43.295, 'BL'),
    ('Tarbes',             0.078, 43.233, 'BL'),
    ('Lourdes',           -0.046, 43.095, 'TR'),
    ('Saint-Gaudens',      0.723, 43.108, 'BL'),
    ('Foix',               1.607, 42.965, 'BL'),
    ('Jaca',              -0.549, 42.570, 'BL'),
    ('Vielha',             0.795, 42.702, 'BL'),
    ('Andorra la Vella',   1.522, 42.507, 'TR'),
]


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def generate_figure(args):
    events = _read_catalog(EventMapsReportParams(
        fileBulletin = args.bulletin,
        figSave      = args.output,
        no_filter    = not args.filter,
        depth_column = args.column,
        position     = args.position,
    ))
    region = args.map_region or list(DEFAULT_REGION)
    lon_min, lon_max, lat_min, lat_max = region
    in_frame = events[events.Longitude.between(lon_min, lon_max)
                      & events.Latitude.between(lat_min, lat_max)]
    print(f'{len(in_frame)} events inside the map frame')

    profiles = DEFAULT_PROFILES
    if args.profile:
        profiles = [(p[0], *map(float, p[1:])) for p in args.profile]

    fig = pg.Figure()
    proj = f'M{args.width_cm}c'

    # -- Main map ------------------------------------------------------------
    with pg.config(MAP_FRAME_TYPE='fancy+', FONT_ANNOT_PRIMARY='8p'):
        fig.basemap(region=region, projection=proj, frame=['WSne', 'af'])
    if args.relief:
        pg.makecpt(cmap='gray', series=[-3000, 3500], output=_tmp('relief.cpt'))
        fig.grdimage(grid=args.relief, cmap=_tmp('relief.cpt'), shading=True,
                     transparency=40)
        fig.coast(water='skyblue', resolution='i', area_thresh='0/0/1',
                  borders='1/0.75p,black')
    else:
        fig.coast(water='skyblue', land='#9A9A9A', resolution='i',
                  area_thresh='0/0/1', borders='1/0.75p,black')

    depth_cpt = _tmp('depth.cpt')
    pg.makecpt(cmap='viridis', series=[0, args.max_depth, 1], reverse=True,
               output=depth_cpt)
    fig.plot(x=in_frame.Longitude, y=in_frame.Latitude, style=f'c{args.dot_size}c',
             fill=in_frame.Depth, cmap=depth_cpt, transparency=15)

    # Towns
    if not args.no_towns:
        for name, lon, lat, just in TOWNS:
            if lon_min < lon < lon_max and lat_min < lat < lat_max:
                fig.plot(x=[lon], y=[lat], style='s0.14c', fill='white', pen='0.5p,black')
                off = '0.12c/0.08c' if just.startswith('B') else '-0.12c/-0.08c'
                fig.text(x=lon, y=lat, text=name, font='7p,Helvetica,black',
                         justify=just, offset=off, fill='white@35', clearance='1p/1p')

    # Profiles, last so they stay legible over the seismicity
    for name, lon0, lat0, az, length, half_w in profiles:
        axis, box = profile_geometry(lon0, lat0, az, length, half_w)
        fig.plot(x=box[0],  y=box[1],  pen=args.box_pen)
        fig.plot(x=axis[0], y=axis[1], pen=args.axis_pen)
        lon_s, lat_s = _dest_point(lon0, lat0, az + 180.0, args.label_gap_km)
        lon_e, lat_e = _dest_point(axis[0][1], axis[1][1], az, args.label_gap_km)
        for txt, x, y in ((name, lon_s, lat_s), (f"{name}'", lon_e, lat_e)):
            fig.text(x=x, y=y, text=txt, font='9p,Helvetica-Bold,red',
                     fill='white', pen='0.4p,red', clearance='1.5p/1.5p')
        lons, lats = box
        print(f"{name}: {lat0:.4f}N {lon0:.4f}E -> {axis[1][1]:.4f}N {axis[0][1]:.4f}E, "
              f"swath lon {min(lons):.3f}-{max(lons):.3f} lat {min(lats):.3f}-{max(lats):.3f}")

    with pg.config(FONT_ANNOT_PRIMARY='7p', FONT_LABEL='7p'):
        fig.basemap(map_scale='jBL+w25k+o0.4c/0.5c+f+lkm')
    label = 'Expectation depth' if args.column == 'expect_z' else 'Depth'
    with pg.config(FONT_ANNOT_PRIMARY='8p', FONT_LABEL='8p'):
        fig.colorbar(cmap=depth_cpt, position='JMR+o0.4c/0c+w6c/0.3c+ef',
                     frame=[f'a5f1+l{label} [km]'])

    # -- Inset: the whole range, frame of the main map in red ----------------
    if not args.no_inset:
        with fig.inset(position=f'j{args.inset_corner}+w{args.inset_width_cm}c+o0.15c',
                       box='+gwhite+p0.6p,black'):
            fig.coast(region=PYRENEES_REGION, projection=f'M{args.inset_width_cm}c',
                      water='skyblue', land='#C8C8C8', resolution='i',
                      area_thresh='0/0/1', borders='1/0.4p,black', frame='+n')
            fig.plot(x=events.Longitude, y=events.Latitude, style='c0.008c',
                     fill='gray25', transparency=30)
            fig.plot(x=[lon_min, lon_max, lon_max, lon_min, lon_min],
                     y=[lat_min, lat_min, lat_max, lat_max, lat_min],
                     pen='1p,red')

    base, _ = os.path.splitext(args.output)
    os.makedirs(os.path.dirname(os.path.abspath(base)), exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{base}.{ext}', dpi=300)
        print(f'Figure saved @ {base}.{ext}')


_TMPDIR = tempfile.mkdtemp(prefix='profile_location_map_')


def _tmp(name):
    return os.path.join(_TMPDIR, name)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description='Location map of cross-section profiles '
                                            'over the relocated seismicity.')
    p.add_argument('--bulletin', required=True, help='RESULT/*.csv or .obs bulletin')
    p.add_argument('--output', required=True,
                   help='Output path; both .png and .pdf are written')
    p.add_argument('--profile', action='append', nargs=6,
                   metavar=('NAME', 'LON0', 'LAT0', 'AZIMUTH', 'LENGTH_KM', 'HALF_WIDTH_KM'),
                   help='Profile to draw; repeatable, replaces the defaults '
                        '(A 0.3 42.7 45 30 5, B 0.1 42.8 45 30 5)')
    p.add_argument('--map-region', nargs=4, type=float,
                   metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'),
                   help=f'Main map extent (default: {" ".join(map(str, DEFAULT_REGION))})')
    p.add_argument('--column', default='expect_z', choices=['depth', 'expect_z'],
                   help='Depth column of a result CSV (default: expect_z)')
    p.add_argument('--position', default='expect', choices=['maxlike', 'expect'],
                   help='Epicentre of a result CSV (default: expect)')
    p.add_argument('--filter', action='store_true',
                   help="Apply event_maps' erh/erv/gap/rms filter (default: off)")
    p.add_argument('--max-depth', type=float, default=15.0,
                   help='Colour-bar saturation depth in km (default: 15)')
    p.add_argument('--dot-size', type=float, default=0.035,
                   help='Epicentre symbol diameter in cm (default: 0.035)')
    p.add_argument('--relief', default=None,
                   help='Shaded relief under the seismicity, e.g. @earth_relief_15s (default: none)')
    p.add_argument('--no-towns', action='store_true', help='Do not draw the orientation towns')
    p.add_argument('--no-inset', action='store_true', help='Do not draw the whole-range inset')
    p.add_argument('--inset-corner', default='TL', choices=['TL', 'TR', 'BL', 'BR'],
                   help='Corner of the inset (default: TL)')
    p.add_argument('--inset-width-cm', type=float, default=5.0)
    p.add_argument('--width-cm', type=float, default=16.0, help='Main map width')
    p.add_argument('--label-gap-km', type=float, default=2.5,
                   help='Distance of the X / X\' labels beyond the axis ends')
    p.add_argument('--axis-pen', default='1p,red')
    p.add_argument('--box-pen', default='0.6p,red,4p_2p',
                   help='Swath outline; dash in absolute points, not "--"')
    generate_figure(p.parse_args())


if __name__ == '__main__':
    main()
