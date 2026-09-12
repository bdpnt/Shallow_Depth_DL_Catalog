"""
cross_section.py
============================
Generate a planimetric map and vertical cross-section for a seismic zone.

Plots seismic events coloured by depth on a PyGMT basemap with fault traces
and station positions, then projects events onto a user-defined cross-section
line and plots the depth profile below the map.

Catalogue formats
-----------------
`--format` selects the reader, and only some of them carry per-event errors.
`--use-err` and the `--uncert-h` / `--uncert-v` filter are meaningful only where
the "errors" column below says so; elsewhere the reader substitutes dummies and
both are silently inert.

    --format  file                                          errors available
    -------------------------------------------------------------------------
    1         parse_nll_output.py output (legacy, zone_Arette/)  erh, erv
    2         RENASS bulletin                                    none
    3         RENASS-like, no gap/phase columns                  none
    4         obs/*.obs (GLOBAL, NLL_result, SSST_result)        none read (*)
    5         Chevrot CSV (#YEAR MONTH ... LAT LON DEPTH MAG)    none
    6         RESULT/NLL_result.csv                              erh, erv
    6         RESULT/SSST_result.csv                             erh, erv, psi, c68z

    (*) obs/NLL_result.obs and obs/SSST_result.obs do carry true_erh/true_erz in
        header columns 13/14, but this reader does not parse them and the
        cross-section is drawn in a flat colour. Use --format 6 on the matching
        RESULT/*.csv instead — same events, full precision.

    psi, c68z and --usable need the location-PDF columns that NLL_run/pdf_metrics.py
    adds, and it annotates SSST_result.csv only.

Usage
-----
    python complem_figures/cross_section.py \\
        --catalog    RESULT/GLOBAL_PR_W.txt \\
        --format     1 \\
        --stations   stations/GTSRCE_W.txt \\
        --output     cross_section/arette_after_erV.pdf

    python complem_figures/cross_section.py \\
        --catalog    RESULT/NLL_result.csv \\
        --format     6 \\
        --stations   stations/GTSRCE_2.txt \\
        --output     cross_section/arette_nll.pdf \\
        --use-err    erv --uncert-h 0.5 --uncert-v 0.5

    python complem_figures/cross_section.py \\
        --catalog    RESULT/SSST_result.csv \\
        --format     6 \\
        --stations   stations/GTSRCE_SSST_2.txt \\
        --output     cross_section/arette_ssst.pdf \\
        --use-err    c68z --usable

    # cross-section panel alone, no map (--stations then not needed)
    python complem_figures/cross_section.py \\
        --catalog    RESULT/SSST_result.csv \\
        --format     6 \\
        --output     cross_section/arette_ssst_section.pdf \\
        --use-err    erv --no-map
"""

import argparse
import os
from dataclasses import dataclass
from math import cos, radians, sin

import numpy as np
import pandas as pd
import pygmt
import xarray as xr
from pygmt.datasets import load_earth_relief
from scipy.ndimage import gaussian_filter

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

# Location-PDF acceptance thresholds — mirror NLL_run/export_quakeml.py, so that
# --usable selects exactly the events exported with pyr:usable="true".
_NOMINAL_COVERAGE = 0.6827   # erf(1/sqrt(2)), one sigma
_C68_Z_MIN        = -2.0


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class CrossSectionParams:
    fichier_seisme:  str
    save_file:       str
    stations_file:   str   = ''     # only read when draw_map is True
    FORMAT_fichier:  int   = 1      # see the module docstring for the format/file mapping
    use_err:         str   = 'erh'  # 'erh', 'erv', or (SSST_result.csv only) 'psi', 'c68z'
    usable_only:     bool  = False  # SSST_result.csv only: keep only pyr:usable events
    draw_map:        bool  = True   # False = cross-section panel alone, no planimetric map
    lon0:            float = -0.6275
    lat0:            float = 43.0
    azimut:          float = 0.0    # degrees from North
    longueur_coupe:  float = 16.0   # km
    largeur_coupe:   float = 8.0    # km
    prof_coupe:      float = 18.0   # km
    prof_min:        float = 0.0    # km
    prof_max:        float = 15.0   # km
    UNCERT_max_H:    float = 1.5    # km
    UNCERT_max_V:    float = 1.5    # km


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _dest_point(lon, lat, azimut, dist_km):
    """
    Compute a destination point from an origin, azimuth, and distance.

    Uses a local flat-Earth approximation (1° ≈ 111 km).

    Parameters
    ----------
    lon, lat   : float — origin coordinates
    azimut     : float — bearing in degrees from North
    dist_km    : float — distance in km

    Returns
    -------
    (float, float) — (longitude, latitude) of the destination
    """
    R   = 111.0  # km per degree
    az  = radians(azimut)
    dlat = (dist_km * cos(az)) / R
    dlon = (dist_km * sin(az)) / (R * cos(radians(lat)))
    return lon + dlon, lat + dlat


def _metric_style(name, parameters):
    """
    Colour-bar styling for the metric colouring the cross-section symbols.

    `higher_is_better` drives the draw order: the best events are plotted last,
    so they end up on top of the worst ones.

    Parameters
    ----------
    name       : str — 'erh', 'erv', 'psi' or 'c68z'
    parameters : CrossSectionParams

    Returns
    -------
    dict with keys: cmap, series, reverse, label, higher_is_better
    """
    if name == 'erv':
        return dict(cmap='magma', series=[0, parameters.UNCERT_max_V],
                    reverse=True,  label='ErV (km)',    higher_is_better=False)
    if name == 'erh':
        return dict(cmap='magma', series=[0, parameters.UNCERT_max_H],
                    reverse=True,  label='ErH (km)',    higher_is_better=False)
    if name == 'psi':
        return dict(cmap='magma', series=[0, 1],
                    reverse=False, label='Psi',         higher_is_better=True)
    return dict(cmap='polar', series=[-4, 4],
                reverse=True, label='C68 z-score', higher_is_better=True)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    """
    Generate the planimetric map and depth cross-section figure.

    Parameters
    ----------
    parameters : CrossSectionParams

    Returns
    -------
    dict with keys: output
    """
    R = 111.0  # km per degree

    # -- Cross-section geometry --
    lon0, lat0       = parameters.lon0, parameters.lat0
    lon1,  lat1      = lon0, lat0
    lon2,  lat2      = _dest_point(lon0, lat0, parameters.azimut, parameters.longueur_coupe)
    # Swath edges sit at +/- largeur_coupe perpendicular to the axis, i.e. along
    # azimut +/- 90 -- not due east-west, which only coincides for a N-S profile
    # and otherwise draws the swath narrower than the width pygmt.project selects.
    lon1a, lat1a     = _dest_point(lon0,  lat0,  parameters.azimut - 90.0, parameters.largeur_coupe)
    lon2a, lat2a     = _dest_point(lon1a, lat1a, parameters.azimut, parameters.longueur_coupe)
    lon1b, lat1b     = _dest_point(lon0,  lat0,  parameters.azimut + 90.0, parameters.largeur_coupe)
    lon2b, lat2b     = _dest_point(lon1b, lat1b, parameters.azimut, parameters.longueur_coupe)
    # Bounding box of the four swath corners, so the frame stays valid whatever
    # the azimuth; a southward profile would otherwise give ymin > ymax.
    corner_lons      = (lon1a, lon2a, lon1b, lon2b)
    corner_lats      = (lat1a, lat2a, lat1b, lat2b)
    Region           = [min(corner_lons) - 0.2, max(corner_lons) + 0.2,
                        min(corner_lats) - 0.2, max(corner_lats) + 0.2]

    # -- Load catalogue --
    quality     = {}     # extra colouring metrics, format 6 only
    usable_mask = True   # broadcasts harmlessly when no PDF-quality filter applies

    if parameters.FORMAT_fichier == 1:
        data    = np.loadtxt(parameters.fichier_seisme)
        lon     = data[:, 7]
        lat     = data[:, 6]
        depth   = data[:, 8]
        erv     = data[:, 13]
        erh     = data[:, 12]
        rms     = data[:, 10]
        gap     = data[:, 14]
        nbphase = data[:, 11]
        year    = [i + 2000 if i < 75 else i + 1900 for i in data[:, 0]]
        month   = data[:, 1]
        day     = data[:, 2]
        hour    = data[:, 3]
        minu    = data[:, 4]

    elif parameters.FORMAT_fichier == 2:
        data    = np.loadtxt(parameters.fichier_seisme)
        lon     = data[:, 7]
        lat     = data[:, 6]
        depth   = data[:, 8]
        erv     = np.zeros_like(depth) + parameters.UNCERT_max_V - 0.01
        erh     = np.zeros_like(depth) + parameters.UNCERT_max_H - 0.01
        rms     = np.zeros_like(depth)
        gap     = data[:, 13]
        nbphase = data[:, 10]
        year    = data[:, 0]
        month   = data[:, 1]
        day     = data[:, 2]
        hour    = data[:, 3]
        minu    = data[:, 4]

    elif parameters.FORMAT_fichier == 3:
        data    = np.loadtxt(parameters.fichier_seisme)
        lon     = data[:, 7]
        lat     = data[:, 6]
        depth   = data[:, 8]
        erv     = np.zeros_like(depth) + parameters.UNCERT_max_V - 0.01
        erh     = np.zeros_like(depth) + parameters.UNCERT_max_H - 0.01
        rms     = np.zeros_like(depth)
        year    = data[:, 0]
        month   = data[:, 1]
        day     = data[:, 2]
        hour    = data[:, 3]
        minu    = data[:, 4]

    elif parameters.FORMAT_fichier == 5:  # CSV: #YEAR MONTH DAY HOUR MINUTE SECOND LAT LON DEPTH MAG
        data    = np.loadtxt(parameters.fichier_seisme, comments='#')
        lon     = data[:, 7]
        lat     = data[:, 6]
        depth   = data[:, 8]
        erv     = np.full(len(depth), -1.0)
        erh     = np.full(len(depth), -1.0)
        rms     = np.full(len(depth), -1.0)
        gap     = np.full(len(depth), -1.0)
        nbphase = np.full(len(depth), -1.0)
        year    = data[:, 0]
        month   = data[:, 1]
        day     = data[:, 2]
        hour    = data[:, 3]
        minu    = data[:, 4]

    elif parameters.FORMAT_fichier == 6:  # RESULT/NLL_result.csv or RESULT/SSST_result.csv
        df      = pd.read_csv(parameters.fichier_seisme)
        lon     = df['longitude'].to_numpy()
        lat     = df['latitude'].to_numpy()
        depth   = df['depth'].to_numpy()
        erh     = df['true_erh'].to_numpy()   # 2-DOF horizontal ellipse, not errH
        erv     = df['true_erz'].to_numpy()   # 1-DOF marginal, not errZ
        rms     = df['RMS'].to_numpy()
        gap     = df['Gap'].to_numpy()
        nbphase = df['Nphs'].to_numpy()

        # The location-PDF columns are written by NLL_run/pdf_metrics.py, which
        # annotates SSST_result.csv only. NLL_result.csv stops at true_erz, so
        # psi/c68z/--usable are unavailable there while erh/erv still are.
        if 'C68' in df.columns:
            c68_z           = ((df['C68'] - _NOMINAL_COVERAGE) / df['C68_sigma_n']).to_numpy()
            quality['psi']  = df['Psi'].to_numpy()
            quality['c68z'] = c68_z

            if parameters.usable_only:
                usable_mask = (df['n_scat'].notna().to_numpy()
                               & (c68_z >= _C68_Z_MIN)
                               & ~df['dip_reject'].to_numpy(dtype=bool))

        elif parameters.usable_only:
            raise ValueError(f"--usable needs the location-PDF columns; "
                             f"{parameters.fichier_seisme} has none. Run "
                             f"NLL_run/pdf_metrics.py on it, or drop --usable.")

    else:  # FORMAT_fichier == 4 (.obs)
        with open(parameters.fichier_seisme, 'r') as f:
            lines = f.readlines()
        num_lines = sum(1 for line in lines if line.startswith('# '))
        year     = np.zeros(num_lines)
        month    = np.zeros(num_lines)
        day      = np.zeros(num_lines)
        hour     = np.zeros(num_lines)
        minu     = np.zeros(num_lines)
        lat      = np.zeros(num_lines)
        lon      = np.zeros(num_lines)
        depth    = np.zeros(num_lines)
        erv      = np.full(num_lines, -1)
        erh      = np.full(num_lines, -1)
        rms      = np.full(num_lines, -1)
        gap      = np.full(num_lines, -1)
        nbphase  = np.full(num_lines, -1)
        idx = 0
        for line in lines:
            if line.startswith('# '):
                d = line.rstrip('\n').lstrip('# ').split()
                year[idx]  = float(d[0])
                month[idx] = float(d[1])
                day[idx]   = float(d[2])
                hour[idx]  = float(d[3])
                minu[idx]  = float(d[4])
                lat[idx]   = float(d[6])
                lon[idx]   = float(d[7])
                depth[idx] = float(d[8])
                if d[12] != 'None':
                    nbphase[idx] = float(d[12])
                idx += 1

    # -- Quality filter --
    mask  = ((erv < parameters.UNCERT_max_V) & (erh < parameters.UNCERT_max_H)
             & (rms < 0.5) & usable_mask)
    lon   = lon[mask]
    lat   = lat[mask]
    depth = depth[mask]
    erv   = erv[mask]
    erh   = erh[mask]

    metrics = {'erh': erh, 'erv': erv, **{k: v[mask] for k, v in quality.items()}}
    if parameters.use_err not in metrics:
        raise ValueError(f"--use-err {parameters.use_err!r} needs a catalogue format "
                         f"providing it; available here: {sorted(metrics)}")
    cvals = metrics[parameters.use_err]

    # -- Planimetric map --
    fig = pygmt.Figure()

    if parameters.draw_map:
        fig.basemap(region=Region, projection='M6i', frame='a')
        fig.coast(shorelines=True, water='lightblue', land='lightgray', resolution='h')

        grid = load_earth_relief('03s', region=Region)
        fig.grdimage(grid=-grid, cmap='gray')

        failles_dir = os.path.join(_PROJECT_ROOT, 'failles')
        fig.plot(os.path.join(failles_dir, 'FNP.dat'),
                 pen='1.25p', style='f1c/0.25c', fill='black')
        fig.plot(os.path.join(failles_dir, 'structures_lacan.dat'),
                 pen='1.25p', style='f1c/0.25c', fill='black')
        fig.plot(os.path.join(failles_dir, 'lacan.thrust'),
                 pen='1.25p', style='f1c/0.25c', fill='blue')
        fig.plot(os.path.join(failles_dir, 'lacan.other'),
                 pen='1.25p', style='f1c/0.25c', fill='blue')

        fig.plot(x=[lon1, lon2],   y=[lat1, lat2],   pen='2p,red')
        fig.plot(x=[lon1a, lon2a], y=[lat1a, lat2a], pen='0.5p,red')
        fig.plot(x=[lon1b, lon2b], y=[lat1b, lat2b], pen='0.5p,red')

        stations = []
        with open(parameters.stations_file, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 7:
                    stations.append((parts[1], float(parts[3]), float(parts[4])))

        for stas, lats, lons in stations:
            fig.plot(x=lons, y=lats, style='t0.5c', fill='red', pen='1p,black')
            fig.text(x=lons + 0.01, y=lats + 0.005, text=stas,
                     font='10p,Helvetica-Bold', justify='LM')

        villes = {
            'Arette':               (-0.717,    43.096),
            'Sarrance':             (-0.6008333, 43.0522),
            'Oloron-Sainte-Marie':  (-0.6056,   43.1947),
        }
        for name, (lon_v, lat_v) in villes.items():
            fig.plot(x=lon_v, y=lat_v, style='s0.4c', fill='yellow', pen='1p,black')
            fig.text(x=lon_v + 0.01, y=lat_v + 0.01, text=name,
                     font='10p,Helvetica-Bold', justify='LM')

        pygmt.makecpt(cmap='viridis', series=[parameters.prof_min, parameters.prof_max], reverse=True)
        fig.plot(x=lon, y=lat, style='c0.3c', fill=depth, cmap=True, pen='black')
        fig.colorbar(frame='af+lProfondeur (km)')

    # -- Cross-section --
    cross_dir  = os.path.dirname(parameters.save_file)
    cross_file = os.path.join(cross_dir, 'cross.dat')
    data_cat   = np.column_stack((lon, lat, depth, cvals))

    pygmt.project(
        data       = data_cat,
        center     = [lon1, lat1],
        endpoint   = [lon2, lat2],
        width      = [-parameters.largeur_coupe, parameters.largeur_coupe],
        convention = 'pz',
        unit       = True,
        outfile    = cross_file,
        output_type = 'file',
    )

    if not os.path.exists(cross_file) or os.path.getsize(cross_file) == 0:
        print('No events projected onto the cross-section (cross.dat empty).')
        plot_coupe = False
    else:
        try:
            data = np.loadtxt(cross_file, ndmin=2)
        except Exception as e:
            print(f'Could not read {cross_file}: {e}')
            plot_coupe = False
        else:
            if data.size == 0:
                print('cross.dat is empty after loading.')
                plot_coupe = False
            else:
                plot_coupe = True
                X    = data[:, 0]
                Z    = data[:, 1]
                cval = data[:, 2]

    if plot_coupe:
        # 1:1 aspect. The panel keeps a fixed plotted width and takes whatever
        # height the depth range needs at the same cm/km, so distance and depth
        # are on one scale. A fixed height instead made the vertical
        # exaggeration a side effect of --length and --depth-max: 0.65x on a
        # 12 x 13 km section, 1.67x on a 50 x 21 km one.
        panel_top    = -1.0   # km, panel top (above sea level)
        panel_width  = 10.0   # cm
        panel_height = panel_width * (parameters.prof_coupe - panel_top) / parameters.longueur_coupe
        if parameters.draw_map:
            # Drop below the map, keeping the same gap whatever the panel height.
            fig.shift_origin(yshift=f'-{panel_height + 3.0}c')
        fig.basemap(
            projection = f'X{panel_width}/-{panel_height}',
            region     = [0, parameters.longueur_coupe, panel_top, parameters.prof_coupe],
            frame      = ['xafg100+lDistance (km)', 'yafg50+lDepth (km)', 'WSen'],
        )

        if parameters.FORMAT_fichier in (4, 5):
            fig.plot(x=X, y=Z, style='c0.15c', fill='#DF2B2B', pen='0.25p,black')
        else:
            style = _metric_style(parameters.use_err, parameters)
            order = np.argsort(cval)
            if not style['higher_is_better']:
                order = order[::-1]   # worst first, so the best land on top
            # Saturate at the ends of the scale: out-of-range values would
            # otherwise take the CPT background/foreground colour, which a
            # reversed CPT turns black. Psi and c68z both overshoot their range.
            shown = np.clip(cval[order], *style['series'])
            pygmt.makecpt(cmap=style['cmap'], series=style['series'],
                          reverse=style['reverse'])
            fig.plot(x=X[order], y=Z[order], style='c0.15c',
                     fill=shown, cmap=True, pen='0.25p,black')
            # Never taller than the panel: at 1:1 a long shallow section can be
            # shorter than the 5 cm the bar would otherwise take.
            fig.colorbar(frame=[f"af+l{style['label']}"],
                         position=f'JMR+w{min(5.0, panel_height)}c/0.5c+o0.5c/0c')

    if not plot_coupe and not parameters.draw_map:
        raise RuntimeError('Nothing to draw: no event projected onto the cross-section '
                           'and --no-map suppresses the only other panel.')

    fig.savefig(parameters.save_file)
    print(f'Figure saved @ {parameters.save_file}')
    return {'output': parameters.save_file}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate a planimetric map and depth cross-section for a seismic zone.'
    )
    parser.add_argument('--catalog',   required=True,
                        help='Seismicity catalogue file')
    parser.add_argument('--format',    type=int, default=1,
                        help='Catalogue format: 1=parse_nll_output output (legacy), '
                             '2/3=RENASS, 4=obs/*.obs, 5=Chevrot CSV, '
                             '6=RESULT/NLL_result.csv or RESULT/SSST_result.csv. '
                             'Only 1 and 6 carry per-event errors; see the module '
                             'docstring for the full mapping (default: 1)')
    parser.add_argument('--stations',
                        help='GTSRCE station file; required unless --no-map')
    parser.add_argument('--output',    required=True,
                        help='Output figure path (PDF or PNG)')
    parser.add_argument('--use-err',   default='erh',
                        choices=['erh', 'erv', 'psi', 'c68z'],
                        help='Metric colouring the cross-section; psi and c68z '
                             'require a CSV annotated by pdf_metrics.py, i.e. '
                             '--format 6 on SSST_result.csv (default: erh)')
    parser.add_argument('--usable',    action='store_true',
                        help='SSST_result.csv only: keep only events flagged '
                             'pyr:usable (C68 z-score >= -2, not dip_reject, '
                             'metrics present)')
    parser.add_argument('--lon0',      type=float, default=-0.6275)
    parser.add_argument('--lat0',      type=float, default=43.0)
    parser.add_argument('--azimut',    type=float, default=0.0)
    parser.add_argument('--length',    type=float, default=16.0,
                        help='Cross-section length in km (default: 16)')
    parser.add_argument('--width',     type=float, default=8.0,
                        help='Cross-section half-width in km (default: 8)')
    parser.add_argument('--depth-max', type=float, default=18.0,
                        help='Maximum depth for cross-section panel in km (default: 18)')
    parser.add_argument('--uncert-h',  type=float, default=1.5,
                        help='Maximum horizontal uncertainty filter in km (default: 1.5)')
    parser.add_argument('--uncert-v',  type=float, default=1.5,
                        help='Maximum vertical uncertainty filter in km (default: 1.5)')
    parser.add_argument('--no-map',    dest='draw_map', action='store_false',
                        help='Output the cross-section panel alone, without the '
                             'planimetric map (--stations then unused)')
    args = parser.parse_args()

    if args.draw_map and not args.stations:
        parser.error('--stations is required when the map is drawn; pass --no-map to omit it')

    generate_figure(CrossSectionParams(
        fichier_seisme = args.catalog,
        save_file      = args.output,
        stations_file  = args.stations or '',
        FORMAT_fichier = args.format,
        use_err        = args.use_err,
        usable_only    = args.usable,
        lon0           = args.lon0,
        lat0           = args.lat0,
        azimut         = args.azimut,
        longueur_coupe = args.length,
        largeur_coupe  = args.width,
        prof_coupe     = args.depth_max,
        UNCERT_max_H   = args.uncert_h,
        UNCERT_max_V   = args.uncert_v,
        draw_map       = args.draw_map,
    ))


if __name__ == '__main__':
    main()
