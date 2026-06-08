"""
cross_section.py
============================
Generate a planimetric map and vertical cross-section for a seismic zone.

Plots seismic events coloured by depth on a PyGMT basemap with fault traces
and station positions, then projects events onto a user-defined cross-section
line and plots the depth profile below the map.

Usage
-----
    python complem_figures/cross_section.py \\
        --catalog    RESULT/GLOBAL_PR_W.txt \\
        --format     1 \\
        --stations   stations/GTSRCE_W.txt \\
        --output     cross_section/arette_after_erV.pdf
"""

import argparse
import os
from dataclasses import dataclass
from math import cos, radians, sin

import numpy as np
import pygmt
import xarray as xr
from pygmt.datasets import load_earth_relief
from scipy.ndimage import gaussian_filter

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class CrossSectionParams:
    fichier_seisme:  str
    save_file:       str
    stations_file:   str
    FORMAT_fichier:  int   = 1      # 1 = NLL output, 2 = RENASS bulletin, 4 = .obs, 5 = CSV (Chevrot)
    use_err:         str   = 'erh'  # 'erh' or 'erv'
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
    lon1a, lat1a     = lon0 - parameters.largeur_coupe / (R * cos(radians(lat0))), lat0
    lon2a, lat2a     = _dest_point(lon1a, lat1a, parameters.azimut, parameters.longueur_coupe)
    lon1b, lat1b     = lon0 + parameters.largeur_coupe / (R * cos(radians(lat0))), lat0
    lon2b, lat2b     = _dest_point(lon1b, lat1b, parameters.azimut, parameters.longueur_coupe)
    Region           = [lon1a - 0.2, lon2b + 0.2, lat1a - 0.2, lat2b + 0.2]

    # -- Load catalogue --
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
    mask  = (erv < parameters.UNCERT_max_V) & (erh < parameters.UNCERT_max_H) & (rms < 0.5)
    lon   = lon[mask]
    lat   = lat[mask]
    depth = depth[mask]
    erv   = erv[mask]
    erh   = erh[mask]

    # -- Planimetric map --
    fig = pygmt.Figure()
    fig.basemap(region=Region, projection='M6i', frame='a')
    fig.coast(shorelines=True, water='lightblue', land='lightgray', resolution='h')

    grid = load_earth_relief('03s', region=Region)
    fig.grdimage(grid=-grid, cmap='gray')

    failles_dir = os.path.join(_PROJECT_ROOT, 'FAILLES')
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
    data_cat   = np.column_stack((lon, lat, depth, erv, erh))

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
                X   = data[:, 0]
                Z   = data[:, 1]
                erv = data[:, 2]
                erh = data[:, 3]

    if plot_coupe:
        if parameters.FORMAT_fichier not in (4, 5):
            if parameters.use_err == 'erv':
                err_col     = erv
                sorted_data = np.column_stack((X, Z, erv))[np.argsort(erv)][::-1]
            else:
                err_col     = erh
                sorted_data = np.column_stack((X, Z, erh))[np.argsort(erh)][::-1]
            X_sorted, Z_sorted, err_sorted = sorted_data.T

        fig.shift_origin(yshift='-10c')
        fig.basemap(
            projection = 'X10/-7',
            region     = [0, parameters.longueur_coupe, -1, parameters.prof_coupe],
            frame      = ['xafg100+lDistance (km)', 'yafg50+lDepth (km)', 'WSen'],
        )

        if parameters.FORMAT_fichier not in (4, 5):
            if parameters.use_err == 'erv':
                pygmt.makecpt(cmap='magma', series=[0, parameters.UNCERT_max_V], reverse=True)
            else:
                pygmt.makecpt(cmap='magma', series=[0, parameters.UNCERT_max_H], reverse=True)
            fig.plot(x=X_sorted, y=Z_sorted, style='c0.15c',
                     fill=err_sorted, cmap=True, pen='0.25p,black')
            label = 'af+lErV' if parameters.use_err == 'erv' else 'af+lErH'
            fig.colorbar(frame=[label], position='JMR+w5c/0.5c+o0.5c/0c')
        else:
            fig.plot(x=X, y=Z, style='c0.15c', fill='#DF2B2B', pen='0.25p,black')

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
                        help='Catalogue format: 1=NLL, 2=RENASS, 4=obs, 5=CSV (default: 1)')
    parser.add_argument('--stations',  required=True,
                        help='GTSRCE station file')
    parser.add_argument('--output',    required=True,
                        help='Output figure path (PDF or PNG)')
    parser.add_argument('--use-err',   default='erh', choices=['erh', 'erv'],
                        help='Error metric to colour the cross-section (default: erh)')
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
    args = parser.parse_args()

    generate_figure(CrossSectionParams(
        fichier_seisme = args.catalog,
        save_file      = args.output,
        stations_file  = args.stations,
        FORMAT_fichier = args.format,
        use_err        = args.use_err,
        lon0           = args.lon0,
        lat0           = args.lat0,
        azimut         = args.azimut,
        longueur_coupe = args.length,
        largeur_coupe  = args.width,
        prof_coupe     = args.depth_max,
        UNCERT_max_H   = args.uncert_h,
        UNCERT_max_V   = args.uncert_v,
    ))


if __name__ == '__main__':
    main()
