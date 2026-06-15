"""
generate_regional_runfiles.py
============================
Generate all input files required for a NonLinLoc run for one geographic zone.

For a given lat/lon bounding box, produces:
  - a child .obs bulletin containing only events within the box
  - a GTSRCE station file with NLL-formatted station coordinates
  - a .in run file with velocity model, grid, and localization parameters

Usage
-----
    python NLL_run/generate_regional_runfiles.py \\
        --bulletin     obs/GLOBAL.obs \\
        --bulletin-in  obs/GLOBAL_1.obs \\
        --inventory    stations/GLOBAL_inventory.xml \\
        --code-map     stations/GLOBAL_code_map.txt \\
        --stations     stations/GTSRCE_1.txt \\
        --run-save     run/run_1.in \\
        --model        model/Pyrenees_1/Pyrenees_1 \\
        --time-file    time/Pyrenees_1/Pyrenees_1 \\
        --bulletin-out loc/GLOBAL_1/GLOBAL_1.obs \\
        --lat-min 42.5 --lat-max 43.5 --lon-min -2.0 --lon-max -0.75
"""

import argparse
import logging
import math
import os
from dataclasses import dataclass, field
from datetime import datetime

from obspy import read_inventory

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

logger = logging.getLogger('generate_regional_runfiles')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GenRunParams:
    fileBulletin:    str
    fileInventory:   str
    fileMap:         str
    fileBulletinIn:  str
    fileStations:    str
    fileRunSave:     str
    latMin_event:    float
    latMax_event:    float
    lonMin_event:    float
    lonMax_event:    float
    fileModel:       str
    fileTime:        str
    fileBulletinOut: str
    VGGRID:          list = field(default_factory=lambda: [9000, 800])


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def _setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    basename  = os.path.splitext(os.path.basename(__file__))[0]
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_path  = os.path.join(log_dir, f"{basename}_{timestamp}.log")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(log_path, encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(handler)
    return log_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_alternate_code_map(file_map):
    """
    Parse the inventory map file and return a mapping from alternate code to
    (network, station) tuple.

    Parameters
    ----------
    file_map : str — path to GLOBAL_code_map.txt

    Returns
    -------
    dict[str, (str, str)]
    """
    with open(file_map, 'r') as f:
        lines = f.readlines()

    code_map = {}
    for i, line in enumerate(lines):
        if line.startswith('Alternate'):
            alternate_code = line.split()[-1]
            j = i + 1
            while j < len(lines):
                code_line = lines[j]
                if code_line.startswith('\n'):
                    break
                if code_line.startswith('  Station'):
                    station_code = code_line.split('.')[-1].rstrip('\n')
                    network_code = code_line.split(':')[-1].split('.')[0].strip()
                j += 1
            code_map[alternate_code] = (network_code, station_code)
    return code_map


def _find_station_info(inventory, alternate_code_map, alternate_code):
    """
    Return (latitude, longitude, elevation_km) for a station identified by
    its alternate code.

    Parameters
    ----------
    inventory          : obspy.Inventory
    alternate_code_map : dict[str, (str, str)]
    alternate_code     : str

    Returns
    -------
    (float, float, float) — lat, lon, elev_km
    """
    codes   = alternate_code_map.get(alternate_code)
    station = inventory.select(network=codes[0], station=codes[1]).networks[0].stations[0]
    elev_km = (station.elevation / 1000
               if hasattr(station, 'elevation') and station.elevation is not None
               else 0.0)
    return station.latitude, station.longitude, elev_km


def _get_station_line(inventory, alternate_code_map, alternate_code):
    """Format and return a NLL GTSRCE line string for a station."""
    lat, lon, elev_km = _find_station_info(inventory, alternate_code_map, alternate_code)
    return f"GTSRCE {alternate_code} LATLON {lat:.6f} {lon:.6f} 0.0 {elev_km:.3f}\n"


def _gen_child_obs(parameters):
    """Extract events within the zone bounding box from the global bulletin."""
    with open(parameters.fileBulletin, 'r') as f:
        lines = f.readlines()

    n_eq = 0
    with open(parameters.fileBulletinIn, 'w') as f:
        f.writelines(lines[:4])
        for i, line in enumerate(lines):
            if line.startswith('# '):
                lat = float(line.split()[7])
                lon = float(line.split()[8])
                if (parameters.latMin_event <= lat <= parameters.latMax_event and
                        parameters.lonMin_event <= lon <= parameters.lonMax_event):
                    n_eq += 1
                    f.write(line)
                    j = i + 1
                    while j < len(lines) and not lines[j].startswith('\n'):
                        f.write(lines[j])
                        j += 1
                    f.write('\n')

    logger.info(f"Child bulletin: {parameters.fileBulletinIn} [{n_eq} events]")


def _gen_gtsrce(parameters):
    """Collect all station codes from the global bulletin and write GTSRCE lines."""
    with open(parameters.fileBulletin, 'r') as f:
        lines = f.readlines()

    unique_stations = {
        line.split()[0]
        for line in lines
        if not (line.startswith('\n') or line.startswith('#') or line.startswith('PUBLIC_ID'))
    }

    inventory         = read_inventory(parameters.fileInventory, format='STATIONXML')
    alternate_code_map = _build_alternate_code_map(parameters.fileMap)

    with open(parameters.fileStations, 'w') as f:
        for code in unique_stations:
            f.write(_get_station_line(inventory, alternate_code_map, code))

    logger.info(f"GTSRCE file  : {parameters.fileStations} [{len(unique_stations)} stations]")


def _verify_folders(parameters):
    """Create all parent directories required by the output file paths."""
    for path in [
        parameters.fileBulletinIn,
        parameters.fileBulletinOut,
        parameters.fileStations,
        parameters.fileRunSave,
        parameters.fileModel,
        parameters.fileTime,
    ]:
        parent = '/'.join(path.split('/')[:-1])
        os.makedirs(parent, exist_ok=True)


def _compute_grid_corners(lat_sw, lon_sw, lat_ne, lon_ne):
    """
    Compute the extended grid corner coordinates and point counts for NLL.

    Extends the zone bounding box by 100 km on each side, then converts back
    to lat/lon and computes the number of grid points at 0.05 km spacing.

    Parameters
    ----------
    lat_sw, lon_sw : float — south-west corner of the event box
    lat_ne, lon_ne : float — north-east corner of the event box

    Returns
    -------
    ((lat_sw_new, lon_sw_new), (lat_ne_new, lon_ne_new), (nx, ny, nz))
    """
    R = 6371.0

    def _latlon_to_km(lat, lon, lat0, lon0):
        lat_r  = math.radians(lat)
        lon_r  = math.radians(lon)
        lat0_r = math.radians(lat0)
        lon0_r = math.radians(lon0)
        x = R * (lon_r - lon0_r) * math.cos((lat0_r + lat_r) / 2.0)
        y = R * (lat_r - lat0_r)
        return x, y

    def _km_to_latlon(x, y, lat0, lon0):
        lat0_r = math.radians(lat0)
        lon0_r = math.radians(lon0)
        lat_r  = lat0_r + y / R
        lon_r  = lon0_r + x / (R * math.cos(lat0_r))
        return math.degrees(lat_r), math.degrees(lon_r)

    x0, y0 = _latlon_to_km(lat_sw, lon_sw, lat_sw, lon_sw)
    x1, y1 = _latlon_to_km(lat_ne, lon_ne, lat_sw, lon_sw)

    x_sw = x0 - 100
    x_ne = x1 + 100
    y_sw = y0 - 100
    y_ne = y1 + 100

    lat_new_sw, lon_new_sw = _km_to_latlon(x_sw, y_sw, lat_sw, lon_sw)
    lat_new_ne, lon_new_ne = _km_to_latlon(x_ne, y_ne, lat_sw, lon_sw)

    dx = x_ne - x_sw
    dy = y_ne - y_sw
    nx = round(dx / 0.05)
    ny = round(dy / 0.05)
    nz = 800

    logger.info(f"Min. VGGRID points: {math.sqrt(nx**2 + ny**2):.0f}")

    return (
        (round(lat_new_sw, 2), round(lon_new_sw, 2)),
        (round(lat_new_ne, 2), round(lon_new_ne, 2)),
        (nx, ny, nz),
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_run(parameters, log_dir=None):
    """
    Generate the child .obs, GTSRCE, and NLL run (.in) files for one zone.

    Parameters
    ----------
    parameters : GenRunParams
    log_dir    : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: output, log
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Zone bounds  : lat [{parameters.latMin_event}, {parameters.latMax_event}]"
                f"  lon [{parameters.lonMin_event}, {parameters.lonMax_event}]")

    _verify_folders(parameters)
    _gen_child_obs(parameters)
    _gen_gtsrce(parameters)

    (lat_sw, lon_sw), (_, _), (dx, dy, dz) = _compute_grid_corners(
        parameters.latMin_event, parameters.lonMin_event,
        parameters.latMax_event, parameters.lonMax_event,
    )

    lines = []

    lines.append('CONTROL 1 54321\n')
    lines.append(f'TRANS  LAMBERT  WGS-84  {lat_sw} {lon_sw}  42 44 0.0\n')
    lines.append('\n')

    lines.append('# Velocity model\n')
    lines.append(f'VGOUT  {parameters.fileModel}\n')
    lines.append('VGTYPE P\n')
    lines.append(f'VGGRID  2 {parameters.VGGRID[0]} {parameters.VGGRID[1]}  0.0 0.0 -3  0.05 0.05 0.05  SLOW_LEN\n')
    lines.append('\n')
    lines.append('LAYER    0.0  5.5 0.0       3.2   0.00   2.72 0.0\n')
    lines.append('LAYER    1    5.6 0.0       3.26  0.00  2.7 0.0\n')
    lines.append('LAYER    4    6.1 0.0       3.55  0.00  2.8 0.0\n')
    lines.append('LAYER    11   6.4 0.0       3.72  0.00  2.8 0.0\n')
    lines.append('LAYER    34   8.0 0.00      4.50  0.00  3.32 0.0\n')
    lines.append('\n')
    lines.append(f'GTFILES  {parameters.fileModel}  {parameters.fileTime} P\n')
    lines.append('GTMODE GRID2D ANGLES_NO\n')
    lines.append('\n')

    lines.append('# Bulletin to read and write\n')
    lines.append(f'LOCFILES {parameters.fileBulletinIn} NLLOC_OBS  {parameters.fileTime}  {parameters.fileBulletinOut}\n')
    lines.append('LOCHYPOUT SAVE_HYPO71_SUM\n')
    lines.append('\n')

    lines.append('# Localization method\n')
    lines.append(f'LOCGRID {dx} {dy} {dz} 0.0 0.0 -3  0.05 0.05 0.05 PROB_DENSITY SAVE\n')
    lines.append('LOCSEARCH  OCT 50 50 5 0.001 50000 500 1 0\n')
    lines.append('LOCMETH EDT_OT_WT 9999 4 -1 -1 1.72 6 -1.0 0\n')
    lines.append('\n')

    lines.append('# Stations coordinates\n')
    lines.append(f'INCLUDE {parameters.fileStations}\n')
    lines.append('\n')

    lines.append('# Localization parameters\n')
    lines.append('GT_PLFD  1.0e-7  0\n')
    lines.append('LOCGAU 0.05 0.0\n')
    lines.append('LOCGAU2 0.01 0.01 2.0\n')
    lines.append('LOCPHASEID  P   P p G PN PG\n')
    lines.append('LOCPHASEID  S   S s G SN SG\n')
    lines.append('LOCQUAL2ERR 0.05 0.15 0.05 0.15 99999.9\n')
    lines.append('LOCPHSTAT 9999.0 -1 9999.0 1.0 1.0 9999.9 -9999.9 9999.9\n')
    lines.append('\n')

    with open(parameters.fileRunSave, 'w', encoding='utf-8') as f:
        f.writelines(lines)

    logger.info(f"Run file     : {parameters.fileRunSave}")

    return {
        'output': parameters.fileRunSave,
        'log':    log_path,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate NLL input files for one geographic zone.'
    )
    parser.add_argument('--bulletin',     required=True, help='Global GLOBAL.obs bulletin')
    parser.add_argument('--bulletin-in',  required=True, help='Output child .obs for this zone')
    parser.add_argument('--inventory',    required=True, help='StationXML inventory file')
    parser.add_argument('--code-map',     required=True, help='GLOBAL_code_map.txt')
    parser.add_argument('--stations',     required=True, help='Output GTSRCE station file')
    parser.add_argument('--run-save',     required=True, help='Output NLL run file (.in)')
    parser.add_argument('--model',        required=True, help='Velocity model base path')
    parser.add_argument('--time-file',    required=True, help='Travel-time grid base path')
    parser.add_argument('--bulletin-out', required=True, help='NLL output loc file path')
    parser.add_argument('--lat-min',      type=float, required=True)
    parser.add_argument('--lat-max',      type=float, required=True)
    parser.add_argument('--lon-min',      type=float, required=True)
    parser.add_argument('--lon-max',      type=float, required=True)
    parser.add_argument('--vggrid',       type=int, nargs=2, default=[9000, 800],
                        help='VGGRID h/v parameters (default: 9000 800)')
    parser.add_argument('--log-dir',      default=None,
                        help='Log directory (default: NLL_run/console_output/)')
    args = parser.parse_args()

    generate_run(
        GenRunParams(
            fileBulletin    = args.bulletin,
            fileInventory   = args.inventory,
            fileMap         = args.code_map,
            fileBulletinIn  = args.bulletin_in,
            fileStations    = args.stations,
            fileRunSave     = args.run_save,
            latMin_event    = args.lat_min,
            latMax_event    = args.lat_max,
            lonMin_event    = args.lon_min,
            lonMax_event    = args.lon_max,
            fileModel       = args.model,
            fileTime        = args.time_file,
            fileBulletinOut = args.bulletin_out,
            VGGRID          = args.vggrid,
        ),
        log_dir = args.log_dir,
    )


if __name__ == '__main__':
    main()
