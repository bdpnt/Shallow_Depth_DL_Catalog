"""
filter_distant_picks.py
============================
Remove picks from stations too far from their associated event.

Reads a GLOBAL.obs bulletin, computes the epicentral distance between each
pick's station and its event, and drops picks beyond a configurable threshold.
The bulletin is overwritten in place.

Usage
-----
    python NLL_run/filter_distant_picks.py \\
        --bulletin   obs/GLOBAL.obs \\
        --inventory  stations/GLOBAL_inventory.xml \\
        --max-dist   80
"""

import argparse
import logging
import math
import os
from dataclasses import dataclass
from datetime import datetime

import pandas as pd
from obspy import read_inventory

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

logger = logging.getLogger('filter_distant_picks')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class RemoveFarPicksParams:
    fileBulletin:  str
    fileInventory: str
    maxDistance:   float  # km


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

def _haversine(lat1, lon1, lat2, lon2):
    """Return the epicentral distance in km between two lat/lon points."""
    R    = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a    = (math.sin(dlat / 2) ** 2 +
            math.cos(math.radians(lat1)) *
            math.cos(math.radians(lat2)) *
            math.sin(dlon / 2) ** 2)
    return 2 * R * math.asin(math.sqrt(a))


def _get_station_coords(inventory_path):
    """
    Build a DataFrame of station alternate codes and coordinates.

    Parameters
    ----------
    inventory_path : str — path to the StationXML inventory

    Returns
    -------
    pd.DataFrame with columns: AlternateCode, Latitude, Longitude
    """
    inv     = read_inventory(inventory_path, format='STATIONXML')
    records = []
    for net in inv.networks:
        for sta in net.stations:
            records.append({
                'AlternateCode': sta.alternate_code,
                'Latitude':      sta.latitude,
                'Longitude':     sta.longitude,
            })
    df = pd.DataFrame(records)
    return df.drop_duplicates(subset='AlternateCode', keep='first')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def remove_far_picks(parameters, log_dir=None):
    """
    Remove picks from stations beyond maxDistance from their event and
    overwrite the bulletin file.

    Parameters
    ----------
    parameters : RemoveFarPicksParams
    log_dir    : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: output, log, n_picks_total, n_picks_removed
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Bulletin    : {parameters.fileBulletin}")
    logger.info(f"Inventory   : {parameters.fileInventory}")
    logger.info(f"Max distance: {parameters.maxDistance} km")

    sta_coords = _get_station_coords(parameters.fileInventory)
    logger.info(f"  {len(sta_coords)} stations loaded.")

    with open(parameters.fileBulletin, 'r') as f:
        lines = f.readlines()

    n_picks   = 0
    remove_id = set()
    event_lat = event_lon = None

    for i, line in enumerate(lines):
        if line.startswith('# '):
            parts     = line.split()
            event_lat = float(parts[7])
            event_lon = float(parts[8])
        elif event_lat is not None and line.strip() and not line.startswith('#') and not line.startswith('PUBLIC_ID'):
            n_picks  += 1
            alt_code  = line.split()[0]
            sta_row   = sta_coords[sta_coords.AlternateCode == alt_code]
            if sta_row.empty:
                continue
            dist = _haversine(
                sta_row.Latitude.iloc[0], sta_row.Longitude.iloc[0],
                event_lat, event_lon,
            )
            if dist > parameters.maxDistance:
                remove_id.add(i)

    lines = [l for i, l in enumerate(lines) if i not in remove_id]
    with open(parameters.fileBulletin, 'w') as f:
        f.writelines(lines)

    logger.info(f"Picks total  : {n_picks}")
    logger.info(f"Picks removed: {len(remove_id)}")
    logger.info(f"Output       : {parameters.fileBulletin}")

    return {
        'output':          parameters.fileBulletin,
        'log':             log_path,
        'n_picks_total':   n_picks,
        'n_picks_removed': len(remove_id),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Remove picks from stations too far from their event.'
    )
    parser.add_argument('--bulletin',  required=True, help='GLOBAL.obs bulletin file (modified in place)')
    parser.add_argument('--inventory', required=True, help='StationXML inventory file')
    parser.add_argument('--max-dist',  type=float, default=80.0,
                        help='Maximum station-event distance in km (default: 80)')
    parser.add_argument('--log-dir',   default=None,
                        help='Log directory (default: NLL_run/console_output/)')
    args = parser.parse_args()

    remove_far_picks(
        RemoveFarPicksParams(
            fileBulletin  = args.bulletin,
            fileInventory = args.inventory,
            maxDistance   = args.max_dist,
        ),
        log_dir = args.log_dir,
    )


if __name__ == '__main__':
    main()
