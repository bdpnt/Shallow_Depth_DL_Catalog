"""
_fill_missing_elevations.py
============================
Fill in missing (zero) elevations in a station CSV by querying the
EU-DEM 25 m model through the OpenTopoData API.

Usage
-----
    python fetch_inventory/_fill_missing_elevations.py \\
        --file  stations/OMP_stations_XML/ADDITIONAL_sta.csv \\
        --save  stations/OMP_stations_XML/ADDITIONAL_sta_filled.csv
"""

import argparse
import logging
import time

import pandas as pd
import requests


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_elevation(lat, lng):
    """
    Query the EU-DEM 25 m model (OpenTopoData) and return the elevation in metres.

    Returns 0 on any error (API failure, network issue, or a point outside the
    model, e.g. at sea). OpenTopoData's public API allows 1 call/s and 1 000 calls/day.
    """
    url = f"https://api.opentopodata.org/v1/eudem25m?locations={lat},{lng}"
    try:
        response = requests.get(url, timeout=10).json()
        if response.get('status') == 'OK' and response['results'][0]['elevation'] is not None:
            return int(response['results'][0]['elevation'])
        logger.warning(f"API error for ({lat}, {lng}): {response}")
        return 0
    except Exception as e:
        logger.error(f"Request failed for ({lat}, {lng}): {e}")
        return 0


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def check_elevation(file, file_save):
    """
    Fill missing elevations (value 0) in a station CSV via EU-DEM.

    Parameters
    ----------
    file      : str — input CSV path
    file_save : str — output CSV path

    Returns
    -------
    dict with key: output
    """
    df = pd.read_csv(file, dtype={'latitude': 'float64', 'longitude': 'float64', 'elevation': 'float32'})

    missing_idx = df.index[df['elevation'] == 0].tolist()
    logger.info(f"Stations with missing elevation (value 0): {len(missing_idx)}")

    n_filled = 0
    for idx in missing_idx:
        lat, lon = df.at[idx, 'latitude'], df.at[idx, 'longitude']
        elev = _get_elevation(lat, lon)
        df.at[idx, 'elevation'] = elev
        if elev != 0:
            logger.info(f"  Filled elevation for station index {idx} ({lat}, {lon}): {elev} m")
            n_filled += 1
        else:
            logger.warning(f"  Elevation remains 0 for station index {idx} ({lat}, {lon})")
        time.sleep(1)  # avoid rate limiting

    logger.info(f"Elevation fill complete: {n_filled}/{len(missing_idx)} filled successfully")
    df.to_csv(file_save, index=False)
    logger.info(f"Saved to {file_save}")
    return {'output': file_save}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Fill missing elevations in a station CSV via the EU-DEM 25 m model (OpenTopoData).'
    )
    parser.add_argument('--file', required=True, help='Input station CSV')
    parser.add_argument('--save', required=True, help='Output CSV path')
    args = parser.parse_args()

    check_elevation(args.file, args.save)


if __name__ == '__main__':
    main()
