"""
fetch_all_bulletins.py
============================
Fetch and convert all source seismic bulletins to .obs format.

Downloads or reads catalogs from RESIF, ICGC, IGN, LDG, and OMP, then
writes each to an individual .obs file in obs/.

Usage
-----
    python fetch_all_bulletins.py
"""

import os

from fetch_obs.ICGC  import ICGCParams
from fetch_obs.IGN   import IGNParams
from fetch_obs.LDG   import LDGParams
from fetch_obs.OMP   import OMPParams
from fetch_obs.RESIF import RESIFParams
import fetch_obs
from obspy import UTCDateTime

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_ORG = os.path.join(_PROJECT_ROOT, 'org_catalogs')
_OBS = os.path.join(_PROJECT_ROOT, 'obs')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Fetch all source bulletins and write them to obs/ as .obs files."""
    # RESIF
    params_resif = RESIFParams(
        file_name   = os.path.join(_ORG, 'RESIF_20-25.xml'),
        save_name   = os.path.join(_OBS, 'RESIF_20-25.obs'),
        client_name = 'RESIF',
        t1          = UTCDateTime('2020-01-01T00:00:00'),
        t2          = UTCDateTime('2026-01-01T00:00:00'),
        lat_min     = 41,
        lat_max     = 44,
        lon_min     = -3,
        lon_max     = 4,
        mag_min     = -5,
        event_type  = 'earthquake',
        mag_type    = 'MLv',
    )

    fetch_obs.RESIF.generate_catalog(params_resif)
    fetch_obs.RESIF.write_catalog_to_obs(params_resif)

    # ICGC
    params_icgc = ICGCParams(
        file_name   = os.path.join(_ORG, 'ICGC_20-25.txt'),
        code_name   = os.path.join(_ORG, 'CODES_ICGC_20-25.txt'),
        error_name  = os.path.join(_ORG, 'ERR_ICGC_20-25.txt'),
        save_name   = os.path.join(_OBS, 'ICGC_20-25.obs'),
        start_year  = 2020,
        start_month = 1,
        end_year    = 2025,
        end_month   = 12,
        mag_min     = 0,
    )

    fetch_obs.ICGC.get_all_codes(params_icgc)
    fetch_obs.ICGC.fetch_catalog(params_icgc)
    fetch_obs.ICGC.write_catalog_to_obs(params_icgc)

    # IGN
    params_ign = IGNParams(
        file_name = os.path.join(_ORG, 'IGN_20-25.txt'),
        save_name = os.path.join(_OBS, 'IGN_20-25.obs'),
    )

    fetch_obs.IGN.write_catalog_to_obs(params_ign)

    # LDG
    params_ldg = LDGParams(
        catalog_file = os.path.join(_ORG, 'LDG_20-25_catalog.txt'),
        arrival_file = os.path.join(_ORG, 'LDG_20-25_arrivals.txt'),
        save_name    = os.path.join(_OBS, 'LDG_20-25.obs'),
    )

    fetch_obs.LDG.write_catalog_to_obs(params_ldg)

    # OMP
    params_omp_1978 = OMPParams(
        file_name = os.path.join(_ORG, 'OMP_78-19.mag'),
        save_name = os.path.join(_OBS, 'OMP_78-19.obs'),
    )

    fetch_obs.OMP.write_catalog_to_obs(params_omp_1978)

    params_omp_2016 = OMPParams(
        file_name = os.path.join(_ORG, 'OMP_2016.mag'),
        save_name = os.path.join(_OBS, 'OMP_2016.obs'),
    )

    fetch_obs.OMP.write_catalog_to_obs(params_omp_2016)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    run_pipeline()


if __name__ == '__main__':
    main()
