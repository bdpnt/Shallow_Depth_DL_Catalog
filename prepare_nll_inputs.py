"""
prepare_nll_inputs.py
============================
Generate NonLinLoc input files for all 6 geographic zones.

For each zone:
  1. Filters obs/GLOBAL.obs to remove picks from stations farther than 80 km.
  2. Generates a zone-specific .obs file, GTSRCE station file, and NLL run file.

NLL runs (Vel2Grid → Grid2Time → NLLoc) are launched automatically after
each zone's run file is generated.

Usage
-----
    python prepare_nll_inputs.py
"""

import os

from NLL_run.filter_distant_picks       import RemoveFarPicksParams
from NLL_run.generate_regional_runfiles import GenRunParams
from NLL_run.run_zone                   import run_zone
import NLL_run

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_OBS      = os.path.join(_PROJECT_ROOT, 'obs')
_STATIONS = os.path.join(_PROJECT_ROOT, 'stations')
_RUN      = os.path.join(_PROJECT_ROOT, 'run')
_MODEL    = os.path.join(_PROJECT_ROOT, 'model')
_TIME     = os.path.join(_PROJECT_ROOT, 'time')
_LOC      = os.path.join(_PROJECT_ROOT, 'loc')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Generate NLL run files and zone bulletins for all 6 geographic zones."""
    # Remove far picks
    params_farpicks = RemoveFarPicksParams(
        fileBulletin  = os.path.join(_OBS, 'GLOBAL.obs'),
        fileInventory = os.path.join(_STATIONS, 'GLOBAL_inventory.xml'),
        maxDistance   = 80,  # max distance between event and station, in kilometers
    )

    NLL_run.filter_distant_picks.remove_far_picks(params_farpicks)

    # Generate run files
    all_runs = {
        "1": ((42.50, -2.00), (43.50, -0.75)),
        "2": ((42.50, -1.00), (43.25,  0.50)),
        "3": ((42.00,  0.25), (43.25,  1.00)),
        "4": ((42.00,  0.75), (43.00,  2.25)),
        "5": ((42.00,  2.00), (43.00,  3.50)),
        "6": ((42.75,  2.25), (43.75,  3.50)),
    }

    for key, item in all_runs.items():
        params_run = GenRunParams(
            fileBulletin   = os.path.join(_OBS, 'GLOBAL.obs'),
            fileInventory  = os.path.join(_STATIONS, 'GLOBAL_inventory.xml'),
            fileMap        = os.path.join(_STATIONS, 'GLOBAL_code_map.txt'),
            fileBulletinIn = os.path.join(_OBS,     f'GLOBAL_{key}.obs'),
            fileStations   = os.path.join(_STATIONS, f'GTSRCE_{key}.txt'),
            fileRunSave    = os.path.join(_RUN,      f'run_{key}.in'),
            latMin_event   = item[0][0],
            latMax_event   = item[1][0],
            lonMin_event   = item[0][1],
            lonMax_event   = item[1][1],
            fileModel      = os.path.join(_MODEL,    f'Pyrenees_{key}', f'Pyrenees_{key}'),
            fileTime       = os.path.join(_TIME,     f'Pyrenees_{key}', f'Pyrenees_{key}'),
            fileBulletinOut = os.path.join(_LOC,     f'GLOBAL_{key}',   f'GLOBAL_{key}.obs'),
            VGGRID         = [9000, 800],
        )

        NLL_run.generate_regional_runfiles.generate_run(params_run)

        run_zone(os.path.join(_RUN, f'run_{key}.in'))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    run_pipeline()


if __name__ == '__main__':
    main()
