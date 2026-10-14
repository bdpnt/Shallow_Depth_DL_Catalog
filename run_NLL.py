"""
run_NLL.py
============================
Generate NLL inputs, run both passes, and finalize the relocated catalog for
all 6 geographic zones.

For each zone (processed with up to 3 zones running concurrently, since zones
are independent and NLL runs are the bottleneck):
  1. Filters obs/GLOBAL.obs to remove picks from stations farther than 80 km
     (done once, globally, before any zone starts).
  2. Generates a zone-specific .obs file, GTSRCE station file, and first-pass
     NLL run file; runs Vel2Grid -> Grid2Time -> NLLoc.
  3. Cleans up .hdr files left by the first-pass NLLoc run.
  4. Derives per-station delay corrections from first-pass residuals,
     generates a second-pass run file, and reruns NLLoc.
  5. Cleans up .hdr files left by the second-pass NLLoc run, and deletes the
     zone's travel-time grids (run/nll_time/Pyrenees_<N>/) — nothing reads
     them after the second pass, and they dominate disk usage.

After all zones complete:
  6. Exports a summary of locdelay corrections.
  7. Merges the 6 zone CSV outputs into RESULT/NLL_result.csv, resolving
     zone-overlap duplicates by keeping the solution with the smallest
     pdfVolume.
  8. Rematches relocated events back to obs/GLOBAL.obs via publicId to
     recover magnitude and picks, and writes obs/NLL_result.obs.

Usage
-----
    python run_NLL.py              # keeps the machine awake (caffeinate)
    python run_NLL.py --decaff     # lets it idle-sleep
"""

import argparse
import glob
import os
import shutil
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

from NLL_run.filter_distant_picks       import RemoveFarPicksParams
from NLL_run.generate_regional_runfiles import GenRunParams
from NLL_run.append_station_delays      import SecondRunParams, append_station_delays
from NLL_run.export_locdelay_info       import export_locdelay_info
from NLL_run.keep_awake                 import keep_awake
from NLL_run.match_pre_post_relocation  import MatchCatalogsParams, save_bulletin
from NLL_run.merge_regional_results     import merge_bulletins
from NLL_run.run_zone                   import run_zone
import NLL_run

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_OBS      = os.path.join(_PROJECT_ROOT, 'obs')
_STATIONS = os.path.join(_PROJECT_ROOT, 'stations')
_RUN      = os.path.join(_PROJECT_ROOT, 'run')
_RUN_NLL  = os.path.join(_RUN, 'nll')
_MODEL    = os.path.join(_RUN, 'nll_model')
_TIME     = os.path.join(_RUN, 'nll_time')
_LOC      = os.path.join(_RUN, 'nll_loc')
_RESULT   = os.path.join(_PROJECT_ROOT, 'RESULT')

_ZONES = {
    "1": ((42.50, -2.00), (43.50, -0.75)),
    "2": ((42.50, -1.00), (43.25,  0.50)),
    "3": ((42.00,  0.25), (43.25,  1.00)),
    "4": ((42.00,  0.75), (43.00,  2.25)),
    "5": ((42.00,  2.00), (43.00,  3.50)),
    "6": ((42.75,  2.25), (43.75,  3.50)),
}

_MAX_PARALLEL_ZONES = 3

# generate_run() and append_station_delays() both reconfigure a shared
# module-level logger's handlers on every call, which is not safe if two
# zones call them at the same time. Serializing just these fast,
# non-NLL calls keeps the slow run_zone() subprocesses fully parallel.
_logger_lock = threading.Lock()


# ---------------------------------------------------------------------------
# Per-zone worker
# ---------------------------------------------------------------------------

def _clean_hdr(key):
    for hdr in glob.glob(os.path.join(_LOC, f'GLOBAL_{key}', '*.hdr')):
        os.remove(hdr)


def _clean_time(key):
    shutil.rmtree(os.path.join(_TIME, f'Pyrenees_{key}'), ignore_errors=True)


def _process_zone(key, item):
    zone_label = f"zone {key}"

    # --- First pass ---
    params_run = GenRunParams(
        fileBulletin   = os.path.join(_OBS, 'GLOBAL.obs'),
        fileInventory  = os.path.join(_STATIONS, 'GLOBAL_inventory.xml'),
        fileMap        = os.path.join(_STATIONS, 'GLOBAL_code_map.txt'),
        fileBulletinIn = os.path.join(_OBS,     f'GLOBAL_{key}.obs'),
        fileStations   = os.path.join(_STATIONS, f'GTSRCE_{key}.txt'),
        fileRunSave    = os.path.join(_RUN_NLL,  f'run_{key}_DELAYS.in'),
        latMin_event   = item[0][0],
        latMax_event   = item[1][0],
        lonMin_event   = item[0][1],
        lonMax_event   = item[1][1],
        fileModel      = os.path.join(_MODEL,    f'Pyrenees_{key}', f'Pyrenees_{key}'),
        fileTime       = os.path.join(_TIME,     f'Pyrenees_{key}', f'Pyrenees_{key}'),
        fileBulletinOut = os.path.join(_LOC,     f'GLOBAL_{key}',   f'GLOBAL_{key}.obs'),
        VGGRID         = [9000, 800],
    )

    with _logger_lock:
        NLL_run.generate_regional_runfiles.generate_run(params_run)

    run_zone(os.path.join(_RUN_NLL, f'run_{key}_DELAYS.in'), zone_label=zone_label)

    # Free disk space before the second pass
    _clean_hdr(key)

    # --- Corrections (second) pass ---
    params_ssst_W = SecondRunParams(
        locFolderName = os.path.join(_LOC, f'GLOBAL_{key}'),
        fileRunName   = os.path.join(_RUN_NLL, f'run_{key}_DELAYS.in'),
        fileRunSave   = os.path.join(_RUN_NLL, f'run_{key}_NLL.in'),
        minPhases     = 100,  # minimal number of phases for the delay to be used
    )

    with _logger_lock:
        append_station_delays(params_ssst_W)

    run_zone(os.path.join(_RUN_NLL, f'run_{key}_NLL.in'), corrections_pass=True, zone_label=zone_label)

    # Free disk space right after this zone's final run; the travel-time
    # grids are no longer read by anything once the second pass is done
    _clean_hdr(key)
    _clean_time(key)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Generate NLL inputs, run both passes for all zones, and finalize the catalog."""
    # Remove far picks (once, applies to the full catalog)
    params_farpicks = RemoveFarPicksParams(
        fileBulletin  = os.path.join(_OBS, 'GLOBAL.obs'),
        fileInventory = os.path.join(_STATIONS, 'GLOBAL_inventory.xml'),
        maxDistance   = 80,  # max distance between event and station, in kilometers
    )
    NLL_run.filter_distant_picks.remove_far_picks(params_farpicks)

    # Process zones, up to _MAX_PARALLEL_ZONES at a time
    with ThreadPoolExecutor(max_workers=_MAX_PARALLEL_ZONES) as executor:
        futures = {executor.submit(_process_zone, key, item): key for key, item in _ZONES.items()}
        for future in as_completed(futures):
            future.result()

    # Export the locdelays (needs all zones' corrected run files)
    export_locdelay_info(
        run_dir      = _RUN_NLL,
        codemap_path = os.path.join(_STATIONS, 'GLOBAL_code_map.txt'),
        output_path  = os.path.join(_RUN_NLL, 'locdelays', 'locdelay_summary.txt'),
    )

    # Merge all zone CSV outputs into RESULT/NLL_result.csv
    csv_files = [
        os.path.join(_LOC, f'GLOBAL_{key}', f'GLOBAL_{key}.obs.sum.grid0.loc.csv')
        for key in _ZONES
    ]
    merge_bulletins(csv_files, os.path.join(_RESULT, 'NLL_result.csv'))

    # Rematch to obs/GLOBAL.obs via publicId and write obs/NLL_result.obs
    save_bulletin(MatchCatalogsParams(
        file_obs   = os.path.join(_OBS, 'GLOBAL.obs'),
        file_final = os.path.join(_RESULT, 'NLL_result.csv'),
        save_file  = os.path.join(_OBS, 'NLL_result.obs'),
    ))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='NLL relocation of obs/GLOBAL.obs over the 6 zones.')
    parser.add_argument('--decaff', action='store_true',
                        help='Do not keep the machine awake with caffeinate')
    args = parser.parse_args()

    if not args.decaff:
        keep_awake()

    run_pipeline()


if __name__ == '__main__':
    main()
