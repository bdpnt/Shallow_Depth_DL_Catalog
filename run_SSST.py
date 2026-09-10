"""
run_SSST.py
============================
Relocate the augmented catalog with iterative SSST corrections, zone by zone.

Third relocation stage of the pipeline, run AFTER run_NLL.py and
add_temp_picks.py:
  1. run_NLL.py          -> RESULT/NLL_result.csv, obs/NLL_result.obs
  2. add_temp_picks.py   -> obs/NLL_result_augmented.obs
  3. run_SSST.py (this)  -> RESULT/SSST_result.csv, obs/SSST_result.obs

For each zone, strictly sequentially (Loc2ssst is memory-bound: zones cannot
run concurrently, but each step inside a zone runs on NLLOC_CORES /
LOC2SSST_CORES parallel processes):
  1. Cuts a zone bulletin from obs/NLL_result_augmented.obs and derives the
     SSST control files from the zone's NLL-stage run file
     (NLL_run/generate_ssst_runfiles.py).
  2. Splits the zone bulletin into per-event .nlloc_obs files
     (NLL_run/reformate_obs.py).
  3. Builds the initial P and S travel-time grids (skipped when present;
     only on a fresh start - a resume never touches them, see 5.).
  4. Runs the iterative NLLoc + Loc2ssst workflow (NLL_run/run_ssst.py):
     len(CHAR_DISTS) SSST iterations + one final NLLoc-only relocation.
  5. Deletes the zone's travel-time grids, the big disk consumers: every
     ssst_corr<i>/ grid set except the last one, which is kept with its
     symlinks materialized (resume input on a partial campaign, reusable
     final SSST model on a finished one). The initial grids
     (run/ssst_time/Pyrenees_<N>/) go too as soon as iteration 0 has run -
     once ssst_corr0 exists and is materialized, resuming at any
     iteration_start > 0 reads it (or a later ssst_corr<i>) instead and
     never touches run/ssst_time/ again.

After all zones complete (all final CSVs present):
  6. Merges the zone CSVs into RESULT/SSST_result.csv (zone-overlap
     duplicates resolved by lowest pdfVolume).
  7. Rematches relocated events to obs/NLL_result_augmented.obs via publicId
     and writes obs/SSST_result.obs.
  8. Annotates RESULT/SSST_result.csv in place with the location-PDF quality
     metrics read from the final iteration's .scat clouds
     (NLL_run/pdf_metrics.py). Non-fatal and idempotent.

Usage
-----
    python run_SSST.py
    python run_SSST.py --zones 1,2       # only these zones
    python run_SSST.py --iteration-stop 2   # partial campaign (iterations 0-2)
    python run_SSST.py --iteration-start 3  # resume (prep + grids are reused)
"""

import argparse
import glob
import os
import shutil

from run_NLL import _ZONES
from NLL_run.generate_ssst_runfiles import GenSSSTRunParams, generate_ssst_run
from NLL_run.reformate_obs          import ReformateObsParams, reformate_obs
from NLL_run.run_ssst               import SSSTRunParams, build_grids, run_ssst
from NLL_run.merge_regional_results import merge_bulletins
from NLL_run.match_pre_post_relocation import MatchCatalogsParams, save_bulletin
from NLL_run.pdf_metrics            import PdfMetricsParams, compute_metrics

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_OBS        = os.path.join(_PROJECT_ROOT, 'obs')
_NLLOC_OBS  = os.path.join(_OBS, 'nlloc_obs')
_STATIONS   = os.path.join(_PROJECT_ROOT, 'stations')
_RUN        = os.path.join(_PROJECT_ROOT, 'run')
_RUN_NLL    = os.path.join(_RUN, 'nll')
_RUN_SSST   = os.path.join(_RUN, 'ssst')
_SSST_MODEL = os.path.join(_RUN, 'ssst_model')
_SSST_TIME  = os.path.join(_RUN, 'ssst_time')
_SSST_LOC   = os.path.join(_RUN, 'ssst_loc')
_RESULT     = os.path.join(_PROJECT_ROOT, 'RESULT')

# ---------------------------------------------------------------------------
# Campaign configuration — the single place to tune an SSST run
# ---------------------------------------------------------------------------

RUN_NAME   = 'ssst_run1'            # identifier of this SSST campaign
CHAR_DISTS = [9999, 50, 15, 5, 1]   # smoothing schedule (km); first huge = static
VPVS       = -9.99                  # -9.99 = use the real S travel-time grids

NLLOC_CORES    = 4                  # parallel NLLoc chunks per iteration
LOC2SSST_CORES = 10                 # parallel Loc2ssst instances (~2 LSGRID buffers each)

# LSPHSTAT: RMSMax NRdgsMin GapMax PResidualMax SResidualMax EllLen3Max.
# NRdgsMin is also the NLLoc min-phases threshold (read back from the
# generated Loc2ssst control file, so the two selections cannot diverge).
LSPHSTAT = [0.15, 6, 200.0, 0.3, 0.5, 10.0]

MIN_PHASES_REFORMATE = 0            # min pick lines per .nlloc_obs event file


# ---------------------------------------------------------------------------
# Per-zone worker
# ---------------------------------------------------------------------------

def _final_csv(key):
    """Path of the zone's final-iteration merged summary CSV."""
    return os.path.join(_SSST_LOC, RUN_NAME, f'Pyrenees_{key}_SSST',
                        f'loc_ssst_corr{len(CHAR_DISTS)}', f'GLOBAL_{key}',
                        f'Pyrenees_{key}.sum.grid0.loc.csv')


def _cleanup_zone_grids(key, iteration_start):
    """Delete the zone's travel-time grids once they are no longer used.

    The LAST ssst_corr set is always kept — it is the resume input of the
    next iteration on a partial campaign, and the reusable final SSST model
    (e.g. to relocate extra events) on a finished one. Its symlinks
    (stations without a correction, resolving into the earlier grid sets)
    are first replaced by real copies, then the earlier ssst_corr sets are
    removed.

    The initial grids (run/ssst_time/Pyrenees_<key>/) are only read by a
    fresh start (iteration_start == 0): build_grids() builds them and
    iteration 0 locates with them, after which ssst_corr0 exists fully
    materialized and every later iteration/resume reads from an ssst_corr<i>
    directory instead. So they are deleted right after a fresh-start call,
    regardless of whether the campaign finished; a resume call
    (iteration_start > 0) never builds or needs them, so there is nothing to
    delete on those calls.
    """
    corr_dirs = sorted(
        glob.glob(os.path.join(_SSST_LOC, RUN_NAME, f'Pyrenees_{key}_SSST',
                               'ssst_corr*')),
        key=lambda path: int(path.rsplit('ssst_corr', 1)[1]))

    if corr_dirs:
        kept = corr_dirs[-1]
        for root, _, files in os.walk(kept):
            for name in files:
                path = os.path.join(root, name)
                if os.path.islink(path):
                    target = os.path.realpath(path)
                    os.remove(path)
                    shutil.copy2(target, path)
        for corr_dir in corr_dirs[:-1]:
            shutil.rmtree(corr_dir)

    if iteration_start == 0:
        shutil.rmtree(os.path.join(_SSST_TIME, f'Pyrenees_{key}'),
                      ignore_errors=True)


def _process_zone(key, item, iteration_start, iteration_stop, rebuild_grids):
    zone_label = f'zone {key}'
    (lat_min, lon_min), (lat_max, lon_max) = item

    # --- zone prep: bulletin, GTSRCE, control files, per-event obs files.
    # Only on a fresh start — regenerating them mid-campaign could change
    # the inputs between iterations of the same run.
    if iteration_start == 0:
        generate_ssst_run(GenSSSTRunParams(
            fileBulletin    = os.path.join(_OBS, 'NLL_result_augmented.obs'),
            fileInventory   = os.path.join(_STATIONS, 'GLOBAL_inventory.xml'),
            fileMap         = os.path.join(_STATIONS, 'GLOBAL_code_map.txt'),
            fileBulletinIn  = os.path.join(_OBS, f'GLOBAL_{key}_SSST.obs'),
            fileStations    = os.path.join(_STATIONS, f'GTSRCE_SSST_{key}.txt'),
            fileRunNLL      = os.path.join(_RUN_NLL,  f'run_{key}_DELAYS.in'),
            fileRunSaveNLL  = os.path.join(_RUN_SSST, f'run_{key}_NLL.in'),
            fileRunSaveSSST = os.path.join(_RUN_SSST, f'run_{key}_SSST.in'),
            fileModel       = os.path.join(_SSST_MODEL, f'Pyrenees_{key}', f'Pyrenees_{key}'),
            fileTime        = os.path.join(_SSST_TIME,  f'Pyrenees_{key}', f'Pyrenees_{key}'),
            latMin_event    = lat_min,
            latMax_event    = lat_max,
            lonMin_event    = lon_min,
            lonMax_event    = lon_max,
            lsphstat        = LSPHSTAT,
        ))

        reformate_obs(ReformateObsParams(
            fileBulletin = os.path.join(_OBS, f'GLOBAL_{key}_SSST.obs'),
            outputDir    = os.path.join(_NLLOC_OBS, f'GLOBAL_{key}'),
            latMin       = lat_min,
            latMax       = lat_max,
            lonMin       = lon_min,
            lonMax       = lon_max,
            minPhases    = MIN_PHASES_REFORMATE,
        ))

    params = SSSTRunParams(
        runName      = RUN_NAME,
        projectName  = f'Pyrenees_{key}',
        model        = f'Pyrenees_{key}',
        catalog      = f'GLOBAL_{key}',
        fileRunNLL   = os.path.join(_RUN_SSST, f'run_{key}_NLL.in'),
        fileRunSSST  = os.path.join(_RUN_SSST, f'run_{key}_SSST.in'),
        locObs       = os.path.join(_NLLOC_OBS, f'GLOBAL_{key}', '*.nlloc_obs'),
        stationsFile = os.path.join(_STATIONS, f'GTSRCE_SSST_{key}.txt'),
        initTimeRoot = os.path.join(_SSST_TIME, f'Pyrenees_{key}', f'Pyrenees_{key}'),
        outRoot      = _SSST_LOC,
        tmpDir       = os.path.join(_RUN_SSST, f'tmp_{key}'),
        logFile      = os.path.join(_RUN_SSST, 'log', f'run_ssst_{key}.log'),
        ssstListFile = os.path.join(_RUN_SSST, 'ssst.list'),
        charDists    = CHAR_DISTS,
        vpvs         = VPVS,
        vpvsIter     = VPVS,
        nllocCores    = NLLOC_CORES,
        loc2ssstCores = LOC2SSST_CORES,
        zoneLabel     = zone_label,
    )

    # initial P+S grids (skipped when already present); only needed on a
    # fresh start - a resume relocates from the previous iteration's
    # corrected grids instead and never touches run/ssst_time/
    if iteration_start == 0:
        build_grids(params, rebuild=rebuild_grids)

    result = run_ssst(params, iteration_start=iteration_start,
                      iteration_stop=iteration_stop)

    # free the zone's travel-time grids before the next zone starts
    _cleanup_zone_grids(key, iteration_start)

    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline(zones=None, iteration_start=0, iteration_stop=None,
                 rebuild_grids=False):
    """Run the SSST stage for the given zones (default: all, sequentially)."""
    bulletin = os.path.join(_OBS, 'NLL_result_augmented.obs')
    if not os.path.isfile(bulletin):
        raise FileNotFoundError(
            f'{bulletin} not found — run run_NLL.py then add_temp_picks.py first')

    zone_keys = zones or list(_ZONES)
    for key in zone_keys:
        _process_zone(key, _ZONES[key], iteration_start, iteration_stop,
                      rebuild_grids)

    # --- finalize only once every zone's final relocation exists (partial
    # campaigns / per-zone runs skip this until the last piece is done)
    csv_files = [_final_csv(key) for key in _ZONES]
    missing = [path for path in csv_files if not os.path.isfile(path)]
    if missing:
        print('Final CSVs still missing for some zones - skipping the merge:')
        for path in missing:
            print(f'  {os.path.relpath(path, _PROJECT_ROOT)}')
        return

    merge_bulletins(csv_files, os.path.join(_RESULT, 'SSST_result.csv'))

    save_bulletin(MatchCatalogsParams(
        file_obs   = bulletin,
        file_final = os.path.join(_RESULT, 'SSST_result.csv'),
        save_file  = os.path.join(_OBS, 'SSST_result.obs'),
    ))

    # --- location-PDF quality metrics, appended in place to the merged CSV.
    # Kept non-fatal: the campaign's results are already written at this
    # point, and a metrics failure must not make a multi-hour run look like
    # it failed. The step is idempotent, so it can simply be rerun.
    try:
        compute_metrics(PdfMetricsParams(
            result_csv = os.path.join(_RESULT, 'SSST_result.csv'),
            ssst_root  = _SSST_LOC,
            run_name   = RUN_NAME,
            zones      = list(_ZONES),
        ))
    except Exception as exc:
        print(f'PDF metrics failed ({exc.__class__.__name__}: {exc})\n'
              f'  RESULT/SSST_result.csv and obs/SSST_result.obs are complete, '
              f'only the quality columns are missing. Rerun with:\n'
              f'    python NLL_run/pdf_metrics.py --run-name {RUN_NAME}')


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='SSST relocation of the augmented catalog, zone by zone.')
    parser.add_argument('--zones', default=None,
                        help='Comma-separated zone keys (default: all), e.g. 1,2')
    parser.add_argument('--iteration-start', type=int, default=0,
                        help='First iteration to run (default 0); > 0 resumes '
                             'from the previous iteration\'s SSST grids')
    parser.add_argument('--iteration-stop', type=int, default=None,
                        help='Last iteration to run, inclusive (default: the '
                             'final NLLoc-only relocation)')
    parser.add_argument('--rebuild-grids', action='store_true',
                        help='Rebuild the P/S grids even when they exist')
    args = parser.parse_args()

    zones = args.zones.split(',') if args.zones else None
    if zones:
        unknown = [z for z in zones if z not in _ZONES]
        if unknown:
            parser.error(f'unknown zone(s): {", ".join(unknown)}')

    run_pipeline(zones=zones,
                 iteration_start=args.iteration_start,
                 iteration_stop=args.iteration_stop,
                 rebuild_grids=args.rebuild_grids)


if __name__ == '__main__':
    main()
