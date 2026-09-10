"""
final_steps.py
============================
Orchestrate the post-SSST finalization stage.

Runs after run_SSST.py, on the artefacts it leaves behind. Each of those holds
half the catalog: obs/SSST_result.obs carries the magnitudes and the phase picks,
RESULT/SSST_result.csv the full-precision hypocentres, the confidence ellipsoid
and the location-PDF quality metrics. This stage merges them on publicId into a
single standard-format bulletin.

Pipeline
--------
1. Export RESULT/FINAL.xml — a QuakeML 1.2 bulletin holding every event with its
   picks, and a boolean usability flag plus the quality metrics behind it.
   Two companions are written alongside it, so that reading the catalog never
   requires loading all of it at once: FINAL_catalog.xml (same events without
   picks) and FINAL_<from>_<to>.xml (the complete bulletin cut into 5-year
   calendar periods).

Usage
-----
    python final_steps.py
"""

import os
import sys

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, _PROJECT_ROOT)

from NLL_run.export_quakeml import ExportQuakeMLParams, write_quakeml

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_OBS       = os.path.join(_PROJECT_ROOT, 'obs',      'SSST_result.obs')
_CSV       = os.path.join(_PROJECT_ROOT, 'RESULT',   'SSST_result.csv')
_INVENTORY = os.path.join(_PROJECT_ROOT, 'stations', 'GLOBAL_inventory.xml')
_FINAL_XML = os.path.join(_PROJECT_ROOT, 'RESULT',   'FINAL.xml')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """
    Run the finalization stage end-to-end.

    Returns
    -------
    dict
        Summary with key 'output' (path to the final QuakeML bulletin).
    """
    for path in (_OBS, _CSV, _INVENTORY):
        if not os.path.exists(path):
            raise FileNotFoundError(f"Required input missing: {path}")

    print("[1/1] Exporting final QuakeML bulletin ...")
    print(f"  {os.path.basename(_OBS)} + {os.path.basename(_CSV)} → {os.path.basename(_FINAL_XML)}")
    summary = write_quakeml(ExportQuakeMLParams(
        file_obs       = _OBS,
        file_csv       = _CSV,
        file_inventory = _INVENTORY,
        save_file      = _FINAL_XML,
    ))

    n_unusable = summary['n_events'] - summary['n_usable']
    print(f"\n  Events   : {summary['n_events']}")
    print(f"  Picks    : {summary['n_picks']}")
    print(f"  Usable   : {summary['n_usable']}")
    print(f"  Unusable : {n_unusable}")
    for reason, count in sorted(summary['reasons'].items()):
        if reason != 'usable':
            print(f"    {reason:<20s} : {count}")
    if summary['unresolved']:
        print(f"  Unresolved station codes: {len(summary['unresolved'])} "
              f"(see {summary['log']})")

    print("\n  Files:")
    for entry in summary['files']:
        print(f"    {os.path.basename(entry['path']):<26s} "
              f"{entry['n_events']:6d} events {entry['size']/1e6:8.1f} MB")

    print(f"\nDone. Final bulletin: {_FINAL_XML}")
    return {'output': _FINAL_XML}


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    run_pipeline()


if __name__ == '__main__':
    main()
