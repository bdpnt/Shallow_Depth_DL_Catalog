"""
generate_nll_corrections.py
============================
Generate second-pass NLL run files with per-station delay corrections.

For each of the 6 geographic zones:
  1. Removes .hdr files left by NLL in the loc/GLOBAL_<key>/ folder.
  2. Derives per-station delay corrections from first-pass residuals and appends
     them to a second-pass run file (run_<key>_PR.in).

Also exports a summary of all locdelay corrections to run/locdelays/.

NLLoc is launched automatically for each zone after its corrections run file
is generated.

Usage
-----
    python generate_nll_corrections.py
"""

import glob
import os

from NLL_run.append_station_delays import SecondRunParams, append_station_delays
from NLL_run.export_locdelay_info  import export_locdelay_info
from NLL_run.run_zone              import run_zone

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_LOC      = os.path.join(_PROJECT_ROOT, 'loc')
_RESULT   = os.path.join(_PROJECT_ROOT, 'RESULT')
_RUN      = os.path.join(_PROJECT_ROOT, 'run')
_STATIONS = os.path.join(_PROJECT_ROOT, 'stations')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Clean first-pass NLL output and generate second-pass run files for all zones."""
    for key in range(1, 7):
        # Clean up .hdr files left by NLL in the loc folder
        for hdr in glob.glob(os.path.join(_LOC, f'GLOBAL_{key}', '*.hdr')):
            os.remove(hdr)

        # Generate the second-pass run file
        params_ssst_W = SecondRunParams(
            locFolderName = os.path.join(_LOC, f'GLOBAL_{key}'),
            fileRunName   = os.path.join(_RUN, f'run_{key}.in'),
            fileRunSave   = os.path.join(_RUN, f'run_{key}_PR.in'),
            minPhases     = 100,  # minimal number of phases for the delay to be used
        )

        append_station_delays(params_ssst_W)

        run_zone(os.path.join(_RUN, f'run_{key}_PR.in'), corrections_pass=True)

    # Export the locdelays
    export_locdelay_info(
        run_dir      = _RUN,
        codemap_path = os.path.join(_STATIONS, 'GLOBAL_code_map.txt'),
        output_path  = os.path.join(_RUN, 'locdelays', 'locdelay_summary.txt'),
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    run_pipeline()


if __name__ == '__main__':
    main()
