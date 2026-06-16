"""
finalize_nll_catalog.py
============================
Finalize the NLL-relocated catalog and produce obs/FINAL.obs.

Steps:
  1. Removes .hdr files left by NLL in each zone's loc/ folder.
  2. Merges the 6 zone CSV outputs into RESULT/FINAL.csv, resolving zone-overlap
     duplicates by keeping the solution with the smallest pdfVolume.
  3. Rematches relocated events back to obs/GLOBAL.obs via publicId to recover
     magnitude and picks, and writes obs/FINAL.obs.

Usage
-----
    python finalize_nll_catalog.py
"""

import glob
import os

from NLL_run.match_pre_post_relocation import MatchCatalogsParams, save_bulletin
from NLL_run.merge_regional_results    import merge_bulletins

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_LOC    = os.path.join(_PROJECT_ROOT, 'loc')
_RESULT = os.path.join(_PROJECT_ROOT, 'RESULT')
_OBS    = os.path.join(_PROJECT_ROOT, 'obs')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Remove .hdr files, merge zone CSVs, and produce obs/FINAL.obs."""
    # Step 1: clean up .hdr files left by NLL in each zone folder
    for key in range(1, 7):
        for hdr in glob.glob(os.path.join(_LOC, f'GLOBAL_{key}', '*.hdr')):
            os.remove(hdr)

    # Step 2: merge all zone CSV outputs into RESULT/FINAL.csv
    csv_files = [
        os.path.join(_LOC, f'GLOBAL_{key}', f'GLOBAL_{key}.obs.sum.grid0.loc.csv')
        for key in range(1, 7)
    ]
    merge_bulletins(csv_files, os.path.join(_RESULT, 'FINAL.csv'))

    # Step 3: rematch to obs/GLOBAL.obs via publicId and write obs/FINAL.obs
    save_bulletin(MatchCatalogsParams(
        file_obs   = os.path.join(_OBS, 'GLOBAL.obs'),
        file_final = os.path.join(_RESULT, 'FINAL.csv'),
        save_file  = os.path.join(_OBS, 'FINAL.obs'),
    ))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    run_pipeline()


if __name__ == '__main__':
    main()
