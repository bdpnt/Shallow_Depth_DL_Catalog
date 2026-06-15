"""
generate_complem_figures.py
============================
Produce matplotlib-based complementary figures for the final catalog.

Requires: seisbench_env

Runs in order:
  1. Depth histograms (GLOBAL.obs and FINAL.obs)
  2. Gutenberg-Richter distribution (GLOBAL.obs — all events, ML and Mw)
  3. Gutenberg-Richter distribution (FINAL.obs — relocated events, ML and Mw)
  4. Depth maps  per 5-year period (RESULT/FINAL.csv)
  5. Error maps  per 5-year period (RESULT/FINAL.csv)

For PyGMT event maps, run generate_complem_maps.py with pygmt_env.

Usage
-----
    python generate_complem_figures.py
"""

import os

from complem_figures.depth_histogram   import DepthHistogramParams,   generate_figure as gen_depth_hist
from complem_figures.depth_maps        import DepthMapsParams,        generate_figure as gen_depth
from complem_figures.error_maps        import ErrorMapsParams,        generate_figure as gen_error
from complem_figures.gutenberg_richter import GutenbergRichterParams, generate_figure as gen_gr

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_OBS          = os.path.join(_PROJECT_ROOT, 'obs')
_RESULT       = os.path.join(_PROJECT_ROOT, 'RESULT')
_FIGS         = os.path.join(_PROJECT_ROOT, 'complem_figures')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Generate depth histograms, Gutenberg-Richter, depth, and error map figures."""

    # 1 & 2. Depth histograms
    for bulletin, tag in [
        (os.path.join(_OBS, 'GLOBAL.obs'), 'GLOBAL'),
        (os.path.join(_OBS, 'FINAL.obs'),  'FINAL'),
    ]:
        gen_depth_hist(DepthHistogramParams(
            file_bulletin = bulletin,
            fig_save      = os.path.join(_FIGS, 'depth_histogram', f'{tag}.png'),
        ))

    # 3 & 4. Gutenberg-Richter
    for bulletin, tag in [
        (os.path.join(_OBS, 'GLOBAL.obs'), 'GLOBAL'),
        (os.path.join(_OBS, 'FINAL.obs'),  'FINAL'),
    ]:
        for mag_type in ('ML', 'Mw'):
            gen_gr(GutenbergRichterParams(
                file_bulletin = bulletin,
                fig_save      = os.path.join(_FIGS, 'gutenberg_richter', f'{tag}_{mag_type}.png'),
                mag_type      = mag_type,
            ))

    # 3. Depth maps
    gen_depth(DepthMapsParams(
        file       = os.path.join(_RESULT, 'FINAL.csv'),
        mapFolder  = os.path.join(_FIGS, 'depth_maps') + os.sep,
        time_range = 5,
    ))

    # 4. Error maps
    gen_error(ErrorMapsParams(
        file       = os.path.join(_RESULT, 'FINAL.csv'),
        mapFolder  = os.path.join(_FIGS, 'error_maps') + os.sep,
        time_range = 5,
    ))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    run_pipeline()


if __name__ == '__main__':
    main()
