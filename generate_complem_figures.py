"""
generate_complem_figures.py
============================
Produce all complementary figures for the final catalog.

Runs in order:
  1. Gutenberg-Richter distribution (GLOBAL.obs — all events, ML and Mw)
  2. Gutenberg-Richter distribution (FINAL.obs — relocated events, ML and Mw)
  3. Depth maps  per 5-year period (RESULT/FINAL.txt)
  4. Error maps  per 5-year period (RESULT/FINAL.txt)
  5. Event map   for each of the 6 NLL zones  (RESULT/GLOBAL_k_PR.txt)
  6. Event map   for the final merged catalog  (obs/FINAL.obs)

Usage
-----
    python generate_complem_figures.py
"""

import os

from complem_figures.depth_maps        import DepthMapsParams,        generate_figure as gen_depth
from complem_figures.error_maps        import ErrorMapsParams,         generate_figure as gen_error
from complem_figures.event_maps        import EventMapsParams,         generate_figure as gen_event
from complem_figures.gutenberg_richter import GutenbergRichterParams,  generate_figure as gen_gr

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_OBS          = os.path.join(_PROJECT_ROOT, 'obs')
_RESULT       = os.path.join(_PROJECT_ROOT, 'RESULT')
_LOC          = os.path.join(_PROJECT_ROOT, 'loc')
_FIGS         = os.path.join(_PROJECT_ROOT, 'complem_figures')

# ---------------------------------------------------------------------------
# Zone configs for event maps (inner box, outer box)
# ---------------------------------------------------------------------------

_ZONE_CONFIGS = {
    1: (((42.50, -2.00), (43.50, -0.75)), ((41.60, -3.22), (44.40,  0.46))),
    2: (((42.50, -1.00), (43.25,  0.50)), ((41.60, -2.22), (44.15,  1.71))),
    3: (((42.00,  0.25), (43.25,  1.00)), ((41.10, -0.96), (44.15,  2.20))),
    4: (((42.00,  0.75), (43.00,  2.25)), ((41.10, -0.46), (43.90,  3.45))),
    5: (((42.00,  2.00), (43.00,  3.50)), ((41.10,  0.79), (43.90,  4.70))),
    6: (((42.75,  2.25), (43.75,  3.50)), ((41.85,  1.03), (44.65,  4.75))),
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Generate all complementary figures for the final catalog."""

    # 1 & 2. Gutenberg-Richter
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
        file       = os.path.join(_RESULT, 'FINAL.txt'),
        mapFolder  = os.path.join(_FIGS, 'depth_maps') + os.sep,
        time_range = 5,
    ))

    # 4. Error maps
    gen_error(ErrorMapsParams(
        file       = os.path.join(_RESULT, 'FINAL.txt'),
        mapFolder  = os.path.join(_FIGS, 'error_maps') + os.sep,
        time_range = 5,
    ))

    # 5. Event maps — per NLL zone
    for key, (region_in, region_out) in _ZONE_CONFIGS.items():
        gen_event(EventMapsParams(
            fileBulletin = os.path.join(_RESULT, f'GLOBAL_{key}_PR.txt'),
            figSave      = os.path.join(_FIGS, 'event_maps', f'GLOBAL_{key}.pdf'),
            fileStations = os.path.join(_LOC,    f'GLOBAL_{key}', 'last.stations'),
            region_in    = region_in,
            region_out   = region_out,
        ))

    # 6. Event map — final merged catalog
    gen_event(EventMapsParams(
        fileBulletin = os.path.join(_OBS, 'FINAL.obs'),
        figSave      = os.path.join(_FIGS, 'event_maps', 'FINAL.pdf'),
    ))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    run_pipeline()


if __name__ == '__main__':
    main()
