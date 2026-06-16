"""
generate_complem_maps.py
============================
Produce PyGMT event maps for all NLL zones and the final catalog.

Requires: pygmt_env

Runs in order:
  1. Event map for each of the 6 NLL zones (loc/GLOBAL_k/GLOBAL_k.obs.sum.grid0.loc.csv)
  2. Event map for the final merged catalog (obs/FINAL.obs)

For matplotlib-based figures (Gutenberg-Richter, depth, error maps),
run generate_complem_figures.py with seisbench_env.

Usage
-----
    python generate_complem_maps.py
"""

import os

from complem_figures.event_maps import EventMapsParams, generate_figure as gen_event

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_OBS          = os.path.join(_PROJECT_ROOT, 'obs')
_RESULT       = os.path.join(_PROJECT_ROOT, 'RESULT')
_LOC          = os.path.join(_PROJECT_ROOT, 'loc')
_FIGS         = os.path.join(_PROJECT_ROOT, 'complem_figures')

# ---------------------------------------------------------------------------
# Zone configs (inner box, outer box)
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
    """Generate PyGMT event maps for each NLL zone and the final catalog."""

    # 1. Per-zone event maps (read second-pass NLL CSV directly)
    for key, (region_in, region_out) in _ZONE_CONFIGS.items():
        gen_event(EventMapsParams(
            fileBulletin = os.path.join(_LOC, f'GLOBAL_{key}',
                                        f'GLOBAL_{key}.obs.sum.grid0.loc.csv'),
            figSave      = os.path.join(_FIGS, 'event_maps', f'GLOBAL_{key}.pdf'),
            fileStations = os.path.join(_LOC,  f'GLOBAL_{key}', 'last.stations'),
            region_in    = region_in,
            region_out   = region_out,
        ))

    # 2. Final merged catalog
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
