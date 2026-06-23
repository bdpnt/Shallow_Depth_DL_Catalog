"""
generate_complem_maps.py
============================
Produce PyGMT event maps for all NLL zones and the final catalog.

Requires: pygmt_env

Runs in order:
  1. Event map for each of the 6 NLL zones
  2. Event map for the final merged catalog

For matplotlib-based figures (Gutenberg-Richter, depth, error maps),
run generate_complem_figures.py with seisbench_env.

Usage
-----
    python generate_complem_maps.py
"""

import glob
import math
import os
import re

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
# Zone configs (inner box, outer box) — derived from loc/GLOBAL_<k>/last.in
# ---------------------------------------------------------------------------

_R = 6371.0
_TRANS_RE  = re.compile(r'^TRANS\s+LAMBERT\s+WGS-84\s+([\d.\-]+)\s+([\d.\-]+)', re.M)
_LOCGRD_RE = re.compile(r'^LOCGRID\s+(\d+)\s+(\d+)', re.M)


def _load_zone_configs(loc_dir):
    configs = {}
    for path in sorted(glob.glob(os.path.join(loc_dir, 'GLOBAL_*', 'last.in'))):
        m = re.search(r'GLOBAL_(\d+)', path)
        if not m:
            continue
        zone_id = int(m.group(1))
        with open(path) as f:
            content = f.read()
        mt = _TRANS_RE.search(content)
        ml = _LOCGRD_RE.search(content)
        if not (mt and ml):
            continue
        lat_sw_out, lon_sw_out = float(mt.group(1)), float(mt.group(2))
        nx, ny = int(ml.group(1)), int(ml.group(2))

        # Recover inner SW: exact inverse of _km_to_latlon(-100, -100, inner_sw)
        lat_sw_in = lat_sw_out + math.degrees(100.0 / _R)
        lon_sw_in = lon_sw_out + math.degrees(100.0 / (_R * math.cos(math.radians(lat_sw_in))))

        # Extents in km from inner SW origin
        x1 = nx * 0.05 - 200
        y1 = ny * 0.05 - 200

        # Inner NE: mirrors _latlon_to_km (average latitude for lon)
        lat_ne_in = lat_sw_in + math.degrees(y1 / _R)
        lon_ne_in = lon_sw_in + math.degrees(x1 / (_R * math.cos(math.radians((lat_sw_in + lat_ne_in) / 2))))

        # Outer NE: mirrors _km_to_latlon (cos of inner SW lat only)
        lat_ne_out = lat_sw_in + math.degrees((y1 + 100) / _R)
        lon_ne_out = lon_sw_in + math.degrees((x1 + 100) / (_R * math.cos(math.radians(lat_sw_in))))

        configs[zone_id] = (
            ((round(lat_sw_in,  2), round(lon_sw_in,  2)), (round(lat_ne_in,  2), round(lon_ne_in,  2))),
            ((round(lat_sw_out, 2), round(lon_sw_out, 2)), (round(lat_ne_out, 2), round(lon_ne_out, 2))),
        )
    return configs


_ZONE_CONFIGS = _load_zone_configs(_LOC)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Generate PyGMT event maps for each NLL zone and the final catalog."""

    # 1. Per-zone event maps (filter deduplicated FINAL.csv by source zone)
    for key, (region_in, region_out) in _ZONE_CONFIGS.items():
        gen_event(EventMapsParams(
            fileBulletin  = os.path.join(_RESULT, 'FINAL.csv'),
            figSave       = os.path.join(_FIGS, 'event_maps', f'GLOBAL_{key}.pdf'),
            fileStations  = os.path.join(_LOC,  f'GLOBAL_{key}', 'last.stations'),
            region_in     = region_in,
            region_out    = region_out,
            source_filter = f'GLOBAL_{key}',
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
