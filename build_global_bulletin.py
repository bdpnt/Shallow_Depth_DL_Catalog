"""
build_global_bulletin.py
============================
Harmonize and fuse all source bulletins into a single GLOBAL.obs catalog.

Steps:
  1. Remap picks to unified station codes.
  2. Generate magnitude conversion models (MLv RESIF, mb_Lg IGN, ML ICGC → ML LDG).
  3. Apply magnitude models to all source bulletins.
  4. Filter each bulletin to the area of interest.
  5. Fuse all bulletins into obs/GLOBAL.obs, then find and merge doubles.

Usage
-----
    python build_global_bulletin.py
"""

import os
import subprocess

from global_obs.apply_magnitude_models     import UpdateMagFilesParams
from global_obs.fuse_bulletins             import FusionParams, MergeDoublesParams
from global_obs.generate_magnitude_models  import MagModelParams
from global_obs.remap_picks_to_unified_codes import AssociatePicksParams
import global_obs

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

_PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
_OBS        = os.path.join(_PROJECT_ROOT, 'obs')
_MAGMODELS  = os.path.join(_PROJECT_ROOT, 'MAGMODELS')
_GLOBAL_OBS = os.path.join(_PROJECT_ROOT, 'global_obs')
_STATIONS   = os.path.join(_PROJECT_ROOT, 'stations')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_pipeline():
    """Run the full catalog harmonization and fusion pipeline."""
    # Associate picks
    params_association = AssociatePicksParams(
        file_inventory  = os.path.join(_STATIONS, 'GLOBAL_inventory.xml'),
        folder_bulletin = os.path.join(_OBS, '*.obs'),
    )

    global_obs.remap_picks_to_unified_codes.remap_picks_to_unified_codes(params_association)

    # Generate magnitude models
    params_magModel_RESIF = MagModelParams(
        file_name1  = os.path.join(_OBS, 'RESIF_20-25.obs'),
        file_name2  = os.path.join(_OBS, 'LDG_20-25.obs'),
        mag_type1   = 'MLv',
        mag_type2   = 'ML',
        mag_name1   = 'MLv RESIF',
        mag_name2   = 'ML LDG',
        dist_thresh = 10.0,
        time_thresh = 2.0,
        save_name   = os.path.join(_MAGMODELS, 'MLv RESIF.joblib'),
        save_figs   = os.path.join(_MAGMODELS, 'FIGURES', ''),
    )

    global_obs.generate_magnitude_models.convert_magnitudes(params_magModel_RESIF, save_figs=True)

    params_magModel_IGN = MagModelParams(
        file_name1  = os.path.join(_OBS, 'IGN_20-25.obs'),
        file_name2  = os.path.join(_OBS, 'LDG_20-25.obs'),
        mag_type1   = 'mb_Lg',
        mag_type2   = 'ML',
        mag_name1   = 'mb_Lg IGN',
        mag_name2   = 'ML LDG',
        dist_thresh = 10.0,
        time_thresh = 2.0,
        save_name   = os.path.join(_MAGMODELS, 'mb_Lg IGN.joblib'),
        save_figs   = os.path.join(_MAGMODELS, 'FIGURES', ''),
    )

    global_obs.generate_magnitude_models.convert_magnitudes(params_magModel_IGN, save_figs=True)

    params_magModel_ICGC = MagModelParams(
        file_name1  = os.path.join(_OBS, 'ICGC_20-25.obs'),
        file_name2  = os.path.join(_OBS, 'LDG_20-25.obs'),
        mag_type1   = 'ML',
        mag_type2   = 'ML',
        mag_name1   = 'ML ICGC',
        mag_name2   = 'ML LDG',
        dist_thresh = 10.0,
        time_thresh = 2.0,
        save_name   = os.path.join(_MAGMODELS, 'ML ICGC.joblib'),
        save_figs   = os.path.join(_MAGMODELS, 'FIGURES', ''),
    )

    global_obs.generate_magnitude_models.convert_magnitudes(params_magModel_ICGC, save_figs=True)

    # Use magnitude models
    parameters_magModels = UpdateMagFilesParams(
        folder_path = os.path.join(_OBS, '*_20-25.obs'),
    )

    global_obs.apply_magnitude_models.apply_magnitude_models(parameters_magModels)

    # Update bulletins AOI
    subprocess.run(
        [
            "conda", "run", "-n", "pygmt_env", "python",
            os.path.join(_GLOBAL_OBS, 'filter_events_by_aoi.py'),
            "--file-names",
                os.path.join(_OBS, 'RESIF_20-25.obs'),
                os.path.join(_OBS, 'IGN_20-25.obs'),
                os.path.join(_OBS, 'ICGC_20-25.obs'),
                os.path.join(_OBS, 'LDG_20-25.obs'),
                os.path.join(_OBS, 'OMP_2016.obs'),
                os.path.join(_OBS, 'OMP_78-19.obs'),
            "--fig-saves",
                os.path.join(_OBS, 'MAPS', 'RESIF_20-25.pdf'),
                os.path.join(_OBS, 'MAPS', 'IGN_20-25.pdf'),
                os.path.join(_OBS, 'MAPS', 'ICGC_20-25.pdf'),
                os.path.join(_OBS, 'MAPS', 'LDG_20-25.pdf'),
                os.path.join(_OBS, 'MAPS', 'OMP_2016.pdf'),
                os.path.join(_OBS, 'MAPS', 'OMP_78-19.pdf'),
        ],
        check=True,
    )

    # Dedup each source catalog before fusion
    _source_catalogs = [
        os.path.join(_OBS, 'RESIF_20-25.obs'),
        os.path.join(_OBS, 'IGN_20-25.obs'),
        os.path.join(_OBS, 'ICGC_20-25.obs'),
        os.path.join(_OBS, 'LDG_20-25.obs'),
        os.path.join(_OBS, 'OMP_2016.obs'),
        os.path.join(_OBS, 'OMP_78-19.obs'),
    ]
    for _catalog_path in _source_catalogs:
        params_dedup = MergeDoublesParams(
            global_bulletin_path = _catalog_path,
            max_dt_seconds       = 1.0,
            max_dist_km          = 50.0,
        )
        global_obs.fuse_bulletins.find_and_merge_doubles(params_dedup)

    # Fusion all bulletins
    params_fusion = FusionParams(
        global_bulletin_path = os.path.join(_OBS, 'GLOBAL.obs'),
        main_bulletin_path   = os.path.join(_OBS, 'RESIF_20-25.obs'),
        folder_path          = os.path.join(_OBS, '*.obs'),
        dist_thresh          = 15,   # km
        loose_dist_thresh    = 50,   # km
        time_thresh          = 2,    # seconds
        loose_time_thresh    = 30,   # seconds
        mag_thresh           = 1.5,  # magnitude units
        sim_pick_thresh      = 2,    # minimum shared P-phases to confirm a loose match
    )

    global_obs.fuse_bulletins.fuse_bulletins(params_fusion)
    subprocess.run(
        [
            "conda", "run", "-n", "pygmt_env", "python",
            os.path.join(_GLOBAL_OBS, 'plot_global_catalog_map.py'),
            "--file-name", os.path.join(_OBS, 'GLOBAL.obs'),
            "--fig-save",  os.path.join(_OBS, 'MAPS', 'GLOBAL.pdf'),
        ],
        check=True,
    )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    run_pipeline()


if __name__ == '__main__':
    main()
