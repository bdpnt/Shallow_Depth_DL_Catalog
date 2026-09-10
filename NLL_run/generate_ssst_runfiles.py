"""
generate_ssst_runfiles.py
============================
Generate the SSST-stage input files for one geographic zone.

The SSST stage relocates the augmented catalog (obs/NLL_result_augmented.obs),
so for a given zone this module produces:
  - a zone bulletin cut from the augmented catalog (obs/GLOBAL_<N>_SSST.obs)
  - a GTSRCE station file for the stations picked in that zone bulletin
  - the NLLoc control file  run/ssst/run_<N>_NLL.in
  - the Loc2ssst control file run/ssst/run_<N>_SSST.in

Both control files are DERIVED from the NLL-stage run file of the same zone
(run/nll/run_<N>_DELAYS.in): TRANS, VGGRID, LAYER, LOCGRID, LOCSEARCH and all
localization parameters are reused verbatim, so the zone information entered
once in run_NLL.py stays the single source of truth (and the TRANS of the two
control files can never diverge). Only the statements that differ for SSST
are rewritten: grid paths (run/ssst_model, run/ssst_time), the GTSRCE include,
and LOCHYPOUT (SAVE_NLLOC_ALL + SAVE_NLLOC_EXPECTATION, without the
SAVE_NLLOC_OCTREE/SAVE_FMAMP flags that run_ssst.py appends itself on the final
iteration only).

The Loc2ssst LSGRID/LSOUTGRID are computed from the parsed LOCGRID extent,
rescaled to a coarse 1.0 km spacing (Loc2ssst holds two full LSGRID buffers
in RAM per instance — never copy the fine LOCGRID dimensions).

Usage
-----
    python NLL_run/generate_ssst_runfiles.py \\
        --bulletin      obs/NLL_result_augmented.obs \\
        --bulletin-out  obs/GLOBAL_1_SSST.obs \\
        --inventory     stations/GLOBAL_inventory.xml \\
        --code-map      stations/GLOBAL_code_map.txt \\
        --stations      stations/GTSRCE_SSST_1.txt \\
        --run-nll       run/nll/run_1_DELAYS.in \\
        --run-save-nll  run/ssst/run_1_NLL.in \\
        --run-save-ssst run/ssst/run_1_SSST.in \\
        --model         run/ssst_model/Pyrenees_1/Pyrenees_1 \\
        --time-file     run/ssst_time/Pyrenees_1/Pyrenees_1 \\
        --lat-min 42.5 --lat-max 43.5 --lon-min -2.0 --lon-max -0.75
"""

import argparse
import logging
import math
import os
import sys
from dataclasses import dataclass, field
from datetime import datetime
from types import SimpleNamespace

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from NLL_run import generate_regional_runfiles as _grf

logger = logging.getLogger('generate_ssst_runfiles')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GenSSSTRunParams:
    fileBulletin:    str    # parent bulletin (obs/NLL_result_augmented.obs)
    fileInventory:   str
    fileMap:         str
    fileBulletinIn:  str    # output: zone bulletin (obs/GLOBAL_<N>_SSST.obs)
    fileStations:    str    # output: GTSRCE file (stations/GTSRCE_SSST_<N>.txt)
    fileRunNLL:      str    # input:  NLL-stage run file (run/nll/run_<N>_DELAYS.in)
    fileRunSaveNLL:  str    # output: NLLoc control (run/ssst/run_<N>_NLL.in)
    fileRunSaveSSST: str    # output: Loc2ssst control (run/ssst/run_<N>_SSST.in)
    fileModel:       str    # run/ssst_model/Pyrenees_<N>/Pyrenees_<N>
    fileTime:        str    # run/ssst_time/Pyrenees_<N>/Pyrenees_<N>
    latMin_event:    float
    latMax_event:    float
    lonMin_event:    float
    lonMax_event:    float
    # LSPHSTAT values: RMSMax NRdgsMin GapMax PResidualMax SResidualMax EllLen3Max.
    # NRdgsMin doubles as the NLLoc min-phases threshold (run_ssst.py reads it
    # back from the generated file, so the two selections cannot drift apart).
    lsphstat:        list = field(default_factory=lambda: [0.15, 6, 200.0, 0.3, 0.5, 10.0])
    lsGridSpacing:   float = 1.0   # km — coarse Loc2ssst correction grid
    lsGridNz:        int   = 40    # 40 nodes from -3 km -> 36 km depth


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

def _setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    basename  = os.path.splitext(os.path.basename(__file__))[0]
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_path  = os.path.join(log_dir, f"{basename}_{timestamp}.log")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(log_path, encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(handler)
    # route the reused generate_regional_runfiles helpers into the same log
    _grf.logger.setLevel(logging.INFO)
    _grf.logger.handlers.clear()
    _grf.logger.addHandler(handler)
    return log_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _verify_folders(parameters):
    """Create all parent directories required by the output file paths."""
    for path in [
        parameters.fileBulletinIn,
        parameters.fileStations,
        parameters.fileRunSaveNLL,
        parameters.fileRunSaveSSST,
        parameters.fileModel,
        parameters.fileTime,
    ]:
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)


def _transform_nll_runfile(nll_lines, parameters):
    """Rewrite the NLL-stage run file into the SSST-stage NLLoc control.

    Everything passes through unchanged except the statements that differ
    for the SSST stage (grid paths, GTSRCE include, LOCHYPOUT) and the two
    statements that run_ssst.py appends per NLLoc chunk anyway (LOCFILES,
    LOCMETH), which are replaced by explanatory comments.
    """
    out = []
    for line in nll_lines:
        stripped = line.lstrip()
        if stripped.startswith('VGOUT'):
            out.append(f'VGOUT  {parameters.fileModel}\n')
        elif stripped.startswith('GTFILES'):
            # P pass; run_ssst.build_grids() derives the S variant itself
            out.append(f'GTFILES  {parameters.fileModel}  {parameters.fileTime} P\n')
        elif stripped.startswith('INCLUDE'):
            out.append(f'INCLUDE {parameters.fileStations}\n')
        elif stripped.startswith('LOCHYPOUT'):
            # must NOT contain SAVE_NLLOC_OCTREE / SAVE_FMAMP: run_ssst.py
            # appends the authoritative LOCHYPOUT (NLL keeps the last one)
            # and adds those flags on the final iteration only
            out.append('LOCHYPOUT SAVE_NLLOC_ALL  NLL_FORMAT_VER_2  '
                       'SAVE_NLLOC_EXPECTATION\n')
        elif stripped.startswith('LOCFILES'):
            out.append('# LOCFILES appended per NLLoc chunk by run_ssst.py\n')
        elif stripped.startswith('LOCMETH'):
            out.append('# LOCMETH appended per NLLoc chunk by run_ssst.py\n')
        else:
            out.append(line)
    return out


def _build_ssst_runfile(nll_lines, parameters):
    """Build the Loc2ssst control file from the parsed NLL run file lines.

    CONTROL, TRANS and the LOCPHASEID mappings are copied verbatim; the
    LSGRID/LSOUTGRID geometry is derived from the LOCGRID extent at the
    coarse lsGridSpacing. The per-iteration statements (LSPARAMS, LSOUT,
    LSLOCFILES, LOCFILES, LOCMETH) and the per-instance LSSTATIONS are
    appended at runtime by run_ssst.py and must not appear here.
    """
    control_line   = None
    trans_line     = None
    locgrid_fields = None
    phaseid_lines  = []

    for line in nll_lines:
        stripped = line.lstrip()
        if stripped.startswith('CONTROL'):
            control_line = stripped
        elif stripped.startswith('TRANS'):
            trans_line = stripped
        elif stripped.startswith('LOCGRID'):
            locgrid_fields = stripped.split()
        elif stripped.startswith('LOCPHASEID'):
            phaseid_lines.append(stripped)

    for name, value in [('CONTROL', control_line), ('TRANS', trans_line),
                        ('LOCGRID', locgrid_fields)]:
        if value is None:
            raise ValueError(f'no {name} statement found in {parameters.fileRunNLL}')

    # LOCGRID nx ny nz x0 y0 z0 dx dy dz ... -> extent in km, re-gridded at
    # lsGridSpacing (+1 closing node) so the correction grid fully covers
    # the search grid horizontally
    nx, ny = int(locgrid_fields[1]), int(locgrid_fields[2])
    dx, dy = float(locgrid_fields[7]), float(locgrid_fields[8])
    x0, y0, z0 = locgrid_fields[4], locgrid_fields[5], locgrid_fields[6]
    spacing = parameters.lsGridSpacing
    nx_ls = math.ceil((nx - 1) * dx / spacing) + 1
    ny_ls = math.ceil((ny - 1) * dy / spacing) + 1
    nz_ls = parameters.lsGridNz

    lsphstat = ' '.join(str(v) for v in parameters.lsphstat)

    lines = []
    lines.append(control_line if control_line.endswith('\n') else control_line + '\n')
    lines.append(trans_line if trans_line.endswith('\n') else trans_line + '\n')
    lines.append('\n')
    lines.append('#LSPARAMS, LSOUT, LSLOCFILES, LOCFILES, LOCMETH and LSSTATIONS\n')
    lines.append('#are appended per iteration / per instance by run_ssst.py\n')
    lines.append('\n')
    lines.append(f'LSGRID  {nx_ls} {ny_ls} {nz_ls} {x0} {y0} {z0}  '
                 f'{spacing} {spacing} {spacing}  SSST_TIMECORR FLOAT\n')
    lines.append(f'LSOUTGRID {nx_ls} {ny_ls} {nz_ls} {x0} {y0} {z0}  '
                 f'{spacing} {spacing} {spacing} TIME FLOAT\n')
    lines.append('\n')
    lines.append('LSMODE ANGLES_NO\n')
    lines.append(f'LSPHSTAT {lsphstat}\n')
    lines.append('\n')
    for phaseid in phaseid_lines:
        lines.append(phaseid if phaseid.endswith('\n') else phaseid + '\n')
    return lines


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_ssst_run(parameters, log_dir=None):
    """
    Generate the zone bulletin, GTSRCE, and both SSST control files.

    Parameters
    ----------
    parameters : GenSSSTRunParams
    log_dir    : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: nll_run, ssst_run, bulletin, stations, log
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Zone bounds  : lat [{parameters.latMin_event}, {parameters.latMax_event}]"
                f"  lon [{parameters.lonMin_event}, {parameters.lonMax_event}]")
    logger.info(f"Source run   : {parameters.fileRunNLL}")

    if not os.path.isfile(parameters.fileRunNLL):
        raise FileNotFoundError(
            f'NLL-stage run file not found: {parameters.fileRunNLL} '
            '(run the NLL pipeline first — the SSST control files are derived from it)')

    _verify_folders(parameters)

    # zone bulletin cut from the augmented catalog (reused NLL-stage helper)
    _grf._gen_child_obs(parameters)

    # GTSRCE from the ZONE bulletin: only the stations actually picked in
    # this zone get travel-time grids (the NLL stage uses the full catalog)
    _grf._gen_gtsrce(SimpleNamespace(
        fileBulletin  = parameters.fileBulletinIn,
        fileInventory = parameters.fileInventory,
        fileMap       = parameters.fileMap,
        fileStations  = parameters.fileStations,
    ))

    with open(parameters.fileRunNLL, 'r') as f:
        nll_lines = f.readlines()

    # grid-containment guard: a LOCGRID reaching the bottom of the 800-node
    # travel-time grid (36.95 km) triggers an interpolation artifact; the
    # project generates 761 nodes (35 km) since 2026-07 — an older value here
    # means the NLL-stage run files predate that change
    for line in nll_lines:
        fields = line.split()
        if fields and fields[0] == 'LOCGRID':
            max_depth = float(fields[6]) + (int(fields[3]) - 1) * float(fields[9])
            if max_depth > 35.01:
                logger.warning(
                    f'LOCGRID of {parameters.fileRunNLL} reaches {max_depth:g} km '
                    '(> 35 km): regenerate the NLL-stage run files (nz=761) '
                    'before running SSST')
            break

    with open(parameters.fileRunSaveNLL, 'w', encoding='utf-8') as f:
        f.writelines(_transform_nll_runfile(nll_lines, parameters))
    logger.info(f"NLLoc control    : {parameters.fileRunSaveNLL}")

    with open(parameters.fileRunSaveSSST, 'w', encoding='utf-8') as f:
        f.writelines(_build_ssst_runfile(nll_lines, parameters))
    logger.info(f"Loc2ssst control : {parameters.fileRunSaveSSST}")

    return {
        'nll_run':  parameters.fileRunSaveNLL,
        'ssst_run': parameters.fileRunSaveSSST,
        'bulletin': parameters.fileBulletinIn,
        'stations': parameters.fileStations,
        'log':      log_path,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate the SSST-stage zone bulletin, GTSRCE and control files.'
    )
    parser.add_argument('--bulletin',      required=True, help='Parent bulletin (obs/NLL_result_augmented.obs)')
    parser.add_argument('--bulletin-out',  required=True, help='Output zone bulletin')
    parser.add_argument('--inventory',     required=True, help='StationXML inventory')
    parser.add_argument('--code-map',      required=True, help='GLOBAL_code_map.txt')
    parser.add_argument('--stations',      required=True, help='Output GTSRCE file')
    parser.add_argument('--run-nll',       required=True, help='Input NLL-stage run file (run/nll/run_<N>_DELAYS.in)')
    parser.add_argument('--run-save-nll',  required=True, help='Output NLLoc control file')
    parser.add_argument('--run-save-ssst', required=True, help='Output Loc2ssst control file')
    parser.add_argument('--model',         required=True, help='SSST velocity-grid file root')
    parser.add_argument('--time-file',     required=True, help='SSST travel-time grid file root')
    parser.add_argument('--lat-min', type=float, required=True)
    parser.add_argument('--lat-max', type=float, required=True)
    parser.add_argument('--lon-min', type=float, required=True)
    parser.add_argument('--lon-max', type=float, required=True)
    parser.add_argument('--lsphstat', default=None,
                        help='LSPHSTAT values, comma-separated '
                             '(default: 0.15,6,200.0,0.3,0.5,10.0)')
    parser.add_argument('--log-dir', default=None,
                        help='Log directory (default: NLL_run/console_output/)')
    args = parser.parse_args()

    kwargs = {}
    if args.lsphstat:
        kwargs['lsphstat'] = [float(v) if '.' in v else int(v)
                              for v in args.lsphstat.split(',')]

    generate_ssst_run(
        GenSSSTRunParams(
            fileBulletin    = args.bulletin,
            fileInventory   = args.inventory,
            fileMap         = args.code_map,
            fileBulletinIn  = args.bulletin_out,
            fileStations    = args.stations,
            fileRunNLL      = args.run_nll,
            fileRunSaveNLL  = args.run_save_nll,
            fileRunSaveSSST = args.run_save_ssst,
            fileModel       = args.model,
            fileTime        = args.time_file,
            latMin_event    = args.lat_min,
            latMax_event    = args.lat_max,
            lonMin_event    = args.lon_min,
            lonMax_event    = args.lon_max,
            **kwargs,
        ),
        log_dir = args.log_dir,
    )


if __name__ == '__main__':
    main()
