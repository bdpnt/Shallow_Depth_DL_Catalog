"""
reformate_obs.py
============================
Convert a multi-event .obs bulletin into one .nlloc_obs file per event.

Each event block of the bulletin becomes a standalone NLLOC_OBS file named
<publicId>.nlloc_obs (the PUBLIC_ID line of the block gives the identity;
NLLoc propagates it into its .hyp/.csv output, which is what allows the
relocated events to be rematched to the source bulletin afterwards). Events
are filtered by a lat/lon bounding box and a minimum phase count; a PriorWt
of 1.00e+00 is appended to every pick line, before any trailing '#' comment.

Port of CODES_SSST/reformate_box.py.

Usage
-----
    python NLL_run/reformate_obs.py \\
        --bulletin  obs/GLOBAL_1_SSST.obs \\
        --output    obs/nlloc_obs/GLOBAL_1 \\
        --lat-min 42.5 --lat-max 43.5 --lon-min -2.0 --lon-max -0.75 \\
        --min-phases 0
"""

import argparse
import glob
import logging
import os
from dataclasses import dataclass
from datetime import datetime

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

logger = logging.getLogger('reformate_obs')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ReformateObsParams:
    fileBulletin: str
    outputDir:    str
    latMin:       float
    latMax:       float
    lonMin:       float
    lonMax:       float
    minPhases:    int = 0


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
    return log_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_QUALITY_LINE = ('QUALITY  Pmax -1 MFmin -1 MFmax -1 RMS -1 Nphs -1 Gap -1 '
                 'Dist -1 Mamp 5.97 4 Mdur -9.90 0\n')
_PHASE_HEADER = ('PHASE ID Ins Cmp On Pha  FM Date     HrMn   Sec     Err  '
                 'ErrMag    Coda      Amp       Per       PriorWt\n')


def _write_event(params, event_id, latitude, longitude, pick_lines):
    """Write one per-event .nlloc_obs file if the event passes the filters.

    Returns True when the file was written, False when the event was skipped
    (outside the bounding box or fewer than minPhases picks).
    """
    if event_id is None:
        return False
    if not (params.latMin <= latitude <= params.latMax and
            params.lonMin <= longitude <= params.lonMax):
        return False
    if len(pick_lines) < params.minPhases:
        return False

    output_file = os.path.join(params.outputDir, f'{event_id}.nlloc_obs')
    with open(output_file, 'w') as f:
        f.write(f'PUBLIC_ID {event_id}\n')
        f.write(_QUALITY_LINE)
        f.write(_PHASE_HEADER)
        f.writelines(pick_lines)
    return True


def _add_prior_wt(line):
    """Append the PriorWt field (1.00e+00) to a pick line, keeping any
    trailing '#' comment after it."""
    stripped = line.rstrip()
    if '#' in stripped:
        fields, comment = stripped.split('#', 1)
        return fields.rstrip() + '  1.00e+00 #' + comment + '\n'
    return stripped + '  1.00e+00\n'


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def reformate_obs(parameters, log_dir=None):
    """
    Split a .obs bulletin into per-event .nlloc_obs files.

    Parameters
    ----------
    parameters : ReformateObsParams
    log_dir    : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: output, log, n_written, n_skipped
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Bulletin     : {parameters.fileBulletin}")
    logger.info(f"Output dir   : {parameters.outputDir}")
    logger.info(f"Box          : lat [{parameters.latMin}, {parameters.latMax}]"
                f"  lon [{parameters.lonMin}, {parameters.lonMax}]")
    logger.info(f"Min phases   : {parameters.minPhases}")

    os.makedirs(parameters.outputDir, exist_ok=True)

    # NLLoc locates every *.nlloc_obs in the folder, so files left by an
    # earlier run would be relocated with this one (and rematched by a
    # publicId that may now name another event)
    stale = glob.glob(os.path.join(parameters.outputDir, '*.nlloc_obs'))
    for path in stale:
        os.remove(path)
    logger.info(f"Stale files removed : {len(stale)}")

    n_written = 0
    n_skipped = 0

    event_id   = None
    latitude   = None
    longitude  = None
    pick_lines = []

    with open(parameters.fileBulletin, 'r') as f:
        for line in f:
            if line.startswith('###'):
                continue

            if line.startswith('#'):
                # new event header: flush the previous event first
                if event_id is not None:
                    if _write_event(parameters, event_id, latitude, longitude, pick_lines):
                        n_written += 1
                    else:
                        n_skipped += 1

                metadata  = line.strip().split()
                latitude  = float(metadata[7])
                longitude = float(metadata[8])
                # fallback identity from the origin time; replaced by the
                # PUBLIC_ID line when the bulletin has one
                date_str = (f"{int(float(metadata[1])):04d}-"
                            f"{int(float(metadata[2])):02d}-"
                            f"{int(float(metadata[3])):02d}")
                time_str = (f"{int(float(metadata[4])):02d}-"
                            f"{int(float(metadata[5])):02d}-"
                            f"{float(metadata[6]):09.6f}")
                event_id   = f'{date_str}T{time_str}'
                pick_lines = []

            elif line.startswith('PUBLIC_ID'):
                event_id = line.split()[1]

            elif line.strip():
                pick_lines.append(_add_prior_wt(line))

    # flush the last event of the file
    if event_id is not None:
        if _write_event(parameters, event_id, latitude, longitude, pick_lines):
            n_written += 1
        else:
            n_skipped += 1

    logger.info(f"Events written : {n_written}")
    logger.info(f"Events skipped : {n_skipped}")

    return {
        'output':    parameters.outputDir,
        'log':       log_path,
        'n_written': n_written,
        'n_skipped': n_skipped,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Split a .obs bulletin into per-event .nlloc_obs files.'
    )
    parser.add_argument('--bulletin',   required=True, help='Input .obs bulletin')
    parser.add_argument('--output',     required=True, help='Output directory for .nlloc_obs files')
    parser.add_argument('--lat-min',    type=float, default=-90.0)
    parser.add_argument('--lat-max',    type=float, default=90.0)
    parser.add_argument('--lon-min',    type=float, default=-180.0)
    parser.add_argument('--lon-max',    type=float, default=180.0)
    parser.add_argument('--min-phases', type=int, default=0,
                        help='Minimum number of pick lines per event (default: 0)')
    parser.add_argument('--log-dir',    default=None,
                        help='Log directory (default: NLL_run/console_output/)')
    args = parser.parse_args()

    reformate_obs(
        ReformateObsParams(
            fileBulletin = args.bulletin,
            outputDir    = args.output,
            latMin       = args.lat_min,
            latMax       = args.lat_max,
            lonMin       = args.lon_min,
            lonMax       = args.lon_max,
            minPhases    = args.min_phases,
        ),
        log_dir = args.log_dir,
    )


if __name__ == '__main__':
    main()
