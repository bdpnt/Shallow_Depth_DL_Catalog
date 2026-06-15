"""
parse_nll_output.py
============================
DEPRECATED — no longer called by finalize_nll_catalog.py.

NLL now outputs a CSV (SAVE_NLLOC_SUM) which is read directly by
merge_regional_results.py. This hypo_71 parser is kept for reference only.

Parse NonLinLoc hypo_71 output and write a plain-text result bulletin.

Removes intermediate .hdr files from the NLL output folder, then reads the
hypo_71 summary file and writes each relocated event as a space-separated
line in the result bulletin.

Usage
-----
    python NLL_run/parse_nll_output.py \\
        --loc-folder loc/GLOBAL_1 \\
        --obs-file   GLOBAL_1.obs \\
        --output     RESULT/GLOBAL_1.txt
"""

import argparse
import logging
import os
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

logger = logging.getLogger('parse_nll_output')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class CleanPostRunParams:
    folderLoc:    str
    obsFile:      str  # base filename only (no path)
    fileBulletin: str  # output result file path


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

def _parse_hypo71(parameters):
    """
    Delete .hdr files from the NLL output folder and parse the hypo_71
    summary file into a DataFrame.

    Parameters
    ----------
    parameters : CleanPostRunParams

    Returns
    -------
    pd.DataFrame — one row per relocated event
    """
    current_folder = os.getcwd()
    os.chdir(parameters.folderLoc)

    for f in os.listdir():
        if f.endswith('.hdr'):
            os.remove(f)

    hypo_file = f'{parameters.obsFile}.sum.grid0.loc.hypo_71'
    data      = []

    with open(hypo_file, 'r') as f:
        f.readline()  # skip header 1
        f.readline()  # skip header 2
        for line in f:
            row = {}
            try:
                row['date']     = line[1:7].strip()
                row['heuremin'] = float(line[7:12])
                row['ssss']     = float(line[12:18])
                row['lat']      = float(line[18:21])
                row['latmin']   = float(line[21:27])
                row['lon']      = float(line[27:31])
                row['lonmin']   = float(line[31:38])
                row['prof']     = float(line[38:45])
                row['mag']      = float(line[47:51])
            except Exception:
                continue
            try:
                row['no']  = float(line[52:55])
                row['dm']  = float(line[55:58])
                row['gap'] = float(line[58:62])
                row['m']   = float(line[62:64])
                row['rms'] = float(line[64:69])
                row['erh'] = float(line[70:74])
                row['erv'] = float(line[75:79])
            except Exception:
                row['no'] = row['dm'] = row['gap'] = row['m'] = np.nan
                row['rms'] = row['erh'] = row['erv'] = np.nan
            data.append(row)

    os.chdir(current_folder)
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def write_events(parameters, log_dir=None):
    """
    Write NLL-relocated events to a plain-text result bulletin.

    Parameters
    ----------
    parameters : CleanPostRunParams
    log_dir    : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: output, log, n_events
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Loc folder   : {parameters.folderLoc}")
    logger.info(f"Hypo file    : {parameters.obsFile}.sum.grid0.loc.hypo_71")
    logger.info(f"Output       : {parameters.fileBulletin}")

    events = _parse_hypo71(parameters)

    with open(parameters.fileBulletin, 'w') as f:
        for _, row in events.iterrows():
            datestr1  = row['date']
            heuremin1 = f"{row['heuremin']:4.0f}"

            if len(heuremin1.strip()) <= 2:
                heure1 = 0
                min1   = int(heuremin1)
            else:
                heure1 = int(heuremin1[:-2])
                min1   = int(heuremin1[-2:])

            if len(datestr1) == 5:
                datestr1 = '0' + datestr1
            if len(datestr1) == 4:
                datestr1 = '00' + datestr1

            year      = int(datestr1[0:2])
            month     = int(datestr1[2:4])
            day       = int(datestr1[4:6])
            latitude  = row['lat'] + row['latmin'] / 60.0
            longitude = row['lon'] + row['lonmin'] / 60.0

            f.write(
                f"{year} {month} {day} "
                f"{heure1} {min1} {row['ssss']} "
                f"{latitude} {longitude} {row['prof']} {row['mag']} "
                f"{row['rms']} {row['no']} {row['erh']} {row['erv']} {row['gap']}\n"
            )

    rms_arr = events['rms'].to_numpy()
    erh_arr = events['erh'].to_numpy()
    erv_arr = events['erv'].to_numpy()

    logger.info(f"Events       : {len(events)}")
    logger.info(f"Mean RMS     : {np.nanmean(rms_arr):.4f}")
    logger.info(f"Median RMS   : {np.nanmedian(rms_arr):.4f}")
    logger.info(f"Mean ERH     : {np.nanmean(erh_arr):.4f}")
    logger.info(f"Mean ERV     : {np.nanmean(erv_arr):.4f}")

    return {
        'output':   parameters.fileBulletin,
        'log':      log_path,
        'n_events': len(events),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Parse NLL hypo_71 output into a plain-text result bulletin.'
    )
    parser.add_argument('--loc-folder', required=True,
                        help='NLL loc output folder (e.g. loc/GLOBAL_1)')
    parser.add_argument('--obs-file',   required=True,
                        help='Base .obs filename, no path (e.g. GLOBAL_1.obs)')
    parser.add_argument('--output',     required=True,
                        help='Output bulletin text file (e.g. RESULT/GLOBAL_1.txt)')
    parser.add_argument('--log-dir',    default=None,
                        help='Log directory (default: NLL_run/console_output/)')
    args = parser.parse_args()

    write_events(
        CleanPostRunParams(
            folderLoc    = args.loc_folder,
            obsFile      = args.obs_file,
            fileBulletin = args.output,
        ),
        log_dir = args.log_dir,
    )


if __name__ == '__main__':
    main()
