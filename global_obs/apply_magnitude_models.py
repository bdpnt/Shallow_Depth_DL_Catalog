"""
apply_magnitude_models.py
============================
Apply pre-computed piecewise regression models to convert all event magnitudes
in a set of .obs bulletins to ML LDG.

Models are loaded from mag_model/<MagType> <Author>.joblib files.

Usage
-----
    python global_obs/apply_magnitude_models.py --folder-path "obs/*.obs"
"""

import argparse
import glob
import logging
import os
import sys
from dataclasses import dataclass
from datetime import datetime as dt

import joblib
import pandas as pd


logger = logging.getLogger('global_obs.apply_magnitude_models')

_DEFAULT_LOG_DIR = 'global_obs/console_output/'


def _setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    basename  = os.path.splitext(os.path.basename(__file__))[0]
    timestamp = dt.now().strftime('%Y%m%d_%H%M%S')
    log_path  = os.path.join(log_dir, f"{basename}_{timestamp}.log")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(log_path, encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(handler)
    return log_path


# ---------------------------------------------------------------------------
# Path constants
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class UpdateMagFilesParams:
    """
    Configuration for applying magnitude conversion models.

    Attributes
    ----------
    folder_path : str — glob pattern for the .obs files to update (e.g. 'obs/*.obs')
    """
    folder_path: str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fetch_events(file):
    """
    Read an .obs file and return all lines together with the indices of event header lines.

    Parameters
    ----------
    file : str — path to the .obs file

    Returns
    -------
    (lines, event_ids) : (list of str, list of int)
    """
    with open(file, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    event_ids = [
        idx for idx, line in enumerate(lines)
        if line.startswith('#') and not line.startswith('###')
    ]
    return lines, event_ids


def _update_magnitudes(lines, org_mag):
    """
    Replace magnitude values in event header lines with converted ML LDG values.

    Parameters
    ----------
    lines   : list of str
    org_mag : pd.DataFrame with columns Mag, MagType, MagAuthor, linesID, ldgMag

    Returns
    -------
    list of str
    """
    for _, row in org_mag.iterrows():
        line_idx = row['linesID']
        columns  = lines[line_idx].split()
        columns[10] = f"{float(row['ldgMag']):.2f}" if row['ldgMag'] is not None else 'None'
        columns[11] = 'ML'
        columns[12] = 'LDG'
        if columns[10] != 'None':
            lines[line_idx] = ' '.join(columns)
    return lines


def _save_magnitudes(lines, file):
    """
    Write updated bulletin lines back to an .obs file.

    Parameters
    ----------
    lines : list of str
    file  : str — output path
    """
    with open(file, 'w') as f:
        for line in lines:
            if not line.endswith('\n'):
                line += '\n'
            f.write(line)

    logger.info(f"Saved: {file}")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def apply_magnitude_models(parameters, log_dir=None):
    """
    Apply all saved magnitude conversion models to every .obs file matching the folder pattern.

    Parameters
    ----------
    parameters : UpdateMagFilesParams
    log_dir    : str, optional — log directory; default: global_obs/console_output/

    Returns
    -------
    dict
        'output'      — glob pattern used
        'n_converted' — total number of magnitudes converted
        'n_total'     — total number of events processed
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file       : {log_path}")
    logger.info(f"Bulletin glob  : {parameters.folder_path}")

    total_converted = 0
    total_events    = 0

    for folder_file in glob.glob(parameters.folder_path):
        folder_author = os.path.basename(folder_file).split('_')[0]

        lines, lines_id = _fetch_events(folder_file)

        org_mag = pd.DataFrame(
            [line.split()[10:13] for line in [lines[idx] for idx in lines_id]],
            columns=['Mag', 'MagType', 'MagAuthor'],
        )
        org_mag['Mag']     = pd.to_numeric(org_mag['Mag'])
        org_mag['linesID'] = lines_id
        org_mag['ldgMag']  = None

        n_total     = len(org_mag)
        n_converted = 0

        logger.info(f"File: {folder_file}  ({n_total} events)")
        for mag_type in org_mag.MagType.unique():
            current_ids     = org_mag[org_mag.MagType == mag_type].index
            model_name_stem = f'{mag_type} {folder_author}'

            try:
                file_model  = f'mag_model/{model_name_stem}.joblib'
                models      = joblib.load(file_model)
                model_ge_2  = list(models.values())[0]
                model_lt_2  = list(models.values())[1]
                logger.info(f"  Model found   : {model_name_stem} ({file_model})")
            except Exception:
                logger.warning(f"  Model missing : {model_name_stem} ({file_model}) — skipped")
                continue

            mask_ge_2 = org_mag.index.isin(current_ids) & (org_mag['Mag'] >= 2)
            mask_lt_2 = org_mag.index.isin(current_ids) & (org_mag['Mag'] < 2)
            org_mag.loc[mask_ge_2, 'ldgMag'] = (
                model_ge_2['slope'] * org_mag.loc[mask_ge_2, 'Mag'] + model_ge_2['intercept']
            )
            org_mag.loc[mask_lt_2, 'ldgMag'] = (
                model_lt_2['slope'] * org_mag.loc[mask_lt_2, 'Mag'] + model_lt_2['intercept']
            )

            n_converted += len(current_ids)
            logger.info(f"  Converted     : {len(current_ids)} {mag_type} magnitude(s) → ML LDG")

        if not org_mag['ldgMag'].isna().all():
            lines = _update_magnitudes(lines, org_mag)
            _save_magnitudes(lines, folder_file)
            logger.info(f"  Total converted: {n_converted}/{n_total}")
        else:
            logger.warning(f"  No magnitudes converted for {folder_file}")

        total_converted += n_converted
        total_events    += n_total

    logger.info(f"Overall: {total_converted}/{total_events} magnitudes converted across all files")
    return {
        'output':      parameters.folder_path,
        'n_converted': total_converted,
        'n_total':     total_events,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Apply pre-computed magnitude models to convert .obs bulletin magnitudes to ML LDG.'
    )
    parser.add_argument('--folder-path', required=True,
                        help='Glob pattern for .obs files (e.g. "obs/*.obs")')
    args = parser.parse_args()

    params = UpdateMagFilesParams(folder_path=args.folder_path)
    apply_magnitude_models(params)


if __name__ == '__main__':
    main()
