"""
match_pre_post_relocation.py
============================
Match pre-NLL .obs events to their NonLinLoc-relocated counterparts via publicId.

Reads obs/GLOBAL.obs and the merged NLL result CSV (RESULT/FINAL.csv).
For each event in the CSV, the corresponding .obs event is found by its
publicId, its header is updated with NLL location and uncertainty parameters,
and its magnitude and picks are preserved unchanged.

Events with no NLL solution (publicId absent from FINAL.csv) are silently dropped.

Usage
-----
    python NLL_run/match_pre_post_relocation.py \\
        --obs    obs/GLOBAL.obs \\
        --final  RESULT/FINAL.csv \\
        --output obs/FINAL.obs
"""

import argparse
import logging
import os
from dataclasses import dataclass
from datetime import datetime

import pandas as pd

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

logger = logging.getLogger('match_pre_post_relocation')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MatchCatalogsParams:
    file_obs:   str
    file_final: str   # path to RESULT/FINAL.csv
    save_file:  str


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

def _build_public_id_index(lines):
    """
    Map each publicId to the line index of its event header.

    Parameters
    ----------
    lines : list[str] — raw bulletin lines

    Returns
    -------
    dict[str, int] — publicId value → index of the '# ' header line
    """
    idx = {}
    for i, line in enumerate(lines):
        if line.startswith('# ') and i + 1 < len(lines):
            next_line = lines[i + 1]
            if next_line.startswith('PUBLIC_ID'):
                pid = next_line.split(None, 1)[1].strip()
                idx[pid] = i
    return idx


def _update_bulletin(lines, final_df, public_id_index):
    """
    Rebuild the bulletin using NLL-relocated parameters.

    For each event in final_df, the corresponding .obs block is found by
    publicId, the header line is rewritten with NLL coordinates and
    uncertainties, and the PUBLIC_ID line plus all pick lines are propagated
    unchanged. Events not present in final_df are dropped.

    Parameters
    ----------
    lines           : list[str]    — raw .obs bulletin lines
    final_df        : pd.DataFrame — merged NLL CSV (one row per event)
    public_id_index : dict[str, int] — publicId → header line index in lines

    Returns
    -------
    (list[str], int, int) — updated lines, n_matched, n_dropped
    """
    output    = []
    n_matched = 0
    n_dropped = 0

    # Preserve file header (### lines + following blank line)
    i = 0
    while i < len(lines) and lines[i].startswith('###'):
        output.append(lines[i])
        i += 1
    while i < len(lines) and lines[i].strip() == '':
        output.append(lines[i])
        i += 1

    for _, row in final_df.iterrows():
        pid = str(row['publicId'])
        if pid not in public_id_index:
            logger.warning(f"publicId {pid!r} not found in obs bulletin — skipping")
            n_dropped += 1
            continue

        header_idx = public_id_index[pid]

        # Extract Mag / MagType / MagAuthor from the original obs header
        orig_parts = lines[header_idx].lstrip('# ').split()
        mag        = orig_parts[9]
        mag_type   = orig_parts[10]
        mag_author = orig_parts[11]

        # Parse NLL datetime (full ISO — no ambiguity)
        dt       = pd.to_datetime(row['date-time'])
        sec_frac = dt.second + dt.microsecond / 1e6

        # Build updated event header with NLL coords + original magnitude
        new_header = (
            f"# {dt.year} {dt.month} {dt.day} {dt.hour} {dt.minute} {sec_frac:.2f} "
            f"{row['latitude']:.4f} {row['longitude']:.4f} {row['depth']:.2f} "
            f"{mag} {mag_type} {mag_author} "
            f"{int(row['Nphs'])} {row['errH']:.2f} {row['errZ']:.2f} "
            f"{row['Gap']:.1f} {row['RMS']:.4f}\n"
        )
        output.append(new_header)

        # Propagate PUBLIC_ID line and all pick lines unchanged
        j = header_idx + 1
        while j < len(lines) and lines[j].strip() != '':
            output.append(lines[j])
            j += 1
        output.append('\n')
        n_matched += 1

    return output, n_matched, n_dropped


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def save_bulletin(parameters, log_dir=None):
    """
    Update pre-NLL .obs headers with NLL locations and write the result.

    Parameters
    ----------
    parameters : MatchCatalogsParams
    log_dir    : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: output, log, n_matched
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Obs file   : {parameters.file_obs}")
    logger.info(f"Final CSV  : {parameters.file_final}")
    logger.info(f"Output     : {parameters.save_file}")

    with open(parameters.file_obs, 'r') as f:
        lines = f.readlines()

    public_id_index = _build_public_id_index(lines)
    logger.info(f"Events with publicId in obs: {len(public_id_index)}")

    final_df = pd.read_csv(parameters.file_final)
    logger.info(f"Events in NLL CSV          : {len(final_df)}")

    updated, n_matched, n_dropped = _update_bulletin(lines, final_df, public_id_index)

    parent = os.path.dirname(parameters.save_file)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(parameters.save_file, 'w') as f:
        f.writelines(updated)

    logger.info(f"Matched  : {n_matched}")
    logger.info(f"Dropped  : {n_dropped} (no NLL solution)")

    return {
        'output':    parameters.save_file,
        'log':       log_path,
        'n_matched': n_matched,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Update pre-NLL .obs headers with NLL-relocated parameters via publicId.'
    )
    parser.add_argument('--obs',     required=True, help='Pre-relocation GLOBAL.obs bulletin')
    parser.add_argument('--final',   required=True, help='Merged NLL result CSV (RESULT/FINAL.csv)')
    parser.add_argument('--output',  required=True, help='Output updated bulletin file')
    parser.add_argument('--log-dir', default=None,
                        help='Log directory (default: NLL_run/console_output/)')
    args = parser.parse_args()

    save_bulletin(MatchCatalogsParams(
        file_obs   = args.obs,
        file_final = args.final,
        save_file  = args.output,
    ))


if __name__ == '__main__':
    main()
