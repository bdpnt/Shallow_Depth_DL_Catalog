"""
merge_regional_results.py
============================
Merge NonLinLoc CSV output files from all geographic zones into one.

For events relocated in multiple overlapping zones (same publicId in several
zone CSVs), the solution with the smallest pdfVolume is kept — a smaller
probability-density volume indicates a tighter, better-constrained location.

Zones can be supplied in any order; no adjacency constraint is needed.

Usage
-----
    python NLL_run/merge_regional_results.py \\
        loc/GLOBAL_1/GLOBAL_1.obs.sum.grid0.loc.csv \\
        loc/GLOBAL_2/GLOBAL_2.obs.sum.grid0.loc.csv ... \\
        -o RESULT/FINAL.csv
"""

import argparse
import logging
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
# from scipy.stats import chi2 as chi2dist

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

logger = logging.getLogger('merge_regional_results')


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
# Ellipsoid helpers — compute true ERH/ERZ from NLLoc confidence ellipsoid axes
# ---------------------------------------------------------------------------

def _ellipsoid_axis_to_xyz(az_deg, dip_deg, length):
    """Convert a NLLoc ellipsoid semi-axis to a 3D vector [East, North, Down]."""
    az  = np.radians(az_deg)
    dip = np.radians(dip_deg)
    return np.array([
        length * np.cos(dip) * np.sin(az),   # East
        length * np.cos(dip) * np.cos(az),   # North
        length * np.sin(dip),                 # Down
    ])

def _build_covariance(az1, dip1, len1, az2, dip2, len2, len3):
    """
    Build the covariance matrix of the NLLoc 3D confidence ellipsoid (R @ R.T).
    The semi-axes are used as-is; no chi-squared rescaling is applied.
    """
    v1 = _ellipsoid_axis_to_xyz(az1, dip1, len1)
    v2 = _ellipsoid_axis_to_xyz(az2, dip2, len2)
    v3_dir = np.cross(v1 / len1, v2 / len2)
    v3 = v3_dir / np.linalg.norm(v3_dir) * len3
    R = np.column_stack([v1, v2, v3])
    return R @ R.T

def _compute_true_erz(az1, dip1, len1, az2, dip2, len2, len3):
    """Maximum vertical extent of the NLLoc 3D confidence ellipsoid (km)."""
    C = _build_covariance(az1, dip1, len1, az2, dip2, len2, len3)
    return float(np.sqrt(C[2, 2]))

def _compute_true_erh(az1, dip1, len1, az2, dip2, len2, len3):
    """Maximum horizontal extent of the NLLoc 3D confidence ellipsoid (km)."""
    C = _build_covariance(az1, dip1, len1, az2, dip2, len2, len3)
    return float(np.sqrt(np.max(np.linalg.eigvalsh(C[:2, :2]))))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def merge_bulletins(csv_files, output_path, log_dir=None):
    """
    Load, deduplicate by publicId, and merge NLL CSV output files into one.

    For events appearing in multiple zones (same publicId), the solution with
    the smallest pdfVolume is kept (smallest volume = tightest location PDF).

    Parameters
    ----------
    csv_files   : list[str] — paths to NLL CSV summary files, any order
    output_path : str       — path for the merged output CSV
    log_dir     : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: output, log, n_merged, n_duplicates
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Files          : {len(csv_files)}")

    frames = []
    for path in csv_files:
        df = pd.read_csv(path, skipinitialspace=True)
        df['_source'] = os.path.basename(os.path.dirname(path))
        frames.append(df)
        logger.info(f"Loaded {len(df):>5d} events from {path!r}")

    all_events = pd.concat(frames, ignore_index=True)
    n_total    = len(all_events)
    logger.info(f"Total raw events   : {n_total}")

    # Log each duplicate group before resolving
    dup_ids = all_events[all_events.duplicated('publicId', keep=False)]['publicId'].unique()
    for pid in sorted(dup_ids):
        rows = all_events[all_events['publicId'] == pid].sort_values('pdfVolume')
        kept = rows.iloc[0]
        others = rows.iloc[1:]
        logger.info(
            f"DUP {pid}  kept={kept['_source']} (pdfVol={kept['pdfVolume']:.4f})"
            f"  dropped: "
            + ', '.join(
                f"{r['_source']} (pdfVol={r['pdfVolume']:.4f})"
                for _, r in others.iterrows()
            )
        )

    # Keep the row with the smallest pdfVolume per publicId
    best_idx = all_events.groupby('publicId')['pdfVolume'].idxmin()
    merged   = all_events.loc[best_idx].copy()
    n_dup    = n_total - len(merged)

    # Compute true ERH/ERZ from the 3D confidence ellipsoid axes
    _ell_args = ['EllipsoidAz1', 'EllipsoidDip1', 'EllipsoidLen1',
                 'EllipsoidAz2', 'EllipsoidDip2', 'EllipsoidLen2',
                 'EllipsoidLen3']
    merged['true_erh'] = merged.apply(
        lambda r: _compute_true_erh(*r[_ell_args]), axis=1
    )
    merged['true_erz'] = merged.apply(
        lambda r: _compute_true_erz(*r[_ell_args]), axis=1
    )

    merged = merged.sort_values('date-time').drop(columns='_source')

    parent = os.path.dirname(output_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    merged.to_csv(output_path, index=False)

    logger.info(f"Duplicates removed  : {n_dup}")
    logger.info(f"Events in merged CSV: {len(merged)}")
    logger.info(f"Output              : {output_path!r}")

    return {
        'output':       output_path,
        'log':          log_path,
        'n_merged':     len(merged),
        'n_duplicates': n_dup,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    """Parse CLI arguments and merge NLL CSV files."""
    parser = argparse.ArgumentParser(
        description="Merge NLL CSV output files and deduplicate by publicId."
    )
    parser.add_argument("csv_files", nargs="+", metavar="FILE",
                        help="NLL CSV summary files to merge (any order).")
    parser.add_argument("-o", "--output", default="RESULT/FINAL.csv",
                        help="Output CSV file (default: RESULT/FINAL.csv).")
    parser.add_argument("--log-dir", default=None,
                        help="Log directory (default: NLL_run/console_output/).")
    args = parser.parse_args()

    if len(args.csv_files) < 2:
        print("ERROR: Please supply at least 2 CSV files.", file=sys.stderr)
        sys.exit(1)

    merge_bulletins(args.csv_files, args.output, args.log_dir)


if __name__ == "__main__":
    main()
