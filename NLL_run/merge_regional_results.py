"""
merge_regional_results.py
============================
Merge NonLinLoc CSV output files from all geographic zones into one.

For events relocated in multiple overlapping zones (same publicId in several
zone CSVs), the solution with the smallest pdfVolume is kept — a smaller
probability-density volume indicates a tighter, better-constrained location.

Zones can be supplied in any order; no adjacency constraint is needed.

The reported hypocenter is the location-PDF *expectation*, not the maximum-
likelihood point: NLLoc runs with `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`, which
re-solves origin time, RMS, residuals, Gap and Dist at the expectation. That is
what puts the published location and the published error on the same point —
the confidence ellipsoid this module turns into true_erh/true_erz is a second
moment *about* the expectation, so quoting it for the mode was quoting an
ellipsoid that is not centred on the point it belongs to. The maximum-likelihood
solution is preserved in the maxlike_* columns, read from the summary .hyp.

Usage
-----
    python NLL_run/merge_regional_results.py \\
        run/nll_loc/GLOBAL_1/GLOBAL_1.obs.sum.grid0.loc.csv \\
        run/nll_loc/GLOBAL_2/GLOBAL_2.obs.sum.grid0.loc.csv ... \\
        -o RESULT/NLL_result.csv
"""

import argparse
import logging
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import chi2 as chi2dist

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
    The semi-axes are used as-is; no chi-squared rescaling is applied here.
    NLLoc's reported axis lengths already carry its own 3-DOF, 68% chi-square
    scaling (_S3_3DOF) — callers must divide it out before reapplying a
    DOF-appropriate factor.
    """
    v1 = _ellipsoid_axis_to_xyz(az1, dip1, len1)
    v2 = _ellipsoid_axis_to_xyz(az2, dip2, len2)
    v3_dir = np.cross(v1 / len1, v2 / len2)
    v3 = v3_dir / np.linalg.norm(v3_dir) * len3
    R = np.column_stack([v1, v2, v3])
    return R @ R.T

# One-sigma chi-square factors used to convert between NLLoc's 3-DOF ellipsoid
# scaling and the DOF-appropriate scaling for each derived quantity.
#
# The coverage is erf(1/sqrt(2)) = 0.6827, i.e. exactly one sigma, not a flat
# 0.68: the two differ by 1.1% in chi-square at 1 DOF (0.9889 vs 1.0000), and
# 0.68 is the value that makes the 1-DOF factor miss unity.
_NOMINAL_COVERAGE = 0.6827            # erf(1/sqrt(2)), one-sigma coverage
_S1_1DOF = chi2dist.ppf(_NOMINAL_COVERAGE, df=1)   # ERZ: 1-DOF marginal std dev
_S2_2DOF = chi2dist.ppf(_NOMINAL_COVERAGE, df=2)   # ERH: 2-DOF horizontal error ellipse

# NLLoc's own ellipsoid-axis scaling. Not a free choice — this factor exists only
# to be divided back out, so it must be the literal the binary applied:
# `#define DELTA_CHI_SQR_68_3 3.53` (matrix_statistics.h:31, quoted from Numerical
# Recipes 2nd ed. sec 15.6), passed to CalcErrorEllipsoid at NLLocLib.c:404. That
# is the one-sigma table value rounded to 3 figures; chi2.ppf(0.6827, df=3) gives
# 3.5267, which is right in principle but 0.09% off what actually scaled the axes.
_S3_3DOF = 3.53

def _compute_true_erz(az1, dip1, len1, az2, dip2, len2, len3):
    """1-DOF, 68% confidence vertical standard deviation (km)."""
    C = _build_covariance(az1, dip1, len1, az2, dip2, len2, len3)
    sigma_zz = C[2, 2] / _S3_3DOF
    return float(np.sqrt(sigma_zz * _S1_1DOF))

def _compute_true_erh(az1, dip1, len1, az2, dip2, len2, len3):
    """2-DOF horizontal error ellipse, reduced via geometric mean of semi-axes (km)."""
    C = _build_covariance(az1, dip1, len1, az2, dip2, len2, len3)
    cov2d = C[:2, :2] / _S3_3DOF
    eigvals = np.linalg.eigvalsh(cov2d)   # ascending: [minor^2, major^2], unscaled
    b, a = np.sqrt(eigvals * _S2_2DOF)    # b = minor semi-axis, a = major semi-axis
    return float(np.sqrt(a * b))


# ---------------------------------------------------------------------------
# Maximum-likelihood hypocenter — recovered from the summary .hyp
# ---------------------------------------------------------------------------

# Both stages run NLLoc with LOCHYPOUT ... SAVE_NLLOC_EXPECTATION, so the
# reported hypocenter is the PDF expectation and the CSV's latitude/longitude/
# depth are identical to its expect_lat/expect_lon/expect_z. The maximum-
# likelihood point is not in the CSV at all: NLLoc writes it only to the .hyp,
# as `MAXIMUM_LIKELIHOOD  MaxLikeLat .. Long .. Depth .. OT ..` (GridLib.c:3511).
# It is recovered here so the final bulletin can publish both estimates.

def _read_maxlike(csv_path):
    """
    Map publicId -> (lat, lon, depth_km, ot_sec) from the summary .hyp beside a zone CSV.

    The summary .hyp sits next to the summary CSV under the same root and holds
    one block per event with no PHASE lines, so a single sequential pass is enough.

    `ot_sec` is seconds within the minute only — NLLoc calls hypotime2hrminsec
    before overwriting the origin time with the expectation's, so the hour and
    minute of the ML solution are not written out. Callers reconstruct the full
    timestamp against the expectation origin (see export_quakeml._maxlike_time).

    Raises FileNotFoundError / ValueError when the run was not made in expectation
    mode — silence there would let a mis-configured campaign run for days.
    """
    hyp_path = csv_path[:-len('.csv')] + '.hyp'
    if not os.path.isfile(hyp_path):
        raise FileNotFoundError(
            f'{hyp_path} not found — the maximum-likelihood hypocenter is read from '
            f'the summary .hyp beside each zone CSV')

    maxlike   = {}
    public_id = None
    with open(hyp_path) as fh:
        for line in fh:
            if line.startswith('PUBLIC_ID'):
                public_id = line.split()[1]
            elif line.startswith('MAXIMUM_LIKELIHOOD'):
                fields = line.split()
                maxlike[public_id] = (float(fields[2]), float(fields[4]),
                                      float(fields[6]), float(fields[8]))

    if not maxlike:
        raise ValueError(
            f'{hyp_path} holds no MAXIMUM_LIKELIHOOD line — that line is written only '
            f'when NLLoc runs with `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`, which is '
            f'also what makes the reported hypocenter the PDF expectation')

    return maxlike


_MAXLIKE_COLUMNS = ['maxlike_latitude', 'maxlike_longitude',
                    'maxlike_depth', 'maxlike_ot_sec']


def _attach_maxlike(df, csv_path):
    """Add the four maxlike_* columns to one zone frame, in place."""
    maxlike = _read_maxlike(csv_path)
    values  = [maxlike.get(pid, (np.nan,) * 4) for pid in df['publicId']]
    df[_MAXLIKE_COLUMNS] = pd.DataFrame(values, index=df.index)

    n_missing = int(df['maxlike_latitude'].isna().sum())
    if n_missing:
        logger.warning(f'{csv_path}: {n_missing} events with no MAXIMUM_LIKELIHOOD block')


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
        _attach_maxlike(df, path)
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

    # Compute true ERH (2-DOF) / ERZ (1-DOF) from the 3D confidence ellipsoid axes
    _ell_args = ['EllipsoidAz1', 'EllipsoidDip1', 'EllipsoidLen1',
                 'EllipsoidAz2', 'EllipsoidDip2', 'EllipsoidLen2',
                 'EllipsoidLen3']
    merged['true_erh'] = merged.apply(
        lambda r: _compute_true_erh(*r[_ell_args]), axis=1
    )
    merged['true_erz'] = merged.apply(
        lambda r: _compute_true_erz(*r[_ell_args]), axis=1
    )

    # Volume of the confidence ellipsoid, for comparison against pdfVolume (NLLoc's
    # OCT-tree-integrated PDF volume) — the two diverge for non-Gaussian events.
    merged['ellipsoidVolume'] = (
        4 / 3 * np.pi
        * merged['EllipsoidLen1'] * merged['EllipsoidLen2'] * merged['EllipsoidLen3']
    )

    merged = merged.sort_values('date-time').rename(columns={'_source': 'source'})

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
    parser.add_argument("-o", "--output", default="RESULT/NLL_result.csv",
                        help="Output CSV file (default: RESULT/NLL_result.csv).")
    parser.add_argument("--log-dir", default=None,
                        help="Log directory (default: NLL_run/console_output/).")
    args = parser.parse_args()

    if len(args.csv_files) < 2:
        print("ERROR: Please supply at least 2 CSV files.", file=sys.stderr)
        sys.exit(1)

    merge_bulletins(args.csv_files, args.output, args.log_dir)


if __name__ == "__main__":
    main()
