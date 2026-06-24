"""
generate_magnitude_models.py
============================
Build piecewise ODR linear regression models to convert one earthquake
magnitude type to another, using matched events from two .obs catalogs.

A piecewise model is fitted with a breakpoint at M=2:
  - M ≥ 2 : unconstrained orthogonal regression (ODR)
  - M < 2 : constrained ODR, continuous at M=2

The fitted model is serialised with joblib and optionally figures are saved
to the directory specified in MagModelParams.save_figs.

Usage
-----
    python global_obs/generate_magnitude_models.py \\
        --file1      obs/RESIF_20-25.obs \\
        --file2      obs/LDG_20-25.obs  \\
        --mag-type1  MLv               \\
        --mag-type2  ML                \\
        --mag-name1  "MLv RESIF"       \\
        --mag-name2  "ML LDG"          \\
        --save-name  mag_model/MLv_RESIF.joblib \\
        --save-figs  mag_model/FIGURES/

    # --dist-thresh (km, default 10) and --time-thresh (s, default 2) are optional
"""

import argparse
import logging
import math
import os
from dataclasses import dataclass
from datetime import datetime as dt

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.odr import Model, ODR, RealData
from scipy.optimize import minimize
from scipy.spatial import KDTree
from sklearn.metrics import r2_score


logger = logging.getLogger('global_obs.generate_magnitude_models')

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
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MagModelParams:
    """
    Configuration for building a magnitude conversion model.

    Attributes
    ----------
    file_name1   : str   — .obs file whose magnitude type is being converted
    file_name2   : str   — .obs file providing the target magnitude type
    mag_type1    : str   — magnitude type token to search in file_name1 headers
    mag_type2    : str   — magnitude type token to search in file_name2 headers
    mag_name1    : str   — human-readable label for the source magnitude
    mag_name2    : str   — human-readable label for the target magnitude
    dist_thresh  : float — max epicentral distance (km) for event matching
    time_thresh  : float — max origin-time difference (s) for event matching
    save_name    : str   — path for the serialised joblib model
    save_figs    : str   — directory path for output PNG figures
    """
    file_name1:  str
    file_name2:  str
    mag_type1:   str
    mag_type2:   str
    mag_name1:   str
    mag_name2:   str
    dist_thresh: float
    time_thresh: float
    save_name:   str
    save_figs:   str


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def haversine(lat1, lon1, lat2, lon2):
    """Return the great-circle distance in km between two lat/lon points."""
    R    = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a    = (math.sin(dlat / 2) ** 2
            + math.cos(math.radians(lat1))
            * math.cos(math.radians(lat2))
            * math.sin(dlon / 2) ** 2)
    return 2 * R * math.asin(math.sqrt(a))


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def _retrieve_events(file_name, mag_type):
    """
    Read an .obs file and return event header lines that contain mag_type.

    Parameters
    ----------
    file_name : str — path to the .obs bulletin file
    mag_type  : str — magnitude type token (e.g. 'MLv', 'ML', 'mb_Lg')

    Returns
    -------
    list[str]
        Event header lines (stripped of the leading '# ').
    """
    with open(file_name, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    event_lines = [
        line.rstrip('\n').lstrip('# ')
        for line in lines
        if line.startswith('#') and not line.startswith('###') and mag_type in line
    ]
    logger.info(f"{len(event_lines)} event(s) from {file_name!r} (type '{mag_type}')")
    return event_lines


def _get_catalog_frame(event_lines):
    """
    Build a DataFrame of event coordinates, times, and magnitudes.

    Parameters
    ----------
    event_lines : list[str] — stripped event header strings from _retrieve_events

    Returns
    -------
    pd.DataFrame
        Columns: latitude, longitude, time (datetime64), magnitude.
    """
    rows = []
    for line in event_lines:
        p = line.split()
        dt = pd.to_datetime(f"{p[0]}-{p[1]}-{p[2]}T{p[3]}:{p[4]}:{p[5]}Z")
        rows.append({
            'latitude':  float(p[6]),
            'longitude': float(p[7]),
            'time':      dt,
            'magnitude': float(p[9]),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def _match_events(catalog1, catalog2, dist_thresh, time_thresh):
    """
    One-to-one event matching between two catalogs.

    For each event in catalog1, find the nearest event in catalog2 that
    falls within dist_thresh (km) and time_thresh (s). Each event is
    matched at most once (greedy closest-first).

    Parameters
    ----------
    catalog1     : pd.DataFrame — events to match (from the source catalog)
    catalog2     : pd.DataFrame — events to match against (target catalog)
    dist_thresh  : float        — maximum epicentral distance in km
    time_thresh  : float        — maximum origin-time difference in seconds

    Returns
    -------
    matched_frame : pd.DataFrame
        Columns: catalog1_idx, catalog2_idx, distance_km, time_diff_seconds,
                 magnitude1, magnitude2.
    matched_idx1  : set — indices used from catalog1
    matched_idx2  : set — indices used from catalog2
    """
    coords2 = catalog2[['latitude', 'longitude']].values
    tree    = KDTree(coords2)

    matched_pairs  = []
    matched_idx1   = set()
    matched_idx2   = set()

    for idx1, row in catalog1.iterrows():
        if idx1 in matched_idx1:
            continue

        _, candidates = tree.query([row['latitude'], row['longitude']], k=100)

        best_idx  = None
        best_dist = float('inf')
        best_dt   = float('inf')

        for i in candidates:
            if i in matched_idx2:
                continue
            cand = catalog2.iloc[i]
            dist = haversine(row['latitude'], row['longitude'],
                             cand['latitude'], cand['longitude'])
            if dist > dist_thresh:
                continue
            dt = abs((row['time'] - cand['time']).total_seconds())
            if dt > time_thresh:
                continue
            if dist < best_dist or (dist == best_dist and dt < best_dt):
                best_idx  = i
                best_dist = dist
                best_dt   = dt

        if best_idx is not None:
            matched_idx1.add(idx1)
            matched_idx2.add(best_idx)
            matched_pairs.append({
                'catalog1_idx':     idx1,
                'catalog2_idx':     best_idx,
                'distance_km':      best_dist,
                'time_diff_seconds': best_dt,
                'magnitude1':       row['magnitude'],
                'magnitude2':       catalog2.iloc[best_idx]['magnitude'],
            })

    return pd.DataFrame(matched_pairs), matched_idx1, matched_idx2


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------

def _linear_func(p, x):
    """ODR linear model: y = slope * x + intercept."""
    slope, intercept = p
    return slope * x + intercept


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _plot_regression(matched_frame, parameters,
                     slope_geq_2, intercept_geq_2,
                     slope_lt_2,  intercept_lt_2,
                     R2_geq_2, R2_lt_2):
    """Save the scatter + piecewise regression figure."""
    sns.set_theme(style='ticks')
    palette = sns.color_palette('deep')

    fig, ax = plt.subplots(figsize=(9, 7))

    label_geq = f'{parameters.mag_name1} ≥ 2'
    label_lt  = f'{parameters.mag_name1} < 2'
    matched_frame['_group'] = matched_frame['magnitude1'].apply(
        lambda x: label_geq if x >= 2 else label_lt
    )

    # Scatter
    group_colors = {label_geq: palette[0], label_lt: palette[1]}
    for label, grp in matched_frame.groupby('_group'):
        ax.scatter(grp['magnitude1'], grp['magnitude2'],
                   s=15, alpha=0.45, color=group_colors[label],
                   linewidths=0, label=label, zorder=1)

    # Regression lines + confidence bands (±1 residual std)
    x_lt  = np.linspace(matched_frame['magnitude1'].min(), 2, 60)
    x_geq = np.linspace(2, matched_frame['magnitude1'].max(), 60)
    y_lt  = slope_lt_2  * x_lt  + intercept_lt_2
    y_geq = slope_geq_2 * x_geq + intercept_geq_2

    grp_lt  = matched_frame[matched_frame['magnitude1'] < 2]
    grp_geq = matched_frame[matched_frame['magnitude1'] >= 2]
    std_lt  = (grp_lt['magnitude2']  - (slope_lt_2  * grp_lt['magnitude1']  + intercept_lt_2)).std()
    std_geq = (grp_geq['magnitude2'] - (slope_geq_2 * grp_geq['magnitude1'] + intercept_geq_2)).std()

    ax.plot(x_lt, y_lt, color=palette[1], lw=2, zorder=3)
    ax.fill_between(x_lt,  y_lt  - std_lt,  y_lt  + std_lt,
                    color=palette[1], alpha=0.18, zorder=2)

    ax.plot(x_geq, y_geq, color=palette[0], lw=2, zorder=3)
    ax.fill_between(x_geq, y_geq - std_geq, y_geq + std_geq,
                    color=palette[0], alpha=0.18, zorder=2)

    # Breakpoint line
    ax.axvline(x=2, color='grey', lw=1, ls='--', alpha=0.6)

    # Equation annotations — bottom-left corner, stacked, colour-coded
    ax.text(
        0.03, 0.03,
        f'y = {slope_lt_2:.3f}x + {intercept_lt_2:.3f}\n$R^2$ = {R2_lt_2:.3f}',
        transform=ax.transAxes,
        fontsize=9, color=palette[1], va='bottom',
    )
    ax.text(
        0.03, 0.17,
        f'y = {slope_geq_2:.3f}x + {intercept_geq_2:.3f}\n$R^2$ = {R2_geq_2:.3f}',
        transform=ax.transAxes,
        fontsize=9, color=palette[0], va='bottom',
    )

    ax.set_xlabel(f'{parameters.mag_name1}', fontsize=12)
    ax.set_ylabel(f'{parameters.mag_name2}', fontsize=12)
    ax.set_title(
        f'Piecewise ODR regression: {parameters.mag_name1} → {parameters.mag_name2}',
        fontsize=13, pad=10,
    )
    ax.legend(markerscale=2, framealpha=0.6)
    sns.despine(fig)
    plt.tight_layout()

    import os
    os.makedirs(parameters.save_figs, exist_ok=True)
    out = os.path.join(parameters.save_figs,
                       f'{parameters.mag_name1}_2_{parameters.mag_name2}.png')
    plt.savefig(out, dpi=300)
    plt.close(fig)


def _plot_residuals(matched_frame, parameters,
                    slope_geq_2, intercept_geq_2,
                    slope_lt_2,  intercept_lt_2,
                    R2_geq_2, R2_lt_2, BIC):
    """Save the residuals figure with outlier highlighting and a stats box."""
    sns.set_theme(style='ticks')
    palette = sns.color_palette('deep')

    label_geq = f'{parameters.mag_name1} ≥ 2'
    label_lt  = f'{parameters.mag_name1} < 2'
    matched_frame['_group'] = matched_frame['magnitude1'].apply(
        lambda x: label_geq if x >= 2 else label_lt
    )

    fig, ax = plt.subplots(figsize=(9, 6))

    # Inliers
    inliers = matched_frame[~matched_frame['is_outlier_iqr']]
    group_colors = {label_geq: palette[0], label_lt: palette[1]}
    for label, grp in inliers.groupby('_group'):
        ax.scatter(grp['magnitude1'], grp['residual'],
                   s=12, alpha=0.45, color=group_colors[label],
                   linewidths=0, label=label, zorder=1)

    # Outliers (diamond marker, drawn on top)
    outliers = matched_frame[matched_frame['is_outlier_iqr']]
    ax.scatter(outliers['magnitude1'], outliers['residual'],
               s=35, marker='D', color=palette[3], linewidths=0.4,
               edgecolors='white', label=f'Outliers (n={len(outliers)})', zorder=3)

    # Zero line
    ax.axhline(0, color='black', lw=1, ls='--', alpha=0.4, zorder=0)

    # Stats annotation box
    stats_text = (
        f'$R^2$ (M≥2) = {R2_geq_2:.3f}\n'
        f'$R^2$ (M<2) = {R2_lt_2:.3f}\n'
        f'BIC = {BIC:.1f}\n'
        f'n = {len(matched_frame)}'
    )
    ax.text(0.02, 0.97, stats_text,
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor='lightgrey', alpha=0.85))

    ax.set_xlabel(f'{parameters.mag_name1}', fontsize=12)
    ax.set_ylabel('Residual (observed − predicted)', fontsize=12)
    ax.set_title(
        f'Residuals: {parameters.mag_name1} → {parameters.mag_name2}',
        fontsize=13, pad=10,
    )
    ax.legend(markerscale=2, framealpha=0.6)
    sns.despine(fig)
    plt.tight_layout()

    import os
    out = os.path.join(parameters.save_figs,
                       f'Residuals_{parameters.mag_name1}_2_{parameters.mag_name2}.png')
    plt.savefig(out, dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def convert_magnitudes(parameters, save_figs=False, log_dir=None):
    """
    Build and save a piecewise ODR linear regression magnitude conversion model.

    Parameters
    ----------
    parameters : MagModelParams
        Full configuration (files, magnitude types, thresholds, output paths).
    save_figs  : bool, optional
        If True, save scatter and residual figures to parameters.save_figs.
    log_dir    : str, optional — log directory; default: global_obs/console_output/

    Returns
    -------
    dict or None
        Summary dict on success (keys: output, n_matched, r2_geq_2, r2_lt_2,
        slope_geq_2, intercept_geq_2, slope_lt_2, intercept_lt_2), or None if
        matching or data requirements are not met.
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file    : {log_path}")
    logger.info(f"Model       : {parameters.mag_name1} → {parameters.mag_name2}")
    logger.info(f"File 1      : {parameters.file_name1}  (type: {parameters.mag_type1})")
    logger.info(f"File 2      : {parameters.file_name2}  (type: {parameters.mag_type2})")
    logger.info(f"Thresholds  : dist={parameters.dist_thresh} km  time={parameters.time_thresh} s")

    # --- Match events from both catalogs ---
    events1   = _retrieve_events(parameters.file_name1, parameters.mag_type1)
    events2   = _retrieve_events(parameters.file_name2, parameters.mag_type2)
    catalog1  = _get_catalog_frame(events1)
    catalog2  = _get_catalog_frame(events2)

    matched_frame, match_idx1, match_idx2 = _match_events(
        catalog1, catalog2,
        dist_thresh=parameters.dist_thresh,
        time_thresh=parameters.time_thresh,
    )

    if len(matched_frame) == 0:
        logger.warning("No matched events — check thresholds or input files.")
        return None

    logger.info(f"Matched        : {len(matched_frame)} event(s)")
    logger.info(f"Unmatched in 1 : {len(catalog1) - len(match_idx1)}/{len(catalog1)}")
    logger.info(f"Unmatched in 2 : {len(catalog2) - len(match_idx2)}/{len(catalog2)}")

    if len(matched_frame) < 100:
        logger.warning(f"Only {len(matched_frame)} matched events — cannot build a reliable model (min 100).")
        return None

    # --- Split at M=2 ---
    label_geq = f'{parameters.mag_name1} ≥ 2'
    label_lt  = f'{parameters.mag_name1} < 2'
    matched_frame['magnitude_group'] = matched_frame['magnitude1'].apply(
        lambda x: label_geq if x >= 2 else label_lt
    )
    groups    = matched_frame.groupby('magnitude_group')
    grp_geq   = groups.get_group(label_geq)
    grp_lt    = groups.get_group(label_lt)

    # --- ODR for M >= 2 ---
    odr_geq = ODR(
        RealData(grp_geq['magnitude1'].values, grp_geq['magnitude2'].values),
        Model(_linear_func), beta0=[1., 0.],
    )
    out_geq          = odr_geq.run()
    slope_geq_2, intercept_geq_2 = out_geq.beta
    y_at_2           = slope_geq_2 * 2 + intercept_geq_2

    # --- Constrained ODR for M < 2 (continuous at M=2) ---
    X_lt = grp_lt['magnitude1'].values
    y_lt = grp_lt['magnitude2'].values

    def _objective(params):
        s, b = params
        return np.sum((y_lt - (b + s * X_lt)) ** 2)

    def _continuity(params):
        s, b = params
        return b + s * 2 - y_at_2

    result = minimize(
        _objective, x0=[1., 0.],
        constraints={'type': 'eq', 'fun': _continuity},
        method='SLSQP',
    )
    slope_lt_2, intercept_lt_2 = result.x

    # --- Predicted values and residuals ---
    matched_frame['predicted'] = np.where(
        matched_frame['magnitude1'] >= 2,
        slope_geq_2 * matched_frame['magnitude1'] + intercept_geq_2,
        slope_lt_2  * matched_frame['magnitude1'] + intercept_lt_2,
    )
    matched_frame['residual'] = matched_frame['magnitude2'] - matched_frame['predicted']

    # --- Statistics ---
    R2_geq_2 = r2_score(grp_geq['magnitude2'],
                        slope_geq_2 * grp_geq['magnitude1'] + intercept_geq_2)
    R2_lt_2  = r2_score(grp_lt['magnitude2'],
                        slope_lt_2  * grp_lt['magnitude1']  + intercept_lt_2)

    n        = len(matched_frame)
    SSR      = np.sum(matched_frame['residual'] ** 2)
    sigma2   = SSR / n
    BIC      = 4 * np.log(n) + n * np.log(2 * np.pi * sigma2) + n

    Q1, Q3   = matched_frame['residual'].quantile([0.25, 0.75])
    IQR      = Q3 - Q1
    matched_frame['is_outlier_iqr'] = (
        (matched_frame['residual'] < Q1 - 1.5 * IQR) |
        (matched_frame['residual'] > Q3 + 1.5 * IQR)
    )
    n_outliers = matched_frame['is_outlier_iqr'].sum()

    logger.info(f"R² (M≥2)   : {R2_geq_2:.3f}")
    logger.info(f"R² (M<2)   : {R2_lt_2:.3f}  (constrained)")
    logger.info(f"BIC        : {BIC:.3f}")
    logger.info(f"Outliers   : {n_outliers}/{n}  (IQR method)")

    # --- Figures ---
    if save_figs:
        _plot_regression(
            matched_frame, parameters,
            slope_geq_2, intercept_geq_2,
            slope_lt_2,  intercept_lt_2,
            R2_geq_2, R2_lt_2,
        )
        _plot_residuals(
            matched_frame, parameters,
            slope_geq_2, intercept_geq_2,
            slope_lt_2,  intercept_lt_2,
            R2_geq_2, R2_lt_2, BIC,
        )

    # --- Save model ---
    models = {
        label_geq: {'slope': slope_geq_2, 'intercept': intercept_geq_2},
        label_lt:  {'slope': slope_lt_2,  'intercept': intercept_lt_2},
    }
    joblib.dump(models, parameters.save_name)

    logger.info(f"Model saved : {parameters.save_name}")
    logger.info(f"  {label_geq}: y = {slope_geq_2:.3f}x + {intercept_geq_2:.3f}")
    logger.info(f"  {label_lt}:  y = {slope_lt_2:.3f}x + {intercept_lt_2:.3f}")

    return {
        'output':          parameters.save_name,
        'n_matched':       n,
        'r2_geq_2':        R2_geq_2,
        'r2_lt_2':         R2_lt_2,
        'slope_geq_2':     slope_geq_2,
        'intercept_geq_2': intercept_geq_2,
        'slope_lt_2':      slope_lt_2,
        'intercept_lt_2':  intercept_lt_2,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Build piecewise ODR magnitude conversion models from two .obs catalogs.'
    )
    parser.add_argument('--file1',       required=True,      help='.obs file with source magnitude type')
    parser.add_argument('--file2',       required=True,      help='.obs file with target magnitude type')
    parser.add_argument('--mag-type1',   required=True,      help='Magnitude type token in file1 (e.g. MLv)')
    parser.add_argument('--mag-type2',   required=True,      help='Magnitude type token in file2 (e.g. ML)')
    parser.add_argument('--mag-name1',   required=True,      help='Human-readable source label (e.g. "MLv RESIF")')
    parser.add_argument('--mag-name2',   required=True,      help='Human-readable target label (e.g. "ML LDG")')
    parser.add_argument('--dist-thresh', type=float, default=10.0, help='Max matching distance in km (default: 10)')
    parser.add_argument('--time-thresh', type=float, default=2.0,  help='Max matching time diff in s (default: 2)')
    parser.add_argument('--save-name',   required=True,      help='Output joblib model path')
    parser.add_argument('--save-figs',   default='mag_model/FIGURES/', help='Output figures directory')
    args = parser.parse_args()

    params = MagModelParams(
        file_name1  = args.file1,
        file_name2  = args.file2,
        mag_type1   = args.mag_type1,
        mag_type2   = args.mag_type2,
        mag_name1   = args.mag_name1,
        mag_name2   = args.mag_name2,
        dist_thresh = args.dist_thresh,
        time_thresh = args.time_thresh,
        save_name   = args.save_name,
        save_figs   = args.save_figs,
    )
    convert_magnitudes(params, save_figs=True)


if __name__ == '__main__':
    main()
