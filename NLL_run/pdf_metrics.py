"""
pdf_metrics.py
============================
Per-event location-PDF quality metrics, computed from the NonLinLoc `.scat`
scatter clouds of the final SSST relocation and appended to a merged result CSV.

The metrics answer one question: does the Gaussian confidence ellipsoid NLLoc
reports actually represent the posterior it was derived from? They measure
internal *consistency*, not accuracy — velocity-model error is a bias none of
them can see.

  J        negentropy KL(p || N(mu, Sigma)) in nats; 0 for an exactly Gaussian PDF
  Psi      exp(-J), the effective-volume ratio V_eff(PDF) / V_eff(ellipsoid)
  C68      fraction of posterior mass inside the nominal one-sigma ellipsoid;
           > 0.6827 means the reported errors are conservative, < 0.6827
           over-confident
  dip_stat Hartigan dip statistic on the depth marginal (multimodality)
  dip_sep_km  width of Hartigan's modal interval, in km — how far apart the
           competing depth solutions are

`dip_stat` says only *whether* a second mode exists, never how far away it is:
the dip is a vertical sup-distance on the ECDF and is invariant under affine
rescaling of the depth axis, so modes 0.5 km apart and modes 50 km apart give
the same statistic and the same p-value. `dip_sep_km` carries the physical scale
the dip structurally cannot, which is why rejection needs both (see depth_dip).

J and C68 are meaningless without an n-matched Gaussian null: the k-NN entropy
estimator is biased in a sample-size-dependent way, and the C68 null spread is
NOT binomial (mu/Sigma are fitted on the very samples being tested, which
suppresses the scatter). `gaussian_null` simulates both per distinct n; the
resulting `J_null_p95` / `C68_sigma_n` columns are what make a later cut a
one-line comparison.

Once the columns are written, the same run also saves the diagnostic figures that
show how the metrics relate to the classical quality indicators (RMS, Nphs, Gap,
Dist, depth), how the published error (true_erh/true_erz) relates to both, and
how everything varies in space, in complem_figures/pdf_metrics/.

Environment: needs `diptest` (imported lazily, only by depth_dip) -> seisbench_env.
The figures additionally need matplotlib/seaborn/plotly, all imported lazily so a
missing plotting stack can never cost the metrics.

Usage
-----
    python NLL_run/pdf_metrics.py --run-name ssst_run1
    python NLL_run/pdf_metrics.py --run-name ssst_run1 --output /tmp/annotated.csv
    python NLL_run/pdf_metrics.py --run-name ssst_run1 --figures false
"""

import argparse
import glob
import itertools
import logging
import os
import re
import struct
import zlib
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.special import digamma, gammaln
from scipy.stats import chi2 as chi2dist
from scipy.stats import rankdata, spearmanr

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')
_DEFAULT_FIG_DIR = os.path.join(_PROJECT_ROOT, 'complem_figures', 'pdf_metrics')

logger = logging.getLogger('pdf_metrics')

# Columns this module writes. Listed so a re-run can drop them before
# recomputing, which is what keeps the CLI idempotent.
METRIC_COLUMNS = [
    'ellipsoidVolume', 'n_scat', 'J', 'Psi', 'C68',
    'J_null_p95', 'C68_sigma_n', 'dip_stat', 'dip_pval', 'dip_sep_km', 'dip_reject',
]

_KNN_K        = 5      # k for the Kozachenko-Leonenko entropy estimator
_NULL_SIMS    = 500    # Gaussian clouds simulated per distinct n
_DIP_COMMON_N = 400    # common subsample size for the dip test (< observed min n)
_DIP_ALPHA    = 0.05

# The modal interval must be wider than this multiple of the depth error the
# catalog publishes (true_erz) before multimodality is worth rejecting over.
# At 1.0 a mode sits a full sigma from the reported depth, i.e. at the edge of
# the quoted confidence. Raising it to 2.0 makes the criterion vacuous: the
# width/true_erz ratio tops out near 2.24 over this catalog.
_DIP_SEP_ERZ_FACTOR = 1.0

# Coverage the confidence ellipsoid is built for: erf(1/sqrt(2)) = 0.6827, exactly
# one sigma, not a flat 0.68 — the same constant merge_regional_results.py scales
# true_erh/true_erz with. It sets both the ellipsoid C68 is measured against and
# the value C68 is compared to, so the two cannot drift apart.
_NOMINAL_COVERAGE = 0.6827
_CHI2_68_3DOF = chi2dist.ppf(_NOMINAL_COVERAGE, df=3)
_GAUSS_ENTROPY_3D = 0.5 * 3 * np.log(2 * np.pi * np.e)   # 4.2568 nats


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PdfMetricsParams:
    result_csv: str
    ssst_root:  str
    run_name:   str
    zones:      list = field(default_factory=lambda: ['1', '2', '3', '4', '5', '6'])
    output:     str = None      # None -> rewrite result_csv in place
    log_dir:    str = None
    figures:    bool = True
    figure_dir: str = None      # None -> complem_figures/pdf_metrics/


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
# .scat reader
# ---------------------------------------------------------------------------

# NLLoc scatter file: int32 nSamples, 12 unused bytes, then nSamples records of
# 4 float32 (x, y, z, log-likelihood), x/y/z in the local Lambert km frame.
_SCAT_DTYPE = np.dtype([('x', '<f4'), ('y', '<f4'), ('z', '<f4'), ('pdf', '<f4')])
_SCAT_HEADER_BYTES = 16


def read_scat(path):
    """
    Read a NonLinLoc .scat file into a structured array (fields x, y, z, pdf).

    Drop-in replacement for obspy.io.nlloc.util.read_nlloc_scatter, which is
    unusable here: it does `np.fromfile(...)[4:]`, dropping the first four
    16-byte *records* where the header is a single 16-byte record — silently
    discarding the first 3 samples of every file.

    Samples are unweighted: NLLoc draws them with density proportional to the
    posterior, so `pdf` is carried for API compatibility but is not a weight.

    Raises ValueError if the declared sample count and the file size disagree.
    """
    with open(path, 'rb') as fh:
        buf = fh.read()

    if len(buf) < _SCAT_HEADER_BYTES:
        raise ValueError(f'{path}: shorter than the {_SCAT_HEADER_BYTES}-byte header')

    n_samples = struct.unpack('<i', buf[:4])[0]
    body = buf[_SCAT_HEADER_BYTES:]
    if n_samples < 0 or len(body) != n_samples * _SCAT_DTYPE.itemsize:
        raise ValueError(
            f'{path}: header declares {n_samples} samples but the body holds '
            f'{len(body) / _SCAT_DTYPE.itemsize:.2f}')

    return np.frombuffer(body, dtype=_SCAT_DTYPE, count=n_samples)


def scat_xyz(scat):
    """Return the (n, 3) float64 coordinate array of a structured scatter array."""
    return np.column_stack([scat['x'], scat['y'], scat['z']]).astype(np.float64)


# ---------------------------------------------------------------------------
# Iteration-directory resolution (shared with complem_figures/)
# ---------------------------------------------------------------------------

def find_iteration_dirs(ssst_root, run_name, zone):
    """Return sorted (step, dir) pairs for every loc_ssst_corr<N>/GLOBAL_<zone> folder."""
    pattern = os.path.join(ssst_root, run_name, f'Pyrenees_{zone}_SSST',
                           'loc_ssst_corr*', f'GLOBAL_{zone}')
    step_dirs = []
    for path in glob.glob(pattern):
        match = re.search(r'loc_ssst_corr(\d+)', path)
        step_dirs.append((int(match.group(1)), path))
    return sorted(step_dirs)


def final_iteration_dir(ssst_root, run_name, zone):
    """Return the highest-numbered loc_ssst_corr<N>/GLOBAL_<zone> folder, or None."""
    step_dirs = find_iteration_dirs(ssst_root, run_name, zone)
    return step_dirs[-1][1] if step_dirs else None


def sibling(hyp_path, ext):
    """Swap the .hyp extension of a per-event NLLoc output path for `ext`."""
    return hyp_path[:-len('.hyp')] + ext


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def whiten(xyz):
    """
    Centre and whiten samples by their own Sigma^(-1/2), so cov(out) = I.

    The reference Gaussian for both J and C68 is the one matched to the cloud's
    own mean and covariance — exactly the distribution NLLoc's confidence
    ellipsoid draws — so whitening reduces both metrics to their standard-normal
    form. Returns None if the covariance is singular (degenerate cloud).

    bias=True (ddof=0) is deliberate: NLLoc derives the STATISTICS CovXX..ZZ line
    — and hence the published ellipsoid — from these same scatter samples with
    that normalisation. Verified against 20 events' .hyp files, where ddof=1 left
    a constant n/(n-1) offset (0.21% at n=479) and ddof=0 leaves none. C68 must
    test the ellipsoid that was actually published.
    """
    centred = xyz - xyz.mean(axis=0)
    cov = np.cov(centred, rowvar=False, bias=True)
    eigvals, eigvecs = np.linalg.eigh(cov)
    if np.any(eigvals <= 0) or not np.all(np.isfinite(eigvals)):
        return None
    return centred @ (eigvecs * eigvals ** -0.5) @ eigvecs.T


def negentropy(white, k=_KNN_K):
    """
    Negentropy J = KL(p || N(mu, Sigma)) in nats, from pre-whitened samples.

    Kozachenko-Leonenko k-NN differential entropy estimator:
        H = -psi(k) + psi(n) + log(c_d) + (d/n) * sum(log(eps_i))
    with eps_i the distance to the i-th sample's k-th neighbour and c_d the
    volume of the unit d-ball. After whitening the reference is the standard
    normal, so J = 0.5*d*log(2*pi*e) - H.

    Duplicate samples would give log(0); they are jittered far below the NLLoc
    grid spacing rather than dropped, so n (and hence the matching null) is
    unchanged.
    """
    n, d = white.shape
    eps = _knn_radius(white, k)
    if eps is None:
        return np.nan

    log_unit_ball = (d / 2) * np.log(np.pi) - gammaln(d / 2 + 1)
    entropy = (-digamma(k) + digamma(n) + log_unit_ball
               + (d / n) * np.sum(np.log(eps)))
    return float(0.5 * d * np.log(2 * np.pi * np.e) - entropy)


def _knn_radius(white, k, _rng=np.random.default_rng(0)):
    """Distance to each sample's k-th neighbour, jittering duplicates away from 0."""
    for attempt in range(3):
        distances, _ = cKDTree(white).query(white, k=k + 1)
        eps = distances[:, k]
        if np.all(eps > 0):
            return eps
        white = white + _rng.normal(scale=1e-6, size=white.shape)
    return None


def coverage_68(white):
    """
    Fraction of samples inside the nominal one-sigma Gaussian confidence ellipsoid.

    For whitened samples the squared Mahalanobis radius is just |w|^2, so this
    is the fraction with |w|^2 <= chi2.ppf(0.6827, 3) = 3.5267. Equals 0.6827 for
    a Gaussian PDF; above means the reported ERH/ERZ are conservative, below means
    they are over-confident.
    """
    return float(np.mean(np.sum(white ** 2, axis=1) <= _CHI2_68_3DOF))


_null_cache = {}


def gaussian_null(n, n_sim=_NULL_SIMS, seed=0):
    """
    Simulate the n-matched Gaussian null for J and C68.

    Returns (J_p95, C68_sigma): the 95th percentile of J and the standard
    deviation of C68 over `n_sim` genuinely Gaussian clouds of size n, run
    through the identical whitening + estimator pipeline.

    This is not optional. The k-NN entropy estimator's bias depends on n, so a
    fixed J threshold would flag sparsely-sampled events for reasons unrelated
    to their PDFs. And the C68 spread is not binomial: sqrt(.68*.32/n)
    overestimates it (0.021 vs ~0.012 at n=479) because mu and Sigma are
    estimated from the same samples being tested.

    Cached per n — the catalog holds only ~82 distinct sample counts.
    """
    key = (n, n_sim, seed)
    if key in _null_cache:
        return _null_cache[key]

    rng = np.random.default_rng(seed + n)
    j_values = np.empty(n_sim)
    c_values = np.empty(n_sim)
    for i in range(n_sim):
        white = whiten(rng.standard_normal((n, 3)))
        j_values[i] = negentropy(white)
        c_values[i] = coverage_68(white)

    result = (float(np.nanpercentile(j_values, 95)), float(np.nanstd(c_values, ddof=1)))
    _null_cache[key] = result
    return result


def depth_dip(depths, public_id, common_n=_DIP_COMMON_N):
    """
    Hartigan's dip test on the depth marginal, at a catalog-wide common n.

    Returns (statistic, p-value, modal-interval width in km).

    A small p-value means the depth PDF has more than one mode, but that alone
    is not grounds for rejection: the dip is a vertical distance on the ECDF and
    is invariant under affine rescaling of x, so it cannot tell modes 0.5 km
    apart from modes 50 km apart. Where both modes sit inside the quoted depth
    error, `depth +/- true_erz` already covers them and is an honest statement.

    The third return value supplies the scale the statistic is blind to.
    Hartigan's algorithm reports the interval over which the ECDF is being
    flattened — the region carrying the departure from unimodality — and its
    width is in the units of the input, i.e. km. `compute_metrics` combines the
    two into the reject flag; neither one decides on its own.

    The dip test's power grows with n, so p-values from different sample sizes
    are not comparable. Every event is therefore subsampled to `common_n`
    (below the smallest observed scatter count) using an RNG seeded from its
    publicId — deterministic across runs and independent between events.
    """
    from diptest import diptest      # lazy: only this function needs seisbench_env

    if len(depths) < common_n:
        return np.nan, np.nan, np.nan

    rng = np.random.default_rng(zlib.crc32(public_id.encode()))
    subsample = depths[rng.choice(len(depths), size=common_n, replace=False)]
    # `xl`/`xu` are the interval's values; `lo`/`hi` would be indices into
    # diptest's own sorted copy of the input (sort_x defaults to True).
    stat, pval, res = diptest(np.asarray(subsample, dtype=np.float64), full_output=True)
    return float(stat), float(pval), float(res['xu'] - res['xl'])


def _blank_metrics():
    """All-NaN metrics for an event whose cloud could not be evaluated.

    Excludes `ellipsoidVolume` (computed from the CSV, not the cloud) and
    `dip_reject` (derived later, once true_erz is in reach).
    """
    return {col: np.nan for col in METRIC_COLUMNS
            if col not in ('ellipsoidVolume', 'dip_reject')}


def event_metrics(scat_path, public_id):
    """
    Compute every per-event metric for one .scat file.

    `dip_reject` is deliberately not produced here: the rule needs `true_erz`,
    which lives in the result CSV and not in the scatter cloud, so it is derived
    in `compute_metrics` once the two are side by side. This keeps the function a
    pure map from one .scat to its own metrics.

    On failure every metric is NaN, and the NaN in `dip_pval` is what marks the
    event as unevaluated downstream.
    """
    blank = _blank_metrics()

    try:
        scat = read_scat(scat_path)
    except (OSError, ValueError) as exc:
        logger.warning(f'{public_id}: unreadable .scat ({exc})')
        return blank

    xyz = scat_xyz(scat)
    n_samples = len(xyz)
    white = whiten(xyz)
    if white is None:
        logger.warning(f'{public_id}: singular covariance over {n_samples} samples')
        return {**blank, 'n_scat': n_samples}

    j_value = negentropy(white)
    c68 = coverage_68(white)
    j_null_p95, c68_sigma = gaussian_null(n_samples)
    dip_stat, dip_pval, dip_sep = depth_dip(xyz[:, 2], public_id)

    return {
        'n_scat':      n_samples,
        'J':           j_value,
        'Psi':         float(np.exp(-j_value)) if np.isfinite(j_value) else np.nan,
        'C68':         c68,
        'J_null_p95':  j_null_p95,
        'C68_sigma_n': c68_sigma,
        'dip_stat':    dip_stat,
        'dip_pval':    dip_pval,
        'dip_sep_km':  dip_sep,
    }


# ---------------------------------------------------------------------------
# Figures — shared gridding
# ---------------------------------------------------------------------------

# Whole-Pyrenees domain, matching complem_figures/error_maps.py, depth_maps.py
# and event_ranking.py.
_LAT_RANGE = (42.0, 44.0)
_LON_RANGE = (-2.25, 3.5)

_MAP_BIN_SIZE    = 0.02          # degrees, 2-D gridmap cell
_MAP_WINDOW      = 4             # cells of smoothing added on each side
_MAP_MIN_COUNT   = 10            # events below which a cell is left blank
_SLICE_WINDOW    = 8             # depth slices hold ~1/4 of the events: smooth wider
_N_DEPTH_SLICES  = 4
_VOXEL_BIN_SIZE  = (0.1, 0.1, 2.5)   # lat deg, lon deg, depth km
_VOXEL_WINDOW    = 1
_VOXEL_MIN_COUNT = 5
_N_CURVE_BINS    = 20            # equal-count bins behind each median curve
_CLIP_PCT        = (0.5, 99.5)   # x-axis clip of the 1-D panels
_REJECT_RATE_FLOOR = 5.0         # % — keeps a near-zero dip row from scaling up to noise


def windowed_stat_grid(coords, values, ranges, bin_sizes,
                       window_size=_MAP_WINDOW, stat=np.median):
    """
    Windowed per-cell statistic of `values` over a regular N-D grid.

    Each cell's statistic is taken over every sample inside the cell *plus*
    `window_size` cells of padding in each direction (clipped at the domain
    edges), which is what smooths a sparse catalog into a readable map.
    Generic in the number of dimensions: 2 for the maps, 3 for the voxels.

    Parameters
    ----------
    coords    : sequence of per-dimension coordinate arrays (all same length)
    values    : array of the quantity to reduce, aligned with `coords`
    ranges    : sequence of (low, high) per dimension
    bin_sizes : sequence of cell sizes per dimension (same units as `ranges`)
    stat      : reducer applied to each window (np.median, np.mean, ...)

    Returns
    -------
    (edges, statistic, count)
        edges     : list of per-dimension bin-edge arrays
        statistic : N-D array, NaN where the window is empty
        count     : N-D int array of samples per window

    Non-finite samples are dropped up front, so `count` is the number of
    samples actually behind the statistic.
    """
    values = np.asarray(values, dtype=float)
    coords = [np.asarray(c, dtype=float) for c in coords]

    finite = np.isfinite(values)
    for coord in coords:
        finite &= np.isfinite(coord)
    values = values[finite]
    coords = [coord[finite] for coord in coords]

    edges = []
    for (low, high), bin_size in zip(ranges, bin_sizes):
        n_bins = max(int(round((high - low) / bin_size)), 1)
        edges.append(np.linspace(low, high, n_bins + 1))

    # The window of a cell only depends on its index along one axis, so the
    # per-axis membership masks are built once and combined per cell — the
    # alternative (rebuilding them inside the cell loop) repeats the same
    # comparisons thousands of times.
    axis_masks = []
    for axis, edge in enumerate(edges):
        step = edge[1] - edge[0]
        masks = []
        for index in range(len(edge) - 1):
            low  = max(edge[index]     - window_size * step, edge[0])
            high = min(edge[index + 1] + window_size * step, edge[-1])
            masks.append((coords[axis] >= low) & (coords[axis] <= high))
        axis_masks.append(masks)

    shape = tuple(len(edge) - 1 for edge in edges)
    statistic = np.full(shape, np.nan)
    count = np.zeros(shape, dtype=int)

    for cell in itertools.product(*(range(size) for size in shape)):
        mask = axis_masks[0][cell[0]]
        for axis in range(1, len(cell)):
            mask = mask & axis_masks[axis][cell[axis]]
        window = values[mask]
        if len(window):
            statistic[cell] = stat(window)
            count[cell] = len(window)

    return edges, statistic, count


# ---------------------------------------------------------------------------
# Figures — what gets plotted
# ---------------------------------------------------------------------------

# 1-D panels: one row per metric, one column per quality indicator. Latitude and
# longitude are deliberately absent — a 1-D latitude panel marginalizes over
# longitude and hides exactly the structure the maps below are there to show.
_METRIC_ROWS = [
    ('Psi',      r'$\Psi$'),
    ('C68',      r'$C_{68}$'),
    ('dip_stat', 'Dip statistic'),
]

_PREDICTORS = [
    ('depth', 'Depth (km)',                     'linear'),
    ('RMS',   'RMS (s)',                        'log'),     # a 4912 s outlier lives here
    ('Nphs',  'Phase count',                    'linear'),
    ('Gap',   'Azimuthal gap (°)',              'linear'),
    ('Dist',  'Nearest station distance (km)',  'linear'),
]

# The inverse view: the error the catalog actually publishes, against every
# indicator that might predict it — the classical four plus the PDF metrics.
#
# Psi rather than J: Psi = exp(-J) is monotone in J, so the rank correlation is
# identical up to sign, and Psi is bounded (0, 1]. Raw C68 rather than the
# null-normalized C68_z the maps use: here C68 is an x-axis, not a cell median
# mixing sample counts, and the raw fraction carries a readable reference at
# 0.6827 (the two differ by rho = 0.001 over this catalog).
_ERROR_ROWS = [
    ('true_erh', 'ERH (km)'),
    ('true_erz', 'ERZ (km)'),
]

_ERROR_PREDICTORS = [
    ('RMS',  'RMS (s)',                'log'),
    ('Gap',  'Azimuthal gap (°)',      'linear'),
    ('Nphs', 'Phase count',            'linear'),
    ('Dist', 'Nearest station distance (km)', 'linear'),
    ('Psi',  r'$\Psi$',                'linear'),
    ('C68',  r'$C_{68}$',              'linear'),
]

# Held fixed in the partial correlation of every _ERROR_PREDICTORS panel.
_ERROR_CONTROL = 'Nphs'

# Map panels. C68 is mapped as its null-normalized z-score and not as the raw
# fraction: a cell median mixes events with different sample counts, and only
# the z-score puts them on a common scale (see _add_null_columns).
_MAP_METRICS = [
    {'column': 'Psi',        'stat': np.median, 'scale': 1.0,
     'label': r'median $\Psi$',
     'cmap': 'viridis', 'vmin': 0.0, 'vmax': 1.0,
     'plotly_cmap': 'Viridis', 'plotly_reverse': False},
    {'column': 'C68_z',      'stat': np.median, 'scale': 1.0,
     'label': r'median $(C_{68}-0.6827)\,/\,\sigma_n$',
     'cmap': 'RdBu', 'vmin': -10.0, 'vmax': 10.0,
     'plotly_cmap': 'RdBu', 'plotly_reverse': False},
    {'column': 'dip_reject', 'stat': np.mean,   'scale': 100.0,
     'label': r'reject rate (%): bimodal & sep > $ERZ$',
     'cmap': 'magma_r', 'vmin': 0.0, 'vmax': None,
     'plotly_cmap': 'Magma', 'plotly_reverse': True},
]


def _add_null_columns(df):
    """
    Working copy carrying the null-normalized quantities the maps colour by.

    J and C68 are only readable against the Gaussian null simulated at the
    event's own sample count, so the maps never plot them raw:

      J_ratio = J / J_null_p95            > 1 -> non-Gaussian beyond the null
      C68_z   = (C68 - 0.6827) / C68_sigma_n  < 0 -> over-confident, in null sigmas

    Both nulls come from the event's own row, never from a catalog-wide
    constant: n_scat happens to be near-constant in the current run, but a
    different LOCSEARCH setting or catalog will not be that uniform.

    Kept out of the CSV on purpose — METRIC_COLUMNS is a documented contract.
    """
    data = df.copy()
    data['J_ratio'] = data['J'] / data['J_null_p95']
    data['C68_z'] = (data['C68'] - _NOMINAL_COVERAGE) / data['C68_sigma_n']
    return data


def _quantile_bins(x, n_bins=_N_CURVE_BINS):
    """
    Equal-count bin assignment for x, with the bin centres.

    Quantile edges rather than equal-width ones: Nphs, Dist and RMS are all
    strongly skewed, and equal-width bins would put nearly every event in the
    first bin and one event in the last.

    Returns (bin_index, centres), or (None, None) if x cannot be split.
    """
    edges = np.unique(np.nanquantile(x, np.linspace(0, 1, n_bins + 1)))
    if len(edges) < 3:
        return None, None
    index = np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2)
    return index, 0.5 * (edges[:-1] + edges[1:])


def _reject_rate_ceiling(data):
    """
    Shared upper limit (%) for the reject-rate twin axis of the dip row.

    One cap across the whole row, so the five panels stay comparable with each
    other rather than each auto-scaling to its own range. Data-driven because
    the rate follows the reject rule: it is ~1.4% catalog-wide under the
    two-condition rule, which a fixed 0-100 axis would flatten onto the floor.

    Mirrors the dropna/log/clip/bin sequence of _plot_panel so the cap really
    bounds what gets drawn.
    """
    peak = 0.0
    for predictor, _, scale in _PREDICTORS:
        panel = data[[predictor, 'dip_stat', 'dip_reject']].dropna(
            subset=[predictor, 'dip_stat'])
        if scale == 'log':
            panel = panel[panel[predictor] > 0]
        if panel.empty:
            continue
        x = panel[predictor].to_numpy(dtype=float)
        x_lo, x_hi = np.nanpercentile(x, _CLIP_PCT)
        keep = (x >= x_lo) & (x <= x_hi)
        bin_index, _ = _quantile_bins(x[keep])
        if bin_index is None:
            continue
        rates = pd.Series(
            panel['dip_reject'].to_numpy(dtype=float)[keep]).groupby(bin_index).mean()
        peak = max(peak, 100 * float(rates.max()))
    return max(_REJECT_RATE_FLOOR, 1.15 * peak)


def _metric_ylim(metric, y):
    """y-axis limits for a 1-D panel: fixed for the bounded Psi, robust otherwise."""
    if metric == 'Psi':
        return 0.0, 1.05
    if metric == 'C68':
        return float(np.nanpercentile(y, 0.2)), float(min(1.0, np.nanpercentile(y, 99.8)))
    return 0.0, float(np.nanpercentile(y, 99.5))


def _plot_panel(ax, data, metric, predictor, log_x, twin_label, reject_ylim=100.0):
    """
    One metric-vs-indicator panel: density + median curve + the event-wise null.

    Density is a log-count hexbin (46 k events make a raw scatter a solid
    block), the blue curve is the median with its IQR band over equal-count
    bins of the predictor, and the red overlay is the Gaussian null of that
    same bin — computed from the J_null_p95 / C68_sigma_n of the events *in
    the bin*, so it bends by itself if n_scat correlates with the predictor.

    The dip row has no null line (the dip test is already n-controlled by the
    common subsample in depth_dip); it carries on a right-hand axis the per-bin
    rate of the two-condition reject rule — bimodal *and* the modes further
    apart than true_erz — scaled to `reject_ylim`, shared across the row.
    """
    columns = [predictor, metric, 'J_null_p95', 'C68_sigma_n', 'dip_reject']
    panel = data[columns].dropna(subset=[predictor, metric])
    if log_x:
        panel = panel[panel[predictor] > 0]
    if panel.empty:
        ax.text(0.5, 0.5, 'no data', transform=ax.transAxes, ha='center', va='center')
        return

    # Clip the predictor tails (RMS carries a 4912 s outlier) and keep everything
    # else — density, curve, rho — on that one visible window, so the median curve
    # cannot run off the edge of the panel it summarizes.
    x_lo, x_hi = np.nanpercentile(panel[predictor], _CLIP_PCT)
    panel = panel[(panel[predictor] >= x_lo) & (panel[predictor] <= x_hi)]

    x = panel[predictor].to_numpy(dtype=float)
    y = panel[metric].to_numpy(dtype=float)
    y_lo, y_hi = _metric_ylim(metric, y)
    extent = ((np.log10(x_lo), np.log10(x_hi)) if log_x else (x_lo, x_hi)) + (y_lo, y_hi)

    ax.set_facecolor('white')       # sparse hexes are invisible on seaborn grey
    ax.hexbin(x, y, gridsize=45, bins='log', cmap='Greys', mincnt=1,
              xscale='log' if log_x else 'linear', extent=extent, zorder=1)

    bin_index, centres = _quantile_bins(x)
    if bin_index is not None:
        binned = pd.DataFrame({
            'bin':    bin_index,
            'metric': y,
            'null_j': panel['J_null_p95'].to_numpy(dtype=float),
            'null_c': panel['C68_sigma_n'].to_numpy(dtype=float),
            'reject': panel['dip_reject'].to_numpy(dtype=float),
        }).groupby('bin')

        centre = centres[binned['metric'].median().index.to_numpy()]
        ax.fill_between(centre, binned['metric'].quantile(0.25), binned['metric'].quantile(0.75),
                        color='tab:blue', alpha=0.25, lw=0, zorder=3)
        ax.plot(centre, binned['metric'].median(), color='tab:blue', lw=2, zorder=4)

        if metric == 'Psi':
            # Psi of a genuinely Gaussian cloud of the same size: exp(-J_null_p95).
            ax.plot(centre, np.exp(-binned['null_j'].median()),
                    color='crimson', ls='--', lw=1.2, zorder=5)
        elif metric == 'C68':
            sigma = binned['null_c'].median().to_numpy()
            ax.axhline(_NOMINAL_COVERAGE, color='crimson', ls='--', lw=1.2, zorder=5)
            ax.fill_between(centre,
                            _NOMINAL_COVERAGE - 3 * sigma, _NOMINAL_COVERAGE + 3 * sigma,
                            color='crimson', alpha=0.15, lw=0, zorder=2)
        else:
            twin = ax.twinx()
            twin.plot(centre, 100 * binned['reject'].mean(),
                      color='crimson', ls='--', lw=1.2, zorder=5)
            twin.set_ylim(0, reject_ylim)
            twin.grid(False)
            if twin_label:
                twin.set_ylabel(r'reject rate (%): bimodal & sep > $ERZ$', color='crimson')
            else:
                twin.set_yticklabels([])

    # Exploratory screening only — a monotone association, not evidence that the
    # two quantities agree (same caveat as event_ranking.py's correlation matrix).
    rho = spearmanr(x, y).statistic
    ax.text(0.97, 0.97, f'$\\rho$ = {rho:+.2f}\nN = {len(x):,}',
            transform=ax.transAxes, ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.75, lw=0))

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)


def _generate_predictor_figure(data, run_name, output_path):
    """3 metrics x 5 quality indicators, one PDF."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme()
    fig, axes = plt.subplots(len(_METRIC_ROWS), len(_PREDICTORS),
                             figsize=(4.2 * len(_PREDICTORS), 3.4 * len(_METRIC_ROWS)),
                             layout='constrained', squeeze=False)

    reject_ylim = _reject_rate_ceiling(data)

    for row, (metric, ylabel) in enumerate(_METRIC_ROWS):
        for col, (predictor, xlabel, scale) in enumerate(_PREDICTORS):
            ax = axes[row][col]
            _plot_panel(ax, data, metric, predictor, scale == 'log',
                        twin_label=(col == len(_PREDICTORS) - 1),
                        reject_ylim=reject_ylim)
            if row == len(_METRIC_ROWS) - 1:
                ax.set_xlabel(xlabel)
            if col == 0:
                ax.set_ylabel(ylabel)

    fig.suptitle(f'Location-PDF quality metrics vs. classical quality indicators — {run_name}\n'
                 r'grey = event density (log counts) · blue = median + IQR over equal-count bins · '
                 'red = per-bin Gaussian null (dip row: bimodal & sep > ERZ reject rate, right axis)',
                 fontweight='bold')
    plt.savefig(output_path)
    plt.close(fig)


def _partial_spearman(x, y, control):
    """
    Spearman rho between x and y with `control` held fixed.

    Rank all three (ties averaged, as spearmanr does), regress the ranks of x
    and of y on the ranks of the control, and correlate the residuals.

    It exists because the marginal rho of RMS against ERH is +0.00 while the
    Nphs-controlled one is +0.62: RMS grows with the phase count and ERH falls
    with it, and over this catalog the two effects cancel almost exactly. A
    panel reporting the marginal number alone would read as "RMS says nothing
    about the location error", which is the opposite of what the data holds.
    """
    ranks = np.column_stack([rankdata(x), rankdata(y), rankdata(control)])
    design = np.column_stack([np.ones(len(ranks)), ranks[:, 2]])
    fitted = design @ np.linalg.lstsq(design, ranks[:, :2], rcond=None)[0]
    return float(np.corrcoef(ranks[:, :2] - fitted, rowvar=False)[0, 1])


def _error_ylim(values):
    """Log y-range of one error row, clipped to the _CLIP_PCT tails.

    Taken once over the whole column and shared by every panel of the row: the
    point of the row is to compare indicators against each other, which per-panel
    autoscaling would quietly defeat.
    """
    positive = values[np.isfinite(values) & (values > 0)]
    low, high = np.nanpercentile(positive, _CLIP_PCT)
    return float(low), float(high)


def _plot_error_panel(ax, data, error_col, predictor, log_x, y_lim):
    """
    One published-error-vs-indicator panel: density + median curve + rho box.

    Same grammar as _plot_panel — log-count hexbin, blue median with its IQR
    band over equal-count bins of the predictor — with two differences: the
    y-axis is logarithmic (ERH spans 0.05-64 km) and there is no red null
    overlay, because a published error has no Gaussian null to be compared to.

    The annotation carries the marginal rho *and* the Nphs-controlled one; see
    _partial_spearman for why the marginal value alone can be actively
    misleading. Both are exploratory screening — a monotone association, not
    evidence that the two quantities agree.
    """
    # dict.fromkeys: on the Nphs panel the predictor IS the control, and selecting
    # it twice gives a frame with duplicate labels that cannot be masked.
    columns = list(dict.fromkeys([predictor, error_col, _ERROR_CONTROL]))
    panel = data[columns].dropna(subset=[predictor, error_col])
    panel = panel[panel[error_col] > 0]          # log y; the errors are positive anyway
    if log_x:
        panel = panel[panel[predictor] > 0]
    if panel.empty:
        ax.text(0.5, 0.5, 'no data', transform=ax.transAxes, ha='center', va='center')
        return

    x_lo, x_hi = np.nanpercentile(panel[predictor], _CLIP_PCT)
    panel = panel[(panel[predictor] >= x_lo) & (panel[predictor] <= x_hi)]

    x = panel[predictor].to_numpy(dtype=float)
    y = panel[error_col].to_numpy(dtype=float)
    y_lo, y_hi = y_lim
    extent = (((np.log10(x_lo), np.log10(x_hi)) if log_x else (x_lo, x_hi))
              + (np.log10(y_lo), np.log10(y_hi)))

    ax.set_facecolor('white')       # sparse hexes are invisible on seaborn grey
    ax.hexbin(x, y, gridsize=45, bins='log', cmap='Greys', mincnt=1,
              xscale='log' if log_x else 'linear', yscale='log',
              extent=extent, zorder=1)

    bin_index, centres = _quantile_bins(x)
    if bin_index is not None:
        binned = pd.DataFrame({'bin': bin_index, 'error': y}).groupby('bin')['error']
        centre = centres[binned.median().index.to_numpy()]
        ax.fill_between(centre, binned.quantile(0.25), binned.quantile(0.75),
                        color='tab:blue', alpha=0.25, lw=0, zorder=3)
        ax.plot(centre, binned.median(), color='tab:blue', lw=2, zorder=4)

    lines = [f'$\\rho$ = {spearmanr(x, y).statistic:+.2f}']
    if predictor != _ERROR_CONTROL:     # a predictor cannot be controlled for itself
        control = panel[_ERROR_CONTROL].to_numpy(dtype=float)
        finite = np.isfinite(control)
        if finite.sum() > 2:
            rho_partial = _partial_spearman(x[finite], y[finite], control[finite])
            lines.append(f'$\\rho\\,|\\,${_ERROR_CONTROL} = {rho_partial:+.2f}')
    lines.append(f'N = {len(x):,}')
    ax.text(0.97, 0.97, '\n'.join(lines), transform=ax.transAxes,
            ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white', alpha=0.75, lw=0))

    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)


def _plot_error_box_panel(ax, data, error_col, y_lim):
    """
    The dip column of the error figure: published error split by dip_reject.

    dip_reject is binary and heavily unbalanced (661 rejected of 46 224), so a
    hexbin would draw two vertical stripes and a rho would compress a large
    median shift into a small number. Two boxes with their medians and counts
    show the actual effect: ERH 0.58 -> 1.61 km, ERZ 1.04 -> 3.72 km.

    Outliers are hidden — the whiskers already reach the row's clipped y-range,
    and 46 k flier points would black out the panel.
    """
    panel = data[[error_col, 'dip_reject']].dropna(subset=[error_col])
    panel = panel[panel[error_col] > 0]
    rejected = panel['dip_reject'].astype(bool)
    groups = [panel.loc[~rejected, error_col].to_numpy(dtype=float),
              panel.loc[rejected, error_col].to_numpy(dtype=float)]

    ax.set_facecolor('white')
    if not all(len(group) for group in groups):
        ax.text(0.5, 0.5, 'no data', transform=ax.transAxes, ha='center', va='center')
        return

    boxes = ax.boxplot(groups, positions=[0, 1], widths=0.5, showfliers=False,
                       patch_artist=True, medianprops=dict(color='black', lw=1.5),
                       zorder=3)
    for patch, colour in zip(boxes['boxes'], ('tab:blue', 'crimson')):
        patch.set_facecolor(colour)
        patch.set_alpha(0.45)

    for position, group in zip((0, 1), groups):
        ax.text(position + 0.3, np.median(group), f'{np.median(group):.2f} km',
                ha='left', va='center', fontsize=8)

    ax.set_yscale('log')
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f'kept\nN = {len(groups[0]):,}',
                        f'rejected\nN = {len(groups[1]):,}'])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(*y_lim)


def _generate_error_figure(data, run_name, output_path):
    """2 published-error rows x (6 indicators + the dip-reject split), one PDF."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme()
    n_cols = len(_ERROR_PREDICTORS) + 1      # + the dip-reject box panel
    fig, axes = plt.subplots(len(_ERROR_ROWS), n_cols,
                             figsize=(4.2 * n_cols, 3.6 * len(_ERROR_ROWS)),
                             layout='constrained', squeeze=False)

    for row, (error_col, ylabel) in enumerate(_ERROR_ROWS):
        last_row = row == len(_ERROR_ROWS) - 1
        y_lim = _error_ylim(data[error_col].to_numpy(dtype=float))

        for col, (predictor, xlabel, scale) in enumerate(_ERROR_PREDICTORS):
            ax = axes[row][col]
            _plot_error_panel(ax, data, error_col, predictor, scale == 'log', y_lim)
            if last_row:
                ax.set_xlabel(xlabel)
            if col == 0:
                ax.set_ylabel(ylabel)

        _plot_error_box_panel(axes[row][-1], data, error_col, y_lim)
        if last_row:
            axes[row][-1].set_xlabel('Dip-validity test')

    fig.suptitle(f'Published location error vs. quality indicators and PDF metrics — {run_name}\n'
                 r'grey = event density (log counts) · blue = median + IQR over equal-count bins · '
                 r'$\rho$ = Spearman, $\rho\,|\,$Nphs = the same with the phase count held fixed '
                 '(log error axis, shared across each row)',
                 fontweight='bold')
    plt.savefig(output_path)
    plt.close(fig)


def _draw_map_panel(ax, edges, statistic, count, spec, vmax, events, title):
    """One lon/lat panel of a windowed statistic; returns the QuadMesh."""
    masked = np.ma.masked_where(count < _MAP_MIN_COUNT, statistic)
    ax.set_facecolor('white')       # blanked cells read as white, not seaborn grey
    mesh = ax.pcolormesh(edges[1], edges[0], masked, cmap=spec['cmap'],
                         shading='auto', vmin=spec['vmin'], vmax=vmax)
    ax.scatter(events['longitude'], events['latitude'],
               s=0.3, color='black', linewidth=0, alpha=0.35, zorder=2)
    ax.set_xlim(*_LON_RANGE)
    ax.set_ylim(*_LAT_RANGE)
    ax.set_title(title, fontsize=10)
    return mesh


def _map_vmax(spec, statistic, count):
    """Colour-scale top: the spec's own, or the largest displayed cell."""
    if spec['vmax'] is not None:
        return spec['vmax']
    shown = statistic[count >= _MAP_MIN_COUNT]
    top = np.nanmax(shown) if shown.size and np.any(np.isfinite(shown)) else None
    return float(top) if top else 1.0


def _map_grid(data, spec, window_size):
    """Windowed lon/lat grid of one map metric over `data`."""
    return windowed_stat_grid(
        [data['latitude'], data['longitude']],
        data[spec['column']].astype(float) * spec['scale'],
        [_LAT_RANGE, _LON_RANGE], [_MAP_BIN_SIZE, _MAP_BIN_SIZE],
        window_size=window_size, stat=spec['stat'],
    )


def _generate_gridmap_figure(data, run_name, output_path):
    """Whole-depth lon/lat map of the three metrics, one row each."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme()
    fig, axes = plt.subplots(len(_MAP_METRICS), 1, figsize=(12, 5 * len(_MAP_METRICS)),
                             layout='constrained', squeeze=False)

    for row, spec in enumerate(_MAP_METRICS):
        ax = axes[row][0]
        edges, statistic, count = _map_grid(data, spec, _MAP_WINDOW)
        vmax = _map_vmax(spec, statistic, count)
        mesh = _draw_map_panel(ax, edges, statistic, count, spec, vmax, data, spec['label'])
        fig.colorbar(mesh, ax=ax, label=spec['label'], shrink=0.85, pad=0.02)
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')

    fig.suptitle(f'Location-PDF quality metrics across the Pyrenees — {run_name}\n'
                 f'windowed median over {_MAP_BIN_SIZE}° cells '
                 f'(±{_MAP_WINDOW} cells), cells with < {_MAP_MIN_COUNT} events blank',
                 fontweight='bold')
    plt.savefig(output_path)
    plt.close(fig)


def _depth_slices(depths, n_slices=_N_DEPTH_SLICES):
    """Depth-quartile edges, so the slicing adapts to whatever catalog is passed."""
    return np.unique(np.nanquantile(depths, np.linspace(0, 1, n_slices + 1)))


def _generate_depth_slice_figure(data, run_name, output_path):
    """The same three maps, one column per depth slice, colour scale shared per row."""
    import matplotlib.pyplot as plt
    import seaborn as sns

    edges_z = _depth_slices(data['depth'])
    slices = []
    for low, high in zip(edges_z[:-1], edges_z[1:]):
        in_slice = (data['depth'] >= low) & (data['depth'] <= high)
        slices.append((f'{low:.1f} – {high:.1f} km', data[in_slice]))

    sns.set_theme()
    fig, axes = plt.subplots(len(_MAP_METRICS), len(slices),
                             figsize=(5.5 * len(slices), 4.2 * len(_MAP_METRICS)),
                             layout='constrained', squeeze=False)

    for row, spec in enumerate(_MAP_METRICS):
        grids = [_map_grid(subset, spec, _SLICE_WINDOW) for _, subset in slices]
        # One scale per row, taken over every slice: otherwise each panel gets its
        # own top value and the slices stop being comparable, which is the point.
        vmax = max(_map_vmax(spec, statistic, count) for _, statistic, count in grids)
        mesh = None
        for col, ((label, subset), (grid_edges, statistic, count)) in enumerate(zip(slices, grids)):
            ax = axes[row][col]
            title = f'{label} — N = {len(subset):,}' if row == 0 else ''
            mesh = _draw_map_panel(ax, grid_edges, statistic, count, spec, vmax, subset, title)
            if row == len(_MAP_METRICS) - 1:
                ax.set_xlabel('Longitude')
            if col == 0:
                ax.set_ylabel('Latitude')
        fig.colorbar(mesh, ax=axes[row], label=spec['label'], shrink=0.85, pad=0.02)

    fig.suptitle(f'Location-PDF quality metrics by depth slice (depth quartiles) — {run_name}\n'
                 f'windowed median over {_MAP_BIN_SIZE}° cells (±{_SLICE_WINDOW} cells), '
                 f'colour scale shared across the slices of a row',
                 fontweight='bold')
    plt.savefig(output_path)
    plt.close(fig)


def _generate_voxel_html(data, spec, run_name, output_path):
    """Interactive lat/lon/depth voxel view of one metric (Plotly, self-contained)."""
    import plotly.graph_objects as go

    depth_range = (float(np.floor(np.nanpercentile(data['depth'], _CLIP_PCT[0]))),
                   float(np.ceil(np.nanpercentile(data['depth'], _CLIP_PCT[1]))))
    edges, statistic, count = windowed_stat_grid(
        [data['latitude'], data['longitude'], data['depth']],
        data[spec['column']].astype(float) * spec['scale'],
        [_LAT_RANGE, _LON_RANGE, depth_range], _VOXEL_BIN_SIZE,
        window_size=_VOXEL_WINDOW, stat=spec['stat'],
    )

    centres = [0.5 * (edge[:-1] + edge[1:]) for edge in edges]
    lat_i, lon_i, depth_i = np.nonzero(count >= _VOXEL_MIN_COUNT)
    values = statistic[lat_i, lon_i, depth_i]
    counts = count[lat_i, lon_i, depth_i]

    vmax = spec['vmax'] if spec['vmax'] is not None else float(np.nanmax(values))
    figure = go.Figure(go.Scatter3d(
        x=centres[1][lon_i], y=centres[0][lat_i], z=centres[2][depth_i],
        mode='markers',
        marker=dict(
            size=4 + 8 * np.sqrt(counts / counts.max()),
            color=values, colorscale=spec['plotly_cmap'], reversescale=spec['plotly_reverse'],
            cmin=spec['vmin'], cmax=vmax, opacity=0.6,
            colorbar=dict(title=spec['column']),
        ),
        customdata=np.column_stack([values, counts]),
        hovertemplate=('lon %{x:.2f}° · lat %{y:.2f}° · depth %{z:.1f} km<br>'
                       f"{spec['column']} " '%{customdata[0]:.3f} · %{customdata[1]} events'
                       '<extra></extra>'),
    ))
    figure.update_layout(
        title=dict(text=(f'{spec["label"]} in 3-D — {run_name}<br>'
                         f'<sub>{_VOXEL_BIN_SIZE[0]}° × {_VOXEL_BIN_SIZE[1]}° × '
                         f'{_VOXEL_BIN_SIZE[2]} km cells (±{_VOXEL_WINDOW} cell window), '
                         f'≥ {_VOXEL_MIN_COUNT} events, marker size ∝ √count</sub>')),
        scene=dict(
            xaxis_title='Longitude (°)',
            yaxis_title='Latitude (°)',
            zaxis=dict(title='Depth (km)', autorange='reversed'),
            aspectmode='manual', aspectratio=dict(x=2.0, y=1.0, z=0.6),
        ),
    )
    figure.write_html(output_path, include_plotlyjs=True)


def generate_figures(df, run_name, figure_dir=None):
    """
    Save every diagnostic figure of the annotated catalog.

    Parameters
    ----------
    df         : pd.DataFrame — the catalog as written by compute_metrics
    run_name   : str          — campaign name, used in titles and filenames
    figure_dir : str          — output folder (default complem_figures/pdf_metrics/)

    Returns
    -------
    list of the paths written. Takes ~30 s for ~46 k events.
    """
    figure_dir = figure_dir or _DEFAULT_FIG_DIR
    os.makedirs(figure_dir, exist_ok=True)
    data = _add_null_columns(df)

    outputs = []
    for name, builder in (('metrics_vs_quality', _generate_predictor_figure),
                          ('error_vs_quality',   _generate_error_figure),
                          ('gridmap',            _generate_gridmap_figure),
                          ('gridmap_depth',      _generate_depth_slice_figure)):
        path = os.path.join(figure_dir, f'{run_name}_{name}.pdf')
        builder(data, run_name, path)
        outputs.append(path)

    for spec in _MAP_METRICS:
        path = os.path.join(figure_dir, f'{run_name}_voxels_{spec["column"]}.html')
        _generate_voxel_html(data, spec, run_name, path)
        outputs.append(path)

    for path in outputs:
        logger.info(f'Figure           : {path!r}')
    print(f'Figures saved @ {figure_dir} ({len(outputs)} files)')
    return outputs


# ---------------------------------------------------------------------------
# Catalog-scale plumbing
# ---------------------------------------------------------------------------

def build_scat_index(ssst_root, run_name, zones):
    """
    Map (zone_source, publicId) -> .scat path over every zone's final iteration.

    The .scat/.hyp filenames carry the *input* event time, not the publicId, and
    the relocated origin time differs from it — so the join key has to come from
    the PUBLIC_ID line (line 2) of each .hyp. Read in bulk here rather than
    grepping per event, which is what plot_pdf_cloud.py does for its single-event
    case and would be hopeless across ~46 k events.
    """
    index = {}
    for zone in zones:
        iter_dir = final_iteration_dir(ssst_root, run_name, zone)
        if iter_dir is None:
            logger.warning(f'zone {zone}: no loc_ssst_corr* directory found')
            continue

        source = os.path.basename(iter_dir)      # GLOBAL_<zone>, matching the CSV
        n_zone = 0
        for hyp_path in glob.glob(os.path.join(iter_dir, '*.hyp')):
            if '.sum.' in os.path.basename(hyp_path):
                continue                         # concatenated summary, no .scat
            with open(hyp_path) as fh:
                fh.readline()                    # NLLOC ... line
                second = fh.readline()
            if not second.startswith('PUBLIC_ID'):
                continue
            scat_path = sibling(hyp_path, '.scat')
            if os.path.exists(scat_path):
                index[(source, second.split()[1])] = scat_path
                n_zone += 1
        logger.info(f'zone {zone}: indexed {n_zone} scatter clouds from {iter_dir!r}')

    return index


def compute_metrics(params):
    """
    Append the PDF-quality columns to a merged result CSV.

    Metrics are computed only for each event's winning zone: `source` already
    records which zone survived the lowest-pdfVolume dedup, so joining on
    (source, publicId) skips the ~7 600 losing solutions.

    Returns dict with keys: output, log, n_events, n_missing.
    """
    log_path = _setup_logger(params.log_dir or _DEFAULT_LOG_DIR)
    logger.info(f'Result CSV : {params.result_csv!r}')
    logger.info(f'Run name   : {params.run_name}')

    df = pd.read_csv(params.result_csv, skipinitialspace=True)
    df = df.drop(columns=[c for c in METRIC_COLUMNS if c in df.columns])

    df['ellipsoidVolume'] = (4 / 3 * np.pi * df['EllipsoidLen1']
                             * df['EllipsoidLen2'] * df['EllipsoidLen3'])

    index = build_scat_index(params.ssst_root, params.run_name, params.zones)
    logger.info(f'Indexed {len(index)} scatter clouds across {len(params.zones)} zones')
    print(f'Indexed {len(index)} scatter clouds')

    records, missing = [], []
    for position, (public_id, source) in enumerate(zip(df['publicId'], df['source']), 1):
        scat_path = index.get((source, public_id))
        if scat_path is None:
            missing.append(public_id)
            records.append(_blank_metrics())
        else:
            records.append(event_metrics(scat_path, public_id))
        if position % 5000 == 0:
            print(f'  {position}/{len(df)} events')

    metrics = pd.DataFrame.from_records(records, index=df.index)
    for column in METRIC_COLUMNS:
        if column not in ('ellipsoidVolume', 'dip_reject'):
            df[column] = metrics[column]
    df['n_scat'] = df['n_scat'].astype('Int64')      # nullable: NaN where no .scat

    # Reject only a multimodality that is both statistically real and wider than
    # the depth error the catalog quotes — the dip statistic alone cannot see km.
    # NaN (unevaluated cloud) compares False throughout, so an unreadable event is
    # never rejected here; it is caught by the no-metrics test downstream instead.
    df['dip_reject'] = (
        (df['dip_pval'] < _DIP_ALPHA)
        & (df['dip_sep_km'] > _DIP_SEP_ERZ_FACTOR * df['true_erz'])
    ).fillna(False).astype(bool)

    output_path = params.output or params.result_csv
    parent = os.path.dirname(output_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    tmp_path = f'{output_path}.tmp'
    df.to_csv(tmp_path, index=False)
    os.replace(tmp_path, output_path)

    logger.info(f'Events            : {len(df)}')
    logger.info(f'Without a .scat   : {len(missing)}')
    for public_id in missing[:50]:
        logger.info(f'  MISSING {public_id}')
    logger.info(f'Output            : {output_path!r}')

    print(f'Metrics saved @ {output_path} ({len(df)} events, {len(missing)} without a .scat)')

    # Non-fatal, exactly like this whole step is inside run_SSST.py: the annotated
    # CSV is already on disk, so a plotting failure must not make it look lost.
    figures = []
    if params.figures:
        try:
            figures = generate_figures(df, params.run_name, params.figure_dir)
        except Exception as exc:
            logger.warning(f'figures failed ({exc.__class__.__name__}: {exc})')
            print(f'Figures failed ({exc.__class__.__name__}: {exc})\n'
                  f'  {output_path} is complete; rerun with:\n'
                  f'    python NLL_run/pdf_metrics.py --run-name {params.run_name}')

    return {
        'output':    output_path,
        'log':       log_path,
        'n_events':  len(df),
        'n_missing': len(missing),
        'figures':   figures,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _str2bool(value):
    if value.lower() in ('true', '1', 'yes'):
        return True
    if value.lower() in ('false', '0', 'no'):
        return False
    raise argparse.ArgumentTypeError(f"expected true/false, got '{value}'")


def main():
    parser = argparse.ArgumentParser(
        description='Append location-PDF quality metrics (J, Psi, C68, dip test) '
                    'to a merged NLLoc result CSV, and save their diagnostic figures.')
    parser.add_argument('--run-name', default='ssst_run1',
                        help='SSST campaign name (default: ssst_run1)')
    parser.add_argument('--result-csv',
                        default=os.path.join(_PROJECT_ROOT, 'RESULT', 'SSST_result.csv'),
                        help='Merged catalog to annotate (default: RESULT/SSST_result.csv)')
    parser.add_argument('--ssst-root', default=os.path.join(_PROJECT_ROOT, 'run', 'ssst_loc'),
                        help='Root holding <run-name>/Pyrenees_<zone>_SSST/... (default: run/ssst_loc)')
    parser.add_argument('--zones', nargs='+', default=['1', '2', '3', '4', '5', '6'],
                        help='Zone keys to index (default: 1 2 3 4 5 6)')
    parser.add_argument('--output', default=None,
                        help='Output CSV path (default: rewrite --result-csv in place)')
    parser.add_argument('--log-dir', default=None,
                        help='Log directory (default: NLL_run/console_output/)')
    parser.add_argument('--figures', type=_str2bool, default=True,
                        help='Save the diagnostic figures (default: true)')
    parser.add_argument('--figure-dir', default=None,
                        help='Figure directory (default: complem_figures/pdf_metrics/)')
    args = parser.parse_args()

    compute_metrics(PdfMetricsParams(
        result_csv = args.result_csv,
        ssst_root  = args.ssst_root,
        run_name   = args.run_name,
        zones      = args.zones,
        output     = args.output,
        log_dir    = args.log_dir,
        figures    = args.figures,
        figure_dir = args.figure_dir,
    ))


if __name__ == '__main__':
    main()
