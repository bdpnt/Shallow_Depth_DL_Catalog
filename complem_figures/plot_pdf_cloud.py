"""
plot_pdf_cloud.py
============================
Visualize the NLLoc PDF location scatter-cloud of a single event across all
SSST iterations (loc_ssst_corr0 ... loc_ssst_corrN, the last being the final
NLLoc-only pass), as an interactive 3D scene: one smooth confidence-ellipsoid
surface per iteration, built directly from the mean/covariance NLLoc already
computes and reports in each .hyp file's STATISTICS line (no re-estimation
from the scattered samples), plus the raw sample cloud and the path traced by
the maximum-likelihood point. Click a legend entry to show/hide that iteration.

Each ellipsoid is centred on the PDF expectation, which is also the hypocentre
the catalog publishes (`LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`). The diamond
markers are the *maximum-likelihood* point of the same PDF, read from the .hyp
MAXIMUM_LIKELIHOOD line — so the gap between a diamond and its ellipsoid centre
is the mode-versus-mean distance for that iteration.

Usage
-----
    python complem_figures/plot_pdf_cloud.py \\
        --run-name ssst_run1 \\
        --event-id PYRENEES_049798 --zone 1

    python complem_figures/plot_pdf_cloud.py \\
        --run-name ssst_run1 \\
        --lat 43.115 --lon -1.500 --date 2023-05-19T03:06:57 \\
        --radius-km 10 --window-days 5
"""

import argparse
import glob
import os
import re
import subprocess
import sys
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import plotly.colors as pcolors
import plotly.graph_objects as go
import pyproj
from scipy.stats import chi2

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
from NLL_run.merge_regional_results import _build_covariance, _S3_3DOF  # noqa: E402
# Not obspy's read_nlloc_scatter: it does `np.fromfile(...)[4:]`, dropping four
# 16-byte records where the header is one, so it silently loses the first 3
# samples of every cloud.
from NLL_run.pdf_metrics import read_scat as read_nlloc_scatter  # noqa: E402

_TRANS_RE = re.compile(
    r'TRANSFORM\s+LAMBERT\s+RefEllipsoid\s+(\S+)\s+'
    r'LatOrig\s+([-\d.]+)\s+LongOrig\s+([-\d.]+)\s+'
    r'FirstStdParal\s+([-\d.]+)\s+SecondStdParal\s+([-\d.]+)\s+'
    r'RotCW\s+([-\d.]+)'
)

_STATISTICS_RE = re.compile(
    r'STATISTICS\s+ExpectX\s+([-\d.eE+]+)\s+Y\s+([-\d.eE+]+)\s+Z\s+([-\d.eE+]+)\s+'
    r'CovXX\s+([-\d.eE+]+)\s+XY\s+([-\d.eE+]+)\s+XZ\s+([-\d.eE+]+)\s+'
    r'YY\s+([-\d.eE+]+)\s+YZ\s+([-\d.eE+]+)\s+ZZ\s+([-\d.eE+]+)'
)

# NLLoc writes this line only under `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`, where
# the reported hypocenter is the PDF expectation and the maximum-likelihood point
# would otherwise be lost. Without it the diamond marker below would sit exactly
# on the ellipsoid centre, which says nothing.
_MAXLIKE_RE = re.compile(
    r'MAXIMUM_LIKELIHOOD\s+MaxLikeLat\s+([-\d.eE+]+)\s+Long\s+([-\d.eE+]+)\s+'
    r'Depth\s+([-\d.eE+]+)'
)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PdfCloudParams:
    ssst_root:      str
    run_name:       str
    event_id:       str
    zone:           str
    nll_result_csv: str = None
    confidence:     float = 0.6827   # erf(1/sqrt(2)), one sigma
    output:       str = None


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _find_iteration_dirs(ssst_root, run_name, zone):
    """Return sorted (step, dir) pairs for every loc_ssst_corr<N>/GLOBAL_<zone> folder."""
    pattern = os.path.join(ssst_root, run_name, f'Pyrenees_{zone}_SSST', 'loc_ssst_corr*', f'GLOBAL_{zone}')
    step_dirs = []
    for d in glob.glob(pattern):
        match = re.search(r'loc_ssst_corr(\d+)', d)
        step_dirs.append((int(match.group(1)), d))
    return sorted(step_dirs)


def _find_hyp_path(iter_dir, event_id):
    """Locate the .hyp file whose PUBLIC_ID line matches event_id, via grep (fast across thousands of files).

    Excludes the per-iteration `*.sum.grid0.loc.hyp` file: it concatenates every event's .hyp
    block for that iteration (so it also matches the grep), but has no companion .scat file.
    """
    proc = subprocess.run(
        ['grep', '-rlF', f'PUBLIC_ID {event_id}', '--exclude=*.sum.*', iter_dir],
        capture_output=True, text=True,
    )
    matches = proc.stdout.strip().splitlines()
    return matches[0] if matches else None


def _sibling(hyp_path, ext):
    return hyp_path[:-len('.hyp')] + ext


def _read_lambert_params(hdr_path):
    """Parse the per-event TRANSFORM LAMBERT line from a .grid0.loc.hdr file."""
    with open(hdr_path) as f:
        text = f.read()
    match = _TRANS_RE.search(text)
    if not match:
        raise ValueError(f'No LAMBERT TRANSFORM found in {hdr_path}')
    _, lat0, lon0, p1, p2, rot = match.groups()
    if float(rot) != 0.0:
        warnings.warn(f'{hdr_path}: non-zero RotCW ({rot}) is not supported by this converter — ignoring.')
    return dict(lat0=float(lat0), lon0=float(lon0), p1=float(p1), p2=float(p2))


def _geographic_to_local_transformer(params):
    """Build a (lon,lat)->(x,y) transformer into the event's local Lambert km frame."""
    crs = pyproj.CRS.from_proj4(
        f"+proj=lcc +lat_1={params['p1']} +lat_2={params['p2']} +lat_0={params['lat0']} "
        f"+lon_0={params['lon0']} +x_0=0 +y_0=0 +ellps=WGS84 +units=km"
    )
    return pyproj.Transformer.from_crs('EPSG:4326', crs, always_xy=True)


def _parse_statistics(hyp_path):
    """Parse the STATISTICS line of a .hyp file: expectation point + covariance, both local km."""
    with open(hyp_path) as f:
        text = f.read()
    match = _STATISTICS_RE.search(text)
    if not match:
        raise ValueError(f'No STATISTICS line found in {hyp_path}')
    ex, ey, ez, cxx, cxy, cxz, cyy, cyz, czz = map(float, match.groups())
    center = np.array([ex, ey, ez])
    cov = np.array([
        [cxx, cxy, cxz],
        [cxy, cyy, cyz],
        [cxz, cyz, czz],
    ])
    return center, cov


def _parse_maxlike(hyp_path, to_local):
    """Parse the MAXIMUM_LIKELIHOOD line of a .hyp file into local km, or None if absent."""
    with open(hyp_path) as f:
        match = _MAXLIKE_RE.search(f.read())
    if not match:
        return None
    lat, lon, depth = map(float, match.groups())
    x, y = to_local.transform(lon, lat)
    return x, y, depth


def _load_iteration(iter_dir, step, event_id):
    """Load one iteration's PDF cloud (x/y/z/pdf), covariance statistics, and the
    maximum-likelihood hypocenter for one event.

    x/y/z are kept in NLLoc's native local Lambert-projected km frame (as recorded in the .scat
    file) rather than converted to lon/lat degrees, matching the frame NLLoc's own STATISTICS
    (expectation + covariance) are reported in.
    """
    hyp_path = _find_hyp_path(iter_dir, event_id)
    if hyp_path is None:
        return None

    scat_path = _sibling(hyp_path, '.scat')
    hdr_path = _sibling(hyp_path, '.hdr')
    if not os.path.exists(scat_path):
        warnings.warn(f'{event_id}: missing .scat at step {step} — skipping this iteration.')
        return None

    cloud = read_nlloc_scatter(scat_path)
    center, cov = _parse_statistics(hyp_path)
    to_local = _geographic_to_local_transformer(_read_lambert_params(hdr_path))

    # The diamond marker is the maximum-likelihood point, read from the .hyp: the
    # CSV's latitude/longitude/depth are the expectation, i.e. `center` above.
    maxlike = _parse_maxlike(hyp_path, to_local)
    hyp_x, hyp_y, hyp_z = maxlike if maxlike else (None, None, None)

    csv_matches = glob.glob(os.path.join(iter_dir, 'Pyrenees_*.sum.grid0.loc.csv'))
    date_str = pdf_volume = ellipsoid_volume = None
    if csv_matches:
        df = pd.read_csv(csv_matches[0], skipinitialspace=True)
        row = df.loc[df['publicId'] == event_id]
        if not row.empty:
            row = row.iloc[0]
            date_str = row['date-time']
            pdf_volume = row['pdfVolume']
            ellipsoid_volume = 4 / 3 * np.pi * row['EllipsoidLen1'] * row['EllipsoidLen2'] * row['EllipsoidLen3']

    return {
        'step': step,
        'x': cloud['x'], 'y': cloud['y'], 'z': cloud['z'], 'pdf': cloud['pdf'],
        'center': center, 'cov': cov, 'to_local': to_local,
        'hyp_x': hyp_x, 'hyp_y': hyp_y, 'hyp_z': hyp_z, 'date': date_str,
        'pdf_volume': pdf_volume, 'ellipsoid_volume': ellipsoid_volume,
    }


def _load_nll_reference(nll_result_csv, event_id, to_local):
    """Load the pre-SSST NLL-only relocation ellipsoid + hypocenter for one event from
    RESULT/NLL_result.csv (the merged, zone-deduplicated catalog).

    The ellipsoid axes there (EllipsoidAz1/Dip1/Len1, Az2/Dip2/Len2, Len3) already carry
    NLLoc's own 3-DOF, 68% chi-square scaling, so the reconstructed covariance is divided by
    that same factor (`_S3_3DOF`, reused from merge_regional_results.py) to get back a raw,
    unscaled covariance — matching what `_ellipsoid_surface` expects (it applies its own
    confidence-level scaling on top).

    `to_local` (an SSST iteration's Lambert-projection transformer) is used to project the
    row's expectation point / hypocenter from lon/lat into the same local km frame the SSST
    iterations are plotted in.

    latitude/longitude/depth in the merged CSV are the PDF expectation — the point the
    ellipsoid is a second moment about, hence the ellipsoid centre. The maximum-likelihood
    point is the separate maxlike_* triplet, drawn as the marker.

    Returns
    -------
    dict with keys: zone, center, cov, hyp_x, hyp_y, hyp_z, pdf_volume, ellipsoid_volume — or None
    if the event isn't in the CSV.
    """
    df = pd.read_csv(nll_result_csv, skipinitialspace=True)
    row = df.loc[df['publicId'] == event_id]
    if row.empty:
        return None
    row = row.iloc[0]

    ell_args = (row['EllipsoidAz1'], row['EllipsoidDip1'], row['EllipsoidLen1'],
                row['EllipsoidAz2'], row['EllipsoidDip2'], row['EllipsoidLen2'],
                row['EllipsoidLen3'])
    if 'maxlike_latitude' not in row:
        raise KeyError(
            f'{nll_result_csv} has no maxlike_* columns — it predates '
            f'`LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`. Rerun the NLL stage, or point '
            f'--nll-result-csv at a merged CSV from a run that carries them.')

    cov = _build_covariance(*ell_args) / _S3_3DOF
    ellipsoid_volume = 4 / 3 * np.pi * row['EllipsoidLen1'] * row['EllipsoidLen2'] * row['EllipsoidLen3']

    center_x, center_y = to_local.transform(row['longitude'], row['latitude'])
    hyp_x, hyp_y = to_local.transform(row['maxlike_longitude'], row['maxlike_latitude'])

    return {
        'zone': row['source'].replace('GLOBAL_', ''),
        'center': np.array([center_x, center_y, row['depth']]),
        'cov': cov,
        'hyp_x': hyp_x, 'hyp_y': hyp_y, 'hyp_z': row['maxlike_depth'],
        'pdf_volume': row['pdfVolume'], 'ellipsoid_volume': ellipsoid_volume,
    }


def _ellipsoid_surface(center, cov, confidence, n=24):
    """
    Parametric confidence-ellipsoid surface from a 3x3 covariance matrix.

    The region {(x-center)^T cov^-1 (x-center) <= chi2.ppf(confidence, df=3)} contains
    `confidence` of the probability mass of a 3D Gaussian with this mean/covariance — the
    standard construction for a location-uncertainty confidence ellipsoid. Always smooth by
    construction (unlike a convex hull of scattered points, which is a faceted polyhedron).

    Returns
    -------
    X, Y, Z : (n, n) arrays suitable for go.Surface
    """
    eigvals, eigvecs = np.linalg.eigh(cov)
    eigvals = np.clip(eigvals, 0, None)
    scale = np.sqrt(chi2.ppf(confidence, df=3) * eigvals)

    u = np.linspace(0, 2 * np.pi, n)
    v = np.linspace(0, np.pi, n)
    unit_sphere = np.stack([
        np.outer(np.cos(u), np.sin(v)),
        np.outer(np.sin(u), np.sin(v)),
        np.outer(np.ones_like(u), np.cos(v)),
    ], axis=-1)

    pts = (unit_sphere * scale) @ eigvecs.T + center
    return pts[..., 0], pts[..., 1], pts[..., 2]


def _build_figure(iterations, confidence, event_id, zone, run_name, nll_ref=None):
    """Build the interactive 3D figure: one confidence ellipsoid + raw cloud + hypocenter per
    iteration, plus the pre-SSST NLL reference ellipsoid (`nll_ref`) if available."""
    steps = [it['step'] for it in iterations]
    n_steps = len(steps)
    colors = pcolors.sample_colorscale('Plasma', [i / max(n_steps - 1, 1) for i in range(n_steps)])

    traces = []
    hyp_xs, hyp_ys, hyp_zs = [], [], []
    vol_rows = []
    rng = np.random.default_rng(0)

    if nll_ref is not None:
        label = 'NLL (pre-SSST)' if nll_ref['zone'] == zone else f"NLL (pre-SSST, Zone {nll_ref['zone']})"
        vol_rows.append((label, nll_ref['pdf_volume'], nll_ref['ellipsoid_volume']))
        nx, ny, nz = _ellipsoid_surface(nll_ref['center'], nll_ref['cov'], confidence)
        traces.append(go.Surface(
            x=nx, y=ny, z=nz,
            colorscale=[[0, 'gray'], [1, 'gray']], showscale=False, opacity=0.35,
            name=label, legendgroup='nll_ref', showlegend=True,
        ))
        traces.append(go.Scatter3d(
            x=[nll_ref['hyp_x']], y=[nll_ref['hyp_y']], z=[nll_ref['hyp_z']],
            mode='markers',
            marker=dict(size=6, color='gray', symbol='diamond', line=dict(color='black', width=1)),
            name=label, legendgroup='nll_ref', showlegend=False,
        ))
        hyp_xs.append(nll_ref['hyp_x'])
        hyp_ys.append(nll_ref['hyp_y'])
        hyp_zs.append(nll_ref['hyp_z'])

    for it, color in zip(iterations, colors):
        label = 'Final' if it['step'] == steps[-1] else f'Iter {it["step"]}'
        group = f'iter{it["step"]}'
        x, y, z = it['x'], it['y'], it['z']
        vol_rows.append((label, it['pdf_volume'], it['ellipsoid_volume']))

        ex, ey, ez = _ellipsoid_surface(it['center'], it['cov'], confidence)
        traces.append(go.Surface(
            x=ex, y=ey, z=ez,
            colorscale=[[0, color], [1, color]], showscale=False, opacity=0.35,
            name=label, legendgroup=group, showlegend=True,
        ))

        n_show = min(len(x), 3000)
        idx = rng.choice(len(x), size=n_show, replace=False) if len(x) > n_show else np.arange(len(x))
        traces.append(go.Scatter3d(
            x=x[idx], y=y[idx], z=z[idx],
            mode='markers', marker=dict(size=1.5, color=color, opacity=0.3),
            name=label, legendgroup=group, showlegend=False,
        ))

        if it['hyp_x'] is not None:
            traces.append(go.Scatter3d(
                x=[it['hyp_x']], y=[it['hyp_y']], z=[it['hyp_z']],
                mode='markers',
                marker=dict(size=6, color=color, symbol='diamond', line=dict(color='black', width=1)),
                name=label, legendgroup=group, showlegend=False,
            ))
            hyp_xs.append(it['hyp_x'])
            hyp_ys.append(it['hyp_y'])
            hyp_zs.append(it['hyp_z'])

    if len(hyp_xs) > 1:
        traces.append(go.Scatter3d(
            x=hyp_xs, y=hyp_ys, z=hyp_zs,
            mode='lines', line=dict(color='black', width=3),
            name='Max-likelihood path', showlegend=True,
        ))

    date_str = next((it['date'] for it in iterations if it['date']), 'unknown date')
    subtitle = f'{confidence:.0%} confidence ellipsoid per iteration'
    if nll_ref is not None:
        subtitle += ' + pre-SSST NLL reference (gray)'
    subtitle += ' — click a legend entry to toggle it'
    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text=(
            f'{event_id} — Zone {zone} ({run_name}) — {date_str}<br>'
            f'<sub>{subtitle}</sub>'
        )),
        scene=dict(
            xaxis_title='Easting (km, local)', yaxis_title='Northing (km, local)', zaxis_title='Depth (km)',
            zaxis=dict(autorange='reversed'),
            aspectmode='data',
        ),
        legend=dict(itemsizing='constant'),
    )

    def _fmt(v):
        return f'{v:.3g}' if v is not None else 'n/a'

    table_lines = [f"{'Step':<24}{'pdfVol (km³)':>14}{'ellVol (km³)':>14}"]
    table_lines += [f'{label:<24}{_fmt(pdf_vol):>14}{_fmt(ell_vol):>14}' for label, pdf_vol, ell_vol in vol_rows]
    fig.add_annotation(
        text='<br>'.join(table_lines), xref='paper', yref='paper', x=0.01, y=0.99,
        xanchor='left', yanchor='top', align='left', showarrow=False,
        font=dict(family='Courier New, monospace', size=11),
        bgcolor='rgba(255,255,255,0.75)', bordercolor='black', borderwidth=1,
    )
    return fig


def _resolve_event_by_location(result_csv, lat, lon, date, radius_km, window_days):
    """Search RESULT/SSST_result.csv for the event nearest (lat, lon, date) within tolerance."""
    df = pd.read_csv(result_csv, skipinitialspace=True)
    df['date-time'] = pd.to_datetime(df['date-time'])
    target_date = pd.to_datetime(date)

    r_earth = 6371.0
    lat1, lon1 = np.radians(df['latitude']), np.radians(df['longitude'])
    lat2, lon2 = np.radians(lat), np.radians(lon)
    dlat, dlon = lat2 - lat1, lon2 - lon1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    dist_km = 2 * r_earth * np.arcsin(np.sqrt(a))
    dt_days = (df['date-time'] - target_date).abs().dt.total_seconds() / 86400.0

    candidates = df[(dist_km <= radius_km) & (dt_days <= window_days)]
    if candidates.empty:
        raise ValueError(f'No event found within {radius_km} km / {window_days} days of ({lat}, {lon}, {date})')
    if len(candidates) > 1:
        print('Multiple candidates found:')
        print(candidates[['publicId', 'source', 'latitude', 'longitude', 'date-time']].to_string(index=False))
        raise ValueError('Ambiguous search — narrow --radius-km / --window-days or use --event-id/--zone directly.')

    row = candidates.iloc[0]
    return row['publicId'], row['source'].replace('GLOBAL_', '')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(params):
    """
    Generate and save the PDF scatter-cloud iteration-overlay figure for one event.

    Parameters
    ----------
    params : PdfCloudParams

    Returns
    -------
    dict with keys: output_path, iterations_found
    """
    iter_dirs = _find_iteration_dirs(params.ssst_root, params.run_name, params.zone)
    if not iter_dirs:
        raise FileNotFoundError(
            f'No SSST iteration folders found under {params.ssst_root}/{params.run_name}/Pyrenees_{params.zone}_SSST'
        )

    iterations = []
    for step, iter_dir in iter_dirs:
        data = _load_iteration(iter_dir, step, params.event_id)
        if data is not None:
            iterations.append(data)

    if not iterations:
        raise ValueError(f'Event {params.event_id} not found in any SSST iteration for zone {params.zone}')

    nll_result_csv = params.nll_result_csv or os.path.join(_PROJECT_ROOT, 'RESULT', 'NLL_result.csv')
    nll_ref = _load_nll_reference(nll_result_csv, params.event_id, iterations[0]['to_local'])
    if nll_ref is None:
        warnings.warn(f'{params.event_id}: not found in {nll_result_csv} — omitting NLL reference ellipsoid.')

    output_path = params.output or os.path.join(_MODULE_DIR, 'pdf_cloud', f'{params.event_id}.html')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig = _build_figure(iterations, params.confidence, params.event_id, params.zone, params.run_name, nll_ref=nll_ref)
    fig.write_html(output_path, include_plotlyjs=True)
    print(f'Figure saved @ {output_path} ({len(iterations)} iteration(s) found)')

    return {'output_path': output_path, 'iterations_found': [it['step'] for it in iterations], 'nll_reference_found': nll_ref is not None}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Visualize the NLLoc PDF scatter-cloud evolution of one event across SSST iterations.'
    )
    parser.add_argument('--ssst-root', default=os.path.join(_PROJECT_ROOT, 'run', 'ssst_loc'),
                        help='Root folder containing <run-name>/Pyrenees_<zone>_SSST/...')
    parser.add_argument('--run-name', default='ssst_run1',
                        help='SSST campaign name (default: ssst_run1)')
    parser.add_argument('--result-csv', default=os.path.join(_PROJECT_ROOT, 'RESULT', 'SSST_result.csv'),
                        help='Merged catalog used for --lat/--lon/--date search (default: RESULT/SSST_result.csv)')
    parser.add_argument('--nll-result-csv', default=os.path.join(_PROJECT_ROOT, 'RESULT', 'NLL_result.csv'),
                        help='Pre-SSST merged catalog, used for the NLL reference ellipsoid (default: RESULT/NLL_result.csv)')
    parser.add_argument('--event-id', help='Event publicId (requires --zone)')
    parser.add_argument('--zone', help='Zone key, e.g. 1..6 (requires --event-id)')
    parser.add_argument('--lat', type=float, help='Approximate latitude for location search')
    parser.add_argument('--lon', type=float, help='Approximate longitude for location search')
    parser.add_argument('--date', help='Approximate ISO date/time for location search')
    parser.add_argument('--radius-km', type=float, default=10.0, help='Search radius for --lat/--lon (default: 10)')
    parser.add_argument('--window-days', type=float, default=5.0, help='Search time window for --date (default: 5)')
    parser.add_argument('--confidence', type=float, default=0.6827,
                        help='Confidence-ellipsoid level per iteration surface (default: 0.6827, one sigma)')
    parser.add_argument('--output', default=None,
                        help='Output HTML path (default: complem_figures/pdf_cloud/<event_id>.html)')
    args = parser.parse_args()

    if args.event_id and args.zone:
        event_id, zone = args.event_id, args.zone
    elif args.lat is not None and args.lon is not None and args.date:
        event_id, zone = _resolve_event_by_location(
            args.result_csv, args.lat, args.lon, args.date, args.radius_km, args.window_days
        )
        print(f'Resolved to event {event_id}, zone {zone}')
    else:
        parser.error('Provide either --event-id and --zone, or --lat, --lon and --date.')

    generate_figure(PdfCloudParams(
        ssst_root      = args.ssst_root,
        run_name       = args.run_name,
        event_id       = event_id,
        zone           = zone,
        nll_result_csv = args.nll_result_csv,
        confidence     = args.confidence,
        output         = args.output,
    ))


if __name__ == '__main__':
    main()
