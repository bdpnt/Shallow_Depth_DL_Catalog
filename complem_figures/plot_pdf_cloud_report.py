"""
plot_pdf_cloud_report.py
============================
Report-ready 2-D counterpart of `plot_pdf_cloud.py`: the NLLoc location PDF of a
single event, one iteration, drawn as two static panels instead of an
interactive 3D scene.

    left   map view          longitude / latitude, always 1:1
    right  vertical section  depth against the horizontal direction in which the
                             confidence ellipsoid is widest (NLLoc's
                             `azMaxHorUnc` direction)

The section fills its panel when the cloud is wider than it is deep, which means
vertical exaggeration; the factor is measured off the drawn axes and written on
the depth axis. Depth is never *compressed* to fill the panel — flattening the
coordinate the figure exists to show would be the wrong trade — so a cloud
deeper than it is wide is drawn at 1:1 and leaves its unused space horizontally.
`--true-scale` forces 1:1 in both cases, at the cost of a mostly empty panel for
badly constrained events.

Each panel shows the full `.scat` sample cloud coloured by its PDF value, the
outline of the confidence ellipsoid *projected* onto that plane, the PDF
expectation (the ellipsoid centre) and the maximum-likelihood point.

Two points of interpretation
----------------------------
1. *Projected*, not sliced. The outline drawn is the shadow of the 3-D
   confidence ellipsoid on the plane of the panel, i.e. the 2x2 sub-block of the
   covariance scaled by chi2.ppf(confidence, df=3). That is the curve that
   contains the same fraction of the plotted (projected) samples as the
   ellipsoid contains of the cloud. It is NOT the 2-D marginal confidence
   ellipse, which is the same shape scaled by chi2.ppf(confidence, df=2) and is
   therefore smaller; `--marginal` adds it as a dashed curve for comparison.

2. Colour is redundant with density, not a weight. NLLoc draws oct-tree scatter
   samples with density proportional to the posterior, so where the points are
   dense the probability is already high; the fourth `.scat` field is carried
   here only to make the shape of the mode easier to read. `--color none`
   disables it.

   That fourth field is the *natural logarithm* of the unnormalised posterior,
   not the posterior itself. This is an inference from the files, not something
   the repository documents: the values of a typical cloud span ~110-119, they
   peak at the maximum-likelihood point, and read as a density they would imply
   a 7 % variation across a 9 km-deep, visibly non-Gaussian cloud, whereas read
   as a log they imply a factor ~7000 — which is what an oct-tree sample set of
   a one-sigma region should show. Confirmed independently: NLLoc draws its
   samples with density proportional to the posterior, and the field correlates
   +0.93 (Spearman +0.95) with the log of a 12-nearest-neighbour density
   estimate over the cloud. Samples are therefore coloured by
   ln(PDF / PDF_max), so 0 is the mode of the sampled posterior and -x is a
   density e^-x times it. Pass `--color linear` to treat the field as a density
   instead, if that inference turns out to be wrong.

Which point is the mode
-----------------------
Under `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION` the GEOGRAPHIC line of the .hyp
holds the expectation and a separate MAXIMUM_LIKELIHOOD line holds the mode.
Without that option — which is the case for the `ssst_run1` files on disk —
NLLoc writes no MAXIMUM_LIKELIHOOD line and the GEOGRAPHIC line *is* the
maximum-likelihood hypocentre, the expectation being reported on the
STATISTICS / STAT_GEOG lines. Both layouts are handled, so the mode marker is
never silently dropped.

Usage
-----
    python complem_figures/plot_pdf_cloud_report.py \\
        --run-name ssst_run1 --event-id PYRENEES_034967 --zone 2

    python complem_figures/plot_pdf_cloud_report.py \\
        --run-name ssst_run1 --event-id PYRENEES_034967 --zone 2 \\
        --iteration 0 --output complem_figures/pdf_cloud/PYRENEES_034967_iter0.png

    python complem_figures/plot_pdf_cloud_report.py \\
        --hyp run/ssst_loc/ssst_run1/Pyrenees_2_SSST/loc_ssst_corr5/GLOBAL_2/Pyrenees_2.20180504.051742.grid0.loc.hyp \\
        --output /tmp/event.png

Writes a 300 dpi PNG and a vector PDF of the same name.
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
import pyproj
from scipy.stats import chi2

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator

_MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
# Not obspy's read_nlloc_scatter: it does `np.fromfile(...)[4:]`, dropping four
# 16-byte records where the header is one, so it silently loses the first 3
# samples of every cloud.
from NLL_run.pdf_metrics import read_scat as read_nlloc_scatter  # noqa: E402


# ---------------------------------------------------------------------------
# .hyp / .hdr parsing (copied from plot_pdf_cloud.py, which stays untouched)
# ---------------------------------------------------------------------------

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

_MAXLIKE_RE = re.compile(
    r'MAXIMUM_LIKELIHOOD\s+MaxLikeLat\s+([-\d.eE+]+)\s+Long\s+([-\d.eE+]+)\s+'
    r'Depth\s+([-\d.eE+]+)'
)

_GEOGRAPHIC_RE = re.compile(
    r'GEOGRAPHIC\s+OT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+'
    r'Lat\s+([-\d.eE+]+)\s+Long\s+([-\d.eE+]+)\s+Depth\s+([-\d.eE+]+)'
)

_HORUNC_RE = re.compile(
    r'QML_OriginUncertainty\s+horUnc\s+([-\d.eE+]+)\s+minHorUnc\s+([-\d.eE+]+)\s+'
    r'maxHorUnc\s+([-\d.eE+]+)\s+azMaxHorUnc\s+([-\d.eE+]+)'
)

_PUBLIC_ID_RE = re.compile(r'PUBLIC_ID\s+(\S+)')


def _find_iteration_dirs(ssst_root, run_name, zone):
    """Return sorted (step, dir) pairs for every loc_ssst_corr<N>/GLOBAL_<zone> folder."""
    pattern = os.path.join(ssst_root, run_name, f'Pyrenees_{zone}_SSST',
                           'loc_ssst_corr*', f'GLOBAL_{zone}')
    step_dirs = []
    for d in glob.glob(pattern):
        match = re.search(r'loc_ssst_corr(\d+)', d)
        step_dirs.append((int(match.group(1)), d))
    return sorted(step_dirs)


def _find_hyp_path(iter_dir, event_id):
    """Locate the .hyp file whose PUBLIC_ID line matches event_id, via grep.

    Excludes the per-iteration `*.sum.grid0.loc.hyp` file: it concatenates every
    event's .hyp block for that iteration, but has no companion .scat file.
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


def _lambert_transformers(params):
    """Build (lon,lat)->(x,y) and (x,y)->(lon,lat) transformers for the local km frame."""
    crs = pyproj.CRS.from_proj4(
        f"+proj=lcc +lat_1={params['p1']} +lat_2={params['p2']} +lat_0={params['lat0']} "
        f"+lon_0={params['lon0']} +x_0=0 +y_0=0 +ellps=WGS84 +units=km"
    )
    to_local = pyproj.Transformer.from_crs('EPSG:4326', crs, always_xy=True)
    to_geo = pyproj.Transformer.from_crs(crs, 'EPSG:4326', always_xy=True)
    return to_local, to_geo


def _parse_statistics(hyp_path):
    """Parse the STATISTICS line: expectation point + raw covariance, both local km."""
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


def _parse_hyp_header(hyp_path, to_local):
    """Origin time, maximum-likelihood point in local km, reported horizontal
    uncertainty, and public id.

    The mode is the MAXIMUM_LIKELIHOOD line when NLLoc was run with
    SAVE_NLLOC_EXPECTATION, and the GEOGRAPHIC line otherwise (see module
    docstring).
    """
    with open(hyp_path) as f:
        text = f.read()

    geo = _GEOGRAPHIC_RE.search(text)
    if not geo:
        raise ValueError(f'No GEOGRAPHIC line found in {hyp_path}')
    y, mo, d, h, mi, s, g_lat, g_lon, g_depth = geo.groups()
    date_str = f'{int(y):04d}-{int(mo):02d}-{int(d):02d}T{int(h):02d}:{int(mi):02d}:{float(s):06.3f}'

    maxlike = _MAXLIKE_RE.search(text)
    if maxlike:
        m_lat, m_lon, m_depth = map(float, maxlike.groups())
        maxlike_source = 'MAXIMUM_LIKELIHOOD line'
    else:
        m_lat, m_lon, m_depth = float(g_lat), float(g_lon), float(g_depth)
        maxlike_source = 'GEOGRAPHIC line'

    mx, my = to_local.transform(m_lon, m_lat)

    hor = _HORUNC_RE.search(text)
    min_hor, max_hor, az_max_hor = (
        (float(hor.group(2)), float(hor.group(3)), float(hor.group(4))) if hor else (None, None, None)
    )

    pid = _PUBLIC_ID_RE.search(text)

    return dict(
        date=date_str,
        maxlike=np.array([mx, my, m_depth]),
        maxlike_source=maxlike_source,
        min_hor_unc=min_hor, max_hor_unc=max_hor, az_max_hor_unc=az_max_hor,
        public_id=pid.group(1) if pid else None,
    )


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def _horizontal_major_axis(cov):
    """Unit vector (east, north) of the widest horizontal direction of the ellipsoid,
    and its azimuth in degrees clockwise from north (0-180).

    This is the major axis of the horizontal 2x2 sub-block of the covariance, i.e.
    of the ellipsoid's horizontal shadow — the same quantity NLLoc reports as
    `azMaxHorUnc`, not the horizontal trace of the 3-D major axis.
    """
    vals, vecs = np.linalg.eigh(cov[:2, :2])
    u = vecs[:, int(np.argmax(vals))]
    if u[1] < 0:                       # fix the sign so the azimuth is in [0, 180)
        u = -u
    az = np.degrees(np.arctan2(u[0], u[1])) % 180.0
    return u, az


def _ellipse_offsets(cov2, k, n=241):
    """(n, 2) offsets tracing {d : d^T cov2^-1 d = k} about the origin."""
    vals, vecs = np.linalg.eigh(cov2)
    vals = np.clip(vals, 0.0, None)
    radii = np.sqrt(k * vals)
    t = np.linspace(0.0, 2.0 * np.pi, n)
    return (np.column_stack([np.cos(t), np.sin(t)]) * radii) @ vecs.T


# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------

_STYLE = {
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.size': 9,
    'axes.titlesize': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'axes.edgecolor': '#444444',
    'axes.linewidth': 0.8,
    'grid.color': '#DDDDDD',
    'grid.linewidth': 0.6,
}

_ELLIPSE_COLOR = '#222222'
_EXPECT_COLOR = '#C44E52'
_SECTION_LINE_COLOR = '#777777'


def _tidy(ax):
    ax.set_axisbelow(True)
    ax.grid(True, linestyle=':', linewidth=0.6, color='#CCCCCC')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)


def _draw_outline(ax, x, y, dashed=False, label=None):
    """Ellipse outline with a white halo, so it stays readable over dark scatter."""
    ax.plot(x, y, color='white', lw=2.6, solid_capstyle='round', zorder=4)
    ax.plot(x, y, color=_ELLIPSE_COLOR, lw=1.2,
            linestyle='--' if dashed else '-', zorder=5, label=label)


def _set_limits(ax, x_groups, y_groups, margin=0.06):
    """Axis limits covering every sample and the whole ellipse outline, plus a margin."""
    xs = np.concatenate([np.ravel(np.asarray(g, dtype=float)) for g in x_groups])
    ys = np.concatenate([np.ravel(np.asarray(g, dtype=float)) for g in y_groups])
    for setter, v in ((ax.set_xlim, xs), (ax.set_ylim, ys)):
        lo, hi = float(v.min()), float(v.max())
        pad = margin * max(hi - lo, 1e-9)
        setter(lo - pad, hi + pad)


def _clip_line_to_axes(ax, point, direction):
    """The two points where the line through `point` along `direction` meets the axes box,
    ordered along `direction`."""
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    ts = []
    for lo, hi, c, d in ((x0, x1, point[0], direction[0]), (y0, y1, point[1], direction[1])):
        if abs(d) > 1e-15:
            ts += [(lo - c) / d, (hi - c) / d]
    eps_x, eps_y = 1e-9 * (x1 - x0), 1e-9 * (y1 - y0)
    hits = []
    for t in sorted(ts):
        px, py = point[0] + t * direction[0], point[1] + t * direction[1]
        if x0 - eps_x <= px <= x1 + eps_x and y0 - eps_y <= py <= y1 + eps_y:
            hits.append((px, py))
    if len(hits) < 2:                       # degenerate: fall back to the axes diagonal
        return (x0, y0), (x1, y1)
    return hits[0], hits[-1]


def _draw_markers(ax, expect_xy, maxlike_xy):
    ax.plot(*expect_xy, marker='o', ms=7, mfc=_EXPECT_COLOR, mec='white', mew=1.2,
            ls='none', zorder=7, label='PDF expectation (ellipsoid centre)')
    ax.plot(*maxlike_xy, marker='D', ms=6.5, mfc='white', mec='black', mew=1.2,
            ls='none', zorder=7, label='Maximum-likelihood hypocentre')


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

@dataclass
class PdfCloudReportParams:
    ssst_root: str = None
    run_name: str = None
    event_id: str = None
    zone: str = None
    hyp: str = None
    iteration: int = None          # None -> last available
    confidence: float = 0.6827     # erf(1/sqrt(2)), one sigma
    color_by: str = 'log'          # 'log' | 'linear' | 'none'
    cmap: str = 'viridis'
    marginal: bool = False
    section_line: bool = True
    true_scale: bool = False       # True -> section at 1:1, at the cost of empty panel
    output: str = None


def _build_figure(cloud, center, cov, info, params, iter_label):
    conf = params.confidence
    k3 = chi2.ppf(conf, df=3)
    k2 = chi2.ppf(conf, df=2)

    x, y, z = cloud['x'].astype(float), cloud['y'].astype(float), cloud['z'].astype(float)
    pdf = cloud['pdf'].astype(float)
    n = len(x)

    to_local, to_geo = info['to_local'], info['to_geo']
    u, az = _horizontal_major_axis(cov)

    # --- panel A: map view, in degrees -------------------------------------
    lon, lat = to_geo.transform(x, y)
    c_lon, c_lat = to_geo.transform(center[0], center[1])
    m_lon, m_lat = to_geo.transform(info['maxlike'][0], info['maxlike'][1])

    ell_h = _ellipse_offsets(cov[:2, :2], k3)
    ell_h_lon, ell_h_lat = to_geo.transform(center[0] + ell_h[:, 0], center[1] + ell_h[:, 1])
    ellm_h = _ellipse_offsets(cov[:2, :2], k2)
    ellm_h_lon, ellm_h_lat = to_geo.transform(center[0] + ellm_h[:, 0], center[1] + ellm_h[:, 1])

    # --- panel B: section along the widest horizontal direction -------------
    s = (x - center[0]) * u[0] + (y - center[1]) * u[1]
    s_max = (info['maxlike'][0] - center[0]) * u[0] + (info['maxlike'][1] - center[1]) * u[1]

    proj = np.array([[u[0], u[1], 0.0], [0.0, 0.0, 1.0]])
    cov_sz = proj @ cov @ proj.T
    ell_v = _ellipse_offsets(cov_sz, k3)
    ellm_v = _ellipse_offsets(cov_sz, k2)

    # --- draw ---------------------------------------------------------------
    with plt.rc_context(_STYLE):
        fig = plt.figure(figsize=(9.2, 5.6))
        has_cbar = params.color_by != 'none'
        gs = fig.add_gridspec(1, 2, left=0.085, right=0.875 if has_cbar else 0.975,
                              bottom=0.155, top=0.845, wspace=0.30)
        ax_map = fig.add_subplot(gs[0])
        ax_sec = fig.add_subplot(gs[1])

        if params.color_by == 'log':
            cvals = pdf - pdf.max()          # ln PDF relative to the mode, in nats
            cbar_label = r'$\ln\,(\mathrm{PDF}\,/\,\mathrm{PDF}_{\max})$'
        elif params.color_by == 'linear':
            cvals = pdf / pdf.max()
            cbar_label = 'PDF value (normalised)'
        else:
            cvals = None

        if cvals is not None:
            order = np.argsort(cvals)        # most likely samples drawn last, on top
            sc_kw = dict(c=cvals[order], cmap=params.cmap)
        else:
            order = np.arange(n)
            sc_kw = dict(color='#4C72B0')

        sc = ax_map.scatter(np.asarray(lon)[order], np.asarray(lat)[order],
                            s=14, alpha=0.85, linewidths=0, zorder=2, **sc_kw)
        ax_sec.scatter(s[order], z[order], s=14, alpha=0.85, linewidths=0, zorder=2, **sc_kw)

        # --- map panel ---
        if params.marginal:
            _draw_outline(ax_map, ellm_h_lon, ellm_h_lat, dashed=True)
        _draw_outline(ax_map, ell_h_lon, ell_h_lat)
        _draw_markers(ax_map, (c_lon, c_lat), (m_lon, m_lat))

        _set_limits(ax_map, [lon, ell_h_lon], [lat, ell_h_lat])
        if params.section_line:
            # direction of the section in lon/lat, taken from the projection itself so
            # that it carries the Lambert meridian convergence exactly
            d_lon, d_lat = to_geo.transform(center[0] + u[0], center[1] + u[1])
            d = (float(d_lon) - float(c_lon), float(d_lat) - float(c_lat))
            (ax_lon, ay_lat), (bx_lon, by_lat) = _clip_line_to_axes(ax_map, (c_lon, c_lat), d)
            ax_map.plot([ax_lon, bx_lon], [ay_lat, by_lat],
                        color=_SECTION_LINE_COLOR, lw=0.9, ls=(0, (5, 4)), zorder=3)
            for lbl, xx, yy, off in (("A", ax_lon, ay_lat, (-12, 3)),
                                     ("A′", bx_lon, by_lat, (5, -10))):
                ax_map.annotate(lbl, (xx, yy), textcoords='offset points', xytext=off,
                                color=_SECTION_LINE_COLOR, fontsize=8, fontweight='bold',
                                zorder=6, annotation_clip=False)

        ax_map.set_xlabel('Longitude (°E)')
        ax_map.set_ylabel('Latitude (°N)')
        ax_map.set_title('Map view', pad=6)
        ax_map.set_aspect(1.0 / np.cos(np.radians(c_lat)), adjustable='box')
        ax_map.set_anchor('N')
        span = max(np.ptp(ax_map.get_xlim()), np.ptp(ax_map.get_ylim()))
        dec = 3 if span < 0.05 else (2 if span < 0.5 else 1)
        ax_map.xaxis.set_major_formatter(FormatStrFormatter(f'%.{dec}f'))
        ax_map.yaxis.set_major_formatter(FormatStrFormatter(f'%.{dec}f'))
        ax_map.xaxis.set_major_locator(MaxNLocator(4))
        ax_map.yaxis.set_major_locator(MaxNLocator(5))
        _tidy(ax_map)

        # --- section panel ---
        ell_v_z = ell_v[:, 1] + center[2]
        if params.marginal:
            _draw_outline(ax_sec, ellm_v[:, 0], ellm_v[:, 1] + center[2], dashed=True,
                          label=f'{conf:.0%} marginal ellipse (2 d.o.f.)')
        _draw_outline(ax_sec, ell_v[:, 0], ell_v_z,
                      label=f'{conf:.0%} confidence ellipsoid, projected')
        _draw_markers(ax_sec, (0.0, center[2]), (s_max, info['maxlike'][2]))

        _set_limits(ax_sec, [s, ell_v[:, 0]], [z, ell_v_z])
        ax_sec.set_xlabel(f"Distance along A–A′, N{az:.0f}°E (km)")
        ax_sec.set_ylabel('Depth (km)')
        ax_sec.set_title(f'Vertical section along A–A′ (N{az:.0f}°E)', pad=6)
        ax_sec.invert_yaxis()
        if params.true_scale:
            ax_sec.set_aspect('equal', adjustable='box')
            ax_sec.set_anchor('N')
        _tidy(ax_sec)

        # --- colour bar, in its own axes so the panels keep their aspect ---
        if has_cbar:
            cax = fig.add_axes([0.905, 0.26, 0.016, 0.44])
            cbar = fig.colorbar(sc, cax=cax)
            cbar.set_label(cbar_label)
            cbar.outline.set_linewidth(0.6)

        # --- titles carry the statistics, so nothing is written over the data ---
        dh = np.hypot(info['maxlike'][0] - center[0], info['maxlike'][1] - center[1])
        dz = info['maxlike'][2] - center[2]
        line2 = f'{n} PDF samples'
        if info['max_hor_unc'] is not None:
            line2 += (f" · horizontal semi-axes {info['min_hor_unc']:.2f}/{info['max_hor_unc']:.2f} km"
                      f" · vertical semi-axis {np.sqrt(k3 * cov[2, 2]):.2f} km")
        line2 += f' · mode–mean offset {dh:.2f} km horizontal, {dz:+.2f} km in depth'
        fig.suptitle(
            f"{info['public_id']} — {info['date']} — zone {info['zone']}"
            f" ({info['run_name']}, {iter_label})\n"
            f"{line2}",
            fontsize=9.5, y=0.975, linespacing=1.5,
        )

        handles, labels = ax_sec.get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=len(labels),
                   frameon=False, bbox_to_anchor=(0.5, 0.005), handletextpad=0.5,
                   columnspacing=1.6)

        # The section fills its panel, so its vertical and horizontal scales differ
        # whenever the cloud is not as deep as it is wide. Measure the ratio off the
        # drawn axes and state it, rather than leaving the reader to assume 1:1.
        # Depth is never compressed: a ratio below 1 would flatten the coordinate the
        # figure exists to show, and the empty space it avoids is horizontal anyway.
        if not params.true_scale:
            fig.canvas.draw()
            box = ax_sec.get_window_extent()
            ve = ((box.height / abs(np.ptp(ax_sec.get_ylim())))
                  / (box.width / abs(np.ptp(ax_sec.get_xlim()))))
            if ve > 1.05:
                ax_sec.set_ylabel(f'Depth (km) — vertical exaggeration ×{ve:.2g}')
            elif ve < 1.0:
                ax_sec.set_aspect('equal', adjustable='box')
                ax_sec.set_anchor('N')

    return fig


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(params):
    """Generate and save the 2-D report figure for one event / one iteration.

    Returns
    -------
    dict with keys: output_paths, hyp_path, n_samples, maxlike_source
    """
    if params.hyp:
        hyp_path = params.hyp
        iter_match = re.search(r'loc_ssst_corr(\d+)', hyp_path)
        step = int(iter_match.group(1)) if iter_match else None
        steps = None
    else:
        iter_dirs = _find_iteration_dirs(params.ssst_root, params.run_name, params.zone)
        if not iter_dirs:
            raise FileNotFoundError(
                f'No SSST iteration folders found under '
                f'{params.ssst_root}/{params.run_name}/Pyrenees_{params.zone}_SSST')
        steps = [s for s, _ in iter_dirs]
        step = steps[-1] if params.iteration is None else params.iteration
        try:
            iter_dir = dict(iter_dirs)[step]
        except KeyError:
            raise ValueError(f'Iteration {step} not found; available: {steps}')
        hyp_path = _find_hyp_path(iter_dir, params.event_id)
        if hyp_path is None:
            raise ValueError(f'Event {params.event_id} not found in iteration {step}, zone {params.zone}')

    scat_path = _sibling(hyp_path, '.scat')
    hdr_path = _sibling(hyp_path, '.hdr')
    if not os.path.exists(scat_path):
        raise FileNotFoundError(f'{hyp_path}: no companion .scat file')

    cloud = read_nlloc_scatter(scat_path)
    center, cov = _parse_statistics(hyp_path)
    to_local, to_geo = _lambert_transformers(_read_lambert_params(hdr_path))
    info = _parse_hyp_header(hyp_path, to_local)
    info.update(to_local=to_local, to_geo=to_geo,
                zone=params.zone if params.zone is not None else '?',
                run_name=params.run_name or '')
    if params.event_id and info['public_id'] and info['public_id'] != params.event_id:
        warnings.warn(f"{hyp_path} carries {info['public_id']}, not {params.event_id}")
    if info['public_id'] is None:
        info['public_id'] = params.event_id or os.path.basename(hyp_path)

    if steps is not None and step == steps[-1]:
        iter_label = f'final pass, iteration {step}'
    elif step is not None:
        iter_label = f'iteration {step}'
    else:
        iter_label = os.path.basename(os.path.dirname(os.path.dirname(hyp_path)))

    fig = _build_figure(cloud, center, cov, info, params, iter_label)

    output = params.output or os.path.join(
        _MODULE_DIR, 'pdf_cloud', f"{info['public_id']}_report.png")
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    base, _ = os.path.splitext(output)
    outputs = [base + '.png', base + '.pdf']
    fig.savefig(outputs[0], dpi=300, bbox_inches='tight')
    fig.savefig(outputs[1], bbox_inches='tight')
    plt.close(fig)

    print(f'Figure saved @ {outputs[0]} and {outputs[1]} '
          f"({len(cloud)} samples, mode from the {info['maxlike_source']})")
    return {'output_paths': outputs, 'hyp_path': hyp_path,
            'n_samples': len(cloud), 'maxlike_source': info['maxlike_source']}


def _resolve_event_by_location(result_csv, lat, lon, date, radius_km, window_days):
    """Search a merged RESULT csv for the event nearest (lat, lon, date) within tolerance."""
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
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Report-ready 2-D view of one event NLLoc location PDF: '
                    'map panel + vertical section along the widest horizontal direction.')
    parser.add_argument('--ssst-root', default=os.path.join(_PROJECT_ROOT, 'run', 'ssst_loc'),
                        help='Root folder containing <run-name>/Pyrenees_<zone>_SSST/...')
    parser.add_argument('--run-name', default='ssst_run1',
                        help='SSST campaign name (default: ssst_run1)')
    parser.add_argument('--result-csv', default=os.path.join(_PROJECT_ROOT, 'RESULT', 'SSST_result.csv'),
                        help='Merged catalog used for --lat/--lon/--date search')
    parser.add_argument('--event-id', help='Event publicId (requires --zone)')
    parser.add_argument('--zone', help='Zone key, e.g. 1..6 (requires --event-id)')
    parser.add_argument('--lat', type=float, help='Approximate latitude for location search')
    parser.add_argument('--lon', type=float, help='Approximate longitude for location search')
    parser.add_argument('--date', help='Approximate ISO date/time for location search')
    parser.add_argument('--radius-km', type=float, default=10.0, help='Search radius for --lat/--lon (default: 10)')
    parser.add_argument('--window-days', type=float, default=5.0, help='Search time window for --date (default: 5)')
    parser.add_argument('--hyp', default=None,
                        help='Path to a .grid0.loc.hyp file, bypassing the event search entirely')
    parser.add_argument('--iteration', type=int, default=None,
                        help='SSST iteration to draw (default: the last one, i.e. the final NLLoc-only pass)')
    parser.add_argument('--confidence', type=float, default=0.6827,
                        help='Confidence level of the ellipsoid (default: 0.6827, one sigma)')
    parser.add_argument('--color', dest='color_by', choices=['log', 'linear', 'none'], default='log',
                        help="Sample colouring: 'log' reads the 4th .scat field as ln(PDF) and colours by "
                             "its offset from the mode (default), 'linear' reads it as a density, "
                             "'none' draws every sample in one colour")
    parser.add_argument('--cmap', default='viridis', help='Colormap for --color pdf (default: viridis)')
    parser.add_argument('--marginal', action='store_true',
                        help='Also draw the 2-D marginal confidence ellipse (2 d.o.f.), dashed')
    parser.add_argument('--no-section-line', dest='section_line', action='store_false',
                        help='Do not draw the A–A′ section trace on the map panel')
    parser.add_argument('--true-scale', action='store_true',
                        help='Draw the section at 1:1 instead of filling the panel; the shape of the '
                             'PDF is then undistorted, but a wide shallow cloud leaves the panel mostly empty')
    parser.add_argument('--output', default=None,
                        help='Output path; the extension is replaced to write both .png and .pdf '
                             '(default: complem_figures/pdf_cloud/<event_id>_report.png)')
    args = parser.parse_args()

    event_id, zone = args.event_id, args.zone
    if not args.hyp:
        if event_id and zone:
            pass
        elif args.lat is not None and args.lon is not None and args.date:
            event_id, zone = _resolve_event_by_location(
                args.result_csv, args.lat, args.lon, args.date, args.radius_km, args.window_days)
            print(f'Resolved to event {event_id}, zone {zone}')
        else:
            parser.error('Provide --hyp, or --event-id and --zone, or --lat, --lon and --date.')

    generate_figure(PdfCloudReportParams(
        ssst_root=args.ssst_root,
        run_name=args.run_name,
        event_id=event_id,
        zone=zone,
        hyp=args.hyp,
        iteration=args.iteration,
        confidence=args.confidence,
        color_by=args.color_by,
        cmap=args.cmap,
        marginal=args.marginal,
        section_line=args.section_line,
        true_scale=args.true_scale,
        output=args.output,
    ))


if __name__ == '__main__':
    main()
