"""
ssst_corrections.py
============================
Reconstruct and map the NonLinLoc SSST travel-time corrections, iteration by
iteration.

Loc2ssst writes an explicit correction grid per station and phase
(`ssst_corr<i>/*.ssst.*`), but run_SSST.py deletes those directories to reclaim
disk once a zone finishes, so nothing of the correction field survives a
completed campaign. It does not need to: the correction is a closed-form
function of the per-event `.hyp` locations, which DO survive under
`run/ssst_loc/<run>/Pyrenees_<N>_SSST/loc_ssst_corr<i>/`. This module recomputes
it with the formula Loc2ssst itself uses (Loc2ssst.c:1091-1102):

    corr_i(sta, pha, x,y,z) = SUM_e res_e * w_e / SUM_e w_e
    w_e                     = exp(-|x_e - x|^2 / L_i^2) + weight_floor

evaluated over the events of iteration i that pass the same LSPHSTAT selection,
with L_i = CHAR_DISTS[i]. The distance is event-to-grid-node: the station's own
position never enters the weight.

`res = obs - pred` and the correction is ADDED to the predicted travel time, so
a positive (red) correction means arrivals are later than the 1-D model
predicts, i.e. the real path is slower than the model.

Each iteration's field is an INCREMENT on top of the grids of the previous one;
the correction actually carried by the final travel-time grids is their sum,
evaluated at the final hypocentre.

Three products, from one reconstruction:

  atlas   Per station/phase: the five increments (9999 -> 50 -> 15 -> 5 -> 1 km)
          plus their cumulative total, as depth-slice maps. This is the SSST
          correction as Loc2ssst defines it - one INDEPENDENT field per station
          and phase, never a shared grid.
  spread  One catalog-wide map: at each event, the spread ACROSS its recording
          stations of the total correction applied to its picks. A correction
          common to every station of an event trades off exactly against origin
          time and cannot move the hypocentre, so the across-station spread -
          not the mean - is the part of the field that CAN relocate. It is not
          a map of where SSST actually did relocate: measured against the true
          displacement, rho = +0.10 (+0.20 given Nphs). See
          event_correction_spread for why.
  displac Where SSST actually moved events, and which way: each event's
          iteration-0 hypocentre against its final one. Same picks, same
          parameters, only the corrected grids differ, so the difference is
          the corrections and nothing else. This is the impact map the spread
          one is often mistaken for.

Extraction is cached as .npz per (zone, iteration) under
run/ssst_corrections_cache/, so only the first run pays for parsing the ~276 k
per-event .hyp files (~100 s).

Validated against the binary, not just against itself
-----------------------------------------------------
Loc2ssst was re-run for real on zone 6 (station FR.0047, P and S) at both ends
of the smoothing schedule, and its `.ssst` grid compared node by node against
this reconstruction over all 303 x 313 x 40 = 3 793 560 nodes:

    iteration 0, L = 9999 km   max |diff| = 0.0000 ms   (exact to float32)
    iteration 4, L = 1 km      max |diff| = 0.0736 ms

against correction amplitudes of +-0.36 s. Loc2ssst independently reported
"652 location files read, 163 accepted" (iteration 0) and 284 accepted
(iteration 4), matching this module's LSPHSTAT selection exactly.

Omitting LOCFILES from the control file is what makes that check cheap:
`ihave_time_input_grids = flag_out_grid * flag_nlloc_outfile` (Loc2ssst.c:590),
so with no LOCFILES statement Loc2ssst writes the correction grid and skips the
travel-time grids entirely - no Grid2Time rebuild needed.

Environment: seisbench_env (numpy/scipy, matplotlib for the figures).

Usage
-----
    python complem_figures/ssst_corrections.py                 # all three
    python complem_figures/ssst_corrections.py --product displacement
    python complem_figures/ssst_corrections.py --product atlas --min-picks 300
    python complem_figures/ssst_corrections.py --extract-only  # fill the cache
    # presentation-quality pages for chosen fields (~55 s each)
    python complem_figures/ssst_corrections.py --product atlas \\
        --stations FR.0041:P,RD.0038:P --map-spacing 0.005

Drawing note
------------
Nothing is interpolated between map nodes: pcolormesh draws one flat cell per
node and every node is an independent evaluation of the formula above. The
smoothness on the page IS the Gaussian kernel. That makes node spacing a real
constraint - the 0.02 deg default is ~2.2 km lat / ~1.6 km lon at 43N, so the
L = 1 km panel is undersampled and part of its speckle is aliasing. The full
atlas keeps the coarse spacing on purpose; --map-spacing 0.005 with --stations
resolves the finest panel for the handful of pages that need it.
"""

import argparse
import glob
import logging
import os
import sys
import time
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_SSST    = os.path.join(_PROJECT_ROOT, 'run', 'ssst_loc')
# under run/ because that path is already git-ignored, unlike *.npz
_DEFAULT_CACHE   = os.path.join(_PROJECT_ROOT, 'run', 'ssst_corrections_cache')
_DEFAULT_FIG_DIR = os.path.join(_MODULE_DIR, 'ssst_corrections')

logger = logging.getLogger('ssst_corrections')

# ---------------------------------------------------------------------------
# Campaign constants - must match what run_SSST.py actually ran with
# ---------------------------------------------------------------------------

_CHAR_DISTS   = [9999.0, 50.0, 15.0, 5.0, 1.0]   # run_SSST.CHAR_DISTS
_WEIGHT_FLOOR = 0.0000001                        # run_ssst.py LSPARAMS 2nd value
_ZONES        = [1, 2, 3, 4, 5, 6]
_N_ITER       = len(_CHAR_DISTS) + 1             # + the final NLLoc-only pass

# LSPHSTAT: RMSMax NRdgsMin GapMax PResidualMax SResidualMax EllLen3Max
_RMS_MAX      = 0.15
_NRDGS_MIN    = 6
_GAP_MAX      = 200.0
_P_RES_MAX    = 0.3
_S_RES_MAX    = 0.5
_ELL_LEN3_MAX = 10.0

_RES_MAX = {'P': _P_RES_MAX, 'S': _S_RES_MAX}

# ---------------------------------------------------------------------------
# Figure constants
# ---------------------------------------------------------------------------

_MAP_SPACING  = 0.02    # degrees, atlas map node spacing (~2.2 km lat,
                        # ~1.6 km lon at 43N). NOTE this is coarser than the
                        # L = 1 km kernel, so the finest panel is undersampled;
                        # --map-spacing renders a chosen station finer.
_SLICE_DEPTH  = None    # km; None = each station's own median event depth.
                        # A fixed slice is the wrong default: the kernel is 3-D,
                        # so a slice sitting dz off the event mass attenuates
                        # every event by exp(-dz^2/L^2). The catalog median is
                        # ~6.6-7.3 km, and slicing at 10 km costs the L = 1 km
                        # panel most of its support (median summed weight 5.7
                        # -> 2.0 for FR.0041) for no reason at all.
_MIN_EFF_N    = 1.0     # grey out nodes whose summed exponential weight is
                        # below this: no event close enough, so the value there
                        # is only the weight-floor static term
_SPREAD_BIN    = 0.05   # degrees, spread-map cell
_SPREAD_WINDOW = 2      # cells of smoothing added on each side of a cell
_SPREAD_MIN_N  = 20     # events in the window below which a cell is left blank


@dataclass
class SsstCorrParams:
    run_name:  str   = 'ssst_run1'
    ssst_root: str   = _DEFAULT_SSST
    cache_dir: str   = _DEFAULT_CACHE
    fig_dir:   str   = _DEFAULT_FIG_DIR
    zones:       list  = field(default_factory=lambda: list(_ZONES))
    min_picks:   int   = 100   # station/phase fields below this are not drawn
    depth:       float = _SLICE_DEPTH   # None = per-station median event depth
    map_spacing: float = _MAP_SPACING
    extent:      tuple = None   # (lon0, lon1, lat0, lat1); None = per-station
    product:     str   = 'all'   # atlas | spread | displacement | all


# ---------------------------------------------------------------------------
# .hyp parsing
# ---------------------------------------------------------------------------

def parse_hyp_file(path):
    """Yield one dict per event block of an NLLoc .hyp file.

    Only what Loc2ssst reads is extracted: the hypocentre in the zone's
    rectangular frame (what latlon2rect hands it), the four event-level LSPHSTAT
    quantities, and every arrival's station / phase / residual.
    """
    with open(path) as f:
        lines = f.read().splitlines()

    ev = None
    in_phase = False
    for line in lines:
        if line.startswith('NLLOC '):
            parts = line.split('"')
            ev = {'status': parts[3] if len(parts) > 3 else '?',
                  'id': '', 'arr': []}
            in_phase = False
        elif ev is None:
            continue
        elif line.startswith('PUBLIC_ID'):
            # the pipeline's canonical event key; the only safe way to follow
            # one event across iterations, since the per-iteration event sets
            # differ slightly as locations flip between LOCATED and REJECTED
            ev['id'] = line.split()[1]
        elif line.startswith('HYPOCENTER'):
            f_ = line.split()
            ev['x'], ev['y'], ev['z'] = float(f_[2]), float(f_[4]), float(f_[6])
        elif line.startswith('GEOGRAPHIC'):
            f_ = line.split()
            ev['lat'] = float(f_[f_.index('Lat') + 1])
            ev['lon'] = float(f_[f_.index('Long') + 1])
        elif line.startswith('QUALITY'):
            f_ = line.split()
            ev['rms']  = float(f_[f_.index('RMS') + 1])
            ev['nphs'] = int(f_[f_.index('Nphs') + 1])
            ev['gap']  = float(f_[f_.index('Gap') + 1])
        elif line.startswith('STATISTICS'):
            f_ = line.split()
            ev['len3'] = float(f_[f_.index('Len3') + 1])
        elif line.startswith('PHASE ID'):
            in_phase = True
        elif line.startswith('END_PHASE'):
            in_phase = False
        elif line.startswith('END_NLLOC'):
            yield ev
            ev = None
        elif in_phase:
            f_ = line.split()
            ev['arr'].append((f_[0], f_[4], float(f_[17])))


def iteration_dir(params, zone, iteration):
    return os.path.join(params.ssst_root, params.run_name,
                        f'Pyrenees_{zone}_SSST', f'loc_ssst_corr{iteration}',
                        f'GLOBAL_{zone}')


def event_hyp_files(directory, zone):
    """Per-event .hyp files of one iteration, excluding the .sum summary.

    The .sum file is a concatenation of the per-event hypocentre blocks WITHOUT
    their phase lines, so it carries no residuals and cannot substitute here.
    """
    pattern = os.path.join(directory, f'Pyrenees_{zone}.*.grid0.loc.hyp')
    return sorted(p for p in glob.glob(pattern)
                  if '.sum.' not in os.path.basename(p))


def extract(params, zone, iteration):
    """Parse one (zone, iteration) into flat arrays, caching the result.

    Events keep their LSPHSTAT verdict in `ev_ok` rather than being dropped, so
    one cache serves both the correction fields (accepted events only) and the
    maps (all events). Station labels are stored as int32 indices into
    `sta_names`: as fixed-width unicode they would cost ~28 bytes per arrival
    over ~7 M arrivals, for no gain.
    """
    os.makedirs(params.cache_dir, exist_ok=True)
    cache = os.path.join(params.cache_dir,
                         f'{params.run_name}_z{zone}_i{iteration}.npz')
    if os.path.exists(cache):
        with np.load(cache, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}

    directory = iteration_dir(params, zone, iteration)
    files = event_hyp_files(directory, zone)
    if not files:
        raise FileNotFoundError(f'no per-event .hyp files in {directory}')

    t0 = time.time()
    ex, ey, ez, elat, elon, eid = [], [], [], [], [], []
    erms, enphs, egap, elen3, eok = [], [], [], [], []
    ar_ev, ar_sta, ar_pha, ar_res = [], [], [], []
    sta_ids = {}

    for path in files:
        for ev in parse_hyp_file(path):
            # Loc2ssst skips ABORTED and, with use_rejected at its default 0,
            # REJECTED locations too
            if ev['status'] != 'LOCATED':
                continue
            idx = len(ex)
            ex.append(ev['x']); ey.append(ev['y']); ez.append(ev['z'])
            elat.append(ev['lat']); elon.append(ev['lon']); eid.append(ev['id'])
            erms.append(ev['rms']); enphs.append(ev['nphs'])
            egap.append(ev['gap']); elen3.append(ev['len3'])
            eok.append(ev['rms'] <= _RMS_MAX and ev['nphs'] >= _NRDGS_MIN
                       and ev['gap'] <= _GAP_MAX and ev['len3'] <= _ELL_LEN3_MAX)
            for sta, pha, res in ev['arr']:
                sid = sta_ids.get(sta)
                if sid is None:
                    sid = sta_ids[sta] = len(sta_ids)
                ar_ev.append(idx); ar_sta.append(sid)
                ar_pha.append(pha == 'P'); ar_res.append(res)

    names = np.empty(len(sta_ids), dtype='<U16')
    for sta, sid in sta_ids.items():
        names[sid] = sta

    out = {
        'ev_x':    np.array(ex, dtype=np.float32),
        'ev_y':    np.array(ey, dtype=np.float32),
        'ev_z':    np.array(ez, dtype=np.float32),
        'ev_lat':  np.array(elat, dtype=np.float64),
        'ev_lon':  np.array(elon, dtype=np.float64),
        'ev_id':   np.array(eid),
        'ev_rms':  np.array(erms, dtype=np.float32),
        'ev_nphs': np.array(enphs, dtype=np.int32),
        'ev_gap':  np.array(egap, dtype=np.float32),
        'ev_len3': np.array(elen3, dtype=np.float32),
        'ev_ok':   np.array(eok, dtype=bool),
        'ar_ev':   np.array(ar_ev, dtype=np.int32),
        'ar_sta':  np.array(ar_sta, dtype=np.int32),
        'ar_is_p': np.array(ar_pha, dtype=bool),
        'ar_res':  np.array(ar_res, dtype=np.float32),
        'sta_names': names,
    }
    np.savez_compressed(cache, **out)
    logger.info('zone %d iter %d: %d events (%d pass LSPHSTAT), %d arrivals, '
                '%d stations [%.1f s]', zone, iteration, len(ex),
                int(out['ev_ok'].sum()), len(ar_ev), len(names),
                time.time() - t0)
    return out


def load_all(params, iterations=range(_N_ITER)):
    """Cache-backed {(zone, iteration): arrays} for every requested iteration."""
    return {(z, i): extract(params, z, i)
            for z in params.zones for i in iterations}


# ---------------------------------------------------------------------------
# Geometry: lat/lon <-> the zone's own rectangular frame
# ---------------------------------------------------------------------------

def _design(lat, lon):
    lat = np.asarray(lat, dtype=np.float64)
    lon = np.asarray(lon, dtype=np.float64)
    return np.column_stack([np.ones_like(lat), lat, lon,
                            lat * lat, lat * lon, lon * lon])


def fit_projection(data):
    """Least-squares quadratic mapping (lat, lon) -> (x, y) for one zone.

    NLLoc prints both the geographic and the rectangular hypocentre of every
    event, so the zone's Lambert projection can be recovered from the events
    themselves instead of re-deriving it from the TRANS line. Over one zone a
    quadratic reproduces it to ~2 m, far below anything that matters here, and
    it cannot drift out of sync with the run's own parameters.
    """
    A = _design(data['ev_lat'], data['ev_lon'])
    cx, *_ = np.linalg.lstsq(A, data['ev_x'].astype(np.float64), rcond=None)
    cy, *_ = np.linalg.lstsq(A, data['ev_y'].astype(np.float64), rcond=None)
    return cx, cy


def project(coef, lat, lon):
    cx, cy = coef
    A = _design(lat, lon)
    return A @ cx, A @ cy


def rect_coords(data):
    """Event coordinates in the zone's rectangular frame, km.

    Deliberately NOT the `HYPOCENTER x y` values NLLoc printed. Loc2ssst
    ignores those too and re-derives the rectangular hypocentre from the
    geographic one with latlon2rect (Loc2ssst.c:734). It matters: the printed
    rect coordinates carry 1 mm of precision, which is 0.7 m of disagreement
    with the re-derived ones, and at char_dist = 1 km that dominates the whole
    reconstruction error (0.26 ms against 0.07 ms - see the validation note in
    the module docstring).
    """
    if '_xyz' not in data:
        x, y = project(fit_projection(data), data['ev_lat'], data['ev_lon'])
        data['_xyz'] = np.column_stack([x, y, data['ev_z'].astype(np.float64)])
    return data['_xyz']


# ---------------------------------------------------------------------------
# The correction field
# ---------------------------------------------------------------------------

def station_index(data, station):
    """Index of a station label in a cache's name table, or None."""
    hit = np.flatnonzero(data['sta_names'] == station)
    return int(hit[0]) if len(hit) else None


def station_arrivals(data, station, phase):
    """Sources feeding one station/phase field: (event indices, xyz, residuals).

    Reproduces Loc2ssst's two-level selection: the event must pass the
    event-level LSPHSTAT cuts (precomputed in `ev_ok`), and the arrival's
    residual must be within its phase's ceiling (Loc2ssst.c:873). Duplicate
    arrivals are kept, as Loc2ssst keeps them (addDuplicates = 1), so the event
    index array can repeat.
    """
    sid = station_index(data, station)
    if sid is None:
        return np.empty(0, dtype=int), np.empty((0, 3)), np.empty(0)

    sel = ((data['ar_sta'] == sid)
           & (data['ar_is_p'] == (phase == 'P'))
           & (np.abs(data['ar_res']) <= _RES_MAX[phase]))
    ev = data['ar_ev'][sel]
    keep = data['ev_ok'][ev]
    ev = ev[keep]
    return ev, rect_coords(data)[ev], data['ar_res'][sel][keep].astype(np.float64)


def correction_field(points, src_xyz, residuals, char_dist,
                     weight_floor=_WEIGHT_FLOOR, return_weight=False):
    """Evaluate one SSST correction field at arbitrary points.

    `points` and `src_xyz` are (n, 3) arrays of km in the same rectangular
    frame; `residuals` holds one value per row of `src_xyz`.

    The weight floor is what makes the field fall back to the station's plain
    static term far from any event instead of going undefined, so it is kept
    even though it is tiny. `return_weight` also hands back the summed
    exponential weight, which is the effective number of events supporting each
    point and the basis for masking.
    """
    points = np.asarray(points, dtype=np.float64)
    src_xyz = np.asarray(src_xyz, dtype=np.float64)
    residuals = np.asarray(residuals, dtype=np.float64)

    n_pts, n_src = len(points), len(src_xyz)
    num = np.zeros(n_pts)          # SUM res * exp
    den = np.zeros(n_pts)          # SUM exp
    if n_src == 0:
        return (num, den) if return_weight else num

    l2 = char_dist * char_dist
    # exp(-d^2/L^2) falls below 1e-15 past 6 L, so for a small characteristic
    # distance only a neighbourhood contributes and a tree pays for itself; for
    # a large one every event is a neighbour and the tree is pure overhead.
    span = float(np.ptp(src_xyz, axis=0).max())
    if 6.0 * char_dist < span:
        tree = cKDTree(src_xyz)
        for i, idx in enumerate(tree.query_ball_point(points, 6.0 * char_dist)):
            if not idx:
                continue
            idx = np.asarray(idx)
            w = np.exp(-np.sum((src_xyz[idx] - points[i]) ** 2, axis=1) / l2)
            num[i] = np.dot(residuals[idx], w)
            den[i] = w.sum()
    else:
        chunk = max(1, int(1e7 // n_src))
        for a in range(0, n_pts, chunk):
            b = min(a + chunk, n_pts)
            d2 = np.sum((points[a:b, None, :] - src_xyz[None, :, :]) ** 2, axis=2)
            w = np.exp(-d2 / l2)
            num[a:b] = w @ residuals
            den[a:b] = w.sum(axis=1)

    value = ((num + weight_floor * residuals.sum())
             / (den + weight_floor * n_src))
    return (value, den) if return_weight else value


# ---------------------------------------------------------------------------
# Station census
# ---------------------------------------------------------------------------

def station_census(cache, zones):
    """Usable-arrival count per (station, phase, zone), at iteration 0."""
    census = {}
    for zone in zones:
        data = cache[(zone, 0)]
        ok = data['ev_ok'][data['ar_ev']]
        names = data['sta_names']
        for phase in ('P', 'S'):
            sel = (ok & (data['ar_is_p'] == (phase == 'P'))
                   & (np.abs(data['ar_res']) <= _RES_MAX[phase]))
            ids, counts = np.unique(data['ar_sta'][sel], return_counts=True)
            for sid, n in zip(ids, counts):
                census[(str(names[sid]), phase, zone)] = int(n)
    return census


def best_zone_per_station(census, min_picks):
    """Keep each (station, phase) in the zone where it has the most arrivals.

    A station inside the overlap of two zones has a different, independent
    correction field in each - separate Loc2ssst runs over different event sets.
    Drawing both would suggest a disagreement that is really two different
    questions; the better-sampled zone is the honest single answer.

    Returns (station, phase, zone, n_arrivals) sorted by decreasing n.
    """
    best = {}
    for (sta, phase, zone), n in census.items():
        key = (sta, phase)
        if n >= min_picks and n > best.get(key, (0, None))[0]:
            best[key] = (n, zone)
    return sorted(((sta, ph, z, n) for (sta, ph), (n, z) in best.items()),
                  key=lambda t: -t[3])


def filter_fields(fields_list, wanted):
    """Restrict the atlas to named stations, as 'STA' or 'STA:PHASE'.

    The full 300-page atlas is deliberately drawn at the coarse default
    spacing; this is what lets a chosen station be re-rendered fine enough to
    resolve the L = 1 km panel without paying for all 300.
    """
    if not wanted:
        return fields_list
    keys = set()
    for item in wanted:
        sta, _, phase = item.partition(':')
        keys.update((sta, p) for p in ((phase,) if phase else ('P', 'S')))
    return [f for f in fields_list if (f[0], f[1]) in keys]


# ---------------------------------------------------------------------------
# Station coordinates
# ---------------------------------------------------------------------------

def read_lsgrid(zone):
    """(x0, x1, y0, y1) km of the zone's Loc2ssst correction grid.

    Read from the control file rather than hardcoded, so it cannot drift away
    from what the campaign actually ran. Outside this box Loc2ssst never
    evaluated a correction at all, which is what a shared map extent has to
    respect: the field is not merely unsupported there, it does not exist.
    """
    path = os.path.join(_PROJECT_ROOT, 'run', 'ssst', f'run_{zone}_SSST.in')
    with open(path) as f:
        for line in f:
            v = line.split()
            if v and v[0] == 'LSGRID':
                nx, ny = int(v[1]), int(v[2])
                x0, y0 = float(v[4]), float(v[5])
                dx, dy = float(v[7]), float(v[8])
                return x0, x0 + (nx - 1) * dx, y0, y0 + (ny - 1) * dy
    raise ValueError(f'no LSGRID statement in {path}')


def read_station_latlon(zone):
    """{station code: (lat, lon)} from the zone's GTSRCE file."""
    path = os.path.join(_PROJECT_ROOT, 'stations', f'GTSRCE_SSST_{zone}.txt')
    coords = {}
    with open(path) as f:
        for line in f:
            f_ = line.split()
            if len(f_) >= 5 and f_[0] == 'GTSRCE' and f_[2] == 'LATLON':
                coords[f_[1]] = (float(f_[3]), float(f_[4]))
    return coords


# ---------------------------------------------------------------------------
# Product 1: the per-station atlas
# ---------------------------------------------------------------------------

def _node_spacing_km(lats, lons):
    """Map node spacing in km, to compare against the smoothing length."""
    deg_lat = float(np.diff(lats).mean()) if len(lats) > 1 else 0.0
    deg_lon = float(np.diff(lons).mean()) if len(lons) > 1 else 0.0
    mid = np.deg2rad(float(np.mean(lats)))
    return deg_lat * 111.32, deg_lon * 111.32 * np.cos(mid)


def _map_grid(src_list, spacing=_MAP_SPACING, pad=0.15):
    """Regular lat/lon grid covering every source cloud of a station, padded."""
    lat = np.concatenate([s[0] for s in src_list])
    lon = np.concatenate([s[1] for s in src_list])
    lat0, lat1 = lat.min() - pad, lat.max() + pad
    lon0, lon1 = lon.min() - pad, lon.max() + pad
    lats = np.arange(lat0, lat1 + spacing, spacing)
    lons = np.arange(lon0, lon1 + spacing, spacing)
    return lats, lons


def station_fields(cache, zone, station, phase, depth=None, spacing=_MAP_SPACING,
                   extent=None):
    """The five increments, their cumulative sum, and the map they live on.

    Every field is evaluated on ONE grid at ONE depth, so the panels are
    directly comparable; each increment comes from its own iteration's
    locations and residuals, which is what makes them increments rather than
    five estimates of the same thing.

    `depth` None slices at the median depth of the events that actually feed
    this station's field. The kernel is 3-D, so a slice dz away from the event
    mass damps every contribution by exp(-dz^2/L^2): harmless at L = 15 km,
    but it empties the L = 1 km panel. Slicing where the events are is the only
    choice that shows the fine panels the support they really have.

    `extent` (lon0, lon1, lat0, lat1) draws every page on one shared frame so
    they can be compared side by side, instead of each page framing its own
    station's events. It is clipped to the zone's LSGRID: outside that box
    Loc2ssst computed nothing, and painting the station's static term across it
    would present a constant as if it were a measurement.
    """
    per_iter = []
    for i in range(len(_CHAR_DISTS)):
        data = cache[(zone, i)]
        ev, src, res = station_arrivals(data, station, phase)
        per_iter.append((data, ev, src, res))

    latlon = [(data['ev_lat'][ev], data['ev_lon'][ev])
              for data, ev, src, _ in per_iter if len(src)]
    if not latlon:
        return None
    if extent is None:
        lats, lons = _map_grid(latlon, spacing=spacing)
    else:
        lon0, lon1, lat0, lat1 = extent
        lats = np.arange(lat0, lat1 + spacing, spacing)
        lons = np.arange(lon0, lon1 + spacing, spacing)

    if depth is None:
        depth = float(np.median(np.concatenate(
            [src[:, 2] for _, _, src, _ in per_iter if len(src)])))

    LON, LAT = np.meshgrid(lons, lats)
    coef = fit_projection(cache[(zone, 0)])
    gx, gy = project(coef, LAT.ravel(), LON.ravel())
    points = np.column_stack([gx, gy, np.full(gx.size, depth)])

    x0, x1, y0, y1 = read_lsgrid(zone)
    outside = ((gx < x0) | (gx > x1) | (gy < y0) | (gy > y1)).reshape(LAT.shape)

    increments, weights = [], []
    for (data, ev, src, res), char_dist in zip(per_iter, _CHAR_DISTS):
        value, weight = correction_field(points, src, res, char_dist,
                                         return_weight=True)
        increments.append(value.reshape(LAT.shape))
        # nodes off the zone's correction grid are marked unsupported, which
        # is what puts them in the grey the panels already use for "no data"
        weights.append(np.where(outside, 0.0, weight.reshape(LAT.shape)))

    return {
        'lats': lats, 'lons': lons, 'depth': depth, 'outside': outside,
        'increments': increments, 'weights': weights,
        'total': np.sum(increments, axis=0),
        'n_arrivals': [len(r) for _, _, _, r in per_iter],
    }


def _draw_map(ax, lats, lons, values, mask, vmax, cmap, station_lat, station_lon):
    import matplotlib.pyplot as plt   # noqa: F401  (kept local, see main)

    plotted = np.ma.masked_where(mask, values)
    # rasterized: 300 pages x 6 vector meshes is a ~90 MB PDF, ~7 MB rasterized,
    # and the meshes carry no detail that survives as vectors anyway
    mesh = ax.pcolormesh(lons, lats, plotted, cmap=cmap, rasterized=True,
                         vmin=-vmax, vmax=vmax, shading='auto')
    ax.set_facecolor('0.88')          # shows through wherever the field is masked
    ax.plot(station_lon, station_lat, marker='v', color='k',
            markersize=9, markeredgecolor='w', markeredgewidth=1.0, zorder=5)
    ax.set_aspect(1.0 / np.cos(np.deg2rad(float(np.mean(lats)))))
    ax.tick_params(labelsize=7)
    return mesh


def atlas_page(fig, fields, station, phase, zone, station_lat, station_lon):
    """One atlas page: five increments, then their cumulative total."""
    increments, weights = fields['increments'], fields['weights']
    lats, lons = fields['lats'], fields['lons']

    masks = [w < _MIN_EFF_N for w in weights]
    stack = np.concatenate([np.asarray(v)[~m].ravel()
                            for v, m in zip(increments, masks) if (~m).any()])
    vmax_inc = float(np.percentile(np.abs(stack), 98)) if stack.size else 1e-3
    inside_total = fields['total'][~fields['outside']]
    vmax_tot = (float(np.percentile(np.abs(inside_total), 98))
                if inside_total.size else 1e-3) or 1e-3

    axes = [fig.add_subplot(2, 3, k + 1) for k in range(6)]
    mesh_inc = None
    for k, (value, mask, char_dist) in enumerate(zip(increments, masks, _CHAR_DISTS)):
        mesh_inc = _draw_map(axes[k], lats, lons, value, mask, vmax_inc,
                             'RdBu_r', station_lat, station_lon)
        label = ('L = 9999 km  (static)' if char_dist > 1000
                 else f'L = {char_dist:g} km')
        # the five panels share one colour scale, so amplitude is comparable
        # between them by eye - but the quiet ones then look empty, and their
        # rms is the only way to read how quiet they actually are
        rms = float(np.sqrt((np.asarray(value)[~mask] ** 2).mean())) if (~mask).any() else 0.0
        axes[k].set_title(f'{label}\niteration {k}, {fields["n_arrivals"][k]} arrivals'
                          f'   rms {rms * 1000:.0f} ms', fontsize=9)

    # the cumulative panel is defined wherever the coarse terms are, i.e.
    # everywhere on the zone's grid - but still nowhere off it
    mesh_tot = _draw_map(axes[5], lats, lons, fields['total'],
                         fields['outside'],
                         vmax_tot, 'RdBu_r', station_lat, station_lon)
    axes[5].set_title('CUMULATIVE  (sum of the five)\ncorrection in the final grids',
                      fontsize=9, fontweight='bold')

    fig.suptitle(f'SSST travel-time correction   {station}  {phase}-phase   '
                 f'zone {zone}   depth slice {fields["depth"]:.1f} km',
                 fontsize=13, fontweight='bold')
    fig.subplots_adjust(left=0.05, right=0.98, top=0.86, bottom=0.17,
                        wspace=0.18, hspace=0.32)

    # colourbar labels go ABOVE the bars, leaving the strip underneath free for
    # the footnote - otherwise the two collide at this figure size
    for cax_rect, mesh, label in (
            ([0.08, 0.095, 0.42, 0.018], mesh_inc,
             'increment (s), shared scale across the five panels'),
            ([0.62, 0.095, 0.30, 0.018], mesh_tot, 'total correction (s)')):
        cax = fig.add_axes(cax_rect)
        cb = fig.colorbar(mesh, cax=cax, orientation='horizontal')
        cb.set_label(label, fontsize=8)
        cb.ax.xaxis.set_label_position('top')
        cb.ax.tick_params(labelsize=7)

    lat_km, lon_km = _node_spacing_km(fields['lats'], fields['lons'])
    fig.text(0.08, 0.010,
             'Axes are longitude / latitude (deg); the triangle is the station.  '
             'red = arrivals LATER than the 1-D model predicts (real path slower); '
             'blue = earlier.\n'
             'grey = no event within reach of the smoothing kernel - the field '
             'there is only the station static term.\n'
             f'Nothing is interpolated: every cell is one independent evaluation, '
             f'{lat_km:.1f} km x {lon_km:.1f} km apart - the smoothness is the '
             f'kernel itself. Panels with L below that spacing are undersampled.',
             fontsize=7.5, color='0.30', linespacing=1.5)


def make_atlas(params, cache, fields_list, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    station_coords = {z: read_station_latlon(z) for z in params.zones}

    drawn = 0
    with PdfPages(path) as pdf:
        for station, phase, zone, n in fields_list:
            fields = station_fields(cache, zone, station, phase, params.depth,
                                    params.map_spacing, params.extent)
            if fields is None:
                continue
            lat, lon = station_coords[zone].get(station, (np.nan, np.nan))
            fig = plt.figure(figsize=(15, 9))
            try:
                atlas_page(fig, fields, station, phase, zone, lat, lon)
                pdf.savefig(fig)
            finally:
                plt.close(fig)
            drawn += 1
            if drawn % 25 == 0:
                logger.info('  ... %d / %d pages', drawn, len(fields_list))
    logger.info('atlas: %d pages -> %s', drawn, path)


# ---------------------------------------------------------------------------
# Product 2: the catalog-wide across-station spread map
# ---------------------------------------------------------------------------

def event_correction_spread(cache, zone, min_picks_per_event=5,
                            return_ids=False):
    """Per-event spread, across recording stations, of the total correction.

    The total correction carried by the final travel-time grids is the sum of
    the five increments; it is evaluated here at each event's FINAL location,
    which is where its picks were actually corrected.

    A correction common to every station of an event is absorbed exactly by the
    origin time and cannot move the hypocentre, so the dispersion across
    stations - not the mean - is the part of the field that relocates.

    P phases only: S residuals are systematically larger, so mixing the two
    would report phase-dependent amplitude as if it were spatial variation.

    NOT a map of where SSST had most impact, however tempting that reading is.
    Measured against the pure SSST displacement (iteration 0 location vs final
    location - same picks, only the grids differ) over 39 080 events:

        rho(spread, displacement)         = +0.10
        rho(spread, displacement | Nphs)  = +0.20
        median displacement, bottom spread quintile -> top: 0.95 km -> 1.27 km

    Disagreement is necessary for an event to move but nowhere near sufficient:
    what actually displaces a hypocentre is the AZIMUTHAL PATTERN of the
    differential correction, and a standard deviation discards all of that
    geometry. Spread evenly around the azimuths largely cancels. (The marginal
    rho understates even this much because Nphs suppresses it - well-recorded
    events carry more spread, rho = +0.42, yet resist moving, rho = -0.18.)
    To show impact, map the displacement itself.
    """
    final = cache[(zone, len(_CHAR_DISTS))]
    final_xyz = rect_coords(final)

    sel = final['ar_is_p'] & (np.abs(final['ar_res']) <= _P_RES_MAX)
    ar_ev = final['ar_ev'][sel]
    ar_sta = final['ar_sta'][sel]
    total = np.zeros(len(ar_ev))

    for sid in np.unique(ar_sta):
        station = str(final['sta_names'][sid])
        rows = np.flatnonzero(ar_sta == sid)
        points = final_xyz[ar_ev[rows]]
        for i, char_dist in enumerate(_CHAR_DISTS):
            _, src, res = station_arrivals(cache[(zone, i)], station, 'P')
            if len(src):
                total[rows] += correction_field(points, src, res, char_dist)

    order = np.argsort(ar_ev, kind='stable')
    ar_ev, total = ar_ev[order], total[order]
    bounds = np.flatnonzero(np.diff(ar_ev)) + 1
    groups = np.split(np.arange(len(ar_ev)), bounds)

    ev_idx, spread = [], []
    for g in groups:
        if len(g) >= min_picks_per_event:
            ev_idx.append(ar_ev[g[0]])
            spread.append(total[g].std())

    ev_idx = np.asarray(ev_idx, dtype=int)
    out = (final['ev_lat'][ev_idx], final['ev_lon'][ev_idx], np.asarray(spread))
    return out + (final['ev_id'][ev_idx],) if return_ids else out


def windowed_median_map(lat, lon, values, bin_size=_SPREAD_BIN,
                        radius=_SPREAD_BIN * _SPREAD_WINDOW,
                        min_count=_SPREAD_MIN_N):
    """Windowed median of a per-event quantity, on a lat/lon grid.

    A single cell holds a handful of events, so a raw per-cell statistic is
    dominated by sampling noise and the map reads as speckle; averaging over a
    window first is what separates real regional structure from that noise.

    Median rather than mean, and a circular window rather than a square one:
    displacement is heavy-tailed (one event moves 99.5 km), so a mean would let
    single events paint whole neighbourhoods, and a degree-square window is
    half as wide in km at this latitude as it is tall.
    """
    lat_bins = np.arange(lat.min(), lat.max() + bin_size, bin_size)
    lon_bins = np.arange(lon.min(), lon.max() + bin_size, bin_size)
    lat_c = lat_bins[:-1] + bin_size / 2
    lon_c = lon_bins[:-1] + bin_size / 2

    # scale longitude so the window is round in km, not in degrees
    scale = np.cos(np.deg2rad(float(np.mean(lat))))
    tree = cKDTree(np.column_stack([lat, lon * scale]))
    LAT, LON = np.meshgrid(lat_c, lon_c, indexing='ij')
    centres = np.column_stack([LAT.ravel(), LON.ravel() * scale])

    grid = np.full(len(centres), np.nan)
    for i, idx in enumerate(tree.query_ball_point(centres, radius)):
        if len(idx) >= min_count:
            grid[i] = np.median(values[idx])
    return lat_bins, lon_bins, grid.reshape(LAT.shape)


def event_displacement(cache, zone):
    """How far, and which way, SSST actually moved each event.

    Iteration 0 located with UNcorrected grids; the final iteration located
    with the fully corrected ones. Same picks, same control parameters, same
    everything else - so the difference between those two hypocentres is the
    SSST corrections and nothing else. That is the honest impact measure, and
    it is not what the spread map shows (see event_correction_spread).

    Matching is by PUBLIC_ID: the per-iteration event sets differ slightly,
    since locations flip between LOCATED and REJECTED as the grids improve.

    dx/dy are the zone's rectangular axes. TRANS sets RotCW 0, so x is east and
    y is north to within the convergence of the Lambert projection - fine for
    drawing arrows, not for measuring an azimuth.
    """
    first, final = cache[(zone, 0)], cache[(zone, len(_CHAR_DISTS))]
    ids, i0, i5 = np.intersect1d(first['ev_id'], final['ev_id'],
                                 return_indices=True)
    delta = rect_coords(final)[i5] - rect_coords(first)[i0]
    return {
        'lat': final['ev_lat'][i5], 'lon': final['ev_lon'][i5], 'ids': ids,
        'dist3d': np.linalg.norm(delta, axis=1),
        'dx': delta[:, 0], 'dy': delta[:, 1], 'dz': delta[:, 2],
    }


def make_displacement_map(params, cache, path):
    """Two panels: how far SSST moved events, and which way in depth."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    parts = [event_displacement(cache, z) for z in params.zones]
    for zone, d in zip(params.zones, parts):
        logger.info('  zone %d: %d events matched across iterations, '
                    'median move %.2f km', zone, len(d['ids']),
                    float(np.median(d['dist3d'])))
    lat = np.concatenate([d['lat'] for d in parts])
    lon = np.concatenate([d['lon'] for d in parts])
    dist = np.concatenate([d['dist3d'] for d in parts])
    dz = np.concatenate([d['dz'] for d in parts])
    dx = np.concatenate([d['dx'] for d in parts])
    dy = np.concatenate([d['dy'] for d in parts])

    lat_bins, lon_bins, g_dist = windowed_median_map(lat, lon, dist)
    _, _, g_dz = windowed_median_map(lat, lon, dz)
    # arrows on a coarser grid than the colour, or they overplot into a smear
    q_lat, q_lon, g_dx = windowed_median_map(lat, lon, dx, bin_size=0.25,
                                             radius=0.20, min_count=_SPREAD_MIN_N)
    _, _, g_dy = windowed_median_map(lat, lon, dy, bin_size=0.25,
                                     radius=0.20, min_count=_SPREAD_MIN_N)

    fig, axes = plt.subplots(2, 1, figsize=(13, 12))
    aspect = 1.0 / np.cos(np.deg2rad(float(np.mean(lat))))

    vmax = float(np.nanpercentile(g_dist, 98))
    mesh = axes[0].pcolormesh(lon_bins, lat_bins, np.ma.masked_invalid(g_dist),
                              cmap='inferno_r', shading='flat',
                              vmin=float(np.nanpercentile(g_dist, 2)), vmax=vmax,
                              rasterized=True)
    cb = fig.colorbar(mesh, ax=axes[0], shrink=0.9, extend='both')
    cb.set_label('median distance moved (km)')
    qy = q_lat[:-1] + 0.125
    qx = q_lon[:-1] + 0.125
    QX, QY = np.meshgrid(qx, qy)
    ok = np.isfinite(g_dx) & np.isfinite(g_dy)
    axes[0].quiver(QX[ok], QY[ok], g_dx[ok], g_dy[ok], color='#00d2ff',
                   angles='xy', scale=12, width=0.0035,
                   edgecolor='k', linewidth=0.4)
    axes[0].set_title('How far SSST actually moved each event\n'
                      '(iteration 0 vs final location - same picks, only the '
                      'corrected grids differ)', fontsize=11, fontweight='bold')
    axes[0].text(0.0, -0.115, 'arrows: median horizontal displacement per '
                 '0.25 deg cell', transform=axes[0].transAxes,
                 fontsize=8, color='0.30')

    vz = float(np.nanpercentile(np.abs(g_dz), 98))
    mesh = axes[1].pcolormesh(lon_bins, lat_bins, np.ma.masked_invalid(g_dz),
                              cmap='RdBu_r', shading='flat', vmin=-vz, vmax=vz,
                              rasterized=True)
    cb = fig.colorbar(mesh, ax=axes[1], shrink=0.9, extend='both')
    cb.set_label('median change in depth (km)')
    axes[1].set_title('Which way in depth\n'
                      'red = SSST pushed events DEEPER, blue = shallower',
                      fontsize=11, fontweight='bold')
    axes[1].text(0.0, -0.19,
                 f'{len(dist)} events. Median move {np.median(dist):.2f} km, '
                 f'p90 {np.percentile(dist, 90):.2f} km, max {dist.max():.1f} km.  '
                 f'{_SPREAD_BIN:g} deg cells, windowed median over '
                 f'{_SPREAD_BIN * _SPREAD_WINDOW:g} deg, blank below '
                 f'{_SPREAD_MIN_N} events.',
                 transform=axes[1].transAxes, fontsize=8, color='0.30')

    for ax in axes:
        ax.set_facecolor('0.92')
        ax.set_aspect(aspect)
        ax.set_ylabel('latitude (deg)')
    axes[1].set_xlabel('longitude (deg)')   # top panel's would land on its note
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logger.info('displacement map: %d events -> %s', len(dist), path)


def make_spread_map(params, cache, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    lat, lon, spread = [], [], []
    for zone in params.zones:
        a, b, c = event_correction_spread(cache, zone)
        lat.append(a); lon.append(b); spread.append(c)
        logger.info('  zone %d: %d events with a spread', zone, len(c))
    lat = np.concatenate(lat); lon = np.concatenate(lon)
    spread = np.concatenate(spread)

    # same helper as the displacement map, so the two bin space identically
    lat_bins, lon_bins, mean = windowed_median_map(lat, lon, spread)

    # the spread never approaches zero (its median is ~0.1 s), so anchoring the
    # scale at 0 would spend the whole colormap on an empty range and flatten
    # every real contrast; clip to the bulk of the distribution instead
    vmin, vmax = (float(v) for v in np.nanpercentile(mean, [2, 98]))

    fig, ax = plt.subplots(figsize=(13, 7))
    mesh = ax.pcolormesh(lon_bins, lat_bins, np.ma.masked_invalid(mean),
                         cmap='viridis', shading='flat', vmin=vmin, vmax=vmax)
    ax.set_facecolor('0.92')
    ax.set_aspect(1.0 / np.cos(np.deg2rad(float(np.mean(lat)))))
    ax.set_xlabel('longitude (deg)')
    ax.set_ylabel('latitude (deg)')
    ax.set_title('SSST correction: across-station spread of the total applied '
                 'correction\n(the component that CAN move a hypocentre - not a '
                 'map of how far one did; P phases)',
                 fontsize=11, fontweight='bold')
    cb = fig.colorbar(mesh, ax=ax, shrink=0.85, extend='both')
    cb.set_label('std across recording stations (s)')
    ax.text(0.0, -0.15,
            f'{len(spread)} events, median spread {np.median(spread) * 1000:.0f} ms.  '
            f'{_SPREAD_BIN:g} deg cells, windowed median over '
            f'{_SPREAD_BIN * _SPREAD_WINDOW:g} deg, blank below '
            f'{_SPREAD_MIN_N} events.  Colour clipped to the 2nd-98th percentile '
            f'({vmin * 1000:.0f}-{vmax * 1000:.0f} ms).\n'
            f'This is what the corrections CAN do, not what they did: against '
            f'the real displacement, rho = +0.10 (+0.20 given Nphs). '
            f'See the displacement map.',
            transform=ax.transAxes, fontsize=8, color='0.30', linespacing=1.6)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    logger.info('spread map: %d events -> %s', len(spread), path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _setup_logger(verbose=True):
    logging.basicConfig(level=logging.INFO if verbose else logging.WARNING,
                        format='%(message)s', stream=sys.stdout)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run-name', default='ssst_run1')
    p.add_argument('--zones', default=None,
                   help='comma-separated zone list (default: all six)')
    p.add_argument('--min-picks', type=int, default=100,
                   help='skip station/phase fields with fewer usable arrivals')
    p.add_argument('--depth', type=float, default=_SLICE_DEPTH,
                   help='fixed depth (km) for the atlas slice; default is each '
                        "station's own median event depth")
    p.add_argument('--map-spacing', type=float, default=_MAP_SPACING,
                   help='atlas map node spacing in degrees (default %(default)s, '
                        '~2.2 km lat); lower it to resolve the L = 1 km panel')
    p.add_argument('--extent', default=None,
                   help='shared map frame for every atlas page, as '
                        'lon0,lon1,lat0,lat1. Needs the = form when lon0 is '
                        'negative: --extent=-2.25,3.5,41.75,43.75. '
                        'Default frames each page on its own station\'s events; '
                        'a shared frame makes pages comparable but is clipped '
                        "to each zone's LSGRID, off which no correction exists")
    p.add_argument('--stations', default=None,
                   help="comma-separated 'STA' or 'STA:PHASE' to draw only those "
                        'fields, e.g. FR.0041:P,RD.0038 - use with a small '
                        '--map-spacing for presentation-quality pages')
    p.add_argument('--product',
                   choices=('atlas', 'spread', 'displacement', 'all'),
                   default='all')
    p.add_argument('--extract-only', action='store_true',
                   help='fill the parse cache and stop')
    args = p.parse_args()

    _setup_logger()
    params = SsstCorrParams(
        run_name=args.run_name,
        zones=[int(z) for z in args.zones.split(',')] if args.zones else list(_ZONES),
        min_picks=args.min_picks,
        depth=args.depth,
        map_spacing=args.map_spacing,
        extent=(tuple(float(v) for v in args.extent.split(','))
                if args.extent else None),
        product=args.product,
    )
    if params.extent is not None and len(params.extent) != 4:
        p.error('--extent needs exactly four numbers: lon0,lon1,lat0,lat1')

    t0 = time.time()
    cache = load_all(params)
    logger.info('extraction ready in %.1f s', time.time() - t0)

    if args.extract_only:
        return

    os.makedirs(params.fig_dir, exist_ok=True)

    if params.product in ('atlas', 'all'):
        census = station_census(cache, params.zones)
        fields = best_zone_per_station(census, params.min_picks)
        logger.info('atlas: %d station/phase fields with >= %d usable arrivals',
                    len(fields), params.min_picks)
        wanted = [s for s in args.stations.split(',')] if args.stations else None
        if wanted:
            fields = filter_fields(fields, wanted)
            logger.info('atlas: restricted to %d field(s)', len(fields))
        name = (f'{params.run_name}_station_atlas'
                + ('_selection' if wanted else '') + '.pdf')
        t = time.time()
        make_atlas(params, cache, fields, os.path.join(params.fig_dir, name))
        logger.info('atlas done in %.1f s', time.time() - t)

    if params.product in ('spread', 'all'):
        t = time.time()
        make_spread_map(params, cache,
                        os.path.join(params.fig_dir,
                                     f'{params.run_name}_spread_map.pdf'))
        logger.info('spread map done in %.1f s', time.time() - t)

    if params.product in ('displacement', 'all'):
        t = time.time()
        make_displacement_map(params, cache,
                              os.path.join(params.fig_dir,
                                           f'{params.run_name}_displacement_map.pdf'))
        logger.info('displacement map done in %.1f s', time.time() - t)


if __name__ == '__main__':
    main()
