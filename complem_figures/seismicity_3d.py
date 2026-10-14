"""
seismicity_3d.py
================
Interactive 3-D map of the seismicity at every location stage of the
workflow, as one self-contained HTML page:

  pre-reloc  the merged source bulletin        obs/GLOBAL.obs
  NLL        first NonLinLoc relocation        RESULT/NLL_result.csv
  SSST       final SSST relocation             RESULT/SSST_result.csv

One stage is drawn at a time; the buttons switch stage without moving the
camera.  Each stage shows **all** of its own events (57 856 / 50 918 /
46 224), so part of what changes between stages is events lost along the way
(zones not covered, too few phases), not relocation — the page prints each
stage's count for that reason.

Events are 1.5 px points coloured by depth; national borders
and shorelines (the GMT cache of `basemap_lines.py`) and the stations NLLoc
located with (`stations/GTSRCE_*.txt`) are drawn at z = 0.  Clicking an event
dims the rest of the catalog, fills a side panel with its metadata at all
three stages, and draws its pre → NLL → SSST trajectory; "Auto-zoom" then
frames the event's solutions on half the panel, keeping the viewing direction.  Filters: year
range, minimum magnitude, and "usable only" (SSST — the same `pyr:usable`
test as the QuakeML export, imported from `NLL_run/export_quakeml.classify`).
A slider sets the vertical exaggeration.  The depth axis stops at
`--max-depth` (40 km); deeper events are hidden like a filter, since Plotly
does not clip 3-D points to the axis range.

Magnitudes are the merged-bulletin ones (`GLOBAL.obs`), joined on publicId
into the two later stages, which carry the same rematched value; the result
CSVs carry none; they feed the hover, the panel and the magnitude filter
(`Mag 0.00` is the OMP placeholder).  Relocated hypocentres are the PDF expectation by default
(`--solution maxlike` for the maximum-likelihood point, see
`event_maps.hypocentre_columns`).

Plotly.js is inlined, so the page works offline; the event columns are
embedded as base64 little-endian typed arrays.

Usage
-----
    python complem_figures/seismicity_3d.py
    python complem_figures/seismicity_3d.py --solution maxlike --exaggeration 3 \\
        --output complem_figures/seismicity_3d/seismicity_3d_maxlike.html
"""

import argparse
import base64
import datetime as dt
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
from complem_figures import basemap_lines                 # noqa: E402
from complem_figures.event_maps import hypocentre_columns  # noqa: E402
from NLL_run.export_quakeml import classify                # noqa: E402

_C68_NOMINAL = 0.6827   # as export_quakeml._NOMINAL_COVERAGE

_DEFAULT_OBS    = os.path.join(_PROJECT_ROOT, 'obs', 'GLOBAL.obs')
_DEFAULT_NLL    = os.path.join(_PROJECT_ROOT, 'RESULT', 'NLL_result.csv')
_DEFAULT_SSST   = os.path.join(_PROJECT_ROOT, 'RESULT', 'SSST_result.csv')
_DEFAULT_OUTPUT = os.path.join(_MODULE_DIR, 'seismicity_3d', 'seismicity_3d.html')
_GTSRCE_GLOB    = os.path.join(_PROJECT_ROOT, 'stations', 'GTSRCE_*.txt')


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def _num(token):
    return np.nan if token == 'None' else float(token)


def read_obs_events(path):
    """
    One row per event of a .obs bulletin: pid, t, lon, lat, z, mag, magtype,
    nphs, erh, erz, gap, rms — the event line, attached to its PUBLIC_ID.
    """
    rows, event = [], None
    with open(path) as f:
        for line in f:
            if line.startswith('# '):
                p = line.split()
                # Seconds can read 60.0 in source bulletins, hence the timedelta
                t = (dt.datetime(int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]),
                                 tzinfo=dt.timezone.utc)
                     + dt.timedelta(seconds=float(p[6]))).timestamp()
                event = (t, float(p[8]), float(p[7]), float(p[9]), _num(p[10]), p[11],
                         _num(p[13]), _num(p[14]), _num(p[15]), _num(p[16]), _num(p[17]))
            elif line.startswith('PUBLIC_ID') and event is not None:
                rows.append((int(line.split()[1].rsplit('_', 1)[1]),) + event)
                event = None
    return pd.DataFrame(rows, columns=['pid', 't', 'lon', 'lat', 'z', 'mag', 'magtype',
                                       'nphs', 'erh', 'erz', 'gap', 'rms'])


def read_result_csv(path, solution):
    """One row per event of a RESULT/*.csv, at the chosen hypocentre."""
    df = pd.read_csv(path, skipinitialspace=True)
    lat_col, lon_col, z_col = hypocentre_columns(df, solution)
    out = pd.DataFrame({
        'pid':  df['publicId'].str.rsplit('_', n=1).str[1].astype(int),
        't':    pd.to_datetime(df['date-time']).dt.tz_localize('UTC').astype('int64') / 1e9,
        'lon':  df[lon_col], 'lat': df[lat_col], 'z': df[z_col],
        'nphs': df['Nphs'], 'erh': df['true_erh'], 'erz': df['true_erz'],
        'gap':  df['Gap'], 'rms': df['RMS'],
    })
    if 'C68' in df.columns:
        out['psi']    = df['Psi']
        out['c68z']   = (df['C68'] - _C68_NOMINAL) / df['C68_sigma_n']
        out['usable'] = [classify(row)[0] for _, row in df.iterrows()]
    return out


def read_stations(pattern=_GTSRCE_GLOB):
    """(code, lon, lat) of every station in the GTSRCE files, one row per code."""
    rows = []
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            for line in f:
                p = line.split()
                if len(p) >= 5 and p[0] == 'GTSRCE':
                    rows.append((p[1], float(p[4]), float(p[3])))
    return pd.DataFrame(rows, columns=['code', 'lon', 'lat']).drop_duplicates('code')


def _basemap_polyline(kind, box):
    """A basemap cache as one NaN-separated polyline, cropped to box."""
    segments = basemap_lines.load(kind)
    if segments is None:
        print(f'Warning: no {kind} cache — run basemap_lines.py --dump in pygmt_env')
        return [], []
    lon0, lon1, lat0, lat1 = box
    xs, ys = [], []
    for seg in segments:
        inside = ((seg[:, 0] >= lon0) & (seg[:, 0] <= lon1)
                  & (seg[:, 1] >= lat0) & (seg[:, 1] <= lat1))
        x = np.where(inside, seg[:, 0], np.nan)
        y = np.where(inside, seg[:, 1], np.nan)
        xs += [*x.round(4), None]
        ys += [*y.round(4), None]
    return [None if v != v else v for v in xs], [None if v != v else v for v in ys]


# ---------------------------------------------------------------------------
# Payload
# ---------------------------------------------------------------------------

_DTYPES = {'pid': '<i4', 't': '<f8', 'usable': '<u1', 'magtype': '<u1'}


def _b64(values, dtype):
    return base64.b64encode(np.ascontiguousarray(values, dtype=dtype).tobytes()).decode('ascii')


def _encode_stage(df):
    return {'n': len(df),
            'cols': {c: [_DTYPES.get(c, '<f4')[1:], _b64(df[c].to_numpy(), _DTYPES.get(c, '<f4'))]
                     for c in df.columns}}


def build_payload(obs_path, nll_path, ssst_path, solution, exaggeration, max_depth):
    pre = read_obs_events(obs_path)
    print(f'pre-reloc  @ {obs_path}: {len(pre)} events')
    nll = read_result_csv(nll_path, solution)
    print(f'NLL        @ {nll_path}: {len(nll)} events')
    ssst = read_result_csv(ssst_path, solution)
    print(f'SSST       @ {ssst_path}: {len(ssst)} events'
          + (f', {int(ssst["usable"].sum())} usable' if 'usable' in ssst else ''))

    # Magnitudes live in the merged bulletin only
    magtypes = sorted(pre['magtype'].unique())
    pre['magtype'] = pre['magtype'].map({m: i for i, m in enumerate(magtypes)})
    mags = pre[['pid', 'mag', 'magtype']]
    nll  = nll.merge(mags, on='pid', how='left')
    ssst = ssst.merge(mags, on='pid', how='left')

    stages = {'pre': pre, 'nll': nll, 'ssst': ssst}
    allz   = pd.concat([s[['lon', 'lat', 'z']] for s in stages.values()])
    box    = [float(np.floor(allz.lon.min() * 10) / 10), float(np.ceil(allz.lon.max() * 10) / 10),
              float(np.floor(allz.lat.min() * 10) / 10), float(np.ceil(allz.lat.max() * 10) / 10)]
    zrange = [float(min(-3.0, np.floor(allz.z.min()))), float(max_depth)]
    print(f'Hidden below {max_depth:g} km: '
          + ', '.join(f'{k} {int((v.z > max_depth).sum())}' for k, v in stages.items()))

    stations = read_stations()
    stations = stations[(stations.lon.between(box[0], box[1])) & (stations.lat.between(box[2], box[3]))]
    bx, by = _basemap_polyline('borders', box)
    sx, sy = _basemap_polyline('shorelines', box)

    years = pd.to_datetime(pd.concat([s.t for s in stages.values()]), unit='s').dt.year
    return {
        'stages':   {k: _encode_stage(v) for k, v in stages.items()},
        'magtypes': magtypes,
        'solution': solution,
        'box':      box,
        'zrange':   zrange,
        'years':    [int(years.min()), int(years.max())],
        'exag':     exaggeration,
        'borders':  [bx, by],
        'coast':    [sx, sy],
        'stations': {'code': stations.code.tolist(),
                     'lon': stations.lon.round(5).tolist(), 'lat': stations.lat.round(5).tolist()},
    }


# ---------------------------------------------------------------------------
# Page
# ---------------------------------------------------------------------------

_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pyrenees seismicity 3-D</title>
<script>__PLOTLYJS__</script>
<style>
  :root { --ink:#1f2328; --muted:#636c76; --line:#d0d7de; --bg:#ffffff; --panel:#f6f8fa; --accent:#0969da; }
  * { box-sizing: border-box; }
  html, body { margin:0; height:100%; background:var(--bg); color:var(--ink);
               font: 13px/1.4 -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; }
  #app { display:flex; height:100vh; }
  #plot { flex:1; min-width:0; }
  #side { width:340px; border-left:1px solid var(--line); background:var(--panel);
          overflow-y:auto; padding:12px 14px; }
  h2 { font-size:13px; margin:14px 0 6px; text-transform:uppercase; letter-spacing:.04em; color:var(--muted); }
  h2:first-child { margin-top:0; }
  .seg { display:flex; gap:0; }
  .seg button { flex:1; padding:6px 0; border:1px solid var(--line); background:#fff; cursor:pointer; font:inherit; }
  .seg button + button { border-left:none; }
  .seg button.on { background:var(--accent); color:#fff; border-color:var(--accent); }
  .row { display:flex; align-items:center; gap:6px; margin:5px 0; }
  .row label { flex:1; }
  input[type=number] { width:70px; font:inherit; padding:2px 4px; }
  input[type=range] { flex:1; }
  #count { color:var(--muted); margin-top:6px; }
  table { border-collapse:collapse; width:100%; font-variant-numeric:tabular-nums; }
  th, td { text-align:right; padding:2px 4px; border-bottom:1px solid var(--line); }
  th:first-child, td:first-child { text-align:left; color:var(--muted); font-weight:normal; }
  .muted { color:var(--muted); }
  button.clear { margin-top:8px; padding:4px 10px; font:inherit; cursor:pointer; }
  @media (max-width: 760px) { #app { flex-direction:column; } #plot { height:60vh; flex:none; }
                              #side { width:auto; border-left:none; border-top:1px solid var(--line); } }
</style>
</head>
<body>
<div id="app">
  <div id="plot"></div>
  <div id="side">
    <h2>Stage</h2>
    <div class="seg" id="stages">
      <button data-s="pre">Pre-reloc</button><button data-s="nll">NLL</button><button data-s="ssst">SSST</button>
    </div>
    <div id="count"></div>

    <h2>Filters</h2>
    <div class="row"><label>Years</label>
      <input type="number" id="ymin"> – <input type="number" id="ymax"></div>
    <div class="row"><label>Min. magnitude</label><input type="number" id="mmin" step="0.1" placeholder="none"></div>
    <div class="row"><label for="usable">Usable only (SSST)</label><input type="checkbox" id="usable"></div>

    <h2>Display</h2>
    <div class="row"><label>Vertical exaggeration</label><span id="exagv"></span></div>
    <div class="row"><input type="range" id="exag" min="1" max="20" step="0.5"></div>
    <div class="row"><label>Other events when one is selected</label><span id="fadev"></span></div>
    <div class="row"><input type="range" id="fade" min="0" max="0.85" step="0.01" value="0.08"></div>

    <h2>Selected event</h2>
    <div id="info" class="muted">Click an event. Esc or “Clear” to deselect.</div>
  </div>
</div>
<script>
const D = __PAYLOAD__;
const STAGES = ['pre', 'nll', 'ssst'];
const LABEL  = {pre: 'Pre-reloc', nll: 'NLL', ssst: 'SSST'};
const SCOL   = {pre: '#6e7781', nll: '#0969da', ssst: '#cf222e'};
const TYPED  = {f4: Float32Array, f8: Float64Array, i4: Int32Array, u1: Uint8Array};

function decode([kind, b64]) {
  const bin = atob(b64), bytes = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
  return new TYPED[kind](bytes.buffer);
}
const S = {};
for (const s of STAGES) {
  const st = {n: D.stages[s].n, c: {}, idx: new Map()};
  for (const [k, v] of Object.entries(D.stages[s].cols)) st.c[k] = decode(v);
  st.year = new Int16Array(st.n);
  for (let i = 0; i < st.n; i++) {
    st.idx.set(st.c.pid[i], i);
    st.year[i] = new Date(st.c.t[i] * 1000).getUTCFullYear();
  }
  S[s] = st;
}

const state = {stage: 'ssst', ymin: D.years[0], ymax: D.years[1], mmin: null,
               usable: false, sel: null, exag: D.exag, fade: 0.08, vis: []};

const fmt  = (v, d = 2) => (v === undefined || Number.isNaN(v)) ? '—' : v.toFixed(d);
const pidS = p => 'PYRENEES_' + String(p).padStart(6, '0');
const date = t => new Date(t * 1000).toISOString().replace('T', ' ').slice(0, 21);
const magS = (st, i) => Number.isNaN(st.c.mag[i]) ? '—'
                        : st.c.mag[i].toFixed(2) + ' ' + D.magtypes[st.c.magtype[i]];

// Aspect ratio in km, so that the exaggeration slider is a true vertical exaggeration
const lat0 = (D.box[2] + D.box[3]) / 2;
const kmX  = (D.box[1] - D.box[0]) * 111.195 * Math.cos(lat0 * Math.PI / 180);
const kmY  = (D.box[3] - D.box[2]) * 111.195;
const kmZ  = D.zrange[1] - D.zrange[0];
function aspect() {
  const m = Math.max(kmX, kmY);
  return {x: kmX / m, y: kmY / m, z: kmZ * state.exag / m};
}

function visible() {
  const st = S[state.stage], out = [];
  const useU = state.usable && state.stage === 'ssst' && st.c.usable;
  for (let i = 0; i < st.n; i++) {
    const y = st.year[i];
    if (y < state.ymin || y > state.ymax) continue;
    if (state.mmin !== null && !(st.c.mag[i] >= state.mmin)) continue;
    if (useU && !st.c.usable[i]) continue;
    if (st.c.z[i] > D.zrange[1]) continue;
    out.push(i);
  }
  return out;
}

function eventTrace() {
  const st = S[state.stage], v = state.vis, n = v.length;
  const x = new Float32Array(n), y = new Float32Array(n), z = new Float32Array(n),
        text = new Array(n);
  for (let k = 0; k < n; k++) {
    const i = v[k];
    x[k] = st.c.lon[i]; y[k] = st.c.lat[i]; z[k] = st.c.z[i];
    text[k] = pidS(st.c.pid[i]) + '<br>' + date(st.c.t[i]) + '<br>M ' + magS(st, i)
              + '<br>depth ' + st.c.z[i].toFixed(2) + ' km';
  }
  return {type: 'scatter3d', mode: 'markers', name: LABEL[state.stage] + ' events',
          x, y, z, text, hovertemplate: '%{text}<extra></extra>', showlegend: false,
          marker: {size: 1.5, color: z, colorscale: 'Viridis', reversescale: true,
                   cmin: D.zrange[0], cmax: Math.min(D.zrange[1], 25),
                   line: {width: 0},
                   colorbar: {title: {text: 'Depth (km)'}, thickness: 12, len: 0.6}}};
}

function surfaceTraces() {
  const zeros = a => a.map(v => v === null ? null : 0);
  return [
    {type: 'scatter3d', mode: 'lines', name: 'Borders', hoverinfo: 'skip',
     x: D.borders[0], y: D.borders[1], z: zeros(D.borders[0]), line: {color: '#24292f', width: 3}},
    {type: 'scatter3d', mode: 'lines', name: 'Shorelines', hoverinfo: 'skip',
     x: D.coast[0], y: D.coast[1], z: zeros(D.coast[0]), line: {color: '#57606a', width: 2}},
    {type: 'scatter3d', mode: 'markers', name: 'Stations', visible: 'legendonly',
     x: D.stations.lon, y: D.stations.lat, z: D.stations.lon.map(() => 0), text: D.stations.code,
     hovertemplate: '%{text}<extra>station</extra>',
     marker: {symbol: 'diamond', size: 3, color: '#bc4c00'}},
  ];
}

function selectionTraces() {
  if (state.sel === null) return [];
  const xs = [], ys = [], zs = [], cs = [], ts = [];
  for (const s of STAGES) {
    const i = S[s].idx.get(state.sel);
    if (i === undefined || S[s].c.z[i] > D.zrange[1]) continue;
    xs.push(S[s].c.lon[i]); ys.push(S[s].c.lat[i]); zs.push(S[s].c.z[i]);
    cs.push(SCOL[s]); ts.push(LABEL[s]);
  }
  const out = [{type: 'scatter3d', mode: 'lines+markers', name: 'Trajectory', showlegend: false,
                x: xs, y: ys, z: zs, text: ts, hovertemplate: '%{text}<extra></extra>',
                line: {color: '#24292f', width: 4}, marker: {size: 5, color: cs}}];
  const i = S[state.stage].idx.get(state.sel);
  if (i !== undefined && S[state.stage].c.z[i] <= D.zrange[1]) {
    const st = S[state.stage];
    out.push({type: 'scatter3d', mode: 'markers', name: 'Selected', showlegend: false,
              x: [st.c.lon[i]], y: [st.c.lat[i]], z: [st.c.z[i]], hoverinfo: 'skip',
              marker: {size: 9, color: SCOL[state.stage],
                       opacity: 1, line: {color: '#000', width: 2}}});
  }
  return out;
}

// Axis ranges, in Plotly's [r0, r1] order (depth reversed).
const RANGES = [[D.box[0], D.box[1]], [D.box[2], D.box[3]], [D.zrange[1], D.zrange[0]]];

// The camera is passed explicitly on every Plotly.react, read back from the
// live scene first, so that neither a redraw nor the auto-zoom loses the view.
let camera = {eye: {x: 0.3, y: -1.6, z: 0.9}, center: {x: 0, y: 0, z: 0}, up: {x: 0, y: 0, z: 1}};
function liveCamera() {
  const sc = plot._fullLayout && plot._fullLayout.scene && plot._fullLayout.scene._scene;
  if (sc) { const c = sc.getCamera(); camera = {eye: c.eye, center: c.center, up: c.up}; }
  return camera;
}

function layout() {
  return {
    uirevision: 'keep', margin: {l: 0, r: 0, t: 0, b: 0},
    legend: {x: 0.01, y: 0.99, bgcolor: 'rgba(255,255,255,0.7)'},
    scene: {
      uirevision: 'keep', aspectmode: 'manual', aspectratio: aspect(),
      xaxis: {title: {text: 'Longitude (°)'}, range: RANGES[0]},
      yaxis: {title: {text: 'Latitude (°)'}, range: RANGES[1]},
      zaxis: {title: {text: 'Depth (km)'}, range: RANGES[2]},
      camera: camera,
    },
  };
}

const plot = document.getElementById('plot');
// The event trace is rebuilt only when the stage or a filter changes.  A
// selection reuses the same arrays, so Plotly.react sees the data unchanged
// and only restyles the opacity and the two small selection traces.
let events = null;
function draw(rebuild = true) {
  if (rebuild) { state.vis = visible(); events = eventTrace(); }
  const ev = {...events, marker: {...events.marker, opacity: state.sel === null ? 0.85 : state.fade}};
  liveCamera();
  Plotly.react(plot, [ev, ...surfaceTraces(), ...selectionTraces()], layout(),
               {responsive: true, displaylogo: false});
  const n = S[state.stage].n;
  document.getElementById('count').textContent =
    `${state.vis.length.toLocaleString('en')} of ${n.toLocaleString('en')} ${LABEL[state.stage]} events shown`;
  document.querySelectorAll('#stages button').forEach(b => b.classList.toggle('on', b.dataset.s === state.stage));
  document.getElementById('usable').disabled = state.stage !== 'ssst';
  renderInfo();
}

function dist(a, i, b, j) {
  const kx = 111.195 * Math.cos(a.c.lat[i] * Math.PI / 180);
  const dx = (b.c.lon[j] - a.c.lon[i]) * kx, dy = (b.c.lat[j] - a.c.lat[i]) * 111.195;
  const dz = b.c.z[j] - a.c.z[i], h = Math.hypot(dx, dy);
  return `${h.toFixed(2)} km horiz., ${dz >= 0 ? '+' : ''}${dz.toFixed(2)} km depth, ${Math.hypot(h, dz).toFixed(2)} km 3-D`;
}

function renderInfo() {
  const el = document.getElementById('info');
  if (state.sel === null) {
    el.className = 'muted'; el.innerHTML = 'Click an event. Esc or “Clear” to deselect.'; return;
  }
  el.className = '';
  const ix = Object.fromEntries(STAGES.map(s => [s, S[s].idx.get(state.sel)]));
  const any = STAGES.find(s => ix[s] !== undefined);
  const cell = (s, f) => ix[s] === undefined ? '—' : f(S[s], ix[s]);
  const col  = (k, d) => (st, i) => st.c[k] ? fmt(st.c[k][i], d) : '—';
  const rows = [
    ['Origin time', (st, i) => date(st.c.t[i]).slice(11)],
    ['Latitude (°)', col('lat', 4)], ['Longitude (°)', col('lon', 4)], ['Depth (km)', col('z', 2)],
    ['RMS (s)', col('rms', 3)], ['Gap (°)', col('gap', 0)], ['Phases', col('nphs', 0)],
    ['ERH (km)', col('erh', 2)], ['ERZ (km)', col('erz', 2)],
    ['Ψ', col('psi', 2)], ['C68 z', col('c68z', 2)],
    ['Usable', (st, i) => st.c.usable ? (st.c.usable[i] ? 'yes' : 'no') : '—'],
  ];
  let h = `<b>${pidS(state.sel)}</b><br>${date(S[any].c.t[ix[any]]).slice(0, 10)} · M ${magS(S[any], ix[any])}`;
  if (ix[state.stage] === undefined) h += `<br><span class="muted">absent from ${LABEL[state.stage]}</span>`;
  h += '<table><tr><th></th>' + STAGES.map(s => `<th style="color:${SCOL[s]}">${LABEL[s]}</th>`).join('') + '</tr>';
  for (const [name, f] of rows) h += `<tr><td>${name}</td>` + STAGES.map(s => `<td>${cell(s, f)}</td>`).join('') + '</tr>';
  h += '</table>';
  const miss = {nll: 'not relocated by NLL', ssst: 'dropped at SSST'};
  for (const s of ['nll', 'ssst']) if (ix[s] === undefined) h += `<div class="muted">${miss[s]}</div>`;
  h += '<h2>Displacement</h2>';
  for (const [a, b] of [['pre', 'nll'], ['nll', 'ssst'], ['pre', 'ssst']]) {
    if (ix[a] === undefined || ix[b] === undefined) continue;
    h += `<div>${LABEL[a]} → ${LABEL[b]}: ${dist(S[a], ix[a], S[b], ix[b])}</div>`;
  }
  h += '<button class="clear" id="zoom">Auto-zoom</button> <button class="clear" id="clear">Clear</button>';
  el.innerHTML = h;
  document.getElementById('clear').onclick = () => { state.sel = null; draw(false); };
  document.getElementById('zoom').onclick = autoZoom;
}

// Auto-zoom.  Plotly places a value v of an axis with range [r0, r1] at
// aspect·(v − (r0 + r1)/2)/(r1 − r0) in scene units — the box spans ±aspect/2
// — and camera eye/center are in those units.  Its perspective has a 45°
// vertical field of view, so from distance L the view is 2·L·tan(22.5°) high
// (× width/height wide).  The camera is aimed at the centroid of the event's
// solutions and keeps its current viewing direction; its distance is set so
// that their projected spread covers half the panel, on whichever screen axis
// it is wider.  Never closer than a 1 km sphere would allow, so that
// solutions that coincide or line up with the view do not zoom without limit.
function toScene(v) {
  const a = aspect(), A = [a.x, a.y, a.z];
  return v.map((u, k) => A[k] * (u - (RANGES[k][0] + RANGES[k][1]) / 2) / (RANGES[k][1] - RANGES[k][0]));
}
function autoZoom() {
  const pts = [];
  for (const s of STAGES) {
    const i = S[s].idx.get(state.sel);
    if (i === undefined || S[s].c.z[i] > D.zrange[1]) continue;
    pts.push(toScene([S[s].c.lon[i], S[s].c.lat[i], S[s].c.z[i]]));
  }
  if (!pts.length) return;
  const c = [0, 1, 2].map(k => pts.reduce((acc, p) => acc + p[k], 0) / pts.length);
  const cam = liveCamera(), T = Math.tan(Math.PI / 8);
  const box = plot.getBoundingClientRect(), wh = box.width / box.height;
  const n = Math.hypot(cam.eye.x - cam.center.x, cam.eye.y - cam.center.y, cam.eye.z - cam.center.z) || 1;
  const f = [(cam.center.x - cam.eye.x) / n, (cam.center.y - cam.eye.y) / n, (cam.center.z - cam.eye.z) / n];
  const u = [cam.up.x, cam.up.y, cam.up.z];
  const right = [f[1] * u[2] - f[2] * u[1], f[2] * u[0] - f[0] * u[2], f[0] * u[1] - f[1] * u[0]];
  const rn = Math.hypot(...right);
  const rx = right.map(v => v / rn);
  const up = [rx[1] * f[2] - rx[2] * f[1], rx[2] * f[0] - rx[0] * f[2], rx[0] * f[1] - rx[1] * f[0]];
  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  // Half-extent of the projected solutions, as a fraction of the panel, from distance L
  const spread = L => {
    const s = pts.map(p => {
      const q = [p[0] - c[0] + f[0] * L, p[1] - c[1] + f[1] * L, p[2] - c[2] + f[2] * L];
      const zc = dot(q, f);
      return [dot(q, rx) / (zc * T * wh), dot(q, up) / (zc * T)];
    });
    const ext = k => (Math.max(...s.map(a => a[k])) - Math.min(...s.map(a => a[k]))) / 2;
    return Math.max(ext(0), ext(1));
  };
  const rMin = aspect().x / kmX;
  const Lmin = 2 * rMin / (Math.min(1, wh) * T);
  let L = Lmin;
  for (let it = 0; it < 4; it++) { const sp = spread(L); if (sp > 0) L *= sp / 0.5; }
  L = Math.max(L, Lmin);
  camera = {up: cam.up, center: {x: c[0], y: c[1], z: c[2]},
            eye: {x: c[0] - f[0] * L, y: c[1] - f[1] * L, z: c[2] - f[2] * L}};
  Plotly.relayout(plot, {'scene.camera': camera});
}

// Controls
document.querySelectorAll('#stages button').forEach(b => b.onclick = () => { state.stage = b.dataset.s; draw(); });
const ymin = document.getElementById('ymin'), ymax = document.getElementById('ymax'),
      mmin = document.getElementById('mmin'), usable = document.getElementById('usable'),
      exag = document.getElementById('exag'), exagv = document.getElementById('exagv');
for (const e of [ymin, ymax]) { e.min = D.years[0]; e.max = D.years[1]; }
ymin.value = D.years[0]; ymax.value = D.years[1];
exag.value = D.exag; exagv.textContent = '×' + D.exag;
ymin.onchange = () => { state.ymin = +ymin.value || D.years[0]; draw(); };
ymax.onchange = () => { state.ymax = +ymax.value || D.years[1]; draw(); };
mmin.onchange = () => { state.mmin = mmin.value === '' ? null : +mmin.value; draw(); };
usable.onchange = () => { state.usable = usable.checked; draw(); };
const fade = document.getElementById('fade'), fadev = document.getElementById('fadev');
fadev.textContent = 'opacity ' + state.fade.toFixed(2);
fade.oninput = () => {
  state.fade = +fade.value; fadev.textContent = 'opacity ' + state.fade.toFixed(2);
  if (state.sel !== null) Plotly.restyle(plot, {'marker.opacity': state.fade}, [0]);
};
exag.oninput = () => {
  state.exag = +exag.value; exagv.textContent = '×' + state.exag;
  Plotly.relayout(plot, {'scene.aspectratio': aspect()});
};
document.addEventListener('keydown', e => { if (e.key === 'Escape' && state.sel !== null) { state.sel = null; draw(false); } });

draw();

// Selection.  In a 3-D scene Plotly emits plotly_click from its render loop on
// *every frame* while a button is held over a point — including the whole of
// a drag that rotates the view.  Redrawing from that handler re-triggers the
// render and loops.  So the handler only records the point under the cursor;
// the selection is committed on button release, if the pointer did not move
// (a click, not a rotation), outside Plotly's render loop.
let press = null, candidate = null;
plot.addEventListener('pointerdown', e => { press = [e.clientX, e.clientY]; candidate = null; }, true);
plot.addEventListener('pointerup', e => {
  const still = press && Math.hypot(e.clientX - press[0], e.clientY - press[1]) < 5;
  const pid = candidate;
  press = null; candidate = null;
  if (still && pid !== null && pid !== state.sel) {
    state.sel = pid;
    setTimeout(() => draw(false), 0);
  }
}, true);
plot.on('plotly_click', ev => {
  const p = ev.points[0];
  if (press === null || p.curveNumber !== 0) return;   // only the event trace selects
  candidate = S[state.stage].c.pid[state.vis[p.pointNumber]];
});
</script>
</body>
</html>
"""


def write_page(payload, output):
    import plotly.offline
    html = (_TEMPLATE
            .replace('__PLOTLYJS__', plotly.offline.get_plotlyjs())
            .replace('__PAYLOAD__', json.dumps(payload, separators=(',', ':'))))
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f'Figure saved @ {output} ({os.path.getsize(output) / 1e6:.1f} MB)')


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--obs',  default=_DEFAULT_OBS,  help='pre-relocation bulletin (.obs)')
    parser.add_argument('--nll',  default=_DEFAULT_NLL,  help='NLL result CSV')
    parser.add_argument('--ssst', default=_DEFAULT_SSST, help='SSST result CSV')
    parser.add_argument('--solution', choices=['expect', 'maxlike'], default='expect',
                        help='relocated hypocentre: PDF expectation (default) or maximum likelihood')
    parser.add_argument('--exaggeration', type=float, default=5.0,
                        help='starting value of the vertical-exaggeration slider (default 5)')
    parser.add_argument('--max-depth', type=float, default=40.0,
                        help='bottom of the depth axis in km; deeper events are hidden (default 40)')
    parser.add_argument('--output', default=_DEFAULT_OUTPUT)
    args = parser.parse_args()

    payload = build_payload(args.obs, args.nll, args.ssst, args.solution, args.exaggeration,
                            args.max_depth)
    write_page(payload, args.output)


if __name__ == '__main__':
    main()
