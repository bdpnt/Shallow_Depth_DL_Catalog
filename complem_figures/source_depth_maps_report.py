"""
source_depth_maps_report.py
===========================
Report figure: two source bulletins side by side (IGN left, LDG right by
default), every event coloured by depth and sized by magnitude, so the depth
distributions of the two catalogues can be compared on the same map.

Style follows the maps drawn by `global_obs/filter_events_by_aoi.py` (same
basemap, `viridis` reversed 0-15 km, `0.03 cm x Mag` symbols), with:

* the figure drawn at its printed size (`--width-cm`, default 16 cm =
  `\\linewidth` of an A4 page with 2.5 cm margins) and 8-9 pt fonts, so LaTeX
  includes it at scale 1 and the text stays legible. Symbol sizes are
  rescaled by panel width so the maps look as dense as the 6-inch originals;
* both panels in one figure, one shared colour bar, the right panel's latitude
  labels dropped, and a bold "a" / "b" in each map's top-left corner;
* read-only: bulletins are only read (a plain `.obs`, or a zip member given as
  `archive.zip::member`). Magnitudes still carrying a source type with a model
  in `mag_model/` are converted to ML LDG in memory, as
  `global_obs/apply_magnitude_models.py` does; already harmonised bulletins
  pass through unchanged;
* both a 300 dpi PNG and a vector PDF.

Usage
-----
    conda run -n pygmt_env python complem_figures/source_depth_maps_report.py \\
        --left   obs/IGN_20-25.obs \\
        --right  obs/LDG_20-25.obs \\
        --output complem_figures/source_depth_maps/IGN_LDG_report.png
"""

import argparse
import io
import os
import pickle
import zipfile
from dataclasses import dataclass
from typing import Optional

import pandas as pd
import pygmt as pg

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

# Map extent and symbol law of global_obs/filter_events_by_aoi.py
DEFAULT_REGION   = [-2.25, 3.5, 42.0, 44.0]
_ORIG_MAP_WIDTH  = 6 * 2.54      # cm — projection 'M6i' of the original
_ORIG_SIZE_PER_M = 0.03          # cm per magnitude unit


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SourceDepthMapsParams:
    left:          str
    right:         str
    figSave:       str
    left_tag:      Optional[str] = 'a'
    right_tag:     Optional[str] = 'b'
    map_region:    Optional[list] = None
    width_cm:      float = 16.0
    gap_cm:        float = 0.35
    size_factor:   float = 1.0
    max_depth:     float = 15.0
    convert_mags:  bool  = True
    font_annot:    str   = '8p'
    font_label:    str   = '9p'
    font_tag:      str   = '11p'


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _open_text(path):
    """Open a plain file, or a member of a zip given as `archive.zip::member`."""
    if '::' in path:
        archive, member = path.split('::', 1)
        return io.TextIOWrapper(zipfile.ZipFile(archive).open(member), encoding='utf-8'), member
    return open(path, 'r', encoding='utf-8'), os.path.basename(path)


def _load_model(stem):
    """Load mag_model/<stem>.joblib; the files are plain pickled dicts of floats."""
    path = os.path.join(_PROJECT_ROOT, 'mag_model', f'{stem}.joblib')
    if not os.path.exists(path):
        return None
    try:
        import joblib
        return joblib.load(path)
    except ImportError:
        with open(path, 'rb') as f:
            return pickle.load(f)


def _convert_magnitudes(df, author):
    """
    In-memory copy of apply_magnitude_models.py: piecewise-linear conversion to
    ML LDG, split at Mag = 2, looked up as '<MagType> <FileAuthor>'. Types with
    no model are left as they are, as in the pipeline.
    """
    df = df.copy()
    for mag_type in df.MagType.unique():
        models = _load_model(f'{mag_type} {author}')
        if models is None:
            continue
        sel  = df.MagType == mag_type
        ge   = next(models[k] for k in models if '≥' in k)
        lt   = next(models[k] for k in models if '<' in k)
        m_ge, m_lt = sel & (df.Magnitude >= 2), sel & (df.Magnitude < 2)
        df.loc[m_ge, 'Magnitude'] = ge['slope'] * df.loc[m_ge, 'Magnitude'] + ge['intercept']
        df.loc[m_lt, 'Magnitude'] = lt['slope'] * df.loc[m_lt, 'Magnitude'] + lt['intercept']
        df.loc[sel, 'MagType'] = 'ML'
        print(f'  converted {sel.sum()} "{mag_type} {author}" magnitudes to ML LDG')
    return df


def read_bulletin(path, convert_mags=True):
    """Event lines of a .obs → DataFrame(Latitude, Longitude, Depth, Magnitude, MagType)."""
    handle, name = _open_text(path)
    with handle:
        events = [line[2:].split() for line in handle if line.startswith('# ')]
    raw = pd.DataFrame(events)
    df = pd.DataFrame({
        'Latitude':  pd.to_numeric(raw[6], errors='coerce'),
        'Longitude': pd.to_numeric(raw[7], errors='coerce'),
        'Depth':     pd.to_numeric(raw[8], errors='coerce'),
        'Magnitude': pd.to_numeric(raw[9], errors='coerce'),
        'MagType':   raw[10],
    }).dropna(subset=['Latitude', 'Longitude', 'Depth'])
    print(f'Catalog read @ {path} — {len(df)} events')
    if convert_mags:
        df = _convert_magnitudes(df, name.split('_')[0])
    return df


def _sizes(mag, scale):
    return (_ORIG_SIZE_PER_M * mag.clip(lower=0) * scale).to_numpy()


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    """
    Draw the two depth maps side by side at their printed width.

    Returns
    -------
    dict with keys: output, outputs, n_events
    """
    region  = parameters.map_region if parameters.map_region else list(DEFAULT_REGION)
    panel_w = (parameters.width_cm - 1.0 - parameters.gap_cm) / 2   # 1 cm for the W labels
    scale   = panel_w / _ORIG_MAP_WIDTH * parameters.size_factor
    proj    = f'M{panel_w:.4f}c'

    frames  = [read_bulletin(parameters.left,  parameters.convert_mags),
               read_bulletin(parameters.right, parameters.convert_mags)]
    tags    = [parameters.left_tag, parameters.right_tag]

    pg.config(
        FONT_ANNOT_PRIMARY=f'{parameters.font_annot},Helvetica,black',
        FONT_LABEL=f'{parameters.font_label},Helvetica,black',
        MAP_FRAME_TYPE='fancy+',
        MAP_FRAME_WIDTH='2.5p',
        MAP_TICK_LENGTH_PRIMARY='3p',
        MAP_ANNOT_OFFSET_PRIMARY='2p',
        MAP_LABEL_OFFSET='3p',
    )

    cpt = os.path.join(_PROJECT_ROOT, '.source_depth_maps_report.cpt')
    pg.makecpt(cmap='viridis', series=[0, parameters.max_depth, 1], reverse=True, output=cpt)

    fig = pg.Figure()
    fig.shift_origin(xshift='1c', yshift='2c')

    for i, (df, tag) in enumerate(zip(frames, tags)):
        if i == 1:
            fig.shift_origin(xshift=f'{panel_w + parameters.gap_cm:.4f}c')
        fig.basemap(region=region, projection=proj,
                    frame=['WSne' if i == 0 else 'wSne', 'xa1f0.5', 'ya1f0.5'])
        fig.coast(water='skyblue', land='#777777', resolution='i',
                  area_thresh='0/0/1', borders='1/0.5p,black')
        fig.plot(x=df.Longitude, y=df.Latitude, style='cc',
                 size=_sizes(df.Magnitude, scale),
                 fill=df.Depth, cmap=cpt, transparency=30)
        if tag:
            # inside the top-left corner, over the Bay of Biscay
            fig.text(text=tag, position='TL', justify='TL', offset='0.2c/-0.15c',
                     font=f'{parameters.font_tag},Helvetica-Bold,black')

    # one colour bar under both panels; the triangle shows the colour of events deeper than max_depth
    fig.shift_origin(xshift=f'{-(panel_w + parameters.gap_cm):.4f}c')
    centre = (2 * panel_w + parameters.gap_cm) / 2
    fig.colorbar(cmap=cpt,
                 position=f'x{centre:.4f}c/-0.85c+w{0.6 * parameters.width_cm:.3f}c/0.22c+h+jTC+ef',
                 frame=['xa5f1+lDepth (km)'])

    base, _ = os.path.splitext(parameters.figSave)
    outputs = []
    for path in (f'{base}.png', f'{base}.pdf'):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path, dpi=300, crop=True)
        outputs.append(path)
        print(f'Figure saved @ {path}')

    try:
        os.remove(cpt)
    except OSError:
        pass

    return {'output': outputs[0], 'outputs': outputs, 'n_events': [len(d) for d in frames]}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Report-ready two-panel depth map of two source bulletins (IGN / LDG by default).'
    )
    parser.add_argument('--left',  default='obs/IGN_20-25.obs',
                        help='Left bulletin: a .obs, or archive.zip::member (default: obs/IGN_20-25.obs)')
    parser.add_argument('--right', default='obs/LDG_20-25.obs',
                        help='Right bulletin (default: obs/LDG_20-25.obs)')
    parser.add_argument('--output', default='complem_figures/source_depth_maps/IGN_LDG_report.png',
                        help='Output path; both .png and .pdf are written')
    parser.add_argument('--left-tag',  default='a', help='Top-left tag of the left map ("" for none)')
    parser.add_argument('--right-tag', default='b', help='Top-left tag of the right map ("" for none)')
    parser.add_argument('--map-region', nargs=4, type=float,
                        metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'), default=None,
                        help=f'Map extent (default: {" ".join(str(v) for v in DEFAULT_REGION)})')
    parser.add_argument('--width-cm', type=float, default=16.0,
                        help='Printed figure width, i.e. your \\linewidth in cm (default: 16)')
    parser.add_argument('--size-factor', type=float, default=1.0,
                        help='Extra multiplier on symbol sizes (default: 1)')
    parser.add_argument('--max-depth', type=float, default=15.0,
                        help='Colour-bar saturation depth in km (default: 15)')
    parser.add_argument('--no-convert-mags', action='store_true',
                        help='Keep source magnitudes (skip the in-memory ML LDG conversion)')
    args = parser.parse_args()

    generate_figure(SourceDepthMapsParams(
        left=args.left, right=args.right, figSave=args.output,
        left_tag=args.left_tag or None, right_tag=args.right_tag or None,
        map_region=args.map_region, width_cm=args.width_cm,
        size_factor=args.size_factor, max_depth=args.max_depth,
        convert_mags=not args.no_convert_mags,
    ))


if __name__ == '__main__':
    main()
