"""
cross_section_report.py
=======================
Report-ready version of `cross_section.py --no-map`: three depth sections along
the same profile side by side — by default the Maladeta/Aneto profile before
relocation (a, flat colour), after NonLinLoc (b) and after NonLinLoc-SSST (c),
the last two coloured by ErV.

Differences from `cross_section.py`
-----------------------------------
* The three panels are one figure, drawn at its printed size (`--width-cm`,
  default 16 cm = `\\linewidth` of an A4 page with 2.5 cm margins) with 7-8 pt
  fonts, so LaTeX includes it at scale 1 and the text stays legible.
* Each panel keeps the 1:1 distance/depth scale of the original.
* Depth labels on panel a only; one ErV colour bar, to the right of the last
  panel; a bold tag (a, b, c) in each panel's top-right corner; the sea-level
  line drawn explicitly.
* Symbol size is the original 0.15 cm rescaled by panel width (the original
  panel is 10 cm wide), so the sections look as dense as before.
* Catalogue reading (formats 4, 5 and 6), the quality filter, the swath selection
  by `pygmt.project`, the draw order and the clipping to the colour range are
  copied from `cross_section.py` unchanged. Events kept per panel are printed:
  state them in the caption.
* Writes both a 300 dpi PNG and a vector PDF.

Usage
-----
    conda run -n pygmt_env python complem_figures/cross_section_report.py

    # other panels / profile: --panel CATALOG FORMAT METRIC, repeatable (METRIC 'none' = flat colour)
    conda run -n pygmt_env python complem_figures/cross_section_report.py \\
        --panel obs/GLOBAL.obs 4 none \\
        --panel RESULT/NLL_result.csv 6 erv \\
        --panel RESULT/SSST_result.csv 6 erv \\
        --lon0 0.6852 --lat0 42.6 --azimut 44 --length 12 --width 3 --depth-max 12 \\
        --output complem_figures/cross_section/Aneto_report.png

    # Arette (Chaînons Béarnais): before relocation, NonLinLoc-SSST, Chevrot et al. (2024)
    conda run -n pygmt_env python complem_figures/cross_section_report.py \\
        --panel obs/GLOBAL.obs 4 none \\
        --panel RESULT/SSST_result.csv 6 erv \\
        --panel complem_figures/cross_section/catalogue_arette_3D_Chevrot.csv 5 none \\
        --lon0 -0.6275 --lat0 43.0 --azimut 0 --length 16 --width 8 --depth-max 18 \\
        --tag-position BL --output complem_figures/cross_section/Arette_report.png

    # Andorra: before relocation, NonLinLoc, NonLinLoc-SSST (default panels)
    conda run -n pygmt_env python complem_figures/cross_section_report.py \\
        --lon0 1.4 --lat0 42.48 --azimut 32 --length 10 --width 3 --depth-max 15 \\
        --output complem_figures/cross_section/Andorra_report.png
"""

import argparse
import os
import tempfile
from dataclasses import dataclass, field
from math import cos, radians, sin
from typing import List, Tuple

import numpy as np
import pandas as pd
import pygmt

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

_FLAT_FILL        = '#DF2B2B'   # cross_section.py, formats 4 and 5
_ORIG_PANEL_W     = 10.0        # cm, cross_section.py panel width
_ORIG_SYMBOL      = 0.15        # cm
_PANEL_TOP        = -1.0        # km, as cross_section.py

DEFAULT_PANELS = [
    ('obs/GLOBAL.obs',         4, 'none'),
    ('RESULT/NLL_result.csv',  6, 'erv'),
    ('RESULT/SSST_result.csv', 6, 'erv'),
]


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class CrossSectionReportParams:
    panels:       List[Tuple[str, int, str]] = field(default_factory=lambda: list(DEFAULT_PANELS))
    save_file:    str   = 'complem_figures/cross_section/Aneto_report.png'
    lon0:         float = 0.6852      # Aneto profile of the report
    lat0:         float = 42.6
    azimut:       float = 44.0
    length:       float = 12.0
    half_width:   float = 3.0
    depth_max:    float = 12.0
    uncert_h:     float = 1.0
    uncert_v:     float = 1.0
    width_cm:     float = 16.0
    size_factor:  float = 1.0
    tags:         str   = 'abc'
    tag_position: str   = 'TR'     # GMT justification code of the tag corner


# ---------------------------------------------------------------------------
# Helpers — copied from cross_section.py
# ---------------------------------------------------------------------------

def _dest_point(lon, lat, azimut, dist_km):
    """Flat-Earth destination point (1 deg = 111 km), as cross_section.py."""
    R    = 111.0
    az   = radians(azimut)
    dlat = (dist_km * cos(az)) / R
    dlon = (dist_km * sin(az)) / (R * cos(radians(lat)))
    return lon + dlon, lat + dlat


def _metric_style(name, uncert_h, uncert_v):
    """cross_section._metric_style restricted to the error metrics."""
    if name == 'erv':
        return dict(cmap='magma', series=[0, uncert_v], reverse=True,
                    label='ErV (km)', higher_is_better=False)
    if name == 'erh':
        return dict(cmap='magma', series=[0, uncert_h], reverse=True,
                    label='ErH (km)', higher_is_better=False)
    raise ValueError(f'unsupported metric {name!r} (erv, erh or none)')


def _read_catalog(path, fmt, uncert_h, uncert_v):
    """Formats 4 (.obs), 5 (Chevrot CSV) and 6 (RESULT/*.csv) of cross_section.py, then its quality filter."""
    if fmt == 6:
        df    = pd.read_csv(path)
        lon   = df['longitude'].to_numpy()
        lat   = df['latitude'].to_numpy()
        depth = df['depth'].to_numpy()
        erh   = df['true_erh'].to_numpy()
        erv   = df['true_erz'].to_numpy()
        rms   = df['RMS'].to_numpy()
    elif fmt == 5:   # Chevrot CSV: #YEAR MONTH DAY HOUR MINUTE SECOND LAT LON DEPTH MAG
        data  = np.loadtxt(path, comments='#', ndmin=2)
        lat, lon, depth = data[:, 6], data[:, 7], data[:, 8]
        erh = erv = rms = np.full(len(depth), -1.0)
    elif fmt == 4:
        lat, lon, depth = [], [], []
        with open(path, 'r') as f:
            for line in f:
                if line.startswith('# '):
                    d = line[2:].split()
                    lat.append(float(d[6]))
                    lon.append(float(d[7]))
                    depth.append(float(d[8]))
        lat, lon, depth = np.array(lat), np.array(lon), np.array(depth)
        erh = erv = rms = np.full(len(depth), -1.0)
    else:
        raise ValueError(f'format {fmt} not supported here (4, 5 or 6)')

    mask = (erv < uncert_v) & (erh < uncert_h) & (rms < 0.5)
    return lon[mask], lat[mask], depth[mask], {'erh': erh[mask], 'erv': erv[mask]}


def _project(lon, lat, depth, cval, p0, p1, half_width, workdir):
    """Swath selection exactly as cross_section.py: pygmt.project, pz convention, km."""
    cross_file = os.path.join(workdir, 'cross.dat')
    pygmt.project(
        data=np.column_stack((lon, lat, depth, cval)),
        center=list(p0), endpoint=list(p1),
        width=[-half_width, half_width],
        convention='pz', unit=True,
        outfile=cross_file, output_type='file',
    )
    if not os.path.exists(cross_file) or os.path.getsize(cross_file) == 0:
        return np.empty((0, 3))
    return np.loadtxt(cross_file, ndmin=2)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    n_panels = len(parameters.panels)
    p0 = (parameters.lon0, parameters.lat0)
    p1 = _dest_point(parameters.lon0, parameters.lat0, parameters.azimut, parameters.length)

    left, gap, cbar_block = 1.0, 0.3, 1.9          # cm
    has_cbar = any(m != 'none' for _, _, m in parameters.panels)
    panel_w  = (parameters.width_cm - left - gap * (n_panels - 1)
                - (cbar_block if has_cbar else 0.2)) / n_panels
    panel_h  = panel_w * (parameters.depth_max - _PANEL_TOP) / parameters.length   # 1:1
    # annotate every 2 km when 2 km spans at least 0.6 cm on the page, else every 5 km,
    # so 7 pt labels never crowd (1:1 scale: the same cm/km on both axes)
    cm_per_km = panel_w / parameters.length
    xa = ya = 2 if 2 * cm_per_km >= 0.6 else 5
    symbol   = _ORIG_SYMBOL * panel_w / _ORIG_PANEL_W * parameters.size_factor

    pygmt.config(
        FONT_ANNOT_PRIMARY='7p,Helvetica,black',
        FONT_LABEL='8p,Helvetica,black',
        MAP_FRAME_PEN='0.6p,black',
        MAP_TICK_PEN_PRIMARY='0.5p,black',
        MAP_TICK_LENGTH_PRIMARY='2.5p',
        MAP_ANNOT_OFFSET_PRIMARY='2p',
        MAP_LABEL_OFFSET='3p',
    )

    fig    = pygmt.Figure()
    counts = []
    style  = None
    with tempfile.TemporaryDirectory() as workdir:
        for i, (path, fmt, metric) in enumerate(parameters.panels):
            if i == 0:
                fig.shift_origin(xshift=f'{left}c', yshift='1.2c')
            else:
                fig.shift_origin(xshift=f'{panel_w + gap:.4f}c')

            lon, lat, depth, metrics = _read_catalog(
                os.path.join(_PROJECT_ROOT, path) if not os.path.isabs(path) else path,
                fmt, parameters.uncert_h, parameters.uncert_v)
            cvals = metrics[metric] if metric != 'none' else np.zeros_like(depth)
            data  = _project(lon, lat, depth, cvals, p0, p1, parameters.half_width, workdir)
            counts.append(len(data))

            fig.basemap(
                projection=f'X{panel_w:.4f}c/-{panel_h:.4f}c',
                region=[0, parameters.length, _PANEL_TOP, parameters.depth_max],
                frame=[f'xa{xa}f1+lDistance (km)',
                       f'ya{ya}f1+lDepth (km)' if i == 0 else f'ya{ya}f1',
                       'WSen' if i == 0 else 'wSen'],
            )
            fig.plot(x=[0, parameters.length], y=[0, 0], pen='0.25p,gray50')

            if len(data):
                X, Z, cval = data[:, 0], data[:, 1], data[:, 2]
                if metric == 'none':
                    fig.plot(x=X, y=Z, style=f'c{symbol:.4f}c', fill=_FLAT_FILL, pen='0.15p,black')
                else:
                    style = _metric_style(metric, parameters.uncert_h, parameters.uncert_v)
                    order = np.argsort(cval)
                    if not style['higher_is_better']:
                        order = order[::-1]            # worst first, best on top
                    shown = np.clip(cval[order], *style['series'])
                    pygmt.makecpt(cmap=style['cmap'], series=style['series'],
                                  reverse=style['reverse'])
                    fig.plot(x=X[order], y=Z[order], style=f'c{symbol:.4f}c',
                             fill=shown, cmap=True, pen='0.15p,black')

            tag = parameters.tags[i] if i < len(parameters.tags) else ''
            if tag:
                pos = parameters.tag_position
                dx  = '-0.12c' if pos[1] == 'R' else '0.12c'
                dy  = '-0.08c' if pos[0] == 'T' else '0.08c'
                fig.text(text=tag, position=pos, justify=pos, offset=f'{dx}/{dy}',
                         font='10p,Helvetica-Bold,black')
            print(f'panel {tag or i}: {path} (format {fmt}, {metric}) — '
                  f'{len(lon)} events after filter, {len(data)} in the swath')

    if style is not None:
        pygmt.makecpt(cmap=style['cmap'], series=style['series'], reverse=style['reverse'])
        # full panel height, clear of the frame ticks, label along the bar
        fig.colorbar(frame=[f"xa0.2f0.1+l{style['label']}"],
                     position=f'JMR+w{panel_h:.3f}c/0.3c+o0.45c/0c+v')

    base, _ = os.path.splitext(parameters.save_file)
    outputs = []
    for path in (f'{base}.png', f'{base}.pdf'):
        path = path if os.path.isabs(path) else os.path.join(_PROJECT_ROOT, path)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fig.savefig(path, dpi=300, crop=True)
        outputs.append(path)
        print(f'Figure saved @ {path}')
    return {'output': outputs[0], 'outputs': outputs, 'counts': counts}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Report-ready side-by-side depth sections (default: Aneto, before / NLL / NLL-SSST).')
    parser.add_argument('--panel', action='append', nargs=3,
                        metavar=('CATALOG', 'FORMAT', 'METRIC'), default=None,
                        help="repeatable; FORMAT 4 (.obs), 5 (Chevrot CSV) or 6 (RESULT/*.csv); "
                             "METRIC erv, erh or none (flat colour). Default: "
                             "GLOBAL.obs 4 none, NLL_result.csv 6 erv, SSST_result.csv 6 erv")
    parser.add_argument('--output', default='complem_figures/cross_section/Aneto_report.png',
                        help='both .png and .pdf are written (default: %(default)s)')
    parser.add_argument('--lon0',      type=float, default=0.6852)
    parser.add_argument('--lat0',      type=float, default=42.6)
    parser.add_argument('--azimut',    type=float, default=44.0)
    parser.add_argument('--length',    type=float, default=12.0, help='km (default: 12)')
    parser.add_argument('--width',     type=float, default=3.0, help='half-width, km (default: 3)')
    parser.add_argument('--depth-max', type=float, default=12.0, help='km (default: 12)')
    parser.add_argument('--uncert-h',  type=float, default=1.0,
                        help='ErH filter in km, format 6 only (default: 1.0)')
    parser.add_argument('--uncert-v',  type=float, default=1.0,
                        help='ErV filter and colour-bar maximum in km, format 6 only (default: 1.0)')
    parser.add_argument('--width-cm',  type=float, default=16.0,
                        help='printed width, i.e. your \\linewidth in cm (default: 16)')
    parser.add_argument('--size-factor', type=float, default=1.0,
                        help='extra multiplier on symbol size (default: 1)')
    parser.add_argument('--tags', default='abc', help='panel tags, one character each ("" for none)')
    parser.add_argument('--tag-position', default='TR', choices=['TL', 'TR', 'BL', 'BR'],
                        help='panel corner holding the tag (default: TR)')
    args = parser.parse_args()

    panels = ([(c, int(f), m) for c, f, m in args.panel] if args.panel else list(DEFAULT_PANELS))
    generate_figure(CrossSectionReportParams(
        panels=panels, save_file=args.output,
        lon0=args.lon0, lat0=args.lat0, azimut=args.azimut, length=args.length,
        half_width=args.width, depth_max=args.depth_max,
        uncert_h=args.uncert_h, uncert_v=args.uncert_v,
        width_cm=args.width_cm, size_factor=args.size_factor, tags=args.tags,
        tag_position=args.tag_position,
    ))


if __name__ == '__main__':
    main()
