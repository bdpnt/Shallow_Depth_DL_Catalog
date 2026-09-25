"""
basemap_lines.py
================
National borders and shorelines for the *matplotlib* report figures, taken
from the same GMT database as the PyGMT maps.

Why a cache
-----------
The PyGMT figures draw their borders with `fig.coast(..., borders='1/...')`
(see `aoi_map_report.py`).  The matplotlib report figures cannot call
`fig.coast`, and no coastline database ships with matplotlib.  So the geometry
is dumped **once** with `gmt coast -M` into multi-segment `.xy` files under
`complem_figures/basemap/`, and every matplotlib figure reads it from there.

The dump uses the same database, resolution and area threshold as
`aoi_map_report.py`'s PyGMT backend (`-Di`, `-A0/0/1`, `-N1`), so the lines
drawn on a matplotlib map are the same lines as on the PyGMT maps.

Dump the caches once, in the environment that has GMT:

    conda run -n pygmt_env python complem_figures/basemap_lines.py --dump

The default dump region (-3.5/4.5/41/45) covers every report map of the
deliverable, so this is a one-off; re-run it only for a map outside that box
(`--region LON_MIN LON_MAX LAT_MIN LAT_MAX`) or at another resolution
(`--resolution h`).

If the cache is missing, `draw_basemap_lines` prints a warning and draws
nothing, so the figure scripts still run in an environment without GMT.

Usage in a figure script
------------------------
    from basemap_lines import draw_basemap_lines
    draw_basemap_lines(ax)                       # borders + shorelines, under the data
    draw_basemap_lines(ax, coastlines=False)     # borders only
    draw_basemap_lines(ax, zorder=5, casing=True)  # on top of the data instead
"""

import argparse
import datetime as _dt
import os
import shutil
import subprocess
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Paths and dump parameters
# ---------------------------------------------------------------------------

_MODULE_DIR  = os.path.dirname(os.path.abspath(__file__))
_BASEMAP_DIR = os.path.join(_MODULE_DIR, 'basemap')

# Generous enough for every map of the deliverable, so the dump is a one-off.
DEFAULT_REGION     = [-3.5, 4.5, 41.0, 45.0]
DEFAULT_RESOLUTION = 'i'          # as aoi_map_report.py's fig.coast(resolution='i')
AREA_THRESH        = '0/0/1'      # as aoi_map_report.py's area_thresh

# GMT flag selecting what `gmt coast -M` dumps
_DUMP_FLAGS = {
    'borders':    ['-N1'],                     # level 1 = national borders
    'shorelines': ['-W', f'-A{AREA_THRESH}'],
}

# ---------------------------------------------------------------------------
# Drawing defaults
#
# The lines are geographic reference, not the message of the figure, so they
# are drawn *under* the data (`DEFAULT_ZORDER` sits above the grid at 0.5 and
# below the mesh and the epicentres at 1): they show through wherever the map
# has no data and are covered wherever it has.  Hence a light grey line and no
# white casing — the casing is only useful for a line drawn on top, and
# `casing=True` restores it.
# ---------------------------------------------------------------------------

BORDER_COLOR   = '#6E6E6E'
BORDER_LW      = 0.5
COAST_COLOR    = '#8A8A8A'
COAST_LW       = 0.4
DEFAULT_ZORDER = 0.8
CASING_COLOR   = 'white'
CASING_EXTRA   = 1.0      # points added to the line width for the casing
CASING_ALPHA   = 0.75


# ---------------------------------------------------------------------------
# Cache
# ---------------------------------------------------------------------------

def cache_path(kind, resolution=DEFAULT_RESOLUTION):
    """Path of the `.xy` cache for 'borders' or 'shorelines' at `resolution`."""
    return os.path.join(_BASEMAP_DIR, f'{kind}_{resolution}.xy')


def dump(region=None, resolution=DEFAULT_RESOLUTION, kinds=('borders', 'shorelines')):
    """
    Write the `.xy` caches with `gmt coast -M`.  Needs GMT on PATH.

    Returns the list of files written.
    """
    region = list(region or DEFAULT_REGION)
    gmt = shutil.which('gmt')
    if gmt is None:
        raise RuntimeError(
            'GMT not found on PATH.  Run this dump in the environment that has it, '
            'e.g.  conda run -n pygmt_env python complem_figures/basemap_lines.py --dump')

    os.makedirs(_BASEMAP_DIR, exist_ok=True)
    stamp   = _dt.date.today().isoformat()
    written = []
    for kind in kinds:
        flags = _DUMP_FLAGS[kind]
        args = [f'-R{region[0]}/{region[1]}/{region[2]}/{region[3]}',
                f'-D{resolution}', *flags, '-M']
        # `gmt coast -M` is the GMT 6 spelling; `gmt pscoast -M` the classic one.
        for module in ('coast', 'pscoast'):
            cmd  = [gmt, module, *args]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode == 0 and proc.stdout.strip():
                break
        else:
            raise RuntimeError(f'{" ".join(cmd)}\n{proc.stderr.strip()}')
        body = proc.stdout
        path = cache_path(kind, resolution)
        with open(path, 'w', encoding='utf-8') as f:
            f.write(f'# {kind} dumped by complem_figures/basemap_lines.py on {stamp}\n')
            f.write(f'# {" ".join(cmd)}\n')
            f.write(body if body.endswith('\n') else body + '\n')
        n = len(load(kind, resolution) or [])
        print(f'{path}: {n} segments')
        written.append(path)
    return written


def load(kind, resolution=DEFAULT_RESOLUTION):
    """
    Read a `.xy` cache into a list of (N, 2) lon/lat arrays.

    Returns None if the cache does not exist.  GMT multi-segment format: '>'
    opens a segment, '#' is a comment, every other line is 'lon lat'.
    """
    path = cache_path(kind, resolution)
    if not os.path.exists(path):
        return None
    segments, current = [], []
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            if s.startswith('>'):
                if len(current) > 1:
                    segments.append(np.asarray(current, dtype=float))
                current = []
                continue
            parts = s.split()
            try:
                current.append((float(parts[0]), float(parts[1])))
            except (IndexError, ValueError):
                continue
    if len(current) > 1:
        segments.append(np.asarray(current, dtype=float))
    return segments


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def draw_basemap_lines(ax, borders=True, coastlines=True,
                       resolution=DEFAULT_RESOLUTION,
                       border_color=BORDER_COLOR, border_lw=BORDER_LW,
                       coast_color=COAST_COLOR,   coast_lw=COAST_LW,
                       casing=False, zorder=DEFAULT_ZORDER, verbose=True):
    """
    Draw the cached borders / shorelines on a matplotlib lon-lat axis.

    Drawn under the data by default (`zorder=0.8`, i.e. above an `axisbelow`
    grid and below the mesh and the scatter); pass a `zorder` above 1 and
    `casing=True` to put them on top instead.

    Nothing is drawn (and a warning is printed) when the cache is missing, so
    a figure script never fails for the want of GMT.  Returns the number of
    segments drawn per kind.
    """
    from matplotlib import patheffects as pe
    from matplotlib.collections import LineCollection

    drawn = {}
    wanted = (('shorelines', coastlines, coast_color, coast_lw),
              ('borders',    borders,    border_color, border_lw))
    for kind, on, color, lw in wanted:
        if not on:
            continue
        segments = load(kind, resolution)
        if segments is None:
            if verbose:
                print(f'WARNING: no {kind} cache at {cache_path(kind, resolution)} — '
                      f'{kind} not drawn.  Create it once with:\n'
                      f'  conda run -n pygmt_env python complem_figures/basemap_lines.py --dump')
            drawn[kind] = 0
            continue
        effects = ([pe.withStroke(linewidth=lw + CASING_EXTRA,
                                  foreground=CASING_COLOR, alpha=CASING_ALPHA)]
                   if casing else None)
        lc = LineCollection(segments, colors=color, linewidths=lw,
                            capstyle='round', joinstyle='round',
                            zorder=zorder, path_effects=effects)
        ax.add_collection(lc, autolim=False)
        drawn[kind] = len(segments)
    return drawn


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Dump the GMT national borders / shorelines used by the matplotlib report figures.')
    parser.add_argument('--dump', action='store_true', help='write the .xy caches (needs GMT)')
    parser.add_argument('--region', nargs=4, type=float, default=DEFAULT_REGION,
                        metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'))
    parser.add_argument('--resolution', default=DEFAULT_RESOLUTION,
                        choices=['c', 'l', 'i', 'h', 'f'], help='GMT -D resolution (default: i)')
    parser.add_argument('--kinds', nargs='+', default=['borders', 'shorelines'],
                        choices=['borders', 'shorelines'])
    args = parser.parse_args()

    if args.dump:
        dump(region=args.region, resolution=args.resolution, kinds=tuple(args.kinds))
        return
    for kind in args.kinds:
        segments = load(kind, args.resolution)
        path = cache_path(kind, args.resolution)
        print(f'{path}: ' + ('missing — run with --dump' if segments is None
                             else f'{len(segments)} segments'))


if __name__ == '__main__':
    sys.exit(main())
