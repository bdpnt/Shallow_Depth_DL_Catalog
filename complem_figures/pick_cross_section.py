"""
pick_cross_section.py
=====================
Interactive picker for a new cross-section: click its start and end points on
the seismicity map, set its half-width (and the depth shown) with the sliders,
and see the resulting depth section update live beside the map.

* Map: the catalogue (default RESULT/SSST_result.csv, PDF expectation), all
  events, coloured by depth, over the cached relief / borders of
  `basemap_lines`. Zoom and pan with the toolbar; clicks are ignored while a
  toolbar tool is active.
* Section: the events passing `cross_section_panels.py`'s quality filter
  (ErH < --uncert-h, ErV < --uncert-v, RMS < 0.5 s), coloured by ErV, at 1:1
  scale, worst errors drawn first. Only events between the two profile ends
  are kept.
* "Print command" prints the `cross_section_panels.py` call drawing that
  profile (lon0/lat0 = start point, azimuth, length, half-width, depth max).

Geometry: the flat-Earth helper of `cross_section.py` (1 deg = 111 km, cos of
the start latitude), so the printed azimuth and length reproduce the clicked
end point. The swath selection is done in numpy rather than by
`pygmt.project` (great circle) so it can redraw on every slider move; on the
Aneto, Arette and Andorra profiles the two differ by 1-3 events, all on the
swath edges.

Needs a display (matplotlib's interactive backend), seisbench_env.

Usage
-----
    python complem_figures/pick_cross_section.py

    # another catalogue, looser filter, start from a 5 km half-width
    python complem_figures/pick_cross_section.py --bulletin RESULT/NLL_result.csv \\
        --uncert-h 1.5 --uncert-v 1.5 --width 5
"""

import argparse
import os
import sys
from dataclasses import dataclass, field
from math import atan2, cos, degrees, hypot, radians
from typing import List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.widgets import Button, Slider

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
from complem_figures.basemap_lines import (  # noqa: E402
    ABOVE_SEA_COLOR, draw_basemap_lines, draw_relief)
from complem_figures.event_maps import (  # noqa: E402
    REPORT_REGION, compass_point, hypocentre_columns)

_KM_PER_DEG = 111.0   # cross_section._dest_point
_PANEL_TOP  = -1.0    # km, as cross_section_panels.py


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class PickCrossSectionParams:
    bulletin:   str   = 'RESULT/SSST_result.csv'
    solution:   str   = 'expect'     # 'expect' | 'maxlike'
    uncert_h:   float = 1.0          # cross_section_panels.py defaults
    uncert_v:   float = 1.0
    half_width: float = 3.0          # km, initial slider value
    depth_max:  float = 15.0         # km, initial slider value
    max_depth_colour: float = 15.0   # map colour saturation, as the PyGMT maps
    map_region: List[float] = field(default_factory=lambda: list(REPORT_REGION))


# ---------------------------------------------------------------------------
# Geometry (flat Earth, as cross_section._dest_point)
# ---------------------------------------------------------------------------

def profile_from_points(lon0, lat0, lon1, lat1):
    """Azimuth (deg from N) and length (km) from start to end, flat Earth at lat0."""
    dx = (lon1 - lon0) * _KM_PER_DEG * cos(radians(lat0))
    dy = (lat1 - lat0) * _KM_PER_DEG
    return degrees(atan2(dx, dy)) % 360.0, hypot(dx, dy)


def project_on_profile(lon, lat, lon0, lat0, azimut, length_km, half_width_km):
    """
    Distance along the profile of the events inside its swath.

    Returns
    -------
    (p, keep) — p: along-profile distance (km) of the kept events; keep: boolean
    mask over the input, |cross distance| <= half_width and 0 <= p <= length.
    """
    x  = (lon - lon0) * _KM_PER_DEG * cos(radians(lat0))
    y  = (lat - lat0) * _KM_PER_DEG
    az = radians(azimut)
    p  = x * np.sin(az) + y * np.cos(az)
    q  = x * np.cos(az) - y * np.sin(az)
    keep = (np.abs(q) <= half_width_km) & (p >= 0) & (p <= length_km)
    return p[keep], keep


def swath_outline(lon0, lat0, azimut, length_km, half_width_km):
    """Closed lon/lat outline of the swath (same corners as event_maps.profile_geometry)."""
    az = radians(azimut)
    ux, uy = np.sin(az), np.cos(az)      # along the profile, km east / north
    nx, ny = uy, -ux                     # perpendicular
    corners = [(-half_width_km, 0.0), (-half_width_km, length_km),
               (half_width_km, length_km), (half_width_km, 0.0), (-half_width_km, 0.0)]
    kx = _KM_PER_DEG * cos(radians(lat0))
    lons = [lon0 + (s * nx + t * ux) / kx for s, t in corners]
    lats = [lat0 + (s * ny + t * uy) / _KM_PER_DEG for s, t in corners]
    return lons, lats


# ---------------------------------------------------------------------------
# Catalogue
# ---------------------------------------------------------------------------

def _read_catalog(parameters):
    """All events (map) and the quality mask of cross_section_panels.py (section)."""
    path = parameters.bulletin
    path = path if os.path.isabs(path) else os.path.join(_PROJECT_ROOT, path)
    df   = pd.read_csv(path)
    lat_col, lon_col, depth_col = hypocentre_columns(df, parameters.solution)
    events = pd.DataFrame({
        'lon':   df[lon_col].to_numpy(),
        'lat':   df[lat_col].to_numpy(),
        'depth': df[depth_col].to_numpy(),
        'erh':   df['true_erh'].to_numpy(),
        'erv':   df['true_erz'].to_numpy(),
        'rms':   df['RMS'].to_numpy(),
    }).dropna(subset=['lon', 'lat', 'depth'])
    good = ((events.erv < parameters.uncert_v) & (events.erh < parameters.uncert_h)
            & (events.rms < 0.5)).to_numpy()
    print(f'Catalog read @ {path} — {len(events)} events, {good.sum()} pass the '
          f'section filter (ErH < {parameters.uncert_h}, ErV < {parameters.uncert_v}, RMS < 0.5)')
    return events, good


# ---------------------------------------------------------------------------
# Interactive picker
# ---------------------------------------------------------------------------

class ProfilePicker:
    """Map + live section; state is the start / end points and the two sliders."""

    def __init__(self, parameters):
        self.parameters = parameters
        self.events, good = _read_catalog(parameters)
        self.sec = self.events[good].reset_index(drop=True)
        self.start = None
        self.end   = None

        self.fig = plt.figure(figsize=(15, 7.5))
        self.ax_map = self.fig.add_axes([0.05, 0.22, 0.48, 0.72])
        self.ax_sec = self.fig.add_axes([0.60, 0.22, 0.33, 0.72])
        self.ax_cb  = self.fig.add_axes([0.945, 0.40, 0.012, 0.36])
        self._draw_map()

        self.s_width = Slider(self.fig.add_axes([0.12, 0.10, 0.32, 0.03]),
                              'Half-width (km)', 0.5, 20.0, valinit=parameters.half_width,
                              valstep=0.5)
        self.s_depth = Slider(self.fig.add_axes([0.12, 0.05, 0.32, 0.03]),
                              'Depth max (km)', 5.0, 40.0, valinit=parameters.depth_max,
                              valstep=1.0)
        self.b_print = Button(self.fig.add_axes([0.60, 0.05, 0.12, 0.05]), 'Print command')
        self.s_width.on_changed(lambda _: self.update())
        self.s_depth.on_changed(lambda _: self.update())
        self.b_print.on_clicked(lambda _: self.print_command())
        self.fig.canvas.mpl_connect('button_press_event', self._on_click)
        self.fig.canvas.mpl_connect('motion_notify_event', self._on_move)
        self.update()

    # -- drawing ------------------------------------------------------------

    def _draw_map(self):
        ax, p = self.ax_map, self.parameters
        lon_min, lon_max, lat_min, lat_max = p.map_region
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_aspect(1.0 / cos(radians(0.5 * (lat_min + lat_max))))
        draw_relief(ax)
        draw_basemap_lines(ax)
        ev   = self.events
        cmap = plt.get_cmap('viridis_r').copy()
        cmap.set_under(ABOVE_SEA_COLOR)
        cmap.set_over('black')
        ax.scatter(ev.lon, ev.lat, c=ev.depth, s=1.5, cmap=cmap, vmin=0,
                   vmax=p.max_depth_colour, linewidths=0, zorder=2, rasterized=True)
        self.swath_line, = ax.plot([], [], color='red', lw=1.0, ls='--', zorder=4)
        self.axis_line,  = ax.plot([], [], color='red', lw=2.0, zorder=5)
        self.end_marks,  = ax.plot([], [], 'o', color='red', ms=5, zorder=6)
        self.labels = [ax.text(0, 0, '', color='red', fontweight='bold', zorder=7,
                               ha='center', va='center') for _ in range(2)]
        ax.set_title('Click the start point, then the end point (a third click restarts)',
                     fontsize=10)

    def _on_click(self, event):
        if (event.inaxes is not self.ax_map or event.button != 1
                or self.fig.canvas.toolbar is not None and self.fig.canvas.toolbar.mode):
            return
        if self.start is None or self.end is not None:
            self.start, self.end = (event.xdata, event.ydata), None
        else:
            self.end = (event.xdata, event.ydata)
        self.update()

    def _on_move(self, event):
        # rubber band from the start point while the end is not chosen yet
        if self.start is None or self.end is not None or event.inaxes is not self.ax_map:
            return
        self.axis_line.set_data([self.start[0], event.xdata], [self.start[1], event.ydata])
        self.fig.canvas.draw_idle()

    def profile(self):
        """(lon0, lat0, azimut, length_km) of the picked profile, or None."""
        if self.start is None or self.end is None:
            return None
        azimut, length = profile_from_points(*self.start, *self.end)
        return self.start[0], self.start[1], azimut, length

    def update(self):
        prof = self.profile()
        for lab in self.labels:
            lab.set_text('')
        if prof is None:
            self.swath_line.set_data([], [])
            self.end_marks.set_data(*([self.start[0]], [self.start[1]]) if self.start else ([], []))
            self.axis_line.set_data([], [])
            self._draw_section(None)
            self.fig.canvas.draw_idle()
            return

        lon0, lat0, azimut, length = prof
        half_w = self.s_width.val
        self.swath_line.set_data(*swath_outline(lon0, lat0, azimut, length, half_w))
        self.axis_line.set_data([self.start[0], self.end[0]], [self.start[1], self.end[1]])
        self.end_marks.set_data([self.start[0], self.end[0]], [self.start[1], self.end[1]])
        for lab, (x, y), az in ((self.labels[0], self.start, azimut + 180.0),
                                (self.labels[1], self.end, azimut)):
            dx = 0.04 * np.sin(radians(az)) / cos(radians(lat0))
            dy = 0.04 * np.cos(radians(az))
            lab.set_position((x + dx, y + dy))
            lab.set_text(compass_point(az))
        self._draw_section(prof)
        self.fig.canvas.draw_idle()

    def _draw_section(self, prof):
        ax = self.ax_sec
        ax.clear()
        self.ax_cb.clear()
        depth_max = self.s_depth.val
        if prof is None:
            ax.text(0.5, 0.5, 'Pick a profile on the map', ha='center', va='center',
                    transform=ax.transAxes, color='0.4')
            ax.set_xticks([])
            ax.set_yticks([])
            self.ax_cb.set_visible(False)
            return

        lon0, lat0, azimut, length = prof
        p, keep = project_on_profile(self.sec.lon.to_numpy(), self.sec.lat.to_numpy(),
                                     lon0, lat0, azimut, length, self.s_width.val)
        z, erv  = self.sec.depth.to_numpy()[keep], self.sec.erv.to_numpy()[keep]
        order   = np.argsort(erv)[::-1]                # worst first, best on top
        visible = int(((z >= _PANEL_TOP) & (z <= depth_max)).sum())

        sc = ax.scatter(p[order], z[order], c=np.clip(erv[order], 0, self.parameters.uncert_v),
                        cmap='magma_r', vmin=0, vmax=self.parameters.uncert_v,
                        s=12, edgecolors='black', linewidths=0.25)
        ax.axhline(0, color='0.5', lw=0.5)
        ax.set_xlim(0, length)
        ax.set_ylim(depth_max, _PANEL_TOP)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel('Distance (km)')
        ax.set_ylabel('Depth (km)')
        ax.text(0, 1.01, compass_point(azimut + 180.0), transform=ax.transAxes,
                ha='left', va='bottom', fontweight='bold')
        ax.text(1, 1.01, compass_point(azimut), transform=ax.transAxes,
                ha='right', va='bottom', fontweight='bold')
        ax.set_title(f'{visible} events shown  ·  az {azimut:.1f}°  ·  {length:.1f} km  '
                     f'·  ±{self.s_width.val:g} km', fontsize=10, pad=16)
        self.ax_cb.set_visible(True)
        self.fig.colorbar(sc, cax=self.ax_cb, label='ErV (km)')

    # -- output ---------------------------------------------------------------

    def print_command(self):
        prof = self.profile()
        if prof is None:
            print('No profile picked yet: click a start and an end point on the map.')
            return None
        lon0, lat0, azimut, length = prof
        cmd = (f'/opt/miniconda3/envs/pygmt_env/bin/python complem_figures/cross_section_panels.py \\\n'
               f'    --lon0 {lon0:.4f} --lat0 {lat0:.4f} --azimut {azimut:.1f} '
               f'--length {length:.1f} --width {self.s_width.val:g} '
               f'--depth-max {self.s_depth.val:g} \\\n'
               f'    --output complem_figures/cross_section/NEW_PROFILE.png')
        print(cmd)
        return cmd


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def pick_cross_section(parameters):
    """
    Open the picker window and block until it is closed.

    Parameters
    ----------
    parameters : PickCrossSectionParams

    Returns
    -------
    dict with keys: profile — (lon0, lat0, azimut, length_km, half_width_km) of
    the last picked profile, or None; command — the matching
    cross_section_panels.py call, or None.
    """
    picker = ProfilePicker(parameters)
    plt.show()
    prof = picker.profile()
    return {
        'profile': None if prof is None else (*prof, picker.s_width.val),
        'command': picker.print_command() if prof is not None else None,
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Pick a cross-section on the seismicity map and preview it live.')
    parser.add_argument('--bulletin', default='RESULT/SSST_result.csv',
                        help='RESULT/*.csv catalogue (default: %(default)s)')
    parser.add_argument('--solution', default='expect', choices=['expect', 'maxlike'],
                        help='PDF expectation (default) or maximum-likelihood hypocentre')
    parser.add_argument('--uncert-h', type=float, default=1.0,
                        help='section filter: ErH below, km (default: 1.0, as cross_section_panels.py)')
    parser.add_argument('--uncert-v', type=float, default=1.0,
                        help='section filter and ErV colour maximum, km (default: 1.0)')
    parser.add_argument('--width', type=float, default=3.0,
                        help='initial half-width, km (default: 3)')
    parser.add_argument('--depth-max', type=float, default=15.0,
                        help='initial depth shown, km (default: 15)')
    parser.add_argument('--map-region', nargs=4, type=float,
                        metavar=('LON_MIN', 'LON_MAX', 'LAT_MIN', 'LAT_MAX'),
                        help=f'initial map extent (default: {" ".join(map(str, REPORT_REGION))})')
    args = parser.parse_args()

    pick_cross_section(PickCrossSectionParams(
        bulletin   = args.bulletin,
        solution   = args.solution,
        uncert_h   = args.uncert_h,
        uncert_v   = args.uncert_v,
        half_width = args.width,
        depth_max  = args.depth_max,
        map_region = args.map_region or list(REPORT_REGION),
    ))


if __name__ == '__main__':
    main()
