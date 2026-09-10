"""
temporary_network_impact.py
============================
Quantify what the temporary-network picks bought the catalog.

`add_temp_picks.py` turns obs/NLL_result.obs into obs/NLL_result_augmented.obs by
folding in picks from eight external deployments. The two bulletins hold the very
same events, with byte-identical hypocentres and origin times -- their event header
lines differ in exactly one field, PhaseCount. That makes them a controlled
before/after pair: anything that moves between them is the added picks alone, with
no relocation mixed in.

Two figures are produced from a single pass over each bulletin:

  phase_count_timeline.pdf     Two chronological strips, 1978 -> 2025, coloured by
                               the monthly median number of phases per event. Upper
                               strip before augmentation, lower strip after, on a
                               shared colour scale so the two are comparable.

  nearest_station_distance.pdf Overlapping histograms of the distance from each
                               event to its closest recording station, before and
                               after. Left panel all events, right panel restricted
                               to the events that actually gained picks.

Station coordinates are read from the NonLinLoc GTSRCE files, i.e. the exact
positions NLLoc itself located with. The resulting minimum distance reproduces
NLLoc's own published `Dist` column to within 0.1 km at the median.

Environment: seisbench_env (numpy/pandas, matplotlib/seaborn for the figures).

Usage
-----
    python complem_figures/temporary_network_impact.py

    # Only one of the two figures
    python complem_figures/temporary_network_impact.py --product timeline
    python complem_figures/temporary_network_impact.py --product distance

    # Widen the time bin, and clip the distance axis harder
    python complem_figures/temporary_network_impact.py \\
        --bin-months 3 --max-distance 40 --distance-bin 0.5
"""

import argparse
import glob
import os
import sys
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
from temp_picks.match_picks import haversine_km  # noqa: E402


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_COLOR_BEFORE = '#2c7fb8'
_COLOR_AFTER  = '#d95f02'
_COLOR_EMPTY  = '#d9d9d9'

_NEAR_KM = 5.0   # "a station within N km", reported in the distance panels


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class TempNetworkImpactParams:
    file_before:   str   = None    # bulletin before augmentation
    file_after:    str   = None    # bulletin after augmentation
    stations_glob: str   = None    # GTSRCE files holding the station coordinates
    output_dir:    str   = None    # where the PDFs land
    product:       str   = 'both'  # both | timeline | distance
    bin_months:    int   = 1       # time bin width of the timeline strips
    vmin:          float = None    # colour scale bounds; None -> data range
    vmax:          float = None
    max_distance:  float = 60.0    # distance axis limit (km)
    distance_bin:  float = 1.0     # distance histogram bin width (km)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_station_coords(stations_glob):
    """
    Read station positions from the NonLinLoc GTSRCE files.

    Lines look like `GTSRCE FR.0058 LATLON 43.255798 -1.261500 0.0 0.190`. The
    per-zone files overlap, so the first occurrence of a code wins.

    Parameters
    ----------
    stations_glob : str  -- glob matching the GTSRCE files

    Returns
    -------
    dict : unified station code -> (lat, lon)
    """
    coords = {}
    files  = sorted(glob.glob(stations_glob))
    if not files:
        raise FileNotFoundError(f'no station file matched {stations_glob!r}')

    for path in files:
        with open(path, 'r') as handle:
            for line in handle:
                parts = line.split()
                if len(parts) >= 5 and parts[0] == 'GTSRCE' and parts[2] == 'LATLON':
                    coords.setdefault(parts[1], (float(parts[3]), float(parts[4])))

    return coords


def _scan_bulletin(path, station_coords):
    """
    Reduce a .obs bulletin to one row per event, in a single streaming pass.

    Only the two quantities the figures need are kept: how many picks the event
    carries, and how far away its closest recording station sits. Station codes are
    deduplicated per event before the distance is computed -- most stations
    contribute both a P and an S pick, so this halves the work.

    Parameters
    ----------
    path           : str   -- .obs bulletin
    station_coords : dict  -- station code -> (lat, lon), from _load_station_coords

    Returns
    -------
    (pd.DataFrame, dict)
        Frame with columns public_id, year, month, n_phases, min_dist_km.
        Dict with the parse tallies: n_events, n_picks, unresolved (code -> count).
    """
    rows        = []
    unresolved  = {}
    n_picks     = 0

    public_id = None
    lat = lon = None
    stations  = set()
    year = month = None
    n_phases  = 0

    def _flush():
        if public_id is None:
            return
        distances = [haversine_km(lat, lon, *station_coords[code])
                     for code in stations]
        rows.append((public_id, year, month, n_phases,
                     min(distances) if distances else np.nan))

    with open(path, 'r') as handle:
        for line in handle:
            if line.startswith('# '):
                _flush()
                tokens    = line[2:].split()
                year      = int(tokens[0])
                month     = int(tokens[1])
                lat       = float(tokens[6])
                lon       = float(tokens[7])
                public_id = None
                stations  = set()
                n_phases  = 0
            elif line.startswith('PUBLIC_ID'):
                public_id = line.split()[1]
            elif line.startswith('###') or not line.strip():
                continue
            else:
                code      = line.split()[0]
                n_phases += 1
                n_picks  += 1
                if code in station_coords:
                    stations.add(code)
                else:
                    unresolved[code] = unresolved.get(code, 0) + 1

    _flush()

    frame = pd.DataFrame(rows, columns=['public_id', 'year', 'month',
                                        'n_phases', 'min_dist_km'])
    return frame, {'n_events': len(frame), 'n_picks': n_picks,
                   'unresolved': unresolved}


def _monthly_medians(frame, bin_months):
    """
    Median phase count per time bin, on a regular grid with no gaps.

    Bins that contain no event stay NaN rather than being dropped, so the strips
    keep a true time axis and quiet months read as absent instead of silently
    shifting everything that follows.

    Parameters
    ----------
    frame      : pd.DataFrame -- needs year, month, n_phases
    bin_months : int          -- bin width in months

    Returns
    -------
    (edges, medians) : both np.ndarray, in fractional years; len(edges) == len(medians) + 1
    """
    absolute = (frame['year'] * 12 + (frame['month'] - 1)).to_numpy()
    start    = absolute.min()
    index    = (absolute - start) // bin_months
    n_bins   = int(index.max()) + 1

    medians = np.full(n_bins, np.nan)
    grouped = pd.Series(frame['n_phases'].to_numpy()).groupby(index).median()
    medians[grouped.index.to_numpy().astype(int)] = grouped.to_numpy()

    edges = (start + np.arange(n_bins + 1) * bin_months) / 12.0
    return edges, medians


def _plot_timeline(before, after, output_path, parameters):
    """Two chronological strips coloured by the median phase count per time bin."""
    sns.set_theme()

    edges, med_before = _monthly_medians(before, parameters.bin_months)
    _,     med_after  = _monthly_medians(after,  parameters.bin_months)

    # A shared colour scale is the whole point: it is what makes the two strips
    # comparable rather than two independently-normalised pictures.
    finite = np.concatenate([med_before[np.isfinite(med_before)],
                             med_after[np.isfinite(med_after)]])
    vmin   = parameters.vmin if parameters.vmin is not None else float(finite.min())
    # A single month out of ~1100 reaches ~52 phases; scaling to it would spend half
    # the colour ramp on one bin and flatten the rest. Clip to p99 and let the
    # colorbar declare the overflow.
    vmax   = parameters.vmax if parameters.vmax is not None else float(np.percentile(finite, 99))

    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad(_COLOR_EMPTY)

    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(13, 4.2), layout='constrained')

    panels = ((axes[0], med_before, 'before'),
              (axes[1], med_after,  'after'))

    for axis, medians, label in panels:
        mesh = axis.pcolormesh(edges, [0, 1],
                               np.ma.masked_invalid(medians)[None, :],
                               cmap=cmap, shading='flat', vmin=vmin, vmax=vmax,
                               rasterized=True)
        axis.set_yticks([])
        axis.set_ylabel(f'{label}\ntemporary picks', rotation=0,
                        ha='right', va='center', fontsize=9)
        axis.grid(False)

    axes[1].set_xlabel('Year')
    axes[1].set_xlim(edges[0], edges[-1])
    axes[1].set_xticks(np.arange(np.ceil(edges[0] / 5) * 5, edges[-1], 5))
    axes[1].xaxis.set_major_formatter(lambda value, _: f'{value:.0f}')

    bin_label = ('monthly' if parameters.bin_months == 1
                 else f'{parameters.bin_months}-month')
    fig.colorbar(mesh, ax=axes, label=f'Median phases per event ({bin_label})',
                 pad=0.01, extend='max' if finite.max() > vmax else 'neither')
    fig.suptitle('Phases available per event, before and after the temporary networks',
                 fontweight='bold')

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    plt.close(fig)

    print(f'Figure saved @ {output_path}')
    return output_path


def _plot_distance(merged, output_path, parameters):
    """Overlapping before/after histograms of the distance to the closest station."""
    sns.set_theme()

    bins = np.arange(0.0, parameters.max_distance + parameters.distance_bin,
                     parameters.distance_bin)

    gained = merged['n_phases_after'] > merged['n_phases_before']
    panels = (('All events', merged),
              ('Events that gained picks', merged[gained]))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False,
                             layout='constrained')

    for axis, (title, subset) in zip(axes, panels):
        # Everything past the axis limit is folded into the last bin rather than
        # dropped, so the counts stay the true event counts (max is ~154 km).
        actual = {
            'before': subset['min_dist_km_before'].to_numpy(),
            'after':  subset['min_dist_km_after'].to_numpy(),
        }

        for label, colour in (('before', _COLOR_BEFORE), ('after', _COLOR_AFTER)):
            drawn  = np.clip(actual[label], None, bins[-1])
            median = np.median(actual[label])
            axis.hist(drawn, bins=bins, histtype='stepfilled',
                      color=colour, alpha=0.55, lw=0, zorder=2)
            axis.hist(drawn, bins=bins, histtype='step',
                      color=colour, lw=1.2, zorder=3)
            axis.axvline(median, color=colour, ls='--', lw=1.1, zorder=4)

        near_before = 100.0 * (actual['before'] < _NEAR_KM).mean()
        near_after  = 100.0 * (actual['after']  < _NEAR_KM).mean()
        closer      = int((actual['after'] < actual['before'] - 1e-3).sum())

        # Everything the figure has to say sits in the title: the panels carry no
        # legend and no annotation box, the colour key belongs to the caption.
        axis.set_xlim(0, parameters.max_distance)
        axis.set_xlabel('Distance to the closest recording station (km)')
        axis.set_ylabel('Events')
        axis.set_title(
            f'{title}   (N = {len(subset)})\n'
            f'median {np.median(actual["before"]):.2f} -> {np.median(actual["after"]):.2f} km   |   '
            f'within {_NEAR_KM:.0f} km {near_before:.1f} % -> {near_after:.1f} %   |   '
            f'closer {closer} ({100.0 * closer / len(subset):.1f} %)',
            fontsize=9)

    fig.suptitle('Distance to the closest recording station, '
                 'before and after the temporary networks', fontweight='bold')

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path)
    plt.close(fig)

    print(f'Figure saved @ {output_path}')
    return output_path


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figures(parameters):
    """
    Build the temporary-network impact figures.

    Parameters
    ----------
    parameters : TempNetworkImpactParams

    Returns
    -------
    dict with keys: output, n_events, n_picks_before, n_picks_after,
                    n_events_gained, median_dist_before, median_dist_after,
                    unresolved
    """
    file_before   = parameters.file_before   or os.path.join(_PROJECT_ROOT, 'obs', 'NLL_result.obs')
    file_after    = parameters.file_after    or os.path.join(_PROJECT_ROOT, 'obs', 'NLL_result_augmented.obs')
    stations_glob = parameters.stations_glob or os.path.join(_PROJECT_ROOT, 'stations', 'GTSRCE_*.txt')
    output_dir    = parameters.output_dir    or os.path.join(_MODULE_DIR, 'temporary_network_impact')

    station_coords = _load_station_coords(stations_glob)
    print(f'Stations         : {len(station_coords)}')

    before, stats_before = _scan_bulletin(file_before, station_coords)
    after,  stats_after  = _scan_bulletin(file_after,  station_coords)
    print(f'Before           : {stats_before["n_events"]} events / {stats_before["n_picks"]} picks')
    print(f'After            : {stats_after["n_events"]} events / {stats_after["n_picks"]} picks')

    unresolved = dict(stats_after['unresolved'])
    for code, count in stats_before['unresolved'].items():
        unresolved[code] = max(unresolved.get(code, 0), count)
    if unresolved:
        detail = ', '.join(f'{code} ({count})' for code, count in sorted(unresolved.items()))
        print(f'Unresolved codes : {len(unresolved)} -- {detail}')

    merged = before.merge(after, on='public_id', suffixes=('_before', '_after'))
    if len(merged) != len(before) or len(merged) != len(after):
        raise ValueError(
            f'the two bulletins are not a matched pair: {len(before)} before, '
            f'{len(after)} after, {len(merged)} joined on publicId'
        )

    outputs = []
    if parameters.product in ('both', 'timeline'):
        outputs.append(_plot_timeline(
            before, after,
            os.path.join(output_dir, 'phase_count_timeline.pdf'), parameters))
    if parameters.product in ('both', 'distance'):
        outputs.append(_plot_distance(
            merged,
            os.path.join(output_dir, 'nearest_station_distance.pdf'), parameters))

    return {
        'output':             outputs,
        'n_events':           len(merged),
        'n_picks_before':     stats_before['n_picks'],
        'n_picks_after':      stats_after['n_picks'],
        'n_events_gained':    int((merged['n_phases_after'] > merged['n_phases_before']).sum()),
        'median_dist_before': float(merged['min_dist_km_before'].median()),
        'median_dist_after':  float(merged['min_dist_km_after'].median()),
        'unresolved':         unresolved,
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Figures on what the temporary-network picks bought the catalog.'
    )
    parser.add_argument('--before',        default=None,
                        help='Bulletin before augmentation (default: obs/NLL_result.obs)')
    parser.add_argument('--after',         default=None,
                        help='Bulletin after augmentation (default: obs/NLL_result_augmented.obs)')
    parser.add_argument('--stations',      default=None,
                        help='Glob of the GTSRCE station files (default: stations/GTSRCE_*.txt)')
    parser.add_argument('--output-dir',    default=None,
                        help='Output folder (default: complem_figures/temporary_network_impact)')
    parser.add_argument('--product',       default='both',
                        choices=('both', 'timeline', 'distance'),
                        help='Which figure(s) to build (default: both)')
    parser.add_argument('--bin-months',    type=int,   default=1,
                        help='Time bin width of the timeline strips, in months (default: 1)')
    parser.add_argument('--vmin',          type=float, default=None,
                        help='Timeline colour scale minimum (default: data range)')
    parser.add_argument('--vmax',          type=float, default=None,
                        help='Timeline colour scale maximum (default: data range)')
    parser.add_argument('--max-distance',  type=float, default=60.0,
                        help='Distance axis limit in km (default: 60.0)')
    parser.add_argument('--distance-bin',  type=float, default=1.0,
                        help='Distance histogram bin width in km (default: 1.0)')
    args = parser.parse_args()

    generate_figures(TempNetworkImpactParams(
        file_before   = args.before,
        file_after    = args.after,
        stations_glob = args.stations,
        output_dir    = args.output_dir,
        product       = args.product,
        bin_months    = args.bin_months,
        vmin          = args.vmin,
        vmax          = args.vmax,
        max_distance  = args.max_distance,
        distance_bin  = args.distance_bin,
    ))


if __name__ == '__main__':
    main()
