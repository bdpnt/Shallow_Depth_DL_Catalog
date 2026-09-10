"""
station_colocation.py
============================
Diagnostic for the 20 m co-location radius of the station inventory fusion
(stage 1, fetch_inventory/merge_station_inventories.py).

Compares what the merge *saw* -- the 2 414 stations of the source inventories --
against what it *did* -- the groups recorded in stations/GLOBAL_inventory.xml --
so that the threshold can be judged against the observed distribution of
station separations rather than taken on trust.

"Merged" is read from the outcome, never assumed from the distance: a pair counts
as merged when both stations end up under the same alternate code, or when one was
absorbed into the other by the within-network duplicate merge. That is what lets
the figure show pairs *beyond* the threshold that merged anyway through the
single-linkage chaining of _combine_close_stations(), and pairs *inside* it that
did not, because check_inventory() removed one of them interactively.

Three panels:
  1. nearest-neighbour distance per station, split merged vs not merged
  2. pairs that would merge as a function of the threshold
  3. per-network station counts, merged share, and median nearest-neighbour distance

Usage
-----
    python complem_figures/station_colocation.py
    python complem_figures/station_colocation.py --threshold 30
"""

import argparse
import glob
import math
import os
import sys
from dataclasses import dataclass

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory
from obspy.core.inventory import Inventory

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from fetch_inventory.merge_station_inventories import haversine  # noqa: E402

_EARTH_RADIUS_KM = 6371.0
_MAX_THRESHOLD_M = 200      # x-range of the sensitivity panel


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class StationColocationParams:
    """
    Configuration for the co-location diagnostic.

    Attributes
    ----------
    source_glob   : str — glob for the source StationXML files (pre-merge)
    file_inventory: str — path to the merged GLOBAL_inventory.xml
    file_output   : str — output PDF path
    threshold     : int — co-location radius in metres, marked on the figure
    """
    source_glob:    str
    file_inventory: str
    file_output:    str
    threshold:      int = 20


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def pairwise_distances(latitudes, longitudes):
    """
    Return the full great-circle distance matrix in metres, diagonal set to inf.

    Vectorized through the chord length rather than looped over haversine(): the
    matrix is n**2 and n is a few thousand. _check_distance_metric() asserts the two
    agree, so the figure provably measures what the merge measured.

    Parameters
    ----------
    latitudes  : array-like of float — degrees
    longitudes : array-like of float — degrees

    Returns
    -------
    numpy.ndarray of shape (n, n) — metres
    """
    lat = np.radians(np.asarray(latitudes, dtype=float))
    lon = np.radians(np.asarray(longitudes, dtype=float))
    xyz = np.stack([np.cos(lat) * np.cos(lon),
                    np.cos(lat) * np.sin(lon),
                    np.sin(lat)], axis=1)

    cosine   = np.clip(xyz @ xyz.T, -1.0, 1.0)
    distance = 2 * _EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.maximum(0.0, (1 - cosine) / 2)))
    distance *= 1000
    np.fill_diagonal(distance, np.inf)
    return distance


def _check_distance_metric(stations, distances, n_samples=200, tolerance_m=1e-2):
    """
    Assert the vectorized matrix reproduces the merge's own haversine().

    The tolerance is a centimetre: the two formulations differ only by float64
    round-off (tens of nanometres at a few hundred metres), and anything that could
    change a decision at a 20 m radius is three orders of magnitude larger.
    """
    rng = np.random.default_rng(0)
    n   = len(stations)
    for i, j in zip(rng.integers(0, n, n_samples), rng.integers(0, n, n_samples)):
        if i == j:
            continue
        expected = haversine(stations[i]['lat'], stations[i]['lon'],
                             stations[j]['lat'], stations[j]['lon']) * 1000
        if abs(expected - distances[i, j]) > tolerance_m:
            raise AssertionError(
                f"distance mismatch for {stations[i]['key']} / {stations[j]['key']}: "
                f"{distances[i, j]} m vectorized vs {expected} m from haversine()"
            )


# ---------------------------------------------------------------------------
# Inventory loading
# ---------------------------------------------------------------------------

def _base_code(station_code):
    """Strip the _N suffix the merge appends to distinct clusters under one code."""
    return station_code.rsplit('_', 1)[0] if '_' in station_code else station_code


def load_source_stations(source_glob):
    """
    Load every station of the source inventories, before any merging.

    Parameters
    ----------
    source_glob : str — glob pattern for the input StationXML files

    Returns
    -------
    list[dict] with keys: network, code, key, lat, lon
    """
    inventory = Inventory()
    for path in sorted(glob.glob(source_glob)):
        inventory.extend(read_inventory(path, format='STATIONXML'))

    return [{'network': network.code,
             'code':    station.code,
             'key':     f'{network.code}.{station.code}',
             'lat':     station.latitude,
             'lon':     station.longitude}
            for network in inventory.networks for station in network.stations]


def load_merge_outcome(file_inventory):
    """
    Read what the merge produced, indexed the way a source station can be found again.

    A surviving station keeps its network, its code up to the cluster suffix, and its
    position, so those identify it on both sides. Stations absorbed by a same-code
    duplicate keep no entry of their own, and are recovered through the surviving
    station of the same network and code.

    Parameters
    ----------
    file_inventory : str — path to GLOBAL_inventory.xml

    Returns
    -------
    dict[(network, base_code), list[(lat, lon, alternate_code)]]
    """
    inventory = read_inventory(file_inventory)
    surviving = {}
    for network in inventory:
        for station in network:
            key = (network.code, _base_code(station.code))
            surviving.setdefault(key, []).append(
                (station.latitude, station.longitude, station.alternate_code))
    return surviving


def resolve_groups(stations, file_inventory):
    """
    Assign every source station the alternate code it ends under, or None if it is gone.

    An absorbed duplicate is given the code of the station that absorbed it -- the
    nearest survivor of the same network and code, which is the cluster it was folded
    into. Without that, a station and the one that swallowed it would read as two
    separate groups and the merge would look like it never happened.

    Returns
    -------
    (list[str or None], dict[str, int])
        Per-station alternate code, and the tallies reported by station_colocation().
    """
    surviving = load_merge_outcome(file_inventory)

    groups  = []
    tallies = {'n_resolved': 0, 'n_removed': 0}
    for station in stations:
        candidates = surviving.get((station['network'], _base_code(station['code'])))
        if not candidates:
            # nothing of that network and code survived: dropped by check_inventory()
            groups.append(None)
            tallies['n_removed'] += 1
            continue

        _, _, alternate_code = min(
            candidates,
            key=lambda c: haversine(station['lat'], station['lon'], c[0], c[1]))
        groups.append(alternate_code)
        tallies['n_resolved'] += 1
    return groups, tallies


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def _panel_nearest_neighbour(axis, nearest_m, merged_with_nearest, threshold):
    """Nearest-neighbour distance per station, split by whether that pair merged."""
    finite = np.isfinite(nearest_m) & (nearest_m > 0)
    edges  = np.logspace(np.log10(max(nearest_m[finite].min(), 1e-2)),
                         np.log10(nearest_m[finite].max()), 60)

    axis.hist([nearest_m[finite & merged_with_nearest],
               nearest_m[finite & ~merged_with_nearest]],
              bins=edges, stacked=True, color=['#2c7fb8', '#d9d9d9'],
              edgecolor='white', linewidth=0.2,
              label=['merged with it', 'not merged'])

    axis.axvline(threshold, color='#d95f02', linestyle='--', linewidth=1.2)
    axis.text(threshold, axis.get_ylim()[1] * 0.95, f'  {threshold} m',
              color='#d95f02', va='top', fontsize=8)
    axis.set_xscale('log')
    axis.set_xlabel('distance to nearest station (m)')
    axis.set_ylabel('stations')
    axis.set_title('Nearest neighbour, and whether that pair was merged', fontsize=10)
    axis.legend(fontsize=8, frameon=False)


def _panel_sensitivity(axis, pair_distances, n_merged_pairs, threshold):
    """Pairs that would merge as a function of the radius."""
    radii  = np.arange(0, _MAX_THRESHOLD_M + 1, 1.0)
    counts = np.searchsorted(np.sort(pair_distances), radii, side='right')

    axis.plot(radii, counts, color='#2c7fb8', linewidth=1.5)
    axis.axvline(threshold, color='#d95f02', linestyle='--', linewidth=1.2)
    axis.axhline(n_merged_pairs, color='#7f7f7f', linestyle=':', linewidth=1.0)

    at_threshold = int(np.searchsorted(np.sort(pair_distances), threshold, side='right'))
    axis.plot([threshold], [at_threshold], 'o', color='#d95f02', markersize=5)
    axis.annotate(f'{at_threshold} pairs at {threshold} m',
                  xy=(threshold, at_threshold), xytext=(threshold + 12, at_threshold * 0.82),
                  fontsize=8, color='#d95f02')
    axis.annotate(f'{n_merged_pairs} actually merged', xy=(_MAX_THRESHOLD_M, n_merged_pairs),
                  xytext=(_MAX_THRESHOLD_M * 0.55, n_merged_pairs * 1.04),
                  fontsize=8, color='#7f7f7f')

    axis.set_xlabel('co-location radius (m)')
    axis.set_ylabel('pairs within the radius')
    axis.set_title('Cost of moving the threshold', fontsize=10)


def _panel_networks(axis, stations, nearest_m, lost_identity, max_networks=46):
    """Per-network station counts, merged share, and median nearest-neighbour distance."""
    per_network = {}
    for index, station in enumerate(stations):
        entry = per_network.setdefault(station['network'], {'total': 0, 'lost': 0, 'nn': []})
        entry['total'] += 1
        entry['lost']  += int(lost_identity[index])
        if np.isfinite(nearest_m[index]):
            entry['nn'].append(nearest_m[index])

    order = sorted(per_network.items(), key=lambda item: item[1]['total'])[-max_networks:]
    names = [name for name, _ in order]
    total = np.array([entry['total'] for _, entry in order], dtype=float)
    lost  = np.array([entry['lost'] for _, entry in order], dtype=float)

    positions = np.arange(len(names))
    axis.barh(positions, total, color='#d9d9d9', height=0.75, label='stations')
    axis.barh(positions, lost, color='#2c7fb8', height=0.75, label='merged away')

    for position, (_, entry) in zip(positions, order):
        median = np.median(entry['nn']) if entry['nn'] else float('nan')
        label  = '—' if math.isnan(median) else (f'{median:.0f} m' if median < 1000
                                                 else f'{median / 1000:.1f} km')
        axis.text(entry['total'] + max(total) * 0.01, position, label,
                  va='center', fontsize=6, color='#555555')

    axis.set_yticks(positions)
    axis.set_yticklabels(names, fontsize=6)
    axis.set_xlabel('stations  (label: median nearest-neighbour distance)')
    axis.set_title('Per network', fontsize=10)
    axis.legend(fontsize=8, frameon=False, loc='lower right')


def build_figure(stations, distances, groups, parameters):
    """
    Draw the three panels and save the PDF.

    Returns
    -------
    dict of the counts annotated on the figure
    """
    n         = len(stations)
    nearest_m = distances.min(axis=1)
    nearest_i = distances.argmin(axis=1)

    def merged(i, j):
        return groups[i] is not None and groups[i] == groups[j]

    merged_with_nearest = np.array([merged(i, nearest_i[i]) for i in range(n)])

    members = {}
    for index, label in enumerate(groups):
        if label is not None:
            members.setdefault(label, []).append(index)
    lost_identity = np.array([label is None or len(members[label]) > 1 for label in groups])

    upper          = np.triu_indices(n, 1)
    pair_distances = distances[upper]

    # Merged pairs, split on the radius: the ones beyond it are the chaining cases,
    # where A-B and B-C are both inside the radius but A-C is not.
    merged_within = merged_beyond = 0
    for indices in members.values():
        for a in range(len(indices)):
            for b in range(a + 1, len(indices)):
                if distances[indices[a], indices[b]] <= parameters.threshold:
                    merged_within += 1
                else:
                    merged_beyond += 1
    n_merged_pairs = merged_within + merged_beyond

    figure = plt.figure(figsize=(11, 12))
    grid   = figure.add_gridspec(2, 2, height_ratios=[1, 1.9], hspace=0.28, wspace=0.22)

    _panel_nearest_neighbour(figure.add_subplot(grid[0, 0]), nearest_m,
                             merged_with_nearest, parameters.threshold)
    _panel_sensitivity(figure.add_subplot(grid[0, 1]), pair_distances,
                       n_merged_pairs, parameters.threshold)
    _panel_networks(figure.add_subplot(grid[1, :]), stations, nearest_m, lost_identity)

    figure.suptitle(
        f'Station co-location — {n} source stations, {parameters.threshold} m radius',
        fontsize=12, y=0.98)

    os.makedirs(os.path.dirname(parameters.file_output), exist_ok=True)
    figure.savefig(parameters.file_output, bbox_inches='tight')
    plt.close(figure)

    n_pairs_within = int((pair_distances <= parameters.threshold).sum())
    return {
        'n_stations':        n,
        'n_pairs_within':    n_pairs_within,
        'n_merged_pairs':    n_merged_pairs,
        'n_merged_beyond':   merged_beyond,
        'n_unmerged_within': n_pairs_within - merged_within,
        'n_lost_identity':   int(lost_identity.sum()),
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def station_colocation(parameters):
    """
    Produce the co-location diagnostic figure.

    Parameters
    ----------
    parameters : StationColocationParams

    Returns
    -------
    dict with key 'output' plus the counts reported on the figure
    """
    stations = load_source_stations(parameters.source_glob)
    print(f"Source stations read: {len(stations)} @ {parameters.source_glob}")

    distances = pairwise_distances([s['lat'] for s in stations],
                                   [s['lon'] for s in stations])
    _check_distance_metric(stations, distances)

    groups, tallies = resolve_groups(stations, parameters.file_inventory)
    print(f"Merge outcome: {tallies['n_resolved']} resolved to a surviving station, "
          f"{tallies['n_removed']} removed outright")

    counts = build_figure(stations, distances, groups, parameters)
    print(f"Pairs within {parameters.threshold} m: {counts['n_pairs_within']}; "
          f"pairs actually merged: {counts['n_merged_pairs']}")
    print(f"Merged pairs beyond the radius (chaining): {counts['n_merged_beyond']}; "
          f"pairs inside it left unmerged: {counts['n_unmerged_within']}")
    print(f"Stations that lost their standalone identity: {counts['n_lost_identity']}")
    print(f"Figure saved: {parameters.file_output}")

    return {'output': parameters.file_output, **tallies, **counts}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Diagnostic for the co-location radius of the station inventory fusion.'
    )
    parser.add_argument('--source-glob',
                        default=os.path.join(_PROJECT_ROOT, 'stations', '*', '*.xml'),
                        help='Glob for the source StationXML files (default: stations/*/*.xml)')
    parser.add_argument('--inventory',
                        default=os.path.join(_PROJECT_ROOT, 'stations', 'GLOBAL_inventory.xml'),
                        help='Merged inventory (default: stations/GLOBAL_inventory.xml)')
    parser.add_argument('--output',
                        default=os.path.join(_MODULE_DIR, 'station_colocation',
                                             'station_colocation.pdf'),
                        help='Output PDF path')
    parser.add_argument('--threshold', type=int, default=20,
                        help='Co-location radius in metres to mark (default: 20)')
    args = parser.parse_args()

    station_colocation(StationColocationParams(
        source_glob    = args.source_glob,
        file_inventory = args.inventory,
        file_output    = args.output,
        threshold      = args.threshold,
    ))


if __name__ == '__main__':
    main()
