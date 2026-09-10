"""
merge_station_inventories.py
============================
Merge all station XML inventories into a single unified inventory with
unique alternate station codes.

Networks are deduplicated, stations that are within accepted_distance metres
of each other receive the same alternate code, and an alternate-code mapping
file is written alongside the unified StationXML.

The code a group of co-located stations keeps is the one of its best-ranked
member under station_priority() -- permanent network first, then RESIF/RENASS,
then a station carrying a real elevation, the XX placeholder last -- with the
oldest start_date breaking ties. Stations left without an elevation inherit one
from their group, then from the Open-Elevation API.

Usage
-----
    python fetch_inventory/merge_station_inventories.py \\
        --folder-path   "stations/*/*.xml" \\
        --save-inventory stations/GLOBAL_inventory.xml \\
        --save-mapping   stations/GLOBAL_code_map.txt \\
        --distance       20
"""

import argparse
import datetime
import logging
import math
import os
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime as dt

from obspy import UTCDateTime, read_inventory
from obspy.core.inventory import Inventory

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)

if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)

from fetch_inventory._fill_missing_elevations import _get_elevation  # noqa: E402


logger = logging.getLogger('fetch_inventory')

_DEFAULT_LOG_DIR = 'fetch_inventory/console_output/'

# Networks operating permanently, from the `type` column of stations/all_networks.xlsx.
# AM (RaspberryShake) is semi-permanent and counted here.
PERMANENT_NETWORKS = frozenset({'FR', 'RA', 'RD', 'ES', 'CA', 'LC', 'AM'})
# RESIF/RENASS networks, preferred among the permanent ones.
RESIF_NETWORKS     = frozenset({'FR', 'RA', 'RD'})
# Placeholder network for uncalled or unknown networks; a real code always wins over it.
UNKNOWN_NETWORK    = 'XX'


def _setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    basename  = os.path.splitext(os.path.basename(__file__))[0]
    timestamp = dt.now().strftime('%Y%m%d_%H%M%S')
    log_path  = os.path.join(log_dir, f"{basename}_{timestamp}.log")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(log_path, encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(handler)
    return log_path


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MergeInventoryParams:
    """
    Configuration for merging station inventories.

    Attributes
    ----------
    folder_path        : str — glob pattern for input StationXML files
    file_save_inventory: str — output path for the unified StationXML
    file_save_mapping  : str — output path for the alternate-code mapping text file
    accepted_distance  : int — max distance (m) for two stations to be considered
                               the same physical site (default: 20)
    fill_elevations    : bool — query Open-Elevation for the stations still left
                               without an elevation (default: True). The only step of
                               the pipeline that reaches the network; disable it to
                               run offline.
    """
    folder_path:         str
    file_save_inventory: str
    file_save_mapping:   str
    accepted_distance:   int = 20
    fill_elevations:     bool = True


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def haversine(lat1, lon1, lat2, lon2):
    """Return the great-circle distance in km between two lat/lon points."""
    R    = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a    = (math.sin(dlat / 2) ** 2
            + math.cos(math.radians(lat1))
            * math.cos(math.radians(lat2))
            * math.sin(dlon / 2) ** 2)
    return 2 * R * math.asin(math.sqrt(a))


# ---------------------------------------------------------------------------
# Merge priority
# ---------------------------------------------------------------------------

def station_priority(net_code, elevation, alternate_code=None):
    """
    Rank a station as the representative of a group of co-located stations.

    Lower sorts first: a permanent network wins over a temporary one, RESIF/RENASS
    wins among the permanent ones, a station carrying a real elevation wins over one
    left at 0, and the XX placeholder loses last. Callers append their own final
    tie-break term (in practice the oldest start_date).

    `alternate_code` adds a last term preferring the member whose network matches the
    code's own prefix, and should be passed by everything that re-derives the
    representative *after* the merge. The merge itself must not pass it: it is what
    assigns the codes, and the term would be circular.

    The term exists because the in-group elevation fill erases the very criterion that
    chose the representative — once a member has inherited the winner's elevation, the
    third term goes neutral and the date tie-break can hand the group to someone else.
    The code's prefix is the durable record of who won, so re-deriving from it
    reproduces the merge's own decision instead of guessing at it again.

    Parameters
    ----------
    net_code       : str            — network code of the station
    elevation      : float or None  — station elevation in metres; 0 and None both
                                      count as missing, as everywhere else here
    alternate_code : str, optional  — the group's code, e.g. 'AM.0043'

    Returns
    -------
    tuple[bool, bool, bool, bool, bool]
    """
    return (net_code not in PERMANENT_NETWORKS,
            net_code not in RESIF_NETWORKS,
            not elevation,
            net_code == UNKNOWN_NETWORK,
            alternate_code is not None and net_code != alternate_code.split('.')[0])


# ---------------------------------------------------------------------------
# Inventory helpers
# ---------------------------------------------------------------------------

def check_inventory(inventory):
    """
    Interactively prompt the user to remove duplicate stations that appear
    in multiple networks.

    Parameters
    ----------
    inventory : obspy.core.inventory.Inventory

    Returns
    -------
    obspy.core.inventory.Inventory
        Inventory with user-selected duplicates removed.
    """
    unique_sta = defaultdict(list)
    for net in inventory.networks:
        for sta in net.stations:
            unique_sta[sta.code.split('_')[0]].append(net.code)

    unique_sta = {
        sta: nets for sta, nets in unique_sta.items()
        if len(nets) >= 2 and any(net != nets[0] for net in nets)
    }

    for station, networks in unique_sta.items():
        print('\n' * 40)
        print(f"\nStation: {station}")
        for network in networks:
            current_inv = inventory.select(network=network, station=station)
            print(f"  In Network: {network}")
            for net in current_inv.networks:
                for sta in net.stations:
                    print(sta)
        print("____________________")

        remove_nets = input('Network(s) to remove, if any, separated by commas (e.g. > FR,RD): ')
        remove_nets = remove_nets.split(',')
        for network in remove_nets:
            inventory = inventory.remove(network=network, station=station)

    return inventory


def _add_alternate_code(inventory):
    """Assign a unique alternate code (NET.XXXX) to every station."""
    for network in inventory.networks:
        net_code = network.code
        sta_id   = 0
        for station in network:
            station.alternate_code = f'{net_code}.{str(sta_id).zfill(4)}'
            sta_id += 1
    return inventory


def _combine_close_stations(inventory, parameters):
    """
    Give the same alternate code to stations from different networks that
    are within accepted_distance metres of each other.

    The code kept is the group's best-ranked station under station_priority(),
    oldest start_date breaking ties. Members left without an elevation then inherit
    the representative's — the group spans at most accepted_distance metres, so the
    value is right to within its own vertical scatter.

    The order matters: back-filling before choosing would erase the elevation term of
    station_priority() and hand some groups to the station that had no elevation.
    """
    all_stations   = [(net.code, sta) for net in inventory.networks for sta in net.stations]
    threshold_km   = parameters.accepted_distance / 1000
    groups         = []

    for net_code, station in all_stations:
        target_group = next(
            (g for g in groups
             if any(haversine(station.latitude, station.longitude,
                              other.latitude, other.longitude) <= threshold_km
                    for _, other in g)),
            None,
        )
        if target_group is not None:
            target_group.append((net_code, station))
        else:
            groups.append([(net_code, station)])

    n_filled = 0
    for group in groups:
        if len(group) > 1:
            best = min(group, key=lambda x: (station_priority(x[0], x[1].elevation),
                                             x[1].start_date or UTCDateTime(datetime.datetime.max)))
            for _, station in group:
                if station is not best[1]:
                    station.alternate_code = best[1].alternate_code
                    if not station.elevation and best[1].elevation:
                        station.elevation = best[1].elevation
                        n_filled += 1

    logger.info(f"Co-located grouping: {len(groups)} group(s), "
                f"{n_filled} elevation(s) inherited within a group")
    return inventory


def _fill_station_elevations(inventory):
    """
    Query Open-Elevation for the stations still left without an elevation.

    Whole nodal deployments are published with the field never populated, and a
    station's elevation enters the NonLinLoc travel-time computation directly through
    its GTSRCE line, so an unfilled 0 is a travel-time error rather than a cosmetic
    gap. Runs after the co-located groups have shared what they had between them, so
    only the stations no neighbour could serve are queried.

    A lookup that fails returns 0 and leaves the station exactly as it was.

    Parameters
    ----------
    inventory : obspy.core.inventory.Inventory

    Returns
    -------
    obspy.core.inventory.Inventory
    """
    missing = [station for network in inventory.networks
               for station in network.stations if not station.elevation]

    logger.info(f"Stations still without an elevation: {len(missing)}")
    if not missing:
        return inventory

    n_filled = 0
    for station in missing:
        elevation = _get_elevation(station.latitude, station.longitude)
        if elevation:
            station.elevation = elevation
            n_filled += 1
        else:
            logger.warning(f"  Elevation lookup failed for {station.code} "
                           f"({station.latitude}, {station.longitude})")
        time.sleep(1)  # avoid rate limiting

    logger.info(f"Elevation fill: {n_filled}/{len(missing)} filled, "
                f"{len(missing) - n_filled} still without one")
    return inventory


def _create_alternate_code_mapping(inventory, parameters):
    """
    Write the alternate-code → station-code mapping file.

    Members of a block are written best-ranked first (station_priority(), oldest
    start_date breaking ties), so the first entry is the station whose code labels
    the block. Consumers parse whole blocks or match on the epochs, so this is for
    readability rather than for them.
    """
    mapping = defaultdict(list)

    for network in inventory.networks:
        for station in network.stations:
            mapping[station.alternate_code].append({
                'station_code': f'{network.code}.{station.code}',
                'start_date':   station.start_date,
                'end_date':     station.end_date,
                'sort_key':     (station_priority(network.code, station.elevation,
                                                  station.alternate_code),
                                 station.start_date or UTCDateTime(datetime.datetime.max)),
            })

    with open(parameters.file_save_mapping, 'w') as f:
        for alt_code, stations in mapping.items():
            f.write(f'Alternate Code: {alt_code}\n')
            for sta in sorted(stations, key=lambda s: s['sort_key']):
                f.write(f"  Station Code: {sta['station_code']}\n")
                f.write(f"  Start Date: {sta['start_date']}\n")
                f.write(f"  End Date: {sta['end_date']}\n")
            f.write('\n')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def merge_inventory(parameters, log_dir=None):
    """
    Merge all station XML inventories into a single unified inventory.

    Steps: load all XMLs → merge duplicate networks → merge duplicate
    stations within each network → interactive duplicate check →
    assign unique alternate codes → combine co-located stations →
    fill the elevations still missing → write unified StationXML and mapping file.

    Parameters
    ----------
    parameters : MergeInventoryParams
    log_dir    : str, optional — log directory; default: fetch_inventory/console_output/

    Returns
    -------
    dict with keys: output (inventory path), mapping (mapping path)
    """
    import glob

    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file             : {log_path}")
    logger.info(f"Input glob           : {parameters.folder_path}")
    logger.info(f"Co-location threshold: {parameters.accepted_distance} m")

    inventory    = Inventory()
    loaded_files = glob.glob(parameters.folder_path)
    for folder_file in loaded_files:
        inventory.extend(read_inventory(folder_file, format='STATIONXML'))

    logger.info(f"Inventory files loaded: {len(loaded_files)}")
    logger.info(f"Networks before merge : {len(inventory.networks)}")

    # --- Merge duplicate networks ---
    all_networks       = [net.code for net in inventory.networks]
    duplicate_networks = [code for code, count in Counter(all_networks).items() if count > 1]

    if duplicate_networks:
        logger.info(f"Duplicate networks to merge: {duplicate_networks}")
    for net_code in duplicate_networks:
        found = [net for net in inventory.networks if net.code == net_code]
        main  = found[0]
        n_merged = 0
        for sub in found[1:]:
            for station in sub:
                main.stations.append(station)
                n_merged += 1
            inventory.networks.remove(sub)
        logger.info(f"  Network {net_code}: merged {n_merged} station(s) from {len(found) - 1} duplicate(s)")

    # --- Merge duplicate stations within each network ---
    total_stations_dropped = 0
    for network in inventory.networks:
        station_groups = defaultdict(list)
        for sta in network.stations:
            station_groups[sta.code].append(sta)

        for _, same_code_stations in station_groups.items():
            if len(same_code_stations) <= 1:
                continue

            same_code_stations.sort(
                key=lambda s: getattr(s, 'start_date', None) or UTCDateTime(datetime.datetime.max)
            )

            n     = len(same_code_stations)
            graph = defaultdict(list)
            for i in range(n):
                for j in range(i + 1, n):
                    s1, s2 = same_code_stations[i], same_code_stations[j]
                    if haversine(s1.latitude, s1.longitude,
                                 s2.latitude, s2.longitude) <= parameters.accepted_distance / 1000:
                        graph[i].append(j)
                        graph[j].append(i)

            visited  = [False] * n
            clusters = []
            for i in range(n):
                if not visited[i]:
                    cluster_idx = []
                    stack = [i]
                    visited[i] = True
                    while stack:
                        node = stack.pop()
                        cluster_idx.append(node)
                        for nb in graph.get(node, []):
                            if not visited[nb]:
                                visited[nb] = True
                                stack.append(nb)
                    clusters.append([same_code_stations[k] for k in cluster_idx])

            for it, cluster in enumerate(clusters):
                if len(cluster) <= 1:
                    cluster[0].code = f'{cluster[0].code}_{it}'
                    continue

                main_sta   = min(cluster, key=lambda s: (
                    station_priority(network.code, s.elevation),
                    getattr(s, 'start_date', None) or UTCDateTime(datetime.datetime.max)))
                main_index = cluster.index(main_sta)

                if len(clusters) > 1:
                    main_sta.code = f'{main_sta.code}_{it}'

                to_remove = [s for i, s in enumerate(cluster) if i != main_index]

                logger.info(
                    f"  {network.code}.{main_sta.code}: merged cluster of {len(cluster)}, "
                    f"dropped {len(to_remove)} duplicate(s)"
                )

                for sta in to_remove:
                    if (main_sta.elevation == 0 or main_sta.elevation is None) and sta.elevation:
                        main_sta.elevation = sta.elevation
                    if hasattr(sta, 'start_date') and sta.start_date is not None:
                        current = getattr(main_sta, 'start_date', None)
                        if current is None or sta.start_date < current:
                            main_sta.start_date = sta.start_date
                    if hasattr(sta, 'end_date'):
                        current = getattr(main_sta, 'end_date', None)
                        if sta.end_date is None:
                            main_sta.end_date = None
                        elif current is not None and sta.end_date > current:
                            main_sta.end_date = sta.end_date
                    for channel in sta.channels:
                        if channel not in main_sta.channels:
                            main_sta.channels.append(channel)

                for sta in to_remove:
                    if sta in network.stations:
                        network.stations.remove(sta)
                total_stations_dropped += len(to_remove)

    logger.info(f"Station deduplication: {total_stations_dropped} station(s) dropped total")

    inventory = check_inventory(inventory)
    inventory = _add_alternate_code(inventory)
    inventory = _combine_close_stations(inventory, parameters)

    if parameters.fill_elevations:
        inventory = _fill_station_elevations(inventory)
    else:
        logger.info("Elevation fill skipped (fill_elevations=False)")

    n_total_stations = sum(len(net.stations) for net in inventory.networks)
    inventory.write(parameters.file_save_inventory, format='STATIONXML')
    logger.info(f"Inventory saved: {parameters.file_save_inventory}")
    logger.info(f"  Networks : {len(inventory.networks)}")
    logger.info(f"  Stations : {n_total_stations}")

    _create_alternate_code_mapping(inventory, parameters)
    logger.info(f"Alternate code mapping saved: {parameters.file_save_mapping}")

    return {
        'output':  parameters.file_save_inventory,
        'mapping': parameters.file_save_mapping,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Merge all station XML inventories into a single unified inventory.'
    )
    parser.add_argument('--folder-path',    required=True, help='Glob pattern for input StationXML files')
    parser.add_argument('--save-inventory', required=True, help='Output unified StationXML path')
    parser.add_argument('--save-mapping',   required=True, help='Output alternate-code mapping path')
    parser.add_argument('--distance',       type=int, default=20,
                        help='Max distance (m) for co-location (default: 20)')
    parser.add_argument('--no-fill-elevations', action='store_true',
                        help='Skip the Open-Elevation lookup for stations left without an elevation')
    args = parser.parse_args()

    params = MergeInventoryParams(
        folder_path         = args.folder_path,
        file_save_inventory = args.save_inventory,
        file_save_mapping   = args.save_mapping,
        accepted_distance   = args.distance,
        fill_elevations     = not args.no_fill_elevations,
    )
    merge_inventory(params)


if __name__ == '__main__':
    main()
