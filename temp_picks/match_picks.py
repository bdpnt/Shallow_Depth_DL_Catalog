"""
match_picks.py
============================
Associate converted pick lines with events in a GLOBAL.obs bulletin.

For each pick in the input file, this script:
  1. Finds all bulletin events whose origin time falls within 60 s before the pick.
  2. For each candidate, computes the epicentral distance and checks whether the
     observed travel time is consistent with the theoretical bounds from the
     travel-time table (±0.1 s for P, ±0.3 s for S), plus ±2.5 s to absorb
     possible origin-time (t0) error.
  3. If exactly one event passes, checks for duplicates and — if new — appends
     the pick line to that event and updates its PhaseCount.

The result is written to a new file; the input bulletin is never modified.

Usage
-----
    python temp_picks/match_picks.py \\
        --picks      temp_picks/pick_files/viehla_final_converted.obs \\
        --bulletin   obs/FINAL.obs \\
        --inventory  stations/GLOBAL_inventory.xml \\
        --tables     temp_picks/tables_Pyr.csv \\
        --output     obs/FINAL_augmented.obs

Adding support for a new input pick format
------------------------------------------
    The script only requires that pick lines follow the GLOBAL.obs column layout
    (output of convert_picks.py). No changes are needed here when a new source
    format is added to convert_picks.py.
"""

import argparse
import bisect
import logging
import math
import os
import sys
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
from obspy import read_inventory

# Ensure the project root is on sys.path so the import works both when the
# script is run directly and when it is imported as part of the package.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from temp_picks.sort_picks import sort_picks

logger = logging.getLogger('match_picks')

T0_TOL = 2.5   # seconds — extra window to absorb origin-time (t0) error


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class Event:
    header_line: str          # full original '# ...' line
    event_dt:    datetime     # parsed origin time (UTC, timezone-aware)
    lat:         float
    lon:         float
    picks:       list = field(default_factory=list)   # list of pick line strings
    pick_keys:   set  = field(default_factory=set)    # (station_code, phase) pairs for O(1) dup check
    public_id:   str  = None                          # PUBLIC_ID value (e.g. 'PYRENEES_000001')


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------

def load_bulletin(path):
    """
    Parse a GLOBAL.obs bulletin file.

    Returns
    -------
    file_header : list[str]
        The '###' comment lines at the top of the file (without trailing newlines).
    events : list[Event]
        All events in the file, in original order.
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    file_header = []
    events      = []
    i           = 0

    # Collect '###' header lines, stop at blank or event line
    while i < len(lines) and lines[i].startswith('###'):
        file_header.append(lines[i].rstrip('\n'))
        i += 1

    # Skip any blank lines between header and first event
    while i < len(lines) and lines[i].strip() == '':
        i += 1

    # Parse event blocks
    while i < len(lines):
        line = lines[i].rstrip('\n')
        if line.startswith('# '):
            event = _parse_event_header(line)
            i += 1
            while i < len(lines) and lines[i].strip() != '':
                pick = lines[i].rstrip('\n')
                if pick.startswith('PUBLIC_ID'):
                    parts = pick.split(None, 1)
                    event.public_id = parts[1] if len(parts) > 1 else None
                else:
                    event.picks.append(pick)
                    parts = pick.split()
                    if len(parts) >= 5:
                        event.pick_keys.add((parts[0].strip(), parts[4].strip()))
                i += 1
            events.append(event)
        i += 1  # advance past blank line (or unrecognised line)

    return file_header, events


def _parse_event_header(line):
    """Parse a '# ...' event header line into an Event object."""
    parts = line[2:].split()  # strip leading '# '
    year   = int(parts[0])
    month  = int(parts[1])
    day    = int(parts[2])
    hour   = int(parts[3])
    minute = int(parts[4])
    second = float(parts[5])
    lat    = float(parts[6])
    lon    = float(parts[7])

    sec_int  = int(second)
    microsec = int(round((second - sec_int) * 1e6))
    event_dt = datetime(year, month, day, hour, minute, sec_int, microsec, tzinfo=timezone.utc)

    return Event(header_line=line, event_dt=event_dt, lat=lat, lon=lon)


def _update_phase_count(header_line, count):
    """Return a copy of the event header line with PhaseCount set to count."""
    parts       = header_line[2:].split()
    parts[12]   = str(count)
    return '# ' + ' '.join(parts)


def load_inventory(path):
    """
    Build a station position lookup from a StationXML file.

    Uses ObsPy's read_inventory(). Stations are keyed by their alternate_code
    (the project's internal code, e.g. 'FR.0035'). When multiple entries share
    the same alternate_code (always within 20 m of each other), the first is used.

    Returns
    -------
    dict[str, (float, float)]
        Mapping from internal station code to (latitude, longitude).
    """
    inv     = read_inventory(path)
    station_map = {}
    for network in inv:
        for station in network:
            code = station.alternate_code
            if code and code not in station_map:
                station_map[code] = (station.latitude, station.longitude)
    return station_map


def load_tables(path):
    """Load the theoretical travel-time table as a DataFrame."""
    return pd.read_csv(path)


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def haversine_km(lat1, lon1, lat2, lon2):
    """Return the epicentral distance in km between two lat/lon points."""
    R    = 6371.0
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a    = math.sin(dphi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlam / 2) ** 2
    return 2.0 * R * math.asin(math.sqrt(a))


# ---------------------------------------------------------------------------
# Travel-time interpolation
# ---------------------------------------------------------------------------

def interpolate_travel_time(tables, dist_km, phase):
    """
    Linearly interpolate theoretical travel-time bounds for a given distance.

    Parameters
    ----------
    tables   : pd.DataFrame  — travel-time table (columns: distance, tp_low, tp_high, ts_low, ts_high)
    dist_km  : float         — epicentral distance in km
    phase    : str           — 'P' or 'S'

    Returns
    -------
    (t_low, t_high) in seconds, or None if dist_km exceeds the table range.
    """
    distances = tables['distance'].values
    if dist_km > distances[-1]:
        return None

    if phase == 'P':
        t_low  = float(np.interp(dist_km, distances, tables['tp_low'].values))
        t_high = float(np.interp(dist_km, distances, tables['tp_high'].values))
    else:
        t_low  = float(np.interp(dist_km, distances, tables['ts_low'].values))
        t_high = float(np.interp(dist_km, distances, tables['ts_high'].values))

    return t_low, t_high


# ---------------------------------------------------------------------------
# Pick line parsing
# ---------------------------------------------------------------------------

def parse_pick_line(line):
    """
    Parse a GLOBAL.obs-format pick line.

    Returns
    -------
    (station_code, phase, arrival_dt) or raises ValueError on malformed lines.
    """
    parts = line.split()
    if len(parts) < 9:
        raise ValueError(f"Too few fields: {line!r}")

    station_code = parts[0].strip()
    phase        = parts[4].strip()

    date_str = parts[6]   # YYYYMMDD
    hhmm_str = parts[7]   # HHMM
    sec_str  = parts[8]   # SS.SSS

    year   = int(date_str[:4])
    month  = int(date_str[4:6])
    day    = int(date_str[6:8])
    hour   = int(hhmm_str[:2])
    minute = int(hhmm_str[2:4])
    second = float(sec_str)

    sec_int  = int(second)
    microsec = int(round((second - sec_int) * 1e6))
    arrival_dt = datetime(year, month, day, hour, minute, sec_int, microsec, tzinfo=timezone.utc)

    return station_code, phase, arrival_dt


# ---------------------------------------------------------------------------
# Logging setup (mirrors convert_picks.py)
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')


def _setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    basename  = os.path.splitext(os.path.basename(__file__))[0]
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_path  = os.path.join(log_dir, f"{basename}_{timestamp}.log")

    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    handler = logging.FileHandler(log_path, encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(handler)

    return log_path


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def match_picks(pick_file, bulletin_file, inventory_file, tables_file,
                output_file=None, log_dir=None):
    """
    Match converted picks to bulletin events and write an augmented bulletin.

    Parameters
    ----------
    pick_file      : str  — path to the converted pick file (GLOBAL.obs pick line format)
    bulletin_file  : str  — path to the event bulletin (obs/GLOBAL.obs)
    inventory_file : str  — path to the StationXML inventory
    tables_file    : str  — path to the theoretical travel-time CSV
    output_file    : str, optional — output path; default: bulletin stem + '_augmented.obs'
    log_dir        : str, optional — log directory; default: temp_picks/console_output/

    Returns
    -------
    dict with keys: output, log, n_picks, n_added, n_skipped_no_event,
                    n_skipped_no_station, n_skipped_no_residual,
                    n_skipped_multi, n_skipped_duplicate
    """
    if output_file is None:
        base, _ = os.path.splitext(bulletin_file)
        output_file = base + '_augmented.obs'

    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)

    # --- Load data ---
    logger.info(f"Loading bulletin     : {bulletin_file}")
    file_header, events = load_bulletin(bulletin_file)
    logger.info(f"  {len(events)} events loaded.")

    logger.info(f"Loading inventory    : {inventory_file}")
    inventory = load_inventory(inventory_file)
    logger.info(f"  {len(inventory)} stations loaded.")

    logger.info(f"Loading travel-time tables: {tables_file}")
    tables = load_tables(tables_file)

    # Sort events by time; build a parallel list of timestamps for bisect
    events.sort(key=lambda e: e.event_dt)
    event_times = [e.event_dt for e in events]

    # --- Process picks ---
    n_picks             = 0
    n_added             = 0
    n_added_near        = 0
    n_skipped_no_event  = 0
    n_skipped_no_sta    = 0
    n_skipped_no_res    = 0
    n_skipped_multi     = 0
    n_skipped_dup       = 0

    with open(pick_file, 'r') as f:
        for raw_line in f:
            line = raw_line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue

            n_picks += 1

            # Parse pick
            try:
                station_code, phase, arrival_dt = parse_pick_line(line)
            except ValueError as exc:
                logger.warning(f"Skipping malformed pick line ({exc}): {line!r}")
                continue

            # 1. Time window: events within [arrival - 60s, arrival]
            window_start = arrival_dt - timedelta(seconds=60)
            idx_lo = bisect.bisect_left(event_times, window_start)
            idx_hi = bisect.bisect_right(event_times, arrival_dt)
            candidates = events[idx_lo:idx_hi]

            if not candidates:
                n_skipped_no_event += 1
                continue

            # 2. Station position
            pos = inventory.get(station_code)
            if pos is None:
                n_skipped_no_sta += 1
                logger.warning(f"Station '{station_code}' not found in inventory. Skipping.")
                continue
            sta_lat, sta_lon = pos

            # 3. Residual filter
            tol     = 0.1 if phase == 'P' else 0.3
            passing = []
            for event in candidates:
                dist_km = haversine_km(event.lat, event.lon, sta_lat, sta_lon)
                tt = interpolate_travel_time(tables, dist_km, phase)
                if tt is None:
                    continue  # station too far from event for this table
                t_low, t_high = tt
                obs_tt = (arrival_dt - event.event_dt).total_seconds()
                if t_low - tol - T0_TOL <= obs_tt <= t_high + tol + T0_TOL:
                    passing.append((event, dist_km))

            if len(passing) == 0:
                n_skipped_no_res += 1
                continue

            if len(passing) > 1:
                n_skipped_multi += 1
                logger.warning(
                    f"Multiple events ({len(passing)}) match pick "
                    f"{station_code} {phase} @ {arrival_dt.isoformat()}. Skipping."
                )
                continue

            # 4. Single match — check for duplicate
            matched, dist_km_matched = passing[0]
            if (station_code, phase) in matched.pick_keys:
                n_skipped_dup += 1
                continue

            # 5. Add pick and update PhaseCount
            matched.picks.append(line)
            matched.pick_keys.add((station_code, phase))
            matched.header_line = _update_phase_count(matched.header_line, len(matched.picks))
            n_added += 1
            if dist_km_matched <= 20.0:
                n_added_near += 1

    # --- Write output ---
    with open(output_file, 'w') as f:
        for h in file_header:
            f.write(h + '\n')
        f.write('\n')
        for event in events:
            f.write(event.header_line + '\n')
            if event.public_id is not None:
                f.write(f'PUBLIC_ID {event.public_id}\n')
            for pick in event.picks:
                f.write(pick + '\n')
            f.write('\n')

    logger.info(f"Input picks                                    : {n_picks}")
    logger.info(f"Added                                          : {n_added}")
    logger.info(f"Added - station within 20 km of event          : {n_added_near}")
    logger.info(f"Skipped - no bulletin event in 60 s window     : {n_skipped_no_event}")
    logger.info(f"Skipped - station not in inventory             : {n_skipped_no_sta}")
    logger.info(f"Skipped - travel time outside theoretical band : {n_skipped_no_res}")
    logger.info(f"Skipped - ambiguous (multiple events matched)  : {n_skipped_multi}")
    logger.info(f"Skipped - duplicate station+phase in event     : {n_skipped_dup}")
    logger.info(f"Output               : {output_file}")

    # Sort picks by arrival time within each event
    logger.info("Sorting picks by arrival time ...")
    sort_result = sort_picks(output_file, output_file)
    logger.info(f"Sorted {sort_result['n_picks']} picks across {sort_result['n_events']} events.")

    return {
        'output':                output_file,
        'log':                   log_path,
        'n_picks':               n_picks,
        'n_added':               n_added,
        'n_added_near':          n_added_near,
        'n_skipped_no_event':    n_skipped_no_event,
        'n_skipped_no_station':  n_skipped_no_sta,
        'n_skipped_no_residual': n_skipped_no_res,
        'n_skipped_multi':       n_skipped_multi,
        'n_skipped_duplicate':   n_skipped_dup,
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Match converted picks to GLOBAL.obs events and write an augmented bulletin.'
    )
    parser.add_argument('--picks',     required=True, help='Converted pick file (GLOBAL.obs pick line format)')
    parser.add_argument('--bulletin',  required=True, help='Input bulletin file (obs/GLOBAL.obs)')
    parser.add_argument('--inventory', required=True, help='StationXML inventory file')
    parser.add_argument('--tables',    required=True, help='Theoretical travel-time CSV')
    parser.add_argument('--output',    default=None,  help='Output bulletin path (default: bulletin stem + _augmented.obs)')
    parser.add_argument('--log-dir',   default=None,  help='Log directory (default: temp_picks/console_output/)')
    args = parser.parse_args()

    match_picks(
        pick_file      = args.picks,
        bulletin_file  = args.bulletin,
        inventory_file = args.inventory,
        tables_file    = args.tables,
        output_file    = args.output,
        log_dir        = args.log_dir,
    )


if __name__ == '__main__':
    main()
