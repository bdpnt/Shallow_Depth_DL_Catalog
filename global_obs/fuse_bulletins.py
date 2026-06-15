"""
fuse_bulletins.py
============================
Merge all source .obs bulletins into a single GLOBAL.obs file by matching
events spatially and temporally, accumulating per-source statistics, and
producing comparison figures.

A second entry point, find_and_merge_doubles, lets the user interactively
resolve potential duplicate events in an existing bulletin.

Usage
-----
    python global_obs/fuse_bulletins.py \
        --global-bulletin-path obs/GLOBAL.obs \
        --main-bulletin-path   obs/LDG_20-25.obs \
        --folder-path          "obs/*.obs" \
        --dist-thresh          15 \
        --loose-dist-thresh    30 \
        --time-thresh          2 \
        --loose-time-thresh    10 \
        --mag-thresh           1.5 \
        --sim-pick-thresh      2

    # Interactive duplicate resolution
    python global_obs/fuse_bulletins.py --mode doubles \
        --global-bulletin-path obs/GLOBAL.obs \
        --max-dt-seconds 1 --max-dist-km 50
"""

import argparse
import glob
import logging
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime as dt
from numpy import mean

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


logger = logging.getLogger('global_obs.fuse_bulletins')

_DEFAULT_LOG_DIR = 'global_obs/console_output/'


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
# Path constants
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class FusionParams:
    """
    Configuration for fusing multiple .obs bulletins.

    Attributes
    ----------
    global_bulletin_path : str   — output path for the merged GLOBAL.obs
    main_bulletin_path   : str   — path to the reference (main) bulletin
    folder_path          : str   — glob pattern for all source bulletins
    dist_thresh          : float — strict distance threshold (km)
    loose_dist_thresh    : float — loose distance threshold (km)
    time_thresh          : float — strict time threshold (s)
    loose_time_thresh    : float — loose time threshold (s)
    mag_thresh           : float — magnitude difference threshold for ML–ML pairs
    """
    global_bulletin_path: str
    main_bulletin_path:   str
    folder_path:          str
    dist_thresh:          float
    loose_dist_thresh:    float
    time_thresh:          float
    loose_time_thresh:    float
    mag_thresh:           float


@dataclass
class MergeDoublesParams:
    """
    Configuration for interactive duplicate resolution.

    Attributes
    ----------
    global_bulletin_path : str   — path to the bulletin to review and clean
    max_dt_seconds       : float — maximum time gap (s) to flag a potential double
    max_dist_km          : float — maximum 3-D distance (km) to flag a potential double
    auto_dt_seconds      : float — Δt threshold (s) below which a pair of 2 is merged automatically
    auto_dist_km         : float — distance threshold (km) below which a pair of 2 is merged automatically
    """
    global_bulletin_path: str
    max_dt_seconds:       float
    max_dist_km:          float
    auto_dt_seconds:      float = 0.15
    auto_dist_km:         float = 10.0


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def haversine(lat1, lon1, lat2, lon2):
    """
    Return the great-circle distance in km between two geographical points.

    Parameters
    ----------
    lat1, lon1, lat2, lon2 : float — coordinates in decimal degrees
    """
    R    = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2) ** 2
    )
    return 2 * R * math.asin(math.sqrt(a))


# ---------------------------------------------------------------------------
# Bulletin I/O helpers
# ---------------------------------------------------------------------------

def _generate_global(parameters):
    """Load and return the lines of the main (reference) bulletin file."""
    with open(parameters.main_bulletin_path, 'r') as f:
        lines = f.readlines()

    event_count = sum(1 for line in lines if line.startswith('# '))
    logger.info(f"Main bulletin loaded: {parameters.main_bulletin_path} ({event_count} events)")
    return lines


def retrieve_events_from_lines(cat_lines):
    """
    Extract event header lines and their line indices from a list of bulletin lines.

    Parameters
    ----------
    cat_lines : list of str

    Returns
    -------
    (event_lines, event_line_ids) : (list of str, list of int)
    """
    event_lines    = []
    event_line_ids = []
    for idx, line in enumerate(cat_lines):
        if line.startswith('# '):
            event_lines.append(line.rstrip('\n').lstrip('# '))
            event_line_ids.append(idx)
    return event_lines, event_line_ids


def retrieve_events_from_file(file_name):
    """
    Read an .obs file and return event header lines, their indices, and all bulletin lines.

    Parameters
    ----------
    file_name : str

    Returns
    -------
    (event_lines, event_line_ids, cat_lines) : (list of str, list of int, list of str)
    """
    with open(file_name, 'r', encoding='utf-8', errors='ignore') as f:
        cat_lines = f.readlines()

    event_lines    = []
    event_line_ids = []
    for idx, line in enumerate(cat_lines):
        if line.startswith('# '):
            event_lines.append(line.rstrip('\n').lstrip('# '))
            event_line_ids.append(idx)

    logger.info(f"Secondary bulletin loaded: {file_name} ({len(event_lines)} events)")
    return event_lines, event_line_ids, cat_lines


def _get_first_non_nan(value):
    """Return the first non-NaN float from a colon-separated statistics string."""
    if isinstance(value, str):
        for part in value.split(':'):
            if part.strip().lower() != 'nan':
                try:
                    return float(part)
                except ValueError:
                    continue
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def get_catalog_frame(event_lines):
    """
    Build a DataFrame with lat, lon, depth, magnitude, and time from event header lines.

    Parameters
    ----------
    event_lines : list of str — event header lines (leading '# ' already stripped)

    Returns
    -------
    pd.DataFrame
    """
    infos = [line.split() for line in event_lines]

    frame = pd.DataFrame({
        'year':       [i[0]  for i in infos],
        'month':      [i[1]  for i in infos],
        'day':        [i[2]  for i in infos],
        'hour':       [i[3]  for i in infos],
        'minute':     [i[4]  for i in infos],
        'second':     [i[5]  for i in infos],
        'latitude':   [_get_first_non_nan(i[6])                   for i in infos],
        'longitude':  [_get_first_non_nan(i[7])                   for i in infos],
        'depth':      [_get_first_non_nan(i[8])                   for i in infos],
        'magnitude':  [_get_first_non_nan(i[9].split(':')[0])     for i in infos],
        'magType':    [i[10] for i in infos],
        'magAuthor':  [i[11] for i in infos],
        'phaseCount': [float(i[12]) if i[12] != 'None' else None  for i in infos],
        'horUncer':   [float(i[13]) if i[13] != 'None' else None  for i in infos],
        'verUncer':   [float(i[14]) if i[14] != 'None' else None  for i in infos],
        'azGap':      [float(i[15]) if i[15] != 'None' else None  for i in infos],
        'rms':        [float(i[16]) if i[16] != 'None' else None  for i in infos],
    })

    frame['second'] = frame['second'].apply(lambda x: f"{float(x):.6f}")
    frame['time']   = pd.to_datetime(
        frame['year'] + '-' + frame['month'] + '-' + frame['day'] + 'T' +
        frame['hour'] + ':' + frame['minute'] + ':' + frame['second'] + 'Z'
    )
    return frame


def find_pick_lines(all_lines, event_id):
    """
    Return all pick lines for the event at position event_id in a bulletin line list.

    Parameters
    ----------
    all_lines : list of str
    event_id  : int — index of the event header line
    """
    picks  = []
    curr   = event_id
    while True:
        curr    += 1
        curr_line = all_lines[curr]
        if curr_line.startswith('\n'):
            break
        if not curr_line.startswith('#') and not curr_line.startswith('PUBLIC_ID'):
            picks.append(curr_line)
    return picks


def check_similar_picks(main_lines, secondary_lines, main_id, secondary_id):
    """
    Count picks shared (same station+phase, time difference ≤ 1 s) between two events.

    Parameters
    ----------
    main_lines      : list of str
    secondary_lines : list of str
    main_id         : int — event header index in main_lines
    secondary_id    : int — event header index in secondary_lines

    Returns
    -------
    int — number of similar P-phase picks
    """
    main_picks      = find_pick_lines(main_lines, main_id)
    secondary_picks = find_pick_lines(secondary_lines, secondary_id)
    all_picks       = main_picks + secondary_picks

    all_phases = {}
    all_times  = []

    for line in all_picks:
        if line[22] == 'S':
            continue
        phase_key  = line[:23]
        phase_time = line[31:51]

        if phase_key not in all_phases:
            all_phases[phase_key] = 1
            all_times.append(phase_time)
        else:
            i       = list(all_phases.keys()).index(phase_key)
            date    = pd.to_datetime(
                f"{all_times[i][:8]} {all_times[i][9:13]}{all_times[i][14:]}",
                format="%Y%m%d %H%M%S.%f",
            )
            date_new = pd.to_datetime(
                f"{phase_time[:8]} {phase_time[9:13]}{phase_time[14:]}",
                format="%Y%m%d %H%M%S.%f",
            )
            if abs(date - date_new) <= pd.Timedelta(seconds=1):
                all_phases[phase_key] += 1

    return sum(1 for v in all_phases.values() if v > 1)


# ---------------------------------------------------------------------------
# Event matching
# ---------------------------------------------------------------------------

_MATCH_COLS = [
    'catalog1_idx', 'catalog2_idx', 'distance_km',
    'time_diff_seconds', 'mag_diff', 'mag_type_ML', 'threshold_used',
]

_ALPHA = 0.95   # weight of normalised time_diff in the assignment cost (1-_ALPHA goes to dist)


def find_match_events(
    catalog1, catalog2,
    dist_thresh, loose_dist_thresh,
    time_thresh, loose_time_thresh,
    mag_thresh,
):
    """
    Find matching event pairs between two catalogs using a two-phase
    best-pair-first greedy algorithm with a weighted cost (α=_ALPHA on time,
    1-_ALPHA on distance, both normalised to their loose thresholds).

    Phase 1 (strict): all candidate pairs within strict thresholds are sorted
    by cost and assigned cheapest-first, one-to-one.  Phase 2 (loose): for
    events unmatched after phase 1, candidate pairs within loose thresholds
    are sorted and assigned the same way.  Sorting globally before any
    assignment makes the algorithm symmetric — it implicitly checks both
    directions (best B for A and best A for B) — bounding errors to
    individual competing pairs rather than propagating chain-shifts across
    a cluster.  Loose pairs are returned separately for pick-based validation.

    Parameters
    ----------
    catalog1, catalog2   : pd.DataFrame — produced by get_catalog_frame()
    dist_thresh          : float — strict distance threshold (km)
    loose_dist_thresh    : float — loose distance threshold (km)
    time_thresh          : float — strict time threshold (s)
    loose_time_thresh    : float — loose time threshold (s)
    mag_thresh           : float — magnitude difference limit for ML–ML pairs

    Returns
    -------
    (strict_matches, possible_matches, unmatched_catalog2)
        strict_matches, possible_matches : pd.DataFrame with columns _MATCH_COLS
        unmatched_catalog2               : list of int (catalog2 index labels)
    """
    # ------------------------------------------------------------------
    # 1. Collect all candidate pairs within loose thresholds
    # ------------------------------------------------------------------
    candidates = []   # (idx1, idx2, time_diff, dist_km, mag_diff, is_ml)

    for idx1, row1 in catalog1.iterrows():
        time_diffs  = abs((catalog2.time - row1.time).dt.total_seconds())
        within_time = time_diffs[time_diffs < loose_time_thresh].index.to_numpy()

        for idx2 in within_time:
            row2    = catalog2.loc[idx2]
            dist_km = haversine(row1.latitude, row1.longitude,
                                row2.latitude, row2.longitude)
            if dist_km > loose_dist_thresh:
                continue
            time_diff = float(time_diffs.loc[idx2])
            mag_diff  = abs(row1.magnitude - row2.magnitude)
            is_ml     = (
                row1.magType == 'ML' and row2.magType == 'ML'
                and row1.magAuthor in ('LDG', 'OMP')
                and row2.magAuthor in ('LDG', 'OMP')
            )
            candidates.append((idx1, idx2, time_diff, dist_km, mag_diff, is_ml))

    _empty = pd.DataFrame(columns=_MATCH_COLS)

    if not candidates:
        logger.info("Strict matches: 0  Possible matches: 0")
        return _empty, _empty, list(catalog2.index)

    # ------------------------------------------------------------------
    # 2. Score and split into strict / loose candidate lists
    # ------------------------------------------------------------------
    strict_candidates = []   # (cost, idx1, idx2, time_diff, dist_km, mag_diff, is_ml)
    loose_candidates  = []

    for idx1, idx2, time_diff, dist_km, mag_diff, is_ml in candidates:
        mag_ok    = not is_ml or mag_diff <= mag_thresh
        is_strict = time_diff <= time_thresh and dist_km <= dist_thresh and mag_ok
        cost      = (
            _ALPHA       * (time_diff / loose_time_thresh)
            + (1 - _ALPHA) * (dist_km   / loose_dist_thresh)
        )
        entry = (cost, idx1, idx2, time_diff, dist_km, mag_diff, is_ml)
        if is_strict:
            strict_candidates.append(entry)
        else:
            loose_candidates.append(entry)

    strict_candidates.sort()
    loose_candidates.sort()

    # ------------------------------------------------------------------
    # 3. Phase 1 — strict best-pair-first greedy (one-to-one)
    # ------------------------------------------------------------------
    matched_id1   = set()
    matched_id2   = set()
    matched_pairs = []

    for cost, idx1, idx2, time_diff, dist_km, mag_diff, is_ml in strict_candidates:
        if idx1 in matched_id1 or idx2 in matched_id2:
            continue
        matched_id1.add(idx1)
        matched_id2.add(idx2)
        matched_pairs.append({
            'catalog1_idx':      idx1,
            'catalog2_idx':      idx2,
            'distance_km':       dist_km,
            'time_diff_seconds': time_diff,
            'mag_diff':          mag_diff,
            'mag_type_ML':       is_ml,
            'threshold_used':    'strict',
        })

    # ------------------------------------------------------------------
    # 4. Phase 2 — loose best-pair-first greedy for remaining events
    # ------------------------------------------------------------------
    # Snapshot strict-only claims before phase 2 so that loose candidates
    # remain in unmatched_catalog2.  _concatenate_bulletin uses that list
    # to gate pick-based validation; removing loose events from it would
    # silently drop them from the output catalog.
    strict_matched_id2 = set(matched_id2)

    possible_match = []

    for cost, idx1, idx2, time_diff, dist_km, mag_diff, is_ml in loose_candidates:
        if idx1 in matched_id1 or idx2 in matched_id2:
            continue
        matched_id1.add(idx1)
        matched_id2.add(idx2)
        possible_match.append({
            'catalog1_idx':      idx1,
            'catalog2_idx':      idx2,
            'distance_km':       dist_km,
            'time_diff_seconds': time_diff,
            'mag_diff':          mag_diff,
            'mag_type_ML':       is_ml,
            'threshold_used':    'loose',
        })

    unmatched_catalog2 = [i for i in catalog2.index if i not in strict_matched_id2]

    logger.info(f"Strict matches: {len(matched_pairs)}  Possible matches: {len(possible_match)}")
    return (
        pd.DataFrame(matched_pairs) if matched_pairs else _empty,
        pd.DataFrame(possible_match) if possible_match else _empty,
        unmatched_catalog2,
    )


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------

def _add_item_for_stats(line_main, line_secondary, field_id, is_nan=False):
    """Append the field at field_id from line_secondary to line_main as a colon-separated entry."""
    to_concat = ':Nan' if is_nan else ':' + line_secondary.split()[field_id]
    new_line  = line_main.split()
    new_line[-1] += '\n'
    new_line[field_id] += to_concat
    return ' '.join(new_line)


def _add_nans_for_stats(line_secondary, field_id, loop_no):
    """Prepend loop_no NaN placeholders to the statistics field at field_id in line_secondary."""
    to_replace = 'Nan:' * loop_no + line_secondary.split()[field_id]
    new_line   = line_secondary.split()
    new_line[-1] += '\n'
    new_line[field_id] = to_replace
    return ' '.join(new_line)


def _add_phases_to_lines(new_lines, old_lines, event_id):
    """Copy all pick lines for the event at event_id from old_lines into new_lines."""
    curr = event_id + 1
    while not old_lines[curr].startswith('\n'):
        new_lines.append(old_lines[curr])
        curr += 1
    return new_lines


# ---------------------------------------------------------------------------
# Bulletin merging
# ---------------------------------------------------------------------------

def _concatenate_bulletin(
    parameters, main_lines, secondary_bulletin_path,
    dist_thresh, loose_dist_thresh, time_thresh, loose_time_thresh, mag_thresh,
    loop_no,
):
    """Merge a secondary bulletin into the main bulletin, matching events and accumulating statistics."""
    main_event_lines, main_ids             = retrieve_events_from_lines(main_lines)
    secondary_event_lines, secondary_ids, secondary_lines = retrieve_events_from_file(secondary_bulletin_path)

    main_bulletin      = get_catalog_frame(main_event_lines)
    secondary_bulletin = get_catalog_frame(secondary_event_lines)

    new_lines = [line for line in main_lines if line.startswith('###')]
    new_lines.append('\n')

    strict_match, possible_match, not_matched_secondary = find_match_events(
        main_bulletin, secondary_bulletin,
        dist_thresh, loose_dist_thresh,
        time_thresh, loose_time_thresh,
        mag_thresh,
    )

    found_possible = []

    for event_idx1 in main_bulletin.index:
        if not (strict_match.empty and possible_match.empty):
            match_row    = strict_match[strict_match.catalog1_idx == event_idx1]
            possible_row = possible_match[possible_match.catalog1_idx == event_idx1]
        else:
            match_row    = pd.DataFrame()
            possible_row = pd.DataFrame()

        event_line_main = main_lines[main_ids[event_idx1]]

        if not match_row.empty:
            event_idx2            = match_row['catalog2_idx'].iloc[0]
            event_line_secondary  = secondary_lines[secondary_ids[event_idx2]]

            if match_row.mag_type_ML.item():
                event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 10)
            else:
                event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 10, is_nan=True)

            event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 7)
            event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 8)
            event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 9)

            new_lines.append(event_line_main)
            new_lines = _add_phases_to_lines(new_lines, main_lines, main_ids[event_idx1])
            new_lines = _add_phases_to_lines(new_lines, secondary_lines, secondary_ids[event_idx2])
            new_lines.append('\n')

        elif not possible_row.empty:
            solution_found = False
            for _, row in possible_row.iterrows():
                event_idx2           = row.catalog2_idx
                event_line_secondary = secondary_lines[secondary_ids[event_idx2]]

                if event_idx2 not in not_matched_secondary:
                    continue

                sim_picks = check_similar_picks(
                    main_lines, secondary_lines,
                    main_ids[event_idx1], secondary_ids[event_idx2],
                )

                if sim_picks >= 1:
                    if row.mag_type_ML:
                        event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 10)
                    else:
                        event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 10, is_nan=True)
                    event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 7)
                    event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 8)
                    event_line_main = _add_item_for_stats(event_line_main, event_line_secondary, 9)
                    new_lines.append(event_line_main)
                    new_lines = _add_phases_to_lines(new_lines, main_lines, main_ids[event_idx1])
                    new_lines = _add_phases_to_lines(new_lines, secondary_lines, secondary_ids[event_idx2])
                    new_lines.append('\n')
                    not_matched_secondary.remove(event_idx2)
                    solution_found = True
                    found_possible.append(possible_row.index[0])
                    break

            if not solution_found:
                for field_id in (7, 8, 9, 10):
                    event_line_main = _add_item_for_stats(event_line_main, '', field_id, is_nan=True)
                new_lines.append(event_line_main)
                new_lines = _add_phases_to_lines(new_lines, main_lines, main_ids[event_idx1])
                new_lines.append('\n')

        else:
            for field_id in (7, 8, 9, 10):
                event_line_main = _add_item_for_stats(event_line_main, '', field_id, is_nan=True)
            new_lines.append(event_line_main)
            new_lines = _add_phases_to_lines(new_lines, main_lines, main_ids[event_idx1])
            new_lines.append('\n')

    for event_idx2 in not_matched_secondary:
        event_line_secondary = secondary_lines[secondary_ids[event_idx2]]
        for field_id in (7, 8, 9, 10):
            event_line_secondary = _add_nans_for_stats(event_line_secondary, field_id, loop_no)
        new_lines.append(event_line_secondary)
        new_lines = _add_phases_to_lines(new_lines, secondary_lines, secondary_ids[event_idx2])
        new_lines.append('\n')

    if found_possible:
        validated    = possible_match.loc[found_possible]
        to_concat    = [df for df in [strict_match, validated] if not df.empty]
        strict_match = pd.concat(to_concat).reset_index(drop=True)
    possible_match = possible_match.drop(found_possible)
    logger.info(f"P-phase pick matching: {len(found_possible)} additional matches "
                f"({len(possible_match)} remaining unmatched)")

    new_lines = _sort_events_chrono(new_lines)
    return new_lines, strict_match, possible_match, main_bulletin, secondary_bulletin


# ---------------------------------------------------------------------------
# Post-processing helpers
# ---------------------------------------------------------------------------

def _sort_events_chrono(lines):
    """Sort all events in a bulletin line list into chronological order."""
    headers     = [line for line in lines if line.startswith('###')]
    headers.append('\n')
    event_lines = [line for line in lines if not line.startswith('###')]

    events      = []
    current     = None

    for line in event_lines:
        if line.startswith('#'):
            if current is not None:
                events.append(current)
            current = [line]
        elif line.strip() == '':
            if current is not None:
                events.append(current)
                current = None
        elif current is not None:
            current.append(line)

    if current is not None:
        events.append(current)

    def _get_timestamp(event):
        parts = event[0].split()
        return tuple(map(float, parts[1:7]))

    events.sort(key=_get_timestamp)

    sorted_lines = list(headers)
    for event in events:
        sorted_lines.extend(event)
        sorted_lines.append('\n')
    return sorted_lines


def _replace_mean_magnitudes(lines):
    """Replace colon-separated multi-source magnitude field with the mean of non-NaN values."""
    for idx, line in enumerate(lines):
        if line.startswith('# '):
            parts  = line.split()
            mags   = [float(m) for m in parts[10].split(':') if m != 'Nan']
            parts[10]  = f'{mean(mags):.2f}'
            parts[-1] += '\n'
            lines[idx] = ' '.join(parts)
    logger.info("Magnitudes replaced by per-event mean across sources")
    return lines


def _remove_stats_values(lines):
    """Strip colon-separated statistics entries from event fields, keeping only the first non-NaN value."""
    for idx, line in enumerate(lines):
        if line.startswith('# '):
            new_line = ''
            for token in line.split():
                items = token.split(':')
                if len(items) > 1:
                    token = str(_get_first_non_nan(token))
                new_line += token + ' '
            lines[idx] = new_line.rstrip() + '\n'
    return lines


def _remove_duplicate_picks(lines):
    """Remove duplicate pick lines (same station and phase) within each event block."""
    to_remove = set()
    for idx, line in enumerate(lines):
        if line.startswith('# '):
            i            = idx + 1
            unique_picks = set()
            while i < len(lines) and not lines[i].startswith('\n'):
                pick = (lines[i][:10], lines[i][22])
                if pick not in unique_picks:
                    unique_picks.add(pick)
                else:
                    to_remove.add(i)
                i += 1

    new_lines = [line for idx, line in enumerate(lines) if idx not in to_remove]
    logger.info(f"Duplicate picks removed: {len(to_remove)}")
    return new_lines


def _remove_magnitudes_under_1(lines):
    """Remove all events with magnitude below ML 1.0 and their associated pick lines."""
    to_remove    = set()
    n_removed    = 0
    for idx, line in enumerate(lines):
        if line.startswith('# '):
            if float(line.split()[10]) < 1:
                to_remove.add(idx)
                n_removed += 1
                i = idx + 1
                while i < len(lines) and not lines[i].startswith('\n'):
                    to_remove.add(i)
                    i += 1

    new_lines = [line for idx, line in enumerate(lines) if idx not in to_remove]
    logger.info(f"Events removed (magnitude < 1.0): {n_removed}")
    return new_lines


def _save_bulletin(lines, parameters):
    """Write the bulletin line list to the global output .obs file, assigning sequential PUBLIC_IDs."""
    output  = []
    counter = 0
    for line in lines:
        output.append(line)
        if line.startswith('# '):
            counter += 1
            output.append(f'PUBLIC_ID PYRENEES_{counter:06d}\n')
    with open(parameters.global_bulletin_path, 'w') as f:
        f.writelines(output)
    logger.info(f"Final bulletin saved: {parameters.global_bulletin_path} ({counter} events)")


# ---------------------------------------------------------------------------
# Statistics figures
# ---------------------------------------------------------------------------

def _stats_figs_versus(main_name, secondary_name, frame):
    """Generate scatter + KDE plots comparing hypocenter parameters between two catalogs."""
    use_cols  = ['latitude', 'longitude', 'depth', 'magnitude']
    use_frame = frame[use_cols]

    _, axs = plt.subplots(nrows=2, ncols=2, figsize=(18, 12))
    axs    = axs.flatten()
    plt.rc('axes', labelsize=13)

    for plot_idx, col in enumerate(use_cols):
        ax   = axs[plot_idx]
        col1, col2 = zip(*[(t[0], t[1]) for t in use_frame[col]])
        data = pd.DataFrame({main_name: col1, secondary_name: col2})

        data_99 = data[
            (data[main_name]      >= data[main_name].quantile(0.005)) &
            (data[main_name]      <= data[main_name].quantile(0.995)) &
            (data[secondary_name] >= data[secondary_name].quantile(0.005)) &
            (data[secondary_name] <= data[secondary_name].quantile(0.995))
        ].copy()

        sns.scatterplot(x=data_99[main_name], y=data_99[secondary_name],
                        ax=ax, color='black', s=2, alpha=0.6, edgecolor=None)
        sns.kdeplot(x=data_99[main_name], y=data_99[secondary_name],
                    ax=ax, cmap='flare', fill=True, alpha=0.65)

        ax.set_xlabel(main_name)
        ax.set_ylabel(secondary_name)
        ax.grid(True)

        r_p, p_p = pearsonr(data_99[main_name], data_99[secondary_name])
        r_s, p_s = spearmanr(data_99[main_name], data_99[secondary_name])
        ax.text(0.5, 1.085,
                f'Pearson: {r_p:.3f} - p-value: {p_p:.3f}\nSpearman: {r_s:.3f} - p-value: {p_s:.3f}',
                transform=ax.transAxes, fontsize=12, ha='center', va='top')
        label = col.capitalize() if col != 'magnitude' else 'Magnitude (ML)'
        ax.text(1.02, 0.5, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', ha='center', va='center', rotation=90)

    plt.subplots_adjust(top=0.88, bottom=0.1, wspace=0.3, hspace=0.25)
    plt.suptitle(f'Correlations and KDE/Distributions ({main_name} vs {secondary_name}) - matched events',
                 fontsize=16, fontweight='bold')
    plt.text(0.5, 0.95, 'Analysis based on the central 99% of the dataset',
             fontsize=14, ha='center', va='center', transform=plt.gcf().transFigure)
    plt.text(0.5, 0.93,
             'Pearson (linear) and Spearman (non-linear) correlation values are statistically significant for p-values under 0.05',
             fontsize=14, ha='center', va='center', transform=plt.gcf().transFigure)

    path = f'obs/STATS/{main_name}_{secondary_name}_versus.pdf'
    plt.savefig(path)
    plt.close()
    print(f'Statistics "versus" figure successfully saved @ {path}')


def _stats_figs_comparison(main_name, secondary_name, frame):
    """Generate difference (residual) plots for hypocenter parameters between two catalogs."""
    use_cols  = ['latitude', 'longitude', 'depth', 'magnitude']
    use_frame = frame[use_cols]

    _, axs = plt.subplots(nrows=2, ncols=2, figsize=(18, 12))
    axs    = axs.flatten()
    plt.rc('axes', labelsize=13)

    for plot_idx, col in enumerate(use_cols):
        ax   = axs[plot_idx]
        col1, col2 = zip(*[(t[0], t[1]) for t in use_frame[col]])
        data = pd.DataFrame({main_name: col1, secondary_name: col2})

        data_99 = data[
            (data[main_name]      >= data[main_name].quantile(0.005)) &
            (data[main_name]      <= data[main_name].quantile(0.995)) &
            (data[secondary_name] >= data[secondary_name].quantile(0.005)) &
            (data[secondary_name] <= data[secondary_name].quantile(0.995))
        ].copy()

        data_99['diff'] = data_99[main_name] - data_99[secondary_name]
        if col == 'latitude':
            data_99['diff'] *= 111.32
        elif col == 'longitude':
            data_99['diff'] *= 111.32 * math.cos(42 * math.pi / 180)

        sns.scatterplot(x=data_99[main_name], y=data_99['diff'],
                        ax=ax, color='black', s=2, alpha=0.6, edgecolor=None)
        sns.kdeplot(x=data_99[main_name], y=data_99['diff'],
                    ax=ax, cmap='flare', fill=True, alpha=0.65)

        ax.set_ylim(-1, 1) if col == 'magnitude' else ax.set_ylim(-20, 20)
        ax.set_xlabel(main_name)
        ax.set_ylabel(f'{main_name} - {secondary_name}')
        ax.grid(True)

        r_p, p_p = pearsonr(data_99[main_name], data_99[secondary_name])
        r_s, p_s = spearmanr(data_99[main_name], data_99[secondary_name])
        ax.text(0.5, 1.085,
                f'Pearson: {r_p:.3f} - p-value: {p_p:.3f}\nSpearman: {r_s:.3f} - p-value: {p_s:.3f}',
                transform=ax.transAxes, fontsize=12, ha='center', va='top')
        label = f'{col} (km)'.capitalize() if col != 'magnitude' else 'Magnitude (ML)'
        ax.text(1.02, 0.5, label, transform=ax.transAxes, fontsize=12,
                fontweight='bold', ha='center', va='center', rotation=90)

    plt.subplots_adjust(top=0.88, bottom=0.1, wspace=0.3, hspace=0.25)
    plt.suptitle(f'Correlations and KDE/Distributions ({main_name} vs {secondary_name}) - matched events',
                 fontsize=16, fontweight='bold')
    plt.text(0.5, 0.95, 'Analysis based on the central 99% of the dataset',
             fontsize=14, ha='center', va='center', transform=plt.gcf().transFigure)
    plt.text(0.5, 0.93,
             'Pearson (linear) and Spearman (non-linear) correlation values are statistically significant for p-values under 0.05',
             fontsize=14, ha='center', va='center', transform=plt.gcf().transFigure)

    path = f'obs/STATS/{main_name}_{secondary_name}_comparison.pdf'
    plt.savefig(path)
    plt.close()
    print(f'Statistics "comparison" figure successfully saved @ {path}')

    path_frame = f'obs/STATS/{main_name}_{secondary_name}.csv'
    use_frame.to_csv(path_frame, index=False)
    print(f'Statistics file successfully saved @ {path_frame}')


def _get_statistics(main_lines, parameters, file_path, file_no):
    """Extract matched parameter pairs from bulletin statistics fields and save comparison figures."""
    lines = [line.lstrip('# ').rstrip('\n').split() for line in main_lines if line.startswith('# ')]

    remove_ids = []
    for idx, line in enumerate(lines):
        for col_idx, token in enumerate(line):
            items = token.split(':')
            if len(items) > 1:
                token = items[0] + ':' + items[file_no]
                if 'Nan' in token:
                    remove_ids.append(idx)
                    continue
                lines[idx][col_idx] = token

    lines = [line for idx, line in enumerate(lines) if idx not in remove_ids]

    if not lines:
        print(f'Not enough matches for a statistical analysis for Bulletin @ {file_path}')
        return

    df = pd.DataFrame(lines).iloc[:, [6, 7, 8, 9]]
    df.columns = ['latitude', 'longitude', 'depth', 'magnitude']

    def _split_to_floats(s):
        try:
            a, b = map(float, s.split(':'))
            return [a, b]
        except Exception:
            return [None, None]

    df = df.map(_split_to_floats)

    main_name      = parameters.main_bulletin_path.split('/')[-1].split('.')[0]
    secondary_name = file_path.split('/')[-1].split('.')[0]

    if len(df) >= 10:
        _stats_figs_versus(main_name, secondary_name, df)
        _stats_figs_comparison(main_name, secondary_name, df)
    else:
        print(f'Not enough matches for a statistical analysis for Bulletin @ {file_path}')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fuse_bulletins(parameters, log_dir=None):
    """
    Merge all source .obs bulletins into a single GLOBAL.obs file.

    Events are matched spatially and temporally across all source files.
    Per-source location statistics are accumulated, then replaced by their mean.
    Events below ML 1.0 are removed.

    Parameters
    ----------
    parameters : FusionParams
    log_dir    : str, optional — log directory; default: global_obs/console_output/

    Returns
    -------
    dict with key: output
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file        : {log_path}")
    logger.info(f"Main bulletin   : {parameters.main_bulletin_path}")
    logger.info(f"Source glob     : {parameters.folder_path}")
    logger.info(f"Output          : {parameters.global_bulletin_path}")
    logger.info(
        f"Thresholds — strict: dist={parameters.dist_thresh} km  time={parameters.time_thresh} s  "
        f"loose: dist={parameters.loose_dist_thresh} km  time={parameters.loose_time_thresh} s  "
        f"mag={parameters.mag_thresh}"
    )

    all_paths  = [
        p for p in glob.glob(parameters.folder_path)
        if p != parameters.main_bulletin_path and 'GLOBAL' not in p
    ]
    logger.info(f"Secondary bulletins to merge: {len(all_paths)}")

    main_lines = _generate_global(parameters)

    for file_no, file_path in enumerate(all_paths):
        logger.info(f"--- Merging [{file_no + 1}/{len(all_paths)}]: {file_path}")
        main_lines, _, _, _, _ = _concatenate_bulletin(
            parameters, main_lines, file_path,
            parameters.dist_thresh, parameters.loose_dist_thresh,
            parameters.time_thresh, parameters.loose_time_thresh,
            parameters.mag_thresh,
            file_no + 1,
        )

    for file_no, file_path in enumerate(all_paths):
        _get_statistics(main_lines, parameters, file_path, file_no + 1)

    main_lines = _replace_mean_magnitudes(main_lines)
    main_lines = _remove_stats_values(main_lines)
    # main_lines = _remove_magnitudes_under_1(main_lines)
    main_lines = _remove_duplicate_picks(main_lines)

    _save_bulletin(main_lines, parameters)
    return {'output': parameters.global_bulletin_path}


def find_and_merge_doubles(parameters, log_dir=None):
    """
    Scan a bulletin for suspiciously close event pairs and let the user resolve them interactively.

    Detection: events are connected when |Δt| ≤ max_dt_seconds AND 3-D distance ≤ max_dist_km.
    Connected components (groups of 2 or more events) are presented one group at a time,
    so triples and larger clusters are handled correctly in a single review step.

    Auto-merge: groups of exactly 2 events where |Δt| ≤ auto_dt_seconds AND distance ≤ auto_dist_km
    are merged automatically (first event kept) without user interaction.

    Interactive choices
    -------------------
    k<n> → keep Event n, drop all others; unique phases from dropped events are merged in
    s    → keep all events in the group (not doubles)
    p    → print phase lines for all events side by side, then re-display the group prompt

    Parameters
    ----------
    parameters : MergeDoublesParams
    log_dir    : str, optional — log directory; default: global_obs/console_output/

    Returns
    -------
    dict with key: output
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file        : {log_path}")
    logger.info(f"Bulletin        : {parameters.global_bulletin_path}")
    logger.info(f"Thresholds      : max_dt={parameters.max_dt_seconds} s  max_dist={parameters.max_dist_km} km  "
                f"auto_dt={parameters.auto_dt_seconds} s  auto_dist={parameters.auto_dist_km} km")

    with open(parameters.global_bulletin_path, 'r') as f:
        bulletin_lines = f.readlines()

    # ------------------------------------------------------------------ #
    # 0. Parse bulletin into event records                                 #
    # ------------------------------------------------------------------ #
    def _parse_header(line):
        tokens = line.split()
        if len(tokens) < 10:
            return None
        try:
            year, month, day = int(tokens[1]), int(tokens[2]), int(tokens[3])
            hour, minute     = int(tokens[4]), int(tokens[5])
            sec_f            = float(tokens[6])
            sec_int          = int(sec_f)
            microsec         = int(round((sec_f - sec_int) * 1e6))
            lat, lon, dep    = float(tokens[7]), float(tokens[8]), float(tokens[9])
            ts = pd.Timestamp(year=year, month=month, day=day,
                              hour=hour, minute=minute, second=sec_int, microsecond=microsec)
            return {'time': ts, 'lat': lat, 'lon': lon, 'dep': dep}
        except (ValueError, IndexError):
            return None

    header_lines = bulletin_lines[:4]
    events       = []
    current      = None

    for abs_idx, line in enumerate(bulletin_lines):
        if line.startswith('# '):
            if current is not None:
                events.append(current)
            parsed = _parse_header(line)
            if parsed is None:
                current = None
                continue
            current = {'bid': abs_idx, 'header': line, 'phases': [], 'blank': [], **parsed}
        elif line.startswith('\n') or line == '':
            if current is not None:
                current['blank'].append(line)
                events.append(current)
                current = None
        else:
            if current is not None:
                current['phases'].append(line.rstrip('\n'))

    if current is not None:
        events.append(current)

    if not events:
        print('  No events found in bulletin — nothing to do.')
        return {'output': parameters.global_bulletin_path}

    # ------------------------------------------------------------------ #
    # 1. Detect groups of potential doubles (connected components)         #
    # ------------------------------------------------------------------ #
    def _to_cartesian(lat_deg, lon_deg, depth_km):
        R   = 6371.0
        lat = np.radians(lat_deg)
        lon = np.radians(lon_deg)
        r   = R - depth_km
        return np.array([r * np.cos(lat) * np.cos(lon),
                         r * np.cos(lat) * np.sin(lon),
                         r * np.sin(lat)])

    events.sort(key=lambda e: e['time'])
    max_dt = pd.Timedelta(seconds=parameters.max_dt_seconds)

    adj = {i: set() for i in range(len(events))}
    for i in range(len(events)):
        for j in range(i + 1, len(events)):
            if abs(events[j]['time'] - events[i]['time']) > max_dt:
                break
            dist_km = float(np.linalg.norm(
                _to_cartesian(events[i]['lat'], events[i]['lon'], events[i]['dep']) -
                _to_cartesian(events[j]['lat'], events[j]['lon'], events[j]['dep'])
            ))
            if dist_km <= parameters.max_dist_km:
                adj[i].add(j)
                adj[j].add(i)

    visited = set()
    groups  = []
    for start in range(len(events)):
        if start in visited or not adj[start]:
            continue
        group = []
        stack = [start]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            group.append(node)
            stack.extend(adj[node] - visited)
        if len(group) > 1:
            groups.append(sorted(group))

    if not groups:
        print('  No potential doubles found.')
        return {'output': parameters.global_bulletin_path}

    print(f'\n  Found {len(groups)} group(s) of potential doubles to review.\n')

    # ------------------------------------------------------------------ #
    # 2. Interactive review                                                #
    # ------------------------------------------------------------------ #
    drop_bids = set()
    merge_map = {}

    def _phase_key(phase_line):
        tokens = phase_line.split()
        return (tokens[0], tokens[4]) if len(tokens) >= 5 else None

    def _unique_phases(kept_phases, donor_phases):
        kept_keys = {_phase_key(p) for p in kept_phases if _phase_key(p) is not None}
        return [p for p in donor_phases if _phase_key(p) not in kept_keys and _phase_key(p) is not None]

    sep = '─' * 72

    for group_num, group in enumerate(groups, start=1):
        group_events = [events[i] for i in group]
        n_ev         = len(group_events)

        if n_ev == 2:
            e1, e2   = group_events[0], group_events[1]
            auto_dt  = abs((e1['time'] - e2['time']).total_seconds())
            auto_d   = float(np.linalg.norm(
                _to_cartesian(e1['lat'], e1['lon'], e1['dep']) -
                _to_cartesian(e2['lat'], e2['lon'], e2['dep'])
            ))
            if auto_dt <= parameters.auto_dt_seconds and auto_d <= parameters.auto_dist_km:
                kept_phases = list(e1['phases'])
                extra       = _unique_phases(kept_phases, e2['phases'])
                kept_phases.extend(extra)
                drop_bids.add(e2['bid'])
                if extra:
                    merge_map.setdefault(e1['bid'], []).extend(extra)
                print(f'  [AUTO] Merged BulletinID={e2["bid"]} into BulletinID={e1["bid"]}  '
                      f'Δt={auto_dt:.3f}s  Δd={auto_d:.2f} km  {len(extra)} unique phase(s) merged.')
                logger.info(f"Auto-merged BulletinID={e2['bid']} into BulletinID={e1['bid']}  "
                            f"dt={auto_dt:.3f}s  dist={auto_d:.2f}km")
                continue

        def _print_group():
            print(f'\n{sep}')
            print(f'  POTENTIAL DOUBLES  Group {group_num}/{len(groups)}  ({n_ev} events)')
            print(sep)
            for k, e in enumerate(group_events, start=1):
                print(f'  Event {k:<4} │ Time: {e["time"]}  '
                      f'Lat: {e["lat"]:.4f}  Lon: {e["lon"]:.4f}  Dep: {e["dep"]:.1f} km')
            print(sep)
            print(f'  k<n>          → keep Event n (1–{n_ev}), merge all others, drop none')
            print(f'  k<n>d<x,y,…>  → keep Event n, drop x,y,… completely, merge the rest')
            print(f'  s             → keep all (not doubles)')
            print(f'  p             → show phases for all events')

        def _print_all_phases():
            for k, e in enumerate(group_events, start=1):
                print(f'\n  Event {k}  (BulletinID={e["bid"]})')
                print(sep)
                for ph in e['phases'] or ['    (no phases)']:
                    print(f'    {ph}')
            print(sep)

        _print_group()

        kept_idx      = None
        discard_set   = set()   # 0-based indices to drop without phase merge

        while True:
            choice = input('  Your choice: ').strip().lower()

            if choice == 's':
                break

            if choice == 'p':
                _print_all_phases()
                _print_group()
                continue

            if choice.startswith('k'):
                rest  = choice[1:]
                k_str, _, d_str = rest.partition('d')

                if not k_str.isdigit():
                    print(f'  Invalid input — please enter k<n>, k<n>d<x,y,…>, s, or p<n>.')
                    continue

                k = int(k_str)
                if not (1 <= k <= n_ev):
                    print(f'  Invalid event number — enter a value between 1 and {n_ev}.')
                    continue

                if d_str:
                    try:
                        drop_nums = [int(x.strip()) for x in d_str.split(',')]
                    except ValueError:
                        print(f'  Invalid drop list — use comma-separated numbers, e.g. k2d1,3.')
                        continue
                    if any(not (1 <= d <= n_ev) for d in drop_nums):
                        print(f'  Drop list contains an invalid event number (must be 1–{n_ev}).')
                        continue
                    if k in drop_nums:
                        print(f'  Cannot discard the kept event (Event {k}).')
                        continue
                    discard_set = {d - 1 for d in drop_nums}

                kept_idx = k - 1
                break

            print(f'  Invalid input — please enter k<n>, k<n>d<x,y,…>, s, or p<n>.')

        if kept_idx is None:   # choice was 's'
            print('  → Kept all.\n')
            logger.info(f"Group {group_num}: kept all events (user skipped)")
            continue

        kept        = group_events[kept_idx]
        kept_phases = list(kept['phases'])

        for i, e in enumerate(group_events):
            if i == kept_idx:
                continue
            drop_bids.add(e['bid'])
            if i not in discard_set:
                extra = _unique_phases(kept_phases, e['phases'])
                kept_phases.extend(extra)
                if extra:
                    merge_map.setdefault(kept['bid'], []).extend(extra)

        n_extra = len(kept_phases) - len(kept['phases'])
        if discard_set:
            discard_label = ', '.join(f'Event {i + 1}' for i in sorted(discard_set))
            discard_bids  = ', '.join(str(group_events[i]['bid']) for i in sorted(discard_set))
            print(f'  → Kept Event {kept_idx + 1} (BulletinID={kept["bid"]}), '
                  f'discarded {discard_label} (no phase merge), '
                  f'merged remaining.  {n_extra} unique phase(s) added.\n')
            logger.info(f"Group {group_num}: kept BulletinID={kept['bid']}, "
                        f"discarded BulletinID(s)={discard_bids} (no phase merge), "
                        f"merged remaining, {n_extra} unique phase(s) added")
        else:
            n_merged = len(group_events) - 1
            merged_bids = ', '.join(str(e['bid']) for i, e in enumerate(group_events) if i != kept_idx)
            print(f'  → Kept Event {kept_idx + 1} (BulletinID={kept["bid"]}), '
                  f'merged {n_merged} other event(s).  {n_extra} unique phase(s) added.\n')
            logger.info(f"Group {group_num}: kept BulletinID={kept['bid']}, "
                        f"merged BulletinID(s)={merged_bids}, {n_extra} unique phase(s) added")

    # ------------------------------------------------------------------ #
    # 3. Rebuild bulletin                                                  #
    # ------------------------------------------------------------------ #
    updated = list(header_lines)
    for e in events:
        if e['bid'] in drop_bids:
            continue
        phase_lines = list(e['phases'])
        phase_lines.extend(merge_map.get(e['bid'], []))
        phase_lines = [p if p.endswith('\n') else p + '\n' for p in phase_lines]
        updated.append(e['header'])
        updated.extend(phase_lines)
        updated.extend(e['blank'])

    n_dropped = len(drop_bids)
    n_merged  = sum(len(v) for v in merge_map.values())
    print(f'  Done.  {n_dropped} event(s) removed, '
          f'{n_merged} phase(s) merged across {len(merge_map)} event(s).')
    logger.info(f"Doubles resolved: {n_dropped} event(s) removed, "
                f"{n_merged} phase(s) merged across {len(merge_map)} event(s)")

    _save_bulletin(updated, parameters)
    return {'output': parameters.global_bulletin_path}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Fuse source .obs bulletins into GLOBAL.obs, or resolve duplicate events interactively.'
    )
    parser.add_argument('--mode', choices=['fuse', 'doubles'], default='fuse',
                        help='"fuse" merges all source bulletins; "doubles" resolves duplicates interactively')

    # --- fuse mode arguments ---
    parser.add_argument('--global-bulletin-path', default=None)
    parser.add_argument('--main-bulletin-path',   default=None)
    parser.add_argument('--folder-path',          default=None,
                        help='Glob pattern for source .obs files')
    parser.add_argument('--dist-thresh',       type=float, default=15.0)
    parser.add_argument('--loose-dist-thresh', type=float, default=30.0)
    parser.add_argument('--time-thresh',       type=float, default=2.0)
    parser.add_argument('--loose-time-thresh', type=float, default=10.0)
    parser.add_argument('--mag-thresh',        type=float, default=1.5)

    # --- doubles mode arguments ---
    parser.add_argument('--max-dt-seconds',  type=float, default=1.0)
    parser.add_argument('--max-dist-km',     type=float, default=50.0)
    parser.add_argument('--auto-dt-seconds', type=float, default=0.15,
                        help='Δt threshold (s) for automatic merging of pairs (default: 0.15)')
    parser.add_argument('--auto-dist-km',    type=float, default=10.0,
                        help='Distance threshold (km) for automatic merging of pairs (default: 10.0)')

    args = parser.parse_args()

    if args.mode == 'fuse':
        params = FusionParams(
            global_bulletin_path = args.global_bulletin_path,
            main_bulletin_path   = args.main_bulletin_path,
            folder_path          = args.folder_path,
            dist_thresh          = args.dist_thresh,
            loose_dist_thresh    = args.loose_dist_thresh,
            time_thresh          = args.time_thresh,
            loose_time_thresh    = args.loose_time_thresh,
            mag_thresh           = args.mag_thresh,
        )
        fuse_bulletins(params)
    else:
        params = MergeDoublesParams(
            global_bulletin_path = args.global_bulletin_path,
            max_dt_seconds       = args.max_dt_seconds,
            max_dist_km          = args.max_dist_km,
            auto_dt_seconds      = args.auto_dt_seconds,
            auto_dist_km         = args.auto_dist_km,
        )
        find_and_merge_doubles(params)


if __name__ == '__main__':
    main()
