"""
match_pre_post_relocation.py
============================
Match pre-NLL .obs events to their NonLinLoc-relocated counterparts.

Reads the global GLOBAL.obs bulletin and the NLL result file, then builds a
bipartite candidate graph, finds connected components (local groups), and
resolves each group — automatically for trivial 1:1 cases, interactively for
ambiguous or multi-event clusters.  Unique picks from duplicate events can be
merged into the kept event.  The matched bulletin is written to a new output
file.

Usage
-----
    python NLL_run/match_pre_post_relocation.py \\
        --obs    obs/GLOBAL.obs \\
        --final  RESULT/FINAL.txt \\
        --output obs/FINAL.obs
"""

import argparse
import logging
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd
from obspy import UTCDateTime as Timing
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

logger = logging.getLogger('match_pre_post_relocation')


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

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
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class MatchCatalogsParams:
    file_obs:   str
    file_final: str
    save_file:  str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_obs(file):
    """
    Parse an .obs bulletin into a DataFrame of event headers plus the raw
    lines list.

    Parameters
    ----------
    file : str — path to the .obs bulletin

    Returns
    -------
    (pd.DataFrame, list[str])
    """
    with open(file, 'r') as f:
        lines = f.readlines()

    obs = [
        [*line.lstrip('# ').rstrip('\n').split(), idx]
        for idx, line in enumerate(lines)
        if line.startswith('# ')
    ]

    df = pd.DataFrame(obs, columns=[
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'Mag', 'MagType', 'MagAuthor',
        'PhaseCount', 'HorUncer', 'VerUncer', 'AzGap', 'RMS', 'BulletinID',
    ])

    cols_to_num = [
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'Mag', 'PhaseCount',
        'HorUncer', 'VerUncer', 'AzGap', 'RMS',
    ]
    df[cols_to_num] = df[cols_to_num].apply(pd.to_numeric, errors='coerce')
    df['Time'] = pd.to_datetime(
        df['Year'].astype(str) + '-'
        + df['Month'].astype(str) + '-'
        + df['Day'].astype(str) + 'T'
        + df['Hour'].astype(str) + ':'
        + df['Min'].astype(str) + ':'
        + df['Sec'].astype(str) + 'Z'
    )
    return df, lines


def _read_final(file):
    """
    Parse a plain-text NLL result file (hypo_71 format) into a DataFrame.

    Parameters
    ----------
    file : str — path to the NLL result bulletin

    Returns
    -------
    pd.DataFrame
    """
    with open(file, 'r') as f:
        lines = f.readlines()

    final = [line.rstrip('\n').split() for line in lines]
    final = [
        [int(line[0]) + 1900 if int(line[0]) > 75 else int(line[0]) + 2000] + line[1:]
        for line in final
    ]

    df = pd.DataFrame(final, columns=[
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'Mag', 'RMS', 'PhaseCount',
        'HorUncer', 'VerUncer', 'AzGap',
    ])

    cols_to_num = [
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'Mag', 'PhaseCount',
        'HorUncer', 'VerUncer', 'AzGap', 'RMS',
    ]
    df[cols_to_num] = df[cols_to_num].apply(pd.to_numeric, errors='coerce')
    df['Time'] = pd.to_datetime(
        df['Year'].astype(str) + '-'
        + df['Month'].astype(str) + '-'
        + df['Day'].astype(str) + 'T'
        + df['Hour'].astype(str) + ':'
        + df['Min'].astype(str) + ':'
        + df['Sec'].astype(str) + 'Z'
    )
    return df


def _phase_key(phase_line: str):
    """Return a (station, phase_type) tuple that identifies a phase, or None."""
    tokens = phase_line.split()
    if len(tokens) < 5:
        return None
    return (tokens[0], tokens[4])


def _diff_phases(kept_bid: int, dropped_bid: int, phase_index: dict) -> list:
    """
    Return phase lines from dropped_bid that are NOT already in kept_bid,
    compared by (station, phase_type) key.
    """
    if not phase_index:
        return []
    kept_phases    = phase_index.get(kept_bid,    [])
    dropped_phases = phase_index.get(dropped_bid, [])
    kept_keys      = {_phase_key(p) for p in kept_phases if _phase_key(p) is not None}
    return [
        ph for ph in dropped_phases
        if _phase_key(ph) is not None and _phase_key(ph) not in kept_keys
    ]


def _update_bulletin(lines, matched_df):
    """
    Rebuild a bulletin from lines, keeping only events present in matched_df
    and updating their header lines with the matched values.

    If a row carries a non-empty '_extra_phases' list, those phase lines are
    appended at the end of the event's phase block.

    Parameters
    ----------
    lines      : list[str]    — raw bulletin lines
    matched_df : pd.DataFrame — output of match_catalogues()

    Returns
    -------
    list[str] — updated bulletin lines ready to be written to disk
    """
    blocks        = []
    current_block = []
    event_indices = []

    for i, line in enumerate(lines):
        if line.startswith('# '):
            if current_block:
                blocks.append(current_block)
            current_block = [line]
            event_indices.append(i)
        else:
            if current_block:
                current_block.append(line)
    if current_block:
        blocks.append(current_block)

    updated_blocks = []
    for block, event_index in zip(blocks, event_indices):
        if not block:
            continue

        row = matched_df[matched_df.BulletinID == event_index]
        if row.empty:
            continue

        r = row.iloc[0]

        extra = r['_extra_phases']
        has_extra = extra is not None and len(extra) > 0
        phase_count = int(r.PhaseCount) + (len(extra) if has_extra else 0)

        updated_event_line = (
            f"# {r.Year} {r.Month} {r.Day} {r.Hour} {r.Min} {r.Sec} "
            f"{r.Lat} {r.Lon} {r.Dep} {r.Mag} {r.MagType} {r.MagAuthor} "
            f"{phase_count} {r.HorUncer} {r.VerUncer} {r.AzGap} {r.RMS}\n"
        )

        phase_lines    = [l for l in block[1:] if not (l.startswith('\n') or l == '')]
        trailing_blank = [l for l in block[1:] if      l.startswith('\n') or l == '']

        if has_extra:
            extra_lines = [ph if ph.endswith('\n') else ph + '\n' for ph in extra]
            phase_lines = phase_lines + extra_lines

        updated_blocks.append([updated_event_line] + phase_lines + trailing_blank)

    lines[0] = f'### Bulletin generated on the {Timing.now()}\n'
    updated_content = lines[0:4]
    for block in updated_blocks:
        updated_content.extend(block)

    return updated_content


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def match_catalogues(
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    max_dt_seconds: float = 10.0,
    max_dist_km: float = 50.0,
    bulletin_lines: list = None,
) -> pd.DataFrame:
    """
    Match seismic events from catalogue 1 (obs) to catalogue 2 (NLL result),
    one-to-one, using a group-based approach.

    Algorithm
    ---------
    1. Build all candidate pairs (df1_i, df2_j) within max_dt_seconds and
       max_dist_km.
    2. Find connected components in the bipartite df1↔df2 candidate graph.
       Each component is a local "group."
    3. Resolve each group:
         1:1  → auto-assign.
         1:N  → auto-pick closest by time if unambiguous (gap ≥ 0.5 s);
                otherwise interactive (candidate number prompt, next df1
                event shown for context).
         M:N  → always interactive (letter/number assignment syntax).
    4. Build the output DataFrame.

    Interactive choices — 1:N prompt
    ---------------------------------
    1…N  → assign NLL candidate N
    s    → skip this obs event (no update)
    p    → print phase lines for the obs event

    Interactive choices — M:N prompt (obs events A/B/C…, NLL candidates 1/2/3…)
    ----------------------------------------------------------------------------
    A1        → assign obs-A to NLL-1; unspecified obs events auto-merge into
                their temporally nearest assigned obs event
    A1,B2     → assign obs-A→NLL-1 and obs-B→NLL-2
    A1,BmA    → assign obs-A→NLL-1, merge obs-B's phases into obs-A
    A1,Bm     → same (single-assignment: target is inferred)
    A1,Bd     → assign obs-A→NLL-1, drop obs-B entirely (no phase merge)
    s         → skip entire group (no update for any event in it)
    p         → print phase lines for all obs events in the group

    Parameters
    ----------
    df1              : pd.DataFrame — primary catalogue (obs); needs 'BulletinID'
                       when bulletin_lines is provided
    df2              : pd.DataFrame — reference catalogue (NLL result)
    max_dt_seconds   : float        — candidate search time window (default: 10)
    max_dist_km      : float        — candidate search distance (default: 50)
    bulletin_lines   : list[str], optional — raw bulletin file lines; enables
                       phase printing and phase merging for duplicates

    Returns
    -------
    pd.DataFrame with columns:
        Year, Month, Day, Hour, Min, Sec, Lat, Lon, Dep, Mag, MagType,
        MagAuthor, PhaseCount, HorUncer, VerUncer, AzGap, RMS,
        BulletinID, _extra_phases
    """
    update_cols = [
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'RMS', 'PhaseCount',
        'HorUncer', 'VerUncer', 'AzGap',
    ]

    for col in update_cols:
        if col not in df1.columns:
            raise ValueError(f"Column '{col}' not found in df1.")
        if col not in df2.columns:
            raise ValueError(f"Column '{col}' not found in df2.")

    if bulletin_lines is not None and 'BulletinID' not in df1.columns:
        raise ValueError("df1 must contain a 'BulletinID' column when bulletin_lines is provided.")

    d1 = df1.copy().sort_values('Time').reset_index(drop=True)
    d2 = df2.copy().sort_values('Time').reset_index(drop=True)
    d1['Time'] = pd.to_datetime(d1['Time'])
    d2['Time'] = pd.to_datetime(d2['Time'])

    d2_times = d2['Time'].values.astype('int64')

    # Pre-index phases from bulletin_lines
    phase_index: dict = {}
    if bulletin_lines is not None:
        current_id = None
        for line_idx, line in enumerate(bulletin_lines):
            if line.startswith('# '):
                current_id = line_idx
                phase_index[current_id] = []
            elif line.startswith('\n') or line == '':
                current_id = None
            elif current_id is not None:
                phase_index[current_id].append(line.rstrip('\n'))

    def to_cartesian(lat_deg, lon_deg, depth_km):
        R   = 6371.0
        lat = np.radians(lat_deg)
        lon = np.radians(lon_deg)
        r   = R - depth_km
        x   = r * np.cos(lat) * np.cos(lon)
        y   = r * np.cos(lat) * np.sin(lon)
        z   = r * np.sin(lat)
        return np.column_stack([x, y, z])

    xyz2 = to_cartesian(d2['Lat'].values, d2['Lon'].values, d2['Dep'].values)
    tree = cKDTree(xyz2)

    max_dt_ns = int(max_dt_seconds * 1e9)

    # ------------------------------------------------------------------
    # Step 1: build bipartite candidate graph
    # ------------------------------------------------------------------
    adj_d1 = defaultdict(set)   # d1_idx → {d2_idx, …}
    adj_d2 = defaultdict(set)   # d2_idx → {d1_idx, …}

    for i, row1 in d1.iterrows():
        t1_ns = row1['Time'].value
        lo = np.searchsorted(d2_times, t1_ns - max_dt_ns, side='left')
        hi = np.searchsorted(d2_times, t1_ns + max_dt_ns, side='right')
        if lo >= hi:
            continue
        xyz1   = to_cartesian(
            np.array([row1['Lat']]), np.array([row1['Lon']]), np.array([row1['Dep']])
        )[0]
        nearby = set(tree.query_ball_point(xyz1, r=max_dist_km))
        for p in range(lo, hi):
            if p in nearby:
                adj_d1[i].add(p)
                adj_d2[p].add(i)

    # ------------------------------------------------------------------
    # Step 2: connected components
    # ------------------------------------------------------------------
    visited_d1 = set()
    visited_d2 = set()
    groups     = []

    for start in sorted(adj_d1.keys()):
        if start in visited_d1:
            continue
        g_d1, g_d2 = set(), set()
        stack = [('d1', start)]
        while stack:
            kind, idx = stack.pop()
            if kind == 'd1':
                if idx in visited_d1:
                    continue
                visited_d1.add(idx)
                g_d1.add(idx)
                for j in adj_d1[idx]:
                    if j not in visited_d2:
                        stack.append(('d2', j))
            else:
                if idx in visited_d2:
                    continue
                visited_d2.add(idx)
                g_d2.add(idx)
                for i_n in adj_d2[idx]:
                    if i_n not in visited_d1:
                        stack.append(('d1', i_n))
        groups.append({'d1': sorted(g_d1), 'd2': sorted(g_d2)})

    # ------------------------------------------------------------------
    # Step 3: resolve each group
    # ------------------------------------------------------------------
    sep = '─' * 72

    def _d1_info(i):
        r = d1.iloc[i]
        return (f"Time: {r['Time']}  Lat: {r['Lat']:.4f}  "
                f"Lon: {r['Lon']:.4f}  Dep: {r['Dep']:.1f} km")

    def _d2_info(j, ref_xyz=None, ref_t_ns=None):
        r  = d2.iloc[j]
        s  = (f"Time: {r['Time']}  Lat: {r['Lat']:.4f}  "
              f"Lon: {r['Lon']:.4f}  Dep: {r['Dep']:.1f} km")
        if ref_xyz is not None and ref_t_ns is not None:
            xyz_j   = to_cartesian(
                np.array([r['Lat']]), np.array([r['Lon']]), np.array([r['Dep']])
            )[0]
            dist_km = float(np.linalg.norm(ref_xyz - xyz_j))
            dt_s    = (d2_times[j] - ref_t_ns) / 1e9
            s += f"  Δt: {dt_s:+.3f}s  Δd: {dist_km:.1f} km"
        return s

    def _print_phases(i, label):
        if not phase_index:
            print("  (no phases — bulletin_lines not provided)")
            return
        bid    = d1.iloc[i].get('BulletinID') if 'BulletinID' in d1.columns else None
        phases = phase_index.get(bid, []) if bid is not None else []
        print(f"\n  {label}  (BulletinID={bid})")
        print(sep)
        for ph in phases or ['    (no phases)']:
            print(f"    {ph}")
        print(sep)

    def _collect_extra_phases(kept_i, donor_is):
        """Return unique phase lines from all donor_is not already in kept_i."""
        all_extra = []
        seen_keys = set()
        kept_bid  = d1.iloc[kept_i].get('BulletinID') if 'BulletinID' in d1.columns else None
        if kept_bid is None or not phase_index:
            return []
        for donor_i in donor_is:
            donor_bid = d1.iloc[donor_i].get('BulletinID') if 'BulletinID' in d1.columns else None
            if donor_bid is None:
                continue
            for p in _diff_phases(kept_bid, donor_bid, phase_index):
                k = _phase_key(p)
                if k is not None and k not in seen_keys:
                    seen_keys.add(k)
                    all_extra.append(p)
        return all_extra

    # ---- inner: parse M:N command ----
    def _parse_mn(cmd_str, n1, n2, letters, letter_to_pos):
        """
        Parse a M:N assignment string.  Returns one of:
            'skip'
            'print_phases'
            {'assignments': {d1_pos: d2_pos}, 'merges': {d1_pos: d1_target|None},
             'drops': set()}
            {'error': str}
        """
        s = cmd_str.strip()
        if s.lower() in ('s', 'skip'):
            return 'skip'
        if s.lower() == 'p':
            return 'print_phases'

        tokens = [t for t in re.split(r'[\s,]+', s) if t]
        assignments = {}
        merges      = {}
        drops       = set()

        for tok in tokens:
            # Assignment: letter + number  e.g. A1, B2
            m = re.fullmatch(r'([A-Za-z])(\d+)', tok)
            if m:
                letter = m.group(1).upper()
                num    = int(m.group(2)) - 1
                if letter not in letter_to_pos:
                    return {'error': f"Unknown obs event '{letter}' in '{tok}'"}
                if not (0 <= num < n2):
                    return {'error': f"NLL candidate number out of range in '{tok}'"}
                d1pos = letter_to_pos[letter]
                if d1pos in merges or d1pos in drops:
                    return {'error': f"obs-{letter} already has a conflicting action"}
                assignments[d1pos] = num
                continue

            # Merge: letter + 'm' + optional target letter  e.g. Bm, BmA
            m = re.fullmatch(r'([A-Za-z])m([A-Za-z])?', tok, re.IGNORECASE)
            if m:
                letter = m.group(1).upper()
                target = m.group(2).upper() if m.group(2) else None
                if letter not in letter_to_pos:
                    return {'error': f"Unknown obs event '{letter}' in '{tok}'"}
                if target is not None and target not in letter_to_pos:
                    return {'error': f"Unknown merge target '{target}' in '{tok}'"}
                d1pos   = letter_to_pos[letter]
                tgt_pos = letter_to_pos[target] if target else None
                if d1pos in assignments or d1pos in drops:
                    return {'error': f"obs-{letter} already has a conflicting action"}
                merges[d1pos] = tgt_pos
                continue

            # Drop: letter + 'd'  e.g. Bd
            m = re.fullmatch(r'([A-Za-z])d', tok, re.IGNORECASE)
            if m:
                letter = m.group(1).upper()
                if letter not in letter_to_pos:
                    return {'error': f"Unknown obs event '{letter}' in '{tok}'"}
                d1pos = letter_to_pos[letter]
                if d1pos in assignments or d1pos in merges:
                    return {'error': f"obs-{letter} already has a conflicting action"}
                drops.add(d1pos)
                continue

            return {'error': f"Unrecognized token '{tok}'"}

        return {'assignments': assignments, 'merges': merges, 'drops': drops}

    def _validate_mn(parsed, n1, n2, d1_idxs, letters):
        """
        Validate and fill defaults for a parsed M:N command.
        Unspecified obs events default to merging into the temporally nearest
        assigned obs event.  Returns filled dict or {'error': str}.
        """
        assignments = dict(parsed['assignments'])
        merges      = dict(parsed['merges'])
        drops       = set(parsed['drops'])

        if not assignments:
            return {'error': "No assignment specified — use e.g. A1 to assign obs-A to NLL-1"}

        # Merge targets must point to an assigned event
        for src, tgt in merges.items():
            if tgt is not None and tgt not in assignments:
                tgt_letter = letters[tgt]
                return {'error': f"Merge target obs-{tgt_letter} is not an assigned event"}

        # Fill unspecified obs events with auto-merge to nearest assigned
        all_explicit = set(assignments) | set(merges) | drops
        for k in range(n1):
            if k in all_explicit:
                continue
            t_k     = d1.iloc[d1_idxs[k]]['Time'].value
            nearest = min(
                assignments.keys(),
                key=lambda a: abs(d1.iloc[d1_idxs[a]]['Time'].value - t_k),
            )
            merges[k] = nearest

        # Fill merge targets that are still None
        for d1pos in list(merges.keys()):
            if merges[d1pos] is not None:
                continue
            if len(assignments) == 1:
                merges[d1pos] = next(iter(assignments))
            else:
                t_k     = d1.iloc[d1_idxs[d1pos]]['Time'].value
                nearest = min(
                    assignments.keys(),
                    key=lambda a: abs(d1.iloc[d1_idxs[a]]['Time'].value - t_k),
                )
                merges[d1pos] = nearest

        return {'assignments': assignments, 'merges': merges, 'drops': drops}

    # ---- resolve loop ----
    match_records = []  # list of {'d1_idx', 'd2_idx', 'extra_phases'}

    for group in groups:
        d1_idxs = group['d1']
        d2_idxs = group['d2']
        n1, n2  = len(d1_idxs), len(d2_idxs)

        # --- 1:1 auto ---
        if n1 == 1 and n2 == 1:
            match_records.append({
                'd1_idx':      d1_idxs[0],
                'd2_idx':      d2_idxs[0],
                'extra_phases': [],
            })
            logger.info(
                f"AUTO 1:1  d1={d1_idxs[0]}  t={d1.iloc[d1_idxs[0]]['Time']}  d2={d2_idxs[0]}"
            )
            continue

        # --- 1:N: pick NLL candidate ---
        if n1 == 1:
            i     = d1_idxs[0]
            row1  = d1.iloc[i]
            t1_ns = row1['Time'].value
            xyz1  = to_cartesian(
                np.array([row1['Lat']]), np.array([row1['Lon']]), np.array([row1['Dep']])
            )[0]

            d2_ns = d2_times[d2_idxs]
            dts   = np.abs(d2_ns - t1_ns)
            order = np.argsort(dts)

            ambiguous_threshold_ns = int(0.5 * 1e9)

            # Auto-resolve if best candidate is clearly separated from the next
            is_ambiguous = (
                len(dts) > 1
                and (dts[order[0]] == dts[order[1]]
                     or (dts[order[1]] - dts[order[0]]) < ambiguous_threshold_ns)
            )

            if not is_ambiguous:
                chosen_j = d2_idxs[order[0]]
                match_records.append({'d1_idx': i, 'd2_idx': chosen_j, 'extra_phases': []})
                logger.info(
                    f"AUTO 1:N  d1={i}  chosen_d2={chosen_j}  "
                    f"dt={dts[order[0]]/1e9:.3f}s"
                )
                continue

            # Interactive 1:N prompt
            next_i = i + 1 if i + 1 < len(d1) else None

            valid_choices = {str(k) for k in range(1, n2 + 1)} | {'s'}
            if phase_index:
                valid_choices.add('p')

            def _print_1n():
                print(f"\n{sep}")
                print(f"  AMBIGUOUS MATCH  ({n2} NLL candidates)")
                print(sep)
                print(f"  obs event │ {_d1_info(i)}")
                if next_i is not None:
                    print(f"  next obs  │ {_d1_info(next_i)}  (context only)")
                print(sep)
                for k, j in enumerate(d2_idxs, start=1):
                    print(f"  Cand {k:<3} │ {_d2_info(j, ref_xyz=xyz1, ref_t_ns=t1_ns)}")
                print(sep)
                for k in range(1, n2 + 1):
                    print(f"  {k}  → assign NLL candidate {k}")
                print(f"  s  → remove (event dropped from output)")
                if phase_index:
                    print(f"  p  → show phases")
                hint = ' / '.join([str(k) for k in range(1, n2 + 1)] + ['s']
                                  + (['p'] if phase_index else []))
                print(f"  [{hint}]")

            _print_1n()

            while True:
                choice = input('  Your choice: ').strip().lower()
                if choice not in valid_choices:
                    print(f"  Invalid — enter: {', '.join(sorted(valid_choices, key=lambda x: (len(x), x)))}")
                    continue
                if choice == 'p':
                    _print_phases(i, 'obs event')
                    if next_i is not None:
                        _print_phases(next_i, 'next obs event (context)')
                    _print_1n()
                    continue
                break

            if choice == 's':
                logger.info(f"SKIPPED 1:N  d1={i}  t={d1.iloc[i]['Time']}")
                continue

            chosen_j = d2_idxs[int(choice) - 1]
            match_records.append({'d1_idx': i, 'd2_idx': chosen_j, 'extra_phases': []})
            logger.info(
                f"RESOLVED 1:N  d1={i}  chosen_d2={chosen_j}  choice={choice!r}"
            )
            continue

        # --- M:N: interactive cluster resolution ---
        letters       = [chr(ord('A') + k) for k in range(n1)]
        letter_to_pos = {ch: k for k, ch in enumerate(letters)}

        def _print_mn():
            print(f"\n{sep}")
            print(f"  CLUSTER  ({n1} obs events  ↔  {n2} NLL candidate{'s' if n2 > 1 else ''})")
            print(sep)
            for k, i_ev in enumerate(d1_idxs):
                print(f"  obs {letters[k]}  │ {_d1_info(i_ev)}")
            print(sep)
            for k, j in enumerate(d2_idxs, start=1):
                nll_t   = d2_times[j]
                nll_xyz = to_cartesian(
                    np.array([d2.iloc[j]['Lat']]),
                    np.array([d2.iloc[j]['Lon']]),
                    np.array([d2.iloc[j]['Dep']]),
                )[0]
                deltas = []
                for lk, i_ev in enumerate(d1_idxs):
                    obs_t   = d1.iloc[i_ev]['Time'].value
                    obs_xyz = to_cartesian(
                        np.array([d1.iloc[i_ev]['Lat']]),
                        np.array([d1.iloc[i_ev]['Lon']]),
                        np.array([d1.iloc[i_ev]['Dep']]),
                    )[0]
                    dt_s  = (nll_t - obs_t) / 1e9
                    dd_km = float(np.linalg.norm(nll_xyz - obs_xyz))
                    deltas.append(f"{letters[lk]}(Δt{dt_s:+.2f}s Δd{dd_km:.1f}km)")
                print(f"  NLL {k}  │ {_d2_info(j)}  [{' | '.join(deltas)}]")
            print(sep)
            print("  Commands (space or comma separated):")
            print("    A1        → assign obs-A to NLL-1")
            print("    A1,B2     → assign obs-A→NLL-1 and obs-B→NLL-2")
            print("    A1,BmA    → assign obs-A→NLL-1, merge obs-B's phases into obs-A")
            print("    A1,Bm     → same (target inferred when only one assignment)")
            print("    A1,Bd     → assign obs-A→NLL-1, drop obs-B (no phase merge)")
            print("  Unspecified obs events auto-merge into their nearest assigned event.")
            print("  s → remove entire group (all obs events dropped)  |  p → show phases")

        _print_mn()

        while True:
            cmd    = input('  Your command: ').strip()
            parsed = _parse_mn(cmd, n1, n2, letters, letter_to_pos)

            if parsed == 'skip':
                logger.info(f"SKIPPED M:N  d1={d1_idxs}  d2={d2_idxs}")
                break

            if parsed == 'print_phases':
                print(f"\n{sep}")
                for k, i_ev in enumerate(d1_idxs):
                    _print_phases(i_ev, f"obs {letters[k]}")
                _print_mn()
                continue

            if 'error' in parsed:
                print(f"  Error: {parsed['error']}")
                continue

            filled = _validate_mn(parsed, n1, n2, d1_idxs, letters)
            if 'error' in filled:
                print(f"  Error: {filled['error']}")
                continue

            assignments = filled['assignments']
            merges      = filled['merges']
            drops       = filled['drops']

            # Group merges by target: target_d1_pos → [source_d1_pos, …]
            merge_into = defaultdict(list)
            for src, tgt in merges.items():
                merge_into[tgt].append(src)

            for d1pos, d2pos in assignments.items():
                kept_i   = d1_idxs[d1pos]
                j        = d2_idxs[d2pos]
                donor_is = [d1_idxs[s] for s in merge_into.get(d1pos, [])]
                extra    = _collect_extra_phases(kept_i, donor_is)
                if extra:
                    print(f"  INFO: {len(extra)} phase(s) merged into obs-{letters[d1pos]}")
                    logger.info(
                        f"MERGED: {len(extra)} phase(s) from d1={donor_is} into d1={kept_i}"
                    )
                match_records.append({
                    'd1_idx':      kept_i,
                    'd2_idx':      j,
                    'extra_phases': extra,
                })

            for d1pos in sorted(drops):
                i_dropped = d1_idxs[d1pos]
                print(
                    f"  DROPPED obs-{letters[d1pos]}  "
                    f"(d1={i_dropped}  Time={d1.iloc[i_dropped]['Time']})"
                )
                logger.info(
                    f"DROPPED d1={i_dropped}  t={d1.iloc[i_dropped]['Time']}"
                )

            logger.info(
                f"RESOLVED M:N  d1={d1_idxs}  d2={d2_idxs}  "
                f"assignments={assignments}  merges={merges}  drops={list(drops)}"
            )
            break

    # ------------------------------------------------------------------
    # Match summary
    # ------------------------------------------------------------------
    matched_d1_set        = {rec['d1_idx'] for rec in match_records}
    matched_d2_set        = {rec['d2_idx'] for rec in match_records}
    n_no_cand             = len(d1) - len(adj_d1)
    n_with_cand_unmatched = len(adj_d1) - len(matched_d1_set)
    n_nll_unmatched       = len(d2) - len(matched_d2_set)

    summary_line = (
        f"Match summary: {len(matched_d1_set)}/{len(d1)} obs matched"
        f" | {n_no_cand} obs with no NLL candidate"
        f" | {n_with_cand_unmatched} obs with candidates but unmatched"
        f" | {n_nll_unmatched}/{len(d2)} NLL events unmatched"
    )
    logger.info(summary_line)
    print(f"\n  {summary_line}")

    # ------------------------------------------------------------------
    # Step 4: build output DataFrame
    # ------------------------------------------------------------------
    if not match_records:
        out_cols = [c for c in d1.columns if c not in update_cols] + update_cols + ['_extra_phases']
        return pd.DataFrame(columns=out_cols)

    rows = []
    for rec in match_records:
        row1 = d1.iloc[rec['d1_idx']]
        row2 = d2.iloc[rec['d2_idx']]
        out  = row1.drop(labels=update_cols).to_dict()
        for col in update_cols:
            out[col] = row2[col]
        out['_extra_phases'] = rec['extra_phases'] if rec['extra_phases'] else None
        rows.append(out)

    result = pd.DataFrame(rows).reset_index(drop=True)
    result = result[[
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'Mag', 'MagType', 'MagAuthor',
        'PhaseCount', 'HorUncer', 'VerUncer', 'AzGap', 'RMS',
        'BulletinID', '_extra_phases',
    ]]
    return result


def save_bulletin(parameters, log_dir=None):
    """
    Match the pre-NLL .obs bulletin to the NLL result file and write the
    updated bulletin to disk.

    Parameters
    ----------
    parameters : MatchCatalogsParams
    log_dir    : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: output, log, n_matched
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Obs file     : {parameters.file_obs}")
    logger.info(f"Final file   : {parameters.file_final}")
    logger.info(f"Output       : {parameters.save_file}")

    obs_df, lines = _read_obs(parameters.file_obs)
    final_df      = _read_final(parameters.file_final)

    matched = match_catalogues(
        obs_df, final_df,
        max_dt_seconds = 5.0,
        max_dist_km    = 50.0,
        bulletin_lines = lines,
    )

    updated_bulletin = _update_bulletin(lines, matched)

    with open(parameters.save_file, 'w') as f:
        f.writelines(updated_bulletin)

    logger.info(f"Matched      : {len(matched)} events")

    return {
        'output':    parameters.save_file,
        'log':       log_path,
        'n_matched': len(matched),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Match pre-NLL .obs events to NLL-relocated events and write the updated bulletin.'
    )
    parser.add_argument('--obs',    required=True, help='Pre-relocation GLOBAL.obs bulletin')
    parser.add_argument('--final',  required=True, help='NLL result file (FINAL.txt)')
    parser.add_argument('--output', required=True, help='Output updated bulletin file')
    args = parser.parse_args()

    save_bulletin(MatchCatalogsParams(
        file_obs   = args.obs,
        file_final = args.final,
        save_file  = args.output,
    ))


if __name__ == '__main__':
    main()
