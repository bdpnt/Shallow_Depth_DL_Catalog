"""
convert_picks.py
============================
Pick file format converter for the Shallow_Depth_DL_Catalog project.

Reads a pick file from temp_picks/pick_files/, strips all non-pick lines (event
headers, blanks), maps station codes to the project's internal codes using
stations/GLOBAL_code_map.txt, and writes pick-only output in the obs/GLOBAL.obs
pick line format.

Usage
-----
    python temp_picks/convert_picks.py --input temp_picks/pick_files/viehla_final.obs --format TEMP_OBS
    python temp_picks/convert_picks.py --input temp_picks/pick_files/merged_pyrenees.txt --format TEMP_RSB
    python temp_picks/convert_picks.py --input temp_picks/pick_files/merged_omp.csv --format TEMP_OMP
    python temp_picks/convert_picks.py --input temp_picks/pick_files/merged_other.txt --format TEMP_OTH
    python temp_picks/convert_picks.py --input temp_picks/all_picks/PICKS_MARC/OMP-picks.pq \
        --format TEMP_STB --output temp_picks/pick_files/OMP-picks_converted.obs

Supported formats
-----------------
    TEMP_OBS : OBS-style .obs file with short station names and floating-point
               year/month/day fields in event headers. Used by Viehla and similar
               generated pick files.
    TEMP_RSB : RaspberryShake/PhaseNet pick files. One pick per line with format:
               NETWORK.STATION.LOCATION PHASE ISO8601_TIMESTAMP prob=PROBABILITY
    TEMP_OMP : OMP/PhaseNet CSV files produced by merge_omp_picks.py. Columns used:
               station_id (fields 0-3), phase_time, phase_type.
    TEMP_OTH : fixed-width quality-coded pick files. One line per P pick, with an
               optional S pick stored as an offset in seconds after the P arrival.
    TEMP_STB : Strasbourg (RENASS/OMP) pick datasets — a directory of hive-partitioned
               parquet files. Columns used: station_id ("NETWORK.STATION"), phase_type,
               phase_time, phase_score. Rows with phase_score below --min-phase-score
               (default DEFAULT_MIN_PHASE_SCORE) are dropped. --input must point to the
               dataset's root directory; --output should be given explicitly, since the
               default naming would otherwise land inside the gitignored source directory
               instead of temp_picks/pick_files/.

Adding a new format
-------------------
    1. Write a converter function: convert_<format>(line, code_map) -> str | list[str] | None
       It receives a single non-header, non-blank pick line and the station code map.
       It must return the converted GLOBAL.obs pick line string, a list of such
       strings (for formats that pack multiple picks into one source line), or
       None/an empty list to skip.
    2. Register it in FORMAT_HANDLERS with a descriptive key string.
    3. Most formats are line-based (see step 1). A format may instead be row-based /
       directory-based (source is a directory of files read into rows, e.g. parquet):
       the handler then takes a row (a DataFrame.itertuples() namedtuple, accessed by
       attribute name — column order from pd.read_parquet(columns=[...]) is not
       guaranteed) instead of a line string, and the format key must also be added to
       ROW_BASED_FORMATS so convert_file() dispatches it through the parquet-glob/
       streaming-write path instead of the per-line path. See TEMP_STB /
       convert_temp_stb / _convert_parquet_dir for the reference implementation.
"""

import argparse
import glob
import logging
import os
import re
from collections import Counter
from datetime import datetime, timedelta, timezone

import pandas as pd

logger = logging.getLogger('convert_picks')

DEFAULT_MIN_PHASE_SCORE = 0.5     # mirrors merge_omp_picks.py / merge_pyrenees_picks.py
ROW_BASED_FORMATS = {'TEMP_STB'}  # formats read from a directory of row-oriented files
                                  # (parquet) instead of opened as a text file line-by-line


# ---------------------------------------------------------------------------
# Station code map
# ---------------------------------------------------------------------------

def load_code_map(codemap_path):
    """
    Parse GLOBAL_code_map.txt into a station lookup dict.

    The file maps project-internal codes (Alternate Code) to canonical
    NETWORK.CODE names (Station Code). This function inverts the lookup:
    given a short station name (the part after the '.' in the canonical code),
    return the list of matching project-internal codes with their validity windows.

    Parameters
    ----------
    codemap_path : str
        Path to stations/GLOBAL_code_map.txt.

    Returns
    -------
    dict[str, list[dict]]
        Keys are short station names (e.g. 'PYLU').
        Values are lists of dicts with keys:
            'internal_code' : str  — project internal code (e.g. 'RA.0012')
            'start'         : datetime or None
            'end'           : datetime or None
    """
    code_map = {}

    with open(codemap_path, 'r') as f:
        lines = f.readlines()

    def _save(internal_code, canonical_code, start_dt, end_dt):
        if not canonical_code:
            return
        short_name = canonical_code.split('.')[-1]
        code_map.setdefault(short_name, []).append({
            'internal_code': internal_code,
            'start':         start_dt,
            'end':           end_dt,
        })

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith('Alternate Code:'):
            i += 1
            continue

        internal_code  = line.split(':', 1)[1].strip()
        canonical_code = None
        start_dt       = None
        end_dt         = None
        i += 1

        while i < len(lines):
            sub = lines[i].strip()
            if sub.startswith('Alternate Code:'):
                break
            if sub.startswith('Station Code:'):
                _save(internal_code, canonical_code, start_dt, end_dt)
                canonical_code = sub.split(':', 1)[1].strip()
                start_dt = None
                end_dt   = None
            elif sub.startswith('Start Date:'):
                start_dt = _parse_map_date(sub.split(':', 1)[1].strip())
            elif sub.startswith('End Date:'):
                end_dt = _parse_map_date(sub.split(':', 1)[1].strip())
            i += 1

        _save(internal_code, canonical_code, start_dt, end_dt)

    return code_map


def _parse_map_date(date_str):
    """Parse ISO date string from GLOBAL_code_map.txt into a timezone-aware datetime."""
    try:
        return datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S.%fZ').replace(tzinfo=timezone.utc)
    except ValueError:
        return None


def resolve_station(short_name, pick_date_str, code_map, fallback_counter=None):
    """
    Return the project internal station code for a short station name.

    Uses the pick date to select the correct entry when multiple time windows
    exist for the same station. Falls back to the last entry if no date-matched
    entry is found; increments fallback_counter[short_name] when provided.

    Parameters
    ----------
    short_name : str
        Short station name from the source pick file (e.g. 'PYLU').
    pick_date_str : str
        Pick date in YYYYMMDD format (e.g. '20201204').
    code_map : dict
        Lookup dict from load_code_map().
    fallback_counter : Counter, optional
        If provided, incremented when the date-window match fails and the
        most-recent-entry fallback is used. Caller logs the summary.

    Returns
    -------
    str or None
        Project internal code (e.g. 'RA.0012'), or None if station is unknown.
    """
    entries = code_map.get(short_name)
    if not entries:
        prefix = short_name + "_"
        entries = [e for k, v in code_map.items() if k.startswith(prefix) for e in v]
    if not entries:
        return None

    try:
        pick_dt = datetime.strptime(pick_date_str, '%Y%m%d').replace(tzinfo=timezone.utc)
    except ValueError:
        pick_dt = None

    if pick_dt is not None:
        # Priority 1: both bounds defined and pick_dt within window
        for entry in entries:
            start, end = entry['start'], entry['end']
            if start is not None and end is not None and start <= pick_dt <= end:
                return entry['internal_code']
        # Priority 2: exactly one bound defined, constraint satisfied
        for entry in entries:
            start, end = entry['start'], entry['end']
            if (start is None) != (end is None):
                if (start is None or pick_dt >= start) and (end is None or pick_dt <= end):
                    return entry['internal_code']
        # Priority 3: no bounds at all — open-ended catch-all
        for entry in entries:
            if entry['start'] is None and entry['end'] is None:
                return entry['internal_code']

    if fallback_counter is not None:
        fallback_counter[short_name] += 1
    return entries[-1]['internal_code']


# ---------------------------------------------------------------------------
# Output formatting helper
# ---------------------------------------------------------------------------

def _format_pick_line(internal_code, phase, date, hhmm, seconds_str, error_str, pick_origin, instrument='?'):
    """
    Format a single pick line in GLOBAL.obs style.

    Field widths match the format used throughout the fetch_obs/ modules.

    Parameters
    ----------
    internal_code : str   — project internal station code (e.g. 'RA.0012')
    phase         : str   — 'P' or 'S'
    date          : str   — YYYYMMDD
    hhmm          : str   — HHMM (4 chars, zero-padded)
    seconds_str   : str   — arrival seconds as string (e.g. '29.8640')
    error_str     : str   — pick uncertainty in seconds (e.g. '0.05e+00')
    pick_origin   : str   — source label appended in trailer (e.g. 'TEMP_OBS')
    instrument    : str   — instrument field; '*' signals NLLoc to use S-P
                            relative timing instead of the absolute pick time
                            (mirrors fetch_obs/OMP.py's quality-9 convention)
    """
    code        = internal_code.ljust(9)
    instrument  = instrument.ljust(4)
    component   = '?'.ljust(4)
    onset       = '?'.ljust(1)
    phase_field = phase.ljust(6)
    direction   = '?'.ljust(1)
    error_type  = 'GAU'.ljust(3)
    seconds     = f"{float(seconds_str):06.3f}"
    error_mag   = f"{float(error_str):.2f}".ljust(9)
    coda_dur    = '-1.00e+00'.ljust(9)
    amp         = '-1.00e+00'.ljust(9)
    period      = '-1.00e+00'.ljust(9)
    real_phase  = phase.ljust(6)
    channel     = 'None'.ljust(4)
    origin      = pick_origin.ljust(9)
    pgv         = 'None'.ljust(4)

    return (
        f"{code} {instrument} {component} {onset} {phase_field} {direction} "
        f"{date} {hhmm} {seconds} {error_type} {error_mag} {coda_dur} {amp} {period}"
        f" # {real_phase} {channel} {origin} {pgv}"
    )


# ---------------------------------------------------------------------------
# Format handlers
# ---------------------------------------------------------------------------

def convert_temp_obs(line, code_map, skipped_stations=None, fallback_counter=None):
    """
    Convert a single pick line from TEMP_OBS format to GLOBAL.obs format.

    TEMP_OBS source format (space-delimited):
        STATION ? ? ? PHASE ? YYYYMMDD HHMM SS.SSSS GAU E.EEe+EE -1.00e+00 -1.00e+00 -1.00e+00

    Returns the converted line string, or None if the station cannot be resolved.
    """
    parts = line.split()
    if len(parts) < 11:
        return None

    short_name  = parts[0]
    phase       = parts[4]
    date        = parts[6]   # YYYYMMDD
    hhmm        = parts[7]   # HHMM
    seconds_str = parts[8]   # SS.SSSS
    error_str   = parts[10]  # E.EEe+EE (e.g. 0.05e+00)

    internal_code = resolve_station(short_name, date, code_map, fallback_counter)
    if internal_code is None:
        if skipped_stations is not None:
            skipped_stations[short_name] += 1
        return None

    return _format_pick_line(internal_code, phase, date, hhmm, seconds_str, error_str, 'TEMP_OBS')


def convert_temp_rsb(line, code_map, skipped_stations=None, fallback_counter=None):
    """
    Convert a single pick line from RaspberryShake/PhaseNet format to GLOBAL.obs format.

    TEMP_RSB source format (space-delimited):
        NETWORK.STATION.LOCATION PHASE ISO8601_TIMESTAMP prob=PROBABILITY

    Pick uncertainty: 0.05 s for P, 0.15 s for S (PhaseNet probability
    is a detection confidence, not a timing error estimate).
    """
    parts = line.split()
    if len(parts) < 3:
        return None

    station_field = parts[0]
    phase         = parts[1]
    timestamp     = parts[2]

    station_parts = station_field.split('.')
    if len(station_parts) < 2:
        return None
    short_name = station_parts[1]

    try:
        dt = datetime.strptime(timestamp, '%Y-%m-%dT%H:%M:%S.%f')
    except ValueError:
        try:
            dt = datetime.strptime(timestamp, '%Y-%m-%dT%H:%M:%S')
        except ValueError:
            logger.warning(f"Cannot parse timestamp '{timestamp}'. Skipping.")
            return None

    date        = dt.strftime('%Y%m%d')
    hhmm        = dt.strftime('%H%M')
    seconds_str = f"{dt.second + dt.microsecond / 1e6:.3f}"
    error_str   = '0.05' if phase == 'P' else '0.15'

    internal_code = resolve_station(short_name, date, code_map, fallback_counter)
    if internal_code is None:
        if skipped_stations is not None:
            skipped_stations[short_name] += 1
        return None

    return _format_pick_line(internal_code, phase, date, hhmm, seconds_str, error_str, 'TEMP_RSB')


def convert_temp_omp(line, code_map, skipped_stations=None, fallback_counter=None):
    """
    Convert a single pick line from OMP/PhaseNet CSV format to GLOBAL.obs format.

    TEMP_OMP source format (CSV columns):
        file_name, begin_time, station_id, phase_index, phase_time,
        phase_score, phase_ampl, phase_type

    Station short name is field [1] of station_id (e.g. 'ARBS' from
    'CA.ARBS.00.HHX.D.2017.100'). Pick uncertainty: 0.05 s for P, 0.15 s for S.
    """
    if line.startswith('file_name,'):
        return None

    parts = line.split(',')
    if len(parts) < 8:
        return None

    station_id = parts[2]
    phase_time = parts[4]
    phase_type = parts[7].strip()

    station_parts = station_id.split('.')
    if len(station_parts) < 2:
        return None
    short_name = station_parts[1]

    try:
        dt = datetime.strptime(phase_time, '%Y-%m-%dT%H:%M:%S.%f')
    except ValueError:
        try:
            dt = datetime.strptime(phase_time, '%Y-%m-%dT%H:%M:%S')
        except ValueError:
            logger.warning(f"Cannot parse timestamp '{phase_time}'. Skipping.")
            return None

    date        = dt.strftime('%Y%m%d')
    hhmm        = dt.strftime('%H%M')
    seconds_str = f"{dt.second + dt.microsecond / 1e6:.3f}"
    error_str   = '0.05' if phase_type == 'P' else '0.15'

    internal_code = resolve_station(short_name, date, code_map, fallback_counter)
    if internal_code is None:
        if skipped_stations is not None:
            skipped_stations[short_name] += 1
        return None

    return _format_pick_line(internal_code, phase_type, date, hhmm, seconds_str, error_str, 'TEMP_OMP')


_TEMP_OTH_RE   = re.compile(r'^(.{5})P (\d) (.{15})(.*)$')
_TEMP_OTH_S_RE = re.compile(r'^\s*(\S+)\s+S\s+(\d)\s*$')


def convert_temp_oth(line, code_map, skipped_stations=None, fallback_counter=None):
    """
    Convert a single pick line from TEMP_OTH format to GLOBAL.obs format(s).

    TEMP_OTH source format (fixed-width, no delimiters):
        STATION(5) "P X "(4) YYMMDDHHMMSS.SS(15) [S_OFFSET " S " QUALITY]

    Station is a 5-char space-padded short name. "P X " holds the P-phase
    quality digit (0-9; >=4 means unusable, per project convention). The
    15-char date/time field packs YY MM DD HH MI SS.SS as six 2-char
    sub-fields, each right-justified with a space instead of a leading zero
    (e.g. " 6" means 06). The optional trailing part gives the S arrival as
    an offset in seconds after the P arrival, followed by " S " and the
    S-phase quality digit.

    P-phase quality 4 lines store an unusable placeholder time (blank or
    "0.00") instead of a real reading; since the S offset is only meaningful
    relative to a real P arrival, such lines are dropped entirely (mirrors
    fetch_obs/OMP.py's handling of missing P seconds). P-phase quality 9
    lines do carry a real timestamp but flag it as absolute-time-unreliable;
    these are kept, with instrument '*' instead of '?' on both the P and S
    lines, so NLLoc uses S-P relative timing instead of the absolute pick
    time (mirrors fetch_obs/OMP.py's quality_p == '9' handling).

    Returns a list of 0, 1, or 2 converted line strings (P pick, S pick,
    both, or neither, depending on quality filtering), or None if the line
    doesn't match the TEMP_OTH shape or the station can't be resolved.
    """
    m = _TEMP_OTH_RE.match(line)
    if not m:
        return None

    station_field, quality_p_str, datefield, rest = m.groups()
    short_name = station_field.strip()
    quality_p  = int(quality_p_str)

    yy = int(datefield[0:2])
    mm = int(datefield[2:4])
    dd = int(datefield[4:6])
    hh = int(datefield[6:8])
    mi = int(datefield[8:10])
    ss_str   = datefield[10:12].strip()
    frac_str = datefield[13:15].strip()
    if not ss_str or not frac_str:
        return None  # P time unusable -> S offset can't be resolved either

    year = 2000 + yy if yy < 78 else 1900 + yy
    p_dt = datetime(year, mm, dd, hh, mi, int(ss_str)) + timedelta(microseconds=int(frac_str) * 10000)

    internal_code = resolve_station(short_name, p_dt.strftime('%Y%m%d'), code_map, fallback_counter)
    if internal_code is None:
        if skipped_stations is not None:
            skipped_stations[short_name] += 1
        return None

    results = []
    marker  = '*' if quality_p == 9 else '?'
    if quality_p < 4 or quality_p == 9:
        err_p = '0.05' if quality_p == 9 else ('0.10' if quality_p > 1 else '0.05')
        results.append(_format_pick_line(
            internal_code, 'P', p_dt.strftime('%Y%m%d'), p_dt.strftime('%H%M'),
            f"{p_dt.second + p_dt.microsecond / 1e6:.3f}", err_p, 'TEMP_OTH', instrument=marker))

    if rest.strip():
        m2 = _TEMP_OTH_S_RE.match(rest)
        if m2:
            offset_str, quality_s_str = m2.groups()
            quality_s = int(quality_s_str)
            if quality_s < 4:
                s_dt  = p_dt + timedelta(seconds=float(offset_str))
                err_s = '0.30' if quality_s > 2 else '0.15'
                results.append(_format_pick_line(
                    internal_code, 'S', s_dt.strftime('%Y%m%d'), s_dt.strftime('%H%M'),
                    f"{s_dt.second + s_dt.microsecond / 1e6:.3f}", err_s, 'TEMP_OTH', instrument=marker))

    return results


def convert_temp_stb(row, code_map, skipped_stations=None, fallback_counter=None, instrument='?'):
    """
    Convert a single pick row from Strasbourg (RENASS/OMP) parquet format to
    GLOBAL.obs format.

    row is a DataFrame.itertuples() namedtuple with attributes station_id
    ("NETWORK.STATION"), phase_type ('P'/'S'), phase_time (Timestamp, already
    filtered by min_phase_score upstream). Access by attribute name, not
    position — pd.read_parquet(columns=[...]) does not guarantee column order
    matches the requested list.

    Pick uncertainty: 0.05 s for P, 0.15 s for S (matches convert_temp_rsb /
    convert_temp_omp — phase_score is a detection confidence, already
    threshold-filtered, not a timing-error estimate).

    instrument is '*' for sources whose absolute timing is unreliable (short-
    period OMP stations), so NLLoc only uses them as S-P pairs.
    """
    station_parts = row.station_id.split('.')
    if len(station_parts) < 2:
        return None
    short_name = station_parts[1]

    dt = row.phase_time
    if pd.isna(dt):
        return None

    phase       = row.phase_type
    date        = dt.strftime('%Y%m%d')
    hhmm        = dt.strftime('%H%M')
    seconds_str = f"{dt.second + dt.microsecond / 1e6:.3f}"
    error_str   = '0.05' if phase == 'P' else '0.15'

    internal_code = resolve_station(short_name, date, code_map, fallback_counter)
    if internal_code is None:
        if skipped_stations is not None:
            skipped_stations[short_name] += 1
        return None

    return _format_pick_line(internal_code, phase, date, hhmm, seconds_str, error_str, 'TEMP_STB',
                             instrument=instrument)


# ---------------------------------------------------------------------------
# Format dispatch table — register new format handlers here
# ---------------------------------------------------------------------------

FORMAT_HANDLERS = {
    'TEMP_OBS': convert_temp_obs,
    'TEMP_RSB': convert_temp_rsb,
    'TEMP_OMP': convert_temp_omp,
    'TEMP_OTH': convert_temp_oth,
    'TEMP_STB': convert_temp_stb,
}


# ---------------------------------------------------------------------------
# Row-based (parquet) reader — used for formats in ROW_BASED_FORMATS
# ---------------------------------------------------------------------------

def _convert_parquet_dir(input_dir, output_path, fmt_handler, code_map,
                          min_phase_score, skipped_stations, fallback_counter, instrument='?'):
    """
    Read every *.parquet file under input_dir, filter by phase_score, convert
    row-by-row, and stream results straight to output_path — one parquet file
    at a time, so peak memory stays bounded regardless of total dataset size
    (pyrope-renass74-bp4_40.pq alone is ~72.5M rows).

    Reads each leaf file individually (never pd.read_parquet(input_dir)):
    some datasets (e.g. OMP-picks.pq) embed year/month as real columns AND as
    hive-partition directories, which pandas can't merge across files. Since
    phase_time already carries the full date, year/month are never requested.
    """
    parquet_files = sorted(glob.glob(os.path.join(input_dir, '**', '*.parquet'), recursive=True))
    n_input = n_converted = n_skipped = n_dropped_low_score = 0

    with open(output_path, 'w') as out_f:
        for fp in parquet_files:
            df = pd.read_parquet(fp, columns=['station_id', 'phase_type', 'phase_time', 'phase_score'])
            n_input += len(df)

            keep = df['phase_score'].notna() & (df['phase_score'] >= min_phase_score)
            n_dropped_low_score += int((~keep).sum())
            df = df[keep]

            for row in df.itertuples(index=False):
                result = fmt_handler(row, code_map, skipped_stations, fallback_counter, instrument=instrument)
                if result is None:
                    n_skipped += 1
                else:
                    out_f.write(result + '\n')
                    n_converted += 1

    return n_input, n_converted, n_skipped, n_dropped_low_score


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_CODEMAP = os.path.join(_MODULE_DIR, '..', 'stations', 'GLOBAL_code_map.txt')
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')


def _setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)

    basename   = os.path.splitext(os.path.basename(__file__))[0]
    timestamp  = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_path   = os.path.join(log_dir, f"{basename}_{timestamp}.log")

    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    handler = logging.FileHandler(log_path, encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(handler)

    return log_path


def convert_file(input_path, fmt, output_path=None, codemap_path=None, log_dir=None, min_phase_score=None,
                 relative_timing=False):
    """
    Convert a pick file to GLOBAL.obs pick line format.

    Parameters
    ----------
    input_path : str
        Path to the source pick file (e.g. 'temp_picks/pick_files/viehla_final.obs'),
        or — for formats in ROW_BASED_FORMATS (e.g. 'TEMP_STB') — a directory of
        parquet files.
    fmt : str
        Source format type. Must be a key in FORMAT_HANDLERS (e.g. 'TEMP_OBS').
    output_path : str, optional
        Destination file path. Defaults to input path with '_converted.obs' suffix.
    codemap_path : str, optional
        Path to GLOBAL_code_map.txt. Defaults to stations/GLOBAL_code_map.txt
        relative to this module's location.
    log_dir : str, optional
        Directory where the log file is written. Defaults to temp_picks/console_output/.
    min_phase_score : float, optional
        Minimum phase_score required to keep a pick. Only used for formats in
        ROW_BASED_FORMATS (e.g. 'TEMP_STB'); ignored otherwise. Defaults to
        DEFAULT_MIN_PHASE_SCORE.
    relative_timing : bool, optional
        Write every pick with instrument '*' so NLLoc uses it only through S-P
        differences at the same station. Only used for formats in
        ROW_BASED_FORMATS; ignored otherwise. Default False.

    Returns
    -------
    dict
        Summary with keys: 'output', 'log', 'n_input', 'n_converted', 'n_skipped'.

    Raises
    ------
    ValueError
        If fmt is not a supported format.

    Examples
    --------
    >>> from temp_picks.convert_picks import convert_file
    >>> result = convert_file('temp_picks/pick_files/viehla_final.obs', 'TEMP_OBS')
    >>> print(result)
    """
    if fmt not in FORMAT_HANDLERS:
        raise ValueError(f"Unknown format '{fmt}'. Supported: {', '.join(FORMAT_HANDLERS)}")

    codemap_path = codemap_path or _DEFAULT_CODEMAP
    log_path     = _setup_logger(log_dir or _DEFAULT_LOG_DIR)

    if output_path is None:
        base, _ = os.path.splitext(input_path)
        output_path = base + '_converted.obs'

    logger.info(f"Converting file: {input_path}")
    logger.info(f"Loading station code map: {codemap_path}")
    code_map = load_code_map(codemap_path)
    n_windows = sum(len(v) for v in code_map.values())
    logger.info(f"Code map loaded: {len(code_map)} station names, {n_windows} validity windows.")

    fmt_handler      = FORMAT_HANDLERS[fmt]
    skipped_stations = Counter()
    fallback_counter = Counter()

    if fmt in ROW_BASED_FORMATS:
        score_threshold = min_phase_score if min_phase_score is not None else DEFAULT_MIN_PHASE_SCORE
        logger.info(f"Min phase score  : {score_threshold}")
        logger.info(f"Relative timing  : {relative_timing}")
        n_input, n_converted, n_skipped, n_dropped_low_score = _convert_parquet_dir(
            input_path, output_path, fmt_handler, code_map,
            score_threshold, skipped_stations, fallback_counter,
            instrument='*' if relative_timing else '?',
        )
    else:
        converted = []
        n_input   = 0
        n_skipped = 0

        with open(input_path, 'r') as f:
            for raw_line in f:
                line = raw_line.rstrip('\n')
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                n_input += 1
                result = fmt_handler(line, code_map, skipped_stations, fallback_counter)
                if result is None:
                    n_skipped += 1
                elif isinstance(result, list):
                    if result:
                        converted.extend(result)
                    else:
                        n_skipped += 1
                else:
                    converted.append(result)

        with open(output_path, 'w') as f:
            for line in converted:
                f.write(line + '\n')

        n_converted         = len(converted)
        n_dropped_low_score = 0

    logger.info(f"Input pick lines : {n_input}")
    logger.info(f"Converted        : {n_converted}")
    logger.info(f"Skipped          : {n_skipped}")
    if fmt in ROW_BASED_FORMATS:
        logger.info(f"Rows dropped (low phase_score) : {n_dropped_low_score}")
    logger.info(f"Output           : {output_path}")
    if skipped_stations:
        summary = ', '.join(f"{s} ({n})" for s, n in sorted(skipped_stations.items()))
        logger.warning(f"Stations not found in code map ({len(skipped_stations)} unique): {summary}")
    if fallback_counter:
        summary = ', '.join(f"{s} ({n})" for s, n in sorted(fallback_counter.items()))
        logger.warning(f"Stations resolved with date fallback ({len(fallback_counter)} unique): {summary}")

    return {
        'output':      output_path,
        'log':         log_path,
        'n_input':     n_input,
        'n_converted': n_converted,
        'n_skipped':   n_skipped,
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert pick files to GLOBAL.obs pick line format.'
    )
    parser.add_argument(
        '--input', required=True,
        help='Path to the input pick file (e.g. Pyrocko/pick_files/viehla_final.obs)'
    )
    parser.add_argument(
        '--format', required=True, choices=list(FORMAT_HANDLERS.keys()),
        help=f'Source format type. Supported: {", ".join(FORMAT_HANDLERS.keys())}'
    )
    parser.add_argument(
        '--output', default=None,
        help='Output file path. Default: input path with _converted.obs suffix.'
    )
    parser.add_argument(
        '--codemap', default=None,
        help='Path to GLOBAL_code_map.txt. Default: stations/GLOBAL_code_map.txt relative to project root.'
    )
    parser.add_argument(
        '--log-dir', default=None,
        help='Directory for log files. Default: temp_picks/console_output/.'
    )
    parser.add_argument(
        '--min-phase-score', type=float, default=None, metavar='SCORE',
        help=f'Minimum phase_score required to keep a pick (TEMP_STB only). Default: {DEFAULT_MIN_PHASE_SCORE}.'
    )
    args = parser.parse_args()
    convert_file(
        args.input, args.format, output_path=args.output, codemap_path=args.codemap,
        log_dir=args.log_dir, min_phase_score=args.min_phase_score,
    )


if __name__ == '__main__':
    main()
