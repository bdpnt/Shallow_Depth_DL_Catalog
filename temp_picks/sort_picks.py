"""
sort_picks.py
============================
Sort pick lines within each event of a GLOBAL.obs bulletin by ascending arrival time.

All event headers, pick associations, and PhaseCount values are preserved unchanged.
Only the order of pick lines within each event block is modified.

Usage
-----
    # In-place
    python temp_picks/sort_picks.py --input obs/GLOBAL.obs

    # To a separate file
    python temp_picks/sort_picks.py --input obs/GLOBAL.obs --output obs/GLOBAL_sorted.obs
"""

import argparse
from datetime import datetime, timezone


# ---------------------------------------------------------------------------
# Bulletin I/O (mirrors load_bulletin in temp_picks/match_picks.py)
# ---------------------------------------------------------------------------

def _parse_arrival_dt(line):
    """
    Return the arrival datetime for a pick line, or None if the line cannot
    be parsed (malformed lines are sorted to the end of the event).
    """
    parts = line.split()
    if len(parts) < 9:
        return None
    try:
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
        return datetime(year, month, day, hour, minute, sec_int, microsec, tzinfo=timezone.utc)
    except (ValueError, IndexError):
        return None


def _sort_key(line):
    dt = _parse_arrival_dt(line)
    # Unparseable lines sort to the end
    return dt if dt is not None else datetime.max.replace(tzinfo=timezone.utc)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def sort_picks(input_path, output_path=None):
    """
    Sort pick lines within each event by ascending arrival time.

    Parameters
    ----------
    input_path  : str  — path to the input bulletin (GLOBAL.obs format)
    output_path : str, optional  — destination path; defaults to input_path (in-place)

    Returns
    -------
    dict with keys: 'output', 'n_events', 'n_picks'
    """
    if output_path is None:
        output_path = input_path

    with open(input_path, 'r') as f:
        lines = f.readlines()

    output    = []
    n_events  = 0
    n_picks   = 0
    i         = 0

    # Preserve file header (### lines + following blank line)
    while i < len(lines) and lines[i].startswith('###'):
        output.append(lines[i])
        i += 1
    while i < len(lines) and lines[i].strip() == '':
        output.append(lines[i])
        i += 1

    # Process event blocks
    while i < len(lines):
        line = lines[i]
        if line.startswith('# '):
            n_events += 1
            output.append(line)  # event header
            i += 1

            # Collect all pick lines for this event; keep PUBLIC_ID separate
            pick_lines    = []
            public_id_line = None
            while i < len(lines) and lines[i].strip() != '':
                raw = lines[i].rstrip('\n')
                if lines[i].startswith('PUBLIC_ID'):
                    public_id_line = raw
                else:
                    pick_lines.append(raw)
                i += 1

            # PUBLIC_ID must stay immediately after the header
            if public_id_line is not None:
                output.append(public_id_line + '\n')

            # Sort by arrival time and write
            pick_lines.sort(key=_sort_key)
            for pick in pick_lines:
                output.append(pick + '\n')
                n_picks += 1

            output.append('\n')  # blank line separator
        i += 1  # advance past blank line (or unrecognised line)

    with open(output_path, 'w') as f:
        f.writelines(output)

    return {'output': output_path, 'n_events': n_events, 'n_picks': n_picks}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Sort pick lines within each bulletin event by ascending arrival time.'
    )
    parser.add_argument('--input',  required=True, help='Input bulletin file (GLOBAL.obs format)')
    parser.add_argument('--output', default=None,  help='Output path (default: same as input)')
    args   = parser.parse_args()
    result = sort_picks(args.input, args.output)
    print(f"Sorted {result['n_picks']} picks across {result['n_events']} events → {result['output']}")


if __name__ == '__main__':
    main()
