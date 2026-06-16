"""
merge_omp_picks.py
============================
Merge OMP/PhaseNet pick CSVs into a single consolidated file.

Reads all PICKS_*.csv files from all yearly subdirectories under picks_OMP/,
skipping stations listed in STATIONS_TO_DROP, and concatenates them into one
merged CSV with a single header row. No format conversion is done here; use
convert_picks.py with --format TEMP_OMP to convert to GLOBAL.obs format.

Usage
-----
    python temp_picks/merge_omp_picks.py

    # Override defaults
    python temp_picks/merge_omp_picks.py --input-dir temp_picks/all_picks/PICKS_PHASENET_TOUS/picks_OMP 
        --output temp_picks/pick_files/merged_omp.csv --drop-years 2026 2027
"""

import argparse
import logging
import os
from datetime import datetime

_MODULE_DIR = os.path.dirname(os.path.abspath(__file__))

_DEFAULT_INPUT_DIR = os.path.join(_MODULE_DIR, 'all_picks', 'PICKS_PHASENET_TOUS', 'picks_OMP')
_DEFAULT_OUTPUT    = os.path.join(_MODULE_DIR, 'pick_files', 'merged_omp.csv')
_DEFAULT_LOG_DIR   = os.path.join(_MODULE_DIR, 'console_output')

# Add station codes here to exclude them from the merged output (permanent exclusions)
STATIONS_TO_DROP = {'SMC'}

# Add years (as strings) here to exclude them from the merged output (permanent exclusions)
YEARS_TO_DROP = {'2026'}

logger = logging.getLogger('merge_omp_picks')


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _setup_logger(log_dir):
    os.makedirs(log_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    basename  = os.path.splitext(os.path.basename(__file__))[0]
    log_path  = os.path.join(log_dir, f"{basename}_{timestamp}.log")

    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    handler = logging.FileHandler(log_path, encoding='utf-8')
    handler.setFormatter(logging.Formatter('%(levelname)s %(message)s'))
    logger.addHandler(handler)

    return log_path


def _station_code_from_filename(filename):
    """Extract station code from PICKS_{CODE}.csv → '{CODE}'."""
    basename = os.path.basename(filename)
    return basename.removeprefix('PICKS_').removesuffix('.csv')


def _year_from_dirname(dirname):
    """Extract year from PICKS_PHASENET_{YEAR}_ALL → '{YEAR}'."""
    parts = dirname.split('_')
    return parts[2] if len(parts) >= 3 else None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def merge_omp(input_dir=None, output_path=None, log_dir=None,
              stations_to_drop=None, years_to_drop=None):
    """
    Merge all OMP PhaseNet pick CSVs into a single file.

    Iterates all yearly subdirectories in input_dir, collects every
    PICKS_*.csv file, skips stations and years in the drop sets, and writes
    one merged CSV with a single header row.

    Parameters
    ----------
    input_dir : str, optional
        Directory containing yearly subdirs (PICKS_PHASENET_*).
        Defaults to temp_picks/all_picks/PICKS_PHASENET_TOUS/picks_OMP/.
    output_path : str, optional
        Destination CSV file. Defaults to temp_picks/pick_files/merged_omp.csv.
    log_dir : str, optional
        Directory for the log file. Defaults to temp_picks/console_output/.
    stations_to_drop : set[str], optional
        Station codes to exclude (e.g. {'SMC'}). Defaults to STATIONS_TO_DROP.
    years_to_drop : set[str], optional
        Years to exclude as strings (e.g. {'2020', '2021'}). Defaults to YEARS_TO_DROP.

    Returns
    -------
    dict
        Summary with keys: 'output', 'log', 'n_files', 'n_rows',
        'n_dropped_files', 'n_dropped_years'.
    """
    input_dir        = input_dir        or _DEFAULT_INPUT_DIR
    output_path      = output_path      or _DEFAULT_OUTPUT
    log_dir          = log_dir          or _DEFAULT_LOG_DIR
    stations_to_drop = stations_to_drop if stations_to_drop is not None else STATIONS_TO_DROP
    years_to_drop    = years_to_drop    if years_to_drop    is not None else YEARS_TO_DROP

    log_path = _setup_logger(log_dir)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    logger.info(f"Input dir  : {input_dir}")
    logger.info(f"Output     : {output_path}")
    if stations_to_drop:
        logger.info(f"Dropping stations : {', '.join(sorted(stations_to_drop))}")
    if years_to_drop:
        logger.info(f"Dropping years    : {', '.join(sorted(years_to_drop))}")

    yearly_dirs = sorted(
        d for d in os.listdir(input_dir)
        if os.path.isdir(os.path.join(input_dir, d))
    )

    header_written  = False
    n_files         = 0
    n_rows          = 0
    n_dropped_files = 0
    n_dropped_years = 0

    with open(output_path, 'w', encoding='utf-8') as out:
        for year_dir in yearly_dirs:
            year = _year_from_dirname(year_dir)
            if year in years_to_drop:
                n_dropped_years += 1
                continue

            year_path = os.path.join(input_dir, year_dir)
            csv_files = sorted(
                f for f in os.listdir(year_path) if f.endswith('.csv')
            )
            for fname in csv_files:
                code = _station_code_from_filename(fname)
                if code in stations_to_drop:
                    n_dropped_files += 1
                    continue

                fpath = os.path.join(year_path, fname)
                with open(fpath, 'r', encoding='utf-8') as f:
                    lines = f.readlines()

                if not lines:
                    continue

                if not header_written:
                    out.write(lines[0])
                    header_written = True

                for line in lines[1:]:
                    out.write(line)
                    n_rows += 1

                n_files += 1

    logger.info(f"Files merged  : {n_files}")
    logger.info(f"Rows written  : {n_rows}")
    logger.info(f"Files dropped : {n_dropped_files}")
    logger.info(f"Years dropped : {n_dropped_years}")
    logger.info(f"Log: {log_path}")

    return {
        'output':          output_path,
        'log':             log_path,
        'n_files':         n_files,
        'n_rows':          n_rows,
        'n_dropped_files': n_dropped_files,
        'n_dropped_years': n_dropped_years,
    }


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Merge OMP PhaseNet pick CSVs into a single consolidated file.'
    )
    parser.add_argument(
        '--input-dir', default=None,
        help='Directory containing yearly PICKS_PHASENET_* subdirectories.'
    )
    parser.add_argument(
        '--output', default=None,
        help='Output CSV file path.'
    )
    parser.add_argument(
        '--log-dir', default=None,
        help='Directory for log files.'
    )
    parser.add_argument(
        '--drop-years', nargs='*', default=None, metavar='YEAR',
        help='Years to exclude (e.g. --drop-years 2020 2021).'
    )
    args = parser.parse_args()
    years_to_drop = set(args.drop_years) if args.drop_years else None
    merge_omp(args.input_dir, args.output, args.log_dir, years_to_drop=years_to_drop)


if __name__ == '__main__':
    main()
