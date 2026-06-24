"""
ICGC.py
============================
Fetch the ICGC (Institut Cartogràfic i Geològic de Catalunya) earthquake
catalog and convert it to the .obs bulletin format.

Step 1 — get_all_codes  : fetch all event codes for the date range.
Step 2 — fetch_catalog  : download the GSE2 bulletin for each code.
Step 3 — write_catalog_to_obs : parse the GSE2 file and write .obs.

Usage
-----
    python fetch_obs/ICGC.py \\
        --file-name   org_catalogs/ICGC_20-25.txt \\
        --code-name   org_catalogs/CODES_ICGC_20-25.txt \\
        --error-name  org_catalogs/ERR_ICGC_20-25.txt \\
        --save-name   obs/ICGC_20-25.obs \\
        --start-year  2020 --start-month 1 \\
        --end-year    2025 --end-month   12
"""

import argparse
import logging
import os
import time
from dataclasses import dataclass
from datetime import datetime as dt

import requests
from obspy import UTCDateTime


logger = logging.getLogger('fetch_obs.ICGC')

_DEFAULT_LOG_DIR = 'fetch_obs/console_output/'


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
class ICGCParams:
    """
    Configuration for fetching and converting the ICGC catalog.

    Attributes
    ----------
    file_name   : str — path for the downloaded GSE2 catalog
    code_name   : str — path for the event codes file
    error_name  : str — path for the error log
    save_name   : str — path for the .obs output file
    start_year  : int — query start year
    start_month : int — query start month
    end_year    : int — query end year
    end_month   : int — query end month
    mag_min     : float — minimum magnitude (default: 0)
    """
    file_name:   str
    code_name:   str
    error_name:  str
    save_name:   str
    start_year:  int
    start_month: int
    end_year:    int
    end_month:   int
    mag_min:     float = 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _iter_months(start_year, start_month, end_year, end_month):
    """Yield (year, month_str) tuples for every month in the date range."""
    year, month = start_year, start_month
    while (year < end_year) or (year == end_year and month <= end_month):
        yield year, f'{month:02d}'
        month += 1
        if month > 12:
            month = 1
            year += 1


def _get_codes(year, month):
    """
    Fetch the list of event codes from the ICGC website for a given month.

    Returns
    -------
    (True, list[str]) on success or (False, error_info) on failure.
    """
    url = (
        f'https://sismocat.icgc.cat/siswebclient/index.php'
        f'?seccio=llistat&area=locals'
        f'&any={str(year).lstrip("0")}&mes={str(month).lstrip("0")}&idioma=ca'
    )
    for _ in range(3):
        try:
            response = requests.get(url, timeout=15)
            break
        except requests.exceptions.RequestException as e:
            print(f'Request failed, retrying... ({e})')
            time.sleep(2)
    else:
        return False, 'Failed after 3 retries'

    if response.status_code != 200:
        return False, response.status_code

    html   = response.text
    lines  = html.split('<a class')[1:]
    codes  = [event.split('>')[1].rstrip('</a') for event in lines]
    return True, codes


def _safe_float(s):
    """Convert a string to float, returning None on failure."""
    try:
        return float(s.strip())
    except Exception:
        return None


def _open_catalog(file_name):
    """Read a catalog file and return its lines."""
    with open(file_name, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    logger.info(f"Catalog loaded: {file_name} ({len(lines)} lines)")
    return lines


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def get_all_codes(parameters):
    """
    Fetch all event codes for the full date range and write them to code_name.

    Parameters
    ----------
    parameters : ICGCParams

    Returns
    -------
    dict with keys: output (code_name), errors (error_name)
    """
    with open(parameters.code_name, 'w') as f, open(parameters.error_name, 'w') as fe:
        fe.write('### ERRORS DURING CODES FETCH\n')
        for year, month in _iter_months(
            parameters.start_year, parameters.start_month,
            parameters.end_year,   parameters.end_month,
        ):
            status, value = _get_codes(year, month)
            if status:
                for code in value:
                    f.write(f'{code}, {year}-{month}\n')
            else:
                fe.write(f'{year}-{month} : error {value}\n')

    print(f'Codes file written  → {parameters.code_name}')
    print(f'Errors file written → {parameters.error_name}')
    return {'output': parameters.code_name, 'errors': parameters.error_name}


def fetch_catalog(parameters):
    """
    Download the GSE2 bulletin for each event code and append to file_name.

    Parameters
    ----------
    parameters : ICGCParams

    Returns
    -------
    dict with key: output
    """
    if os.path.exists(parameters.file_name):
        os.remove(parameters.file_name)

    with open(parameters.code_name, 'r', encoding='utf-8', errors='ignore') as f:
        codes = [line.split(',')[0] for line in f]

    first_error = True

    for code in codes:
        if not code.strip():
            continue

        print(f'Processing code: {code}')
        url = f'http://sismocat.icgc.cat/siswebclient/index.php?seccio=gse&codi={code}'

        for _ in range(3):
            try:
                response = requests.get(url, timeout=10)
                break
            except requests.exceptions.RequestException as e:
                print(f'Request failed, retrying... ({e})')
                time.sleep(2)
        else:
            print(f'Failed to fetch {url} after 3 retries.')
            with open(parameters.error_name, 'a') as fe:
                if first_error:
                    first_error = False
                    fe.write('### ERRORS DURING EVENTS FETCH\n')
                fe.write(f'{code} : Failed after retries\n')
            continue

        if "S'ha produit un error" in response.text:
            print(f'Error page for code: {code}')
            with open(parameters.error_name, 'a') as fe:
                if first_error:
                    first_error = False
                    fe.write('### ERRORS DURING EVENTS FETCH\n')
                fe.write(f'{code} : Error page received\n')
        else:
            with open(parameters.file_name, 'ab') as f:
                f.write(response.content)
            print(f'Written: {code}')

    print(f'Catalog written → {parameters.file_name}')
    print(f'Errors written  → {parameters.error_name}')
    return {'output': parameters.file_name}


def write_catalog_to_obs(parameters, log_dir=None):
    """
    Convert the ICGC GSE2 catalog to the .obs bulletin format.

    Only manual P/S picks are written.

    Parameters
    ----------
    parameters : ICGCParams
    log_dir    : str, optional — log directory; default: fetch_obs/console_output/

    Returns
    -------
    dict with key: output
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file   : {log_path}")
    logger.info(f"Input file : {parameters.file_name}")
    logger.info(f"Output file: {parameters.save_name}")
    logger.info(f"mag_min    : {parameters.mag_min}")

    lines = _open_catalog(parameters.file_name)

    n_events            = 0
    n_skipped_mag       = 0
    n_picks             = 0
    n_picks_skip_phase  = 0
    n_picks_skip_manual = 0

    with open(parameters.save_name, 'w') as f:
        f.write(f'### Catalog generated on the {UTCDateTime()}\n')
        f.write('### Year Month Day Hour Min Sec Lat Lon Dep Mag MagType MagAuthor PhaseCount HorUncer VerUncer AzGap RMS\n')
        f.write('### Code Ins Comp Onset Phase Dir Date HHMM S.MS Err ErrMag CodaDur P2PAmp PeriodAmp # RealPhase Channel PickOrigin PGV\n')
        f.write('\n')

        for ind, line in enumerate(lines):
            if not line.startswith('DATA_TYPE'):
                continue

            # Event header: 3rd line after DATA_TYPE
            event_info = lines[ind + 3].rstrip('\n')
            year       = event_info[0:4].strip()
            month      = event_info[5:7].strip()
            day        = event_info[8:10].strip()
            ev_hour    = event_info[11:13].strip()
            minute     = event_info[14:16].strip()
            second     = event_info[17:22].strip()
            latitude   = _safe_float(event_info[36:44])
            longitude  = _safe_float(event_info[45:54])
            depth      = _safe_float(event_info[71:76])
            az_gap     = _safe_float(event_info[92:97])
            rms        = _safe_float(event_info[30:35])

            # Magnitude: 6th line after DATA_TYPE
            mag_info       = lines[ind + 6].rstrip('\n')
            magnitude      = _safe_float(mag_info[7:10])
            magnitude_type = mag_info[0:6].strip()
            mag_author     = mag_info[20:29].strip()

            phases_count = _safe_float(event_info[89:93])

            if magnitude is None or magnitude < parameters.mag_min:
                n_skipped_mag += 1
                continue

            f.write(
                f"# {year} {month.lstrip('0')} {day.lstrip('0')} "
                f"{ev_hour.lstrip('0') if ev_hour != '00' else '0'} "
                f"{minute.lstrip('0') if minute != '00' else '0'} "
                f"{second[1:] if second.startswith('00') else second.lstrip('0')} "
                f"{latitude} {longitude} {depth} {magnitude} "
                f"{magnitude_type} {mag_author} {phases_count} None None {az_gap} {rms}\n"
            )
            n_events += 1

            # Phases: from 11th line after DATA_TYPE
            phase_ind = ind + 11
            while phase_ind < len(lines) and lines[phase_ind].strip():
                phase_info = lines[phase_ind].rstrip('\n')

                phase_name = phase_info[19:27].strip()
                if (not phase_name.lower().startswith('p') and
                        not phase_name.lower().startswith('s')):
                    n_picks_skip_phase += 1
                    phase_ind += 1
                    continue
                if phase_info[99:102] != 'm__':
                    n_picks_skip_manual += 1
                    phase_ind += 1
                    continue

                network   = phase_info[114:116].strip()
                station   = phase_info[0:7].strip()
                phase     = phase_name
                hr        = phase_info[28:30].strip()
                mn        = phase_info[31:33].strip()
                sc        = phase_info[34:36].strip()
                ms        = phase_info[37:41].strip()
                error_mag = '0.05' if phase.lower().startswith('p') else '0.15'

                code       = (network + '.' + station).ljust(9)
                phase_type = phase[0].ljust(6)
                if ev_hour == '23' and hr == '00':
                    phase_day = str(int(day.lstrip('0')) + 1).ljust(2)
                    date = year + month + phase_day
                else:
                    date = year + month + day
                hours   = hr + mn
                seconds = sc + '.' + ms

                f.write(
                    f"{code} {'?'.ljust(4)} {'?'.ljust(4)} {'?'.ljust(1)} {phase_type} {'?'.ljust(1)} "
                    f"{date} {hours} {seconds} {'GAU'.ljust(3)} {error_mag.ljust(9)} "
                    f"{'-1.00e+00'.ljust(9)} {'-1.00e+00'.ljust(9)} {'-1.00e+00'.ljust(9)}"
                    f" # {phase.ljust(6)} {'None'.ljust(4)} {'ICGC'.ljust(9)} {'None'.ljust(4)}\n"
                )
                n_picks += 1
                phase_ind += 1

            f.write('\n')

    logger.info(f"Events written                   : {n_events}")
    logger.info(f"Events skipped (magnitude filter): {n_skipped_mag}")
    logger.info(f"Picks written                    : {n_picks}")
    logger.info(f"Picks skipped (wrong phase)      : {n_picks_skip_phase}")
    logger.info(f"Picks skipped (non-manual)       : {n_picks_skip_manual}")
    logger.info(f"Output: {parameters.save_name}")
    return {'output': parameters.save_name}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Fetch the ICGC catalog (codes → GSE2 → .obs).'
    )
    parser.add_argument('--file-name',   required=True)
    parser.add_argument('--code-name',   required=True)
    parser.add_argument('--error-name',  required=True)
    parser.add_argument('--save-name',   required=True)
    parser.add_argument('--start-year',  type=int, required=True)
    parser.add_argument('--start-month', type=int, required=True)
    parser.add_argument('--end-year',    type=int, required=True)
    parser.add_argument('--end-month',   type=int, required=True)
    parser.add_argument('--mag-min',     type=float, default=-5)
    args = parser.parse_args()

    params = ICGCParams(
        file_name   = args.file_name,
        code_name   = args.code_name,
        error_name  = args.error_name,
        save_name   = args.save_name,
        start_year  = args.start_year,
        start_month = args.start_month,
        end_year    = args.end_year,
        end_month   = args.end_month,
        mag_min     = args.mag_min,
    )
    get_all_codes(params)
    fetch_catalog(params)
    write_catalog_to_obs(params)


if __name__ == '__main__':
    main()
