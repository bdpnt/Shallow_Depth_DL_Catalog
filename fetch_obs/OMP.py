"""
OMP.py
============================
Convert the OMP (Observatoire Midi-Pyrénées) .mag catalog to the .obs
bulletin format, extracting P and S picks.

Usage
-----
    python fetch_obs/OMP.py \\
        --file-name org_catalogs/OMP_78-19.mag \\
        --save-name obs/OMP_78-19.obs
"""

import argparse
import datetime
import logging
import os
from dataclasses import dataclass
from datetime import datetime as dt

from obspy import UTCDateTime


logger = logging.getLogger('fetch_obs.OMP')

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
class OMPParams:
    """
    Configuration for converting the OMP catalog.

    Attributes
    ----------
    file_name : str — path to the OMP .mag input file
    save_name : str — path for the .obs output file
    """
    file_name: str
    save_name: str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _open_catalog(file_name):
    """Read a catalog file and return its lines."""
    with open(file_name, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
    logger.info(f"Catalog loaded: {file_name} ({len(lines)} lines)")
    return lines


def _format_arrival_datetime(arrival):
    """Return (date, hours, seconds) strings formatted for the .obs bulletin."""
    ms_str  = str(arrival.microsecond)
    ms_str  = ms_str.zfill(3) if len(ms_str) < 3 else ms_str[:3]
    date    = f'{arrival.year:04d}{arrival.month:02d}{arrival.day:02d}'
    hours   = f'{arrival.hour:02d}{arrival.minute:02d}'
    seconds = f'{arrival.second:02d}.{ms_str}'
    return date, hours, seconds


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def write_catalog_to_obs(parameters, log_dir=None):
    """
    Convert the OMP .mag catalog to the .obs bulletin format.

    Parameters
    ----------
    parameters : OMPParams
    log_dir    : str, optional — log directory; default: fetch_obs/console_output/

    Returns
    -------
    dict with key: output
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file   : {log_path}")
    logger.info(f"Input file : {parameters.file_name}")
    logger.info(f"Output file: {parameters.save_name}")

    lines = _open_catalog(parameters.file_name)

    n_events          = 0
    n_skipped_neg     = 0
    n_skipped_mag99   = 0
    n_skipped_invalid = 0
    n_picks_p         = 0
    n_picks_s         = 0
    n_skipped_quality = 0
    n_skipped_miss_p  = 0
    n_skipped_miss_s  = 0

    with open(parameters.save_name, 'w') as f:
        f.write(f'### Catalog generated on the {UTCDateTime()}\n')
        f.write('### Year Month Day Hour Min Sec Lat Lon Dep Mag MagType MagAuthor PhaseCount HorUncer VerUncer AzGap RMS\n')
        f.write('### Code Ins Comp Onset Phase Dir Date HHMM S.MS Err ErrMag CodaDur P2PAmp PeriodAmp # RealPhase Channel PickOrigin PGV\n')
        f.write('\n')

        for ind, line in enumerate(lines):
            if not line.startswith(' Localisation'):
                continue

            event_header = lines[ind + 2]
            inverted_lon = -1 if event_header[35] == 'W' else 1

            event_info = lines[ind + 4]
            year       = int(event_info[0:3])
            month      = int(event_info[3:5])
            day        = int(event_info[5:7])
            hour       = int(event_info[8:10])
            minute     = int(event_info[10:12])
            second     = float(event_info[13:18])
            lat_deg    = int(event_info[19:21])
            lat_sec    = float(event_info[22:27])
            lon_deg    = int(event_info[29:31])
            lon_sec    = float(event_info[32:37])
            depth      = float(event_info[39:44])
            rms        = float(event_info[45:50])
            magnitude  = float(event_info[51:54])

            if magnitude < -5:
                n_skipped_neg += 1
                continue

            latitude  = lat_deg + lat_sec / 60
            longitude = inverted_lon * (lon_deg + lon_sec / 60)

            year = 2000 + year if year < 78 else 1900 + year

            weird_hour = False
            if hour >= 24:
                hour = 23
                weird_hour = True

            if second >= 60:
                second -= 60
                minute += 1

            try:
                UTCDateTime(f'{year}-{month}-{day}T{hour}:{minute}:{second}Z')
            except Exception:
                n_skipped_invalid += 1
                continue

            event_date = datetime.datetime(year, month, day, 0, 0, 0)

            if magnitude == 9.9:
                n_skipped_mag99 += 1
                continue

            f.write(
                f'# {year} {month} {day} {hour} {minute} {second} '
                f'{latitude} {longitude} {depth} {magnitude} '
                f'ML OMP None None None None {rms}\n'
            )
            n_events += 1

            phase_ind = ind + 7
            while phase_ind < len(lines) and lines[phase_ind].strip():
                double_error  = False
                phase_info    = lines[phase_ind]

                station    = phase_info[1:5].strip()
                code       = ('.' + station).ljust(9)
                instrument = '?'.ljust(4)
                component  = '?'.ljust(4)
                onset      = '?'.ljust(1)
                first_mot  = '?'.ljust(1)
                error_type = 'GAU'.ljust(3)
                coda_dur   = '-1.00e+00'.ljust(9)
                amp        = '-1.00e+00'.ljust(9)
                period     = '-1.00e+00'.ljust(9)

                quality_p = phase_info[23:24]
                quality_s = phase_info[102:103]

                if quality_p == '4' or quality_s == '0':
                    n_skipped_quality += 1
                    phase_ind += 1
                    continue
                elif quality_p == '9' or station == 'LARF':
                    instrument = '*'.ljust(4)
                elif int(quality_p) >= 2 or int(quality_s) >= 3:
                    double_error = True

                # --- P phase ---
                hour_p   = 23 if weird_hour else int(phase_info[25:27])
                minute_p = int(phase_info[27:29])
                second_p = phase_info[30:35].strip()

                if not second_p or '*' in second_p:
                    n_skipped_miss_p += 1
                    phase_ind += 1
                    continue

                second_p = float(second_p)
                if second_p < 0:
                    minute_p  = minute - 1
                    second_p += 60

                arrival_p  = event_date + datetime.timedelta(seconds=hour_p * 3600 + minute_p * 60 + second_p)

                try:
                    UTCDateTime(f'{arrival_p.year}-{arrival_p.month}-{arrival_p.day}T'
                                f'{arrival_p.hour}:{arrival_p.minute}:{arrival_p.second}'
                                f'.{arrival_p.microsecond}Z')
                except Exception:
                    phase_ind += 1
                    continue

                date_p, hours_p, secs_p = _format_arrival_datetime(arrival_p)
                err_p = ('0.05' if not double_error else '0.10').ljust(9)

                f.write(
                    f"{code} {instrument} {component} {onset} {'P'.ljust(6)} {first_mot} "
                    f"{date_p} {hours_p} {secs_p} {error_type} {err_p} "
                    f"{coda_dur} {amp} {period}"
                    f" # {'P'.ljust(6)} {'None'.ljust(4)} {'OMP'.ljust(9)} {'None'.ljust(4)}\n"
                )
                n_picks_p += 1

                # --- S phase ---
                hour_s   = 23 if weird_hour else int(phase_info[25:27])
                minute_s = int(phase_info[27:29])
                second_s = phase_info[105:110].strip()

                if not second_s or '*' in second_s:
                    n_skipped_miss_s += 1
                    phase_ind += 1
                    continue

                second_s = float(second_s)
                if second_s < 0:
                    minute_s  = minute - 1
                    second_s += 60

                arrival_s = event_date + datetime.timedelta(seconds=hour_s * 3600 + minute_s * 60 + second_s)

                try:
                    UTCDateTime(f'{arrival_s.year}-{arrival_s.month}-{arrival_s.day}T'
                                f'{arrival_s.hour}:{arrival_s.minute}:{arrival_s.second}'
                                f'.{arrival_s.microsecond}Z')
                except Exception:
                    phase_ind += 1
                    continue

                date_s, hours_s, secs_s = _format_arrival_datetime(arrival_s)
                err_s = ('0.05' if not double_error else '0.10').ljust(9)

                f.write(
                    f"{code} {instrument} {component} {onset} {'S'.ljust(6)} {first_mot} "
                    f"{date_s} {hours_s} {secs_s} {error_type} {err_s} "
                    f"{coda_dur} {amp} {period}"
                    f" # {'S'.ljust(6)} {'None'.ljust(4)} {'OMP'.ljust(9)} {'None'.ljust(4)}\n"
                )
                n_picks_s += 1

                phase_ind += 1

            f.write('\n')

    logger.info(f"Events written                        : {n_events}")
    logger.info(f"Events skipped (magnitude < -5)        : {n_skipped_neg}")
    logger.info(f"Events skipped (magnitude == 9.9)     : {n_skipped_mag99}")
    logger.info(f"Events skipped (invalid datetime)     : {n_skipped_invalid}")
    logger.info(f"P picks written                       : {n_picks_p}")
    logger.info(f"S picks written                       : {n_picks_s}")
    logger.info(f"Picks skipped (quality filter)        : {n_skipped_quality}")
    logger.info(f"Picks skipped (missing P second)      : {n_skipped_miss_p}")
    logger.info(f"Picks skipped (missing S second)      : {n_skipped_miss_s}")
    logger.info(f"Output: {parameters.save_name}")
    return {'output': parameters.save_name}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert the OMP .mag catalog to .obs format.'
    )
    parser.add_argument('--file-name', required=True, help='Input OMP .mag file')
    parser.add_argument('--save-name', required=True, help='Output .obs file')
    args = parser.parse_args()

    params = OMPParams(file_name=args.file_name, save_name=args.save_name)
    write_catalog_to_obs(params)


if __name__ == '__main__':
    main()
