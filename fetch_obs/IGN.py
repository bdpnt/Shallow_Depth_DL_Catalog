"""
IGN.py
============================
Convert the IGN (Instituto Geográfico Nacional) GSE2 catalog to the
.obs bulletin format.

Only manual P/S picks are written.

Usage
-----
    python fetch_obs/IGN.py \\
        --file-name org_catalogs/IGN_20-25.txt \\
        --save-name obs/IGN_20-25.obs
"""

import argparse
import logging
import os
from dataclasses import dataclass
from datetime import datetime as dt

from obspy import UTCDateTime


logger = logging.getLogger('fetch_obs.IGN')

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
class IGNParams:
    """
    Configuration for converting the IGN catalog.

    Attributes
    ----------
    file_name : str — path to the IGN GSE2 input file
    save_name : str — path for the .obs output file
    """
    file_name: str
    save_name: str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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

def write_catalog_to_obs(parameters, log_dir=None):
    """
    Convert the IGN GSE2 catalog to the .obs bulletin format.

    Parameters
    ----------
    parameters : IGNParams
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

            # Event header: 4th line after DATA_TYPE
            event_info = lines[ind + 4].rstrip('\n')
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

            # Magnitude: 11th line after DATA_TYPE
            mag_info       = lines[ind + 11].rstrip('\n')
            magnitude      = _safe_float(mag_info[7:10])
            magnitude_type = mag_info[0:6].strip()
            mag_author     = mag_info[20:29].strip()

            phases_count = _safe_float(event_info[89:93])

            if magnitude is None:
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

            # Phases: from 15th line after DATA_TYPE
            phase_ind = ind + 15
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
                channel   = phase_info[119:123].strip()

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
                    f" # {phase.ljust(6)} {channel.ljust(4)} {'IGN'.ljust(9)} {'None'.ljust(4)}\n"
                )
                n_picks += 1
                phase_ind += 1

            f.write('\n')

    logger.info(f"Events written                   : {n_events}")
    logger.info(f"Events skipped (no magnitude)    : {n_skipped_mag}")
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
        description='Convert the IGN GSE2 catalog to .obs format.'
    )
    parser.add_argument('--file-name', required=True, help='Input IGN GSE2 file')
    parser.add_argument('--save-name', required=True, help='Output .obs file')
    args = parser.parse_args()

    params = IGNParams(file_name=args.file_name, save_name=args.save_name)
    write_catalog_to_obs(params)


if __name__ == '__main__':
    main()
