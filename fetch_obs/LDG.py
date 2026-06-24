"""
LDG.py
============================
Convert the LDG (Laboratoire de Détection et de Géophysique) catalog CSV
files to the .obs bulletin format.

Usage
-----
    python fetch_obs/LDG.py \\
        --catalog-file org_catalogs/LDG_20-25_catalog.txt \\
        --arrival-file org_catalogs/LDG_20-25_arrivals.txt \\
        --save-name    obs/LDG_20-25.obs
"""

import argparse
import logging
import os
from dataclasses import dataclass
from datetime import datetime as dt

import pandas as pd
from obspy import UTCDateTime


logger = logging.getLogger('fetch_obs.LDG')

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
class LDGParams:
    """
    Configuration for converting the LDG catalog.

    Attributes
    ----------
    catalog_file : str — path to the LDG events CSV
    arrival_file : str — path to the LDG arrivals CSV
    save_name    : str — path for the .obs output file
    """
    catalog_file: str
    arrival_file: str
    save_name:    str


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def write_catalog_to_obs(parameters, log_dir=None):
    """
    Convert the LDG catalog and arrivals CSV files to the .obs bulletin format.

    Parameters
    ----------
    parameters : LDGParams
    log_dir    : str, optional — log directory; default: fetch_obs/console_output/

    Returns
    -------
    dict with key: output
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file      : {log_path}")
    logger.info(f"Catalog file  : {parameters.catalog_file}")
    logger.info(f"Arrivals file : {parameters.arrival_file}")
    logger.info(f"Output file   : {parameters.save_name}")

    catalog  = pd.read_csv(parameters.catalog_file, sep=';', header=0)
    arrivals = pd.read_csv(parameters.arrival_file,  sep=';', header=0)

    logger.info(f"Events loaded: {len(catalog)}")
    logger.info(f"Picks loaded : {len(arrivals)}")

    n_events           = 0
    n_skipped_mag      = 0
    n_events_no_picks  = 0
    n_picks            = 0
    n_picks_skip_phase = 0

    with open(parameters.save_name, 'w') as f:
        f.write(f'### Catalog generated on the {UTCDateTime()}\n')
        f.write('### Year Month Day Hour Min Sec Lat Lon Dep Mag MagType MagAuthor PhaseCount HorUncer VerUncer AzGap RMS\n')
        f.write('### Code Ins Comp Onset Phase Dir Date HHMM S.MS Err ErrMag CodaDur P2PAmp PeriodAmp # RealPhase Channel PickOrigin PGV\n')
        f.write('\n')

        arrivals_by_orid = arrivals.groupby('orid')

        for row in catalog.itertuples():
            year   = row.datetime[6:10]
            month  = row.datetime[3:5]
            day    = row.datetime[0:2]
            hour   = row.datetime[11:13]
            minute = row.datetime[14:16]
            second = row.datetime[17:19]

            magnitude = row.ML if row.ML != -999 else (row.MD if row.MD != -999 else None)
            if magnitude is None:
                n_skipped_mag += 1
                continue

            magnitude_type = 'ML' if row.ML != -999 else 'MD'
            phases_count   = row.nbmagML if row.ML != -999 else row.nbmagMD

            f.write(
                f"# {year} {month.lstrip('0')} {day.lstrip('0')} "
                f"{hour.lstrip('0') if hour != '00' else '0'} "
                f"{minute.lstrip('0') if minute != '00' else '0'} "
                f"{second[1:] if second.startswith('00') else second.lstrip('0')} "
                f"{row.lat} {row.lon} {row.depth} {magnitude} "
                f"{magnitude_type} LDG {phases_count} None None {row.gap1} {row.rms}\n"
            )
            n_events += 1

            event_id = row.orid
            if event_id not in arrivals_by_orid.groups:
                n_events_no_picks += 1
                f.write('\n')
                continue

            for row_p in arrivals_by_orid.get_group(event_id).itertuples():
                phase = row_p.phase
                if (not phase.strip().lower().startswith('p') and
                        not phase.strip().lower().startswith('s')):
                    n_picks_skip_phase += 1
                    continue

                station   = row_p.sta
                yr        = row_p.arrtime[6:10]
                mo        = row_p.arrtime[3:5]
                dy        = row_p.arrtime[0:2]
                hr        = row_p.arrtime[11:13]
                mn        = row_p.arrtime[14:16]
                sc        = row_p.arrtime[17:22].ljust(6, '0')
                error_mag = '0.05' if phase.lower().startswith('p') else '0.15'

                code       = ('.' + station).ljust(9)
                phase_type = phase[0].ljust(6)
                date       = yr + mo + dy
                hours      = hr + mn

                f.write(
                    f"{code} {'?'.ljust(4)} {'?'.ljust(4)} {'?'.ljust(1)} {phase_type} {'?'.ljust(1)} "
                    f"{date} {hours} {sc} {'GAU'.ljust(3)} {error_mag.ljust(9)} "
                    f"{'-1.00e+00'.ljust(9)} {'-1.00e+00'.ljust(9)} {'-1.00e+00'.ljust(9)}"
                    f" # {phase.ljust(6)} {'None'.ljust(4)} {'LDG'.ljust(9)} {'None'.ljust(4)}\n"
                )
                n_picks += 1

            f.write('\n')

    logger.info(f"Events written                        : {n_events}")
    logger.info(f"Events skipped (no valid magnitude)   : {n_skipped_mag}")
    logger.info(f"Events written with no picks          : {n_events_no_picks}")
    logger.info(f"Picks written                         : {n_picks}")
    logger.info(f"Picks skipped (wrong phase)           : {n_picks_skip_phase}")
    logger.info(f"Output: {parameters.save_name}")
    return {'output': parameters.save_name}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert the LDG catalog CSV files to .obs format.'
    )
    parser.add_argument('--catalog-file', required=True, help='LDG events CSV')
    parser.add_argument('--arrival-file', required=True, help='LDG arrivals CSV')
    parser.add_argument('--save-name',    required=True, help='Output .obs file')
    args = parser.parse_args()

    params = LDGParams(
        catalog_file = args.catalog_file,
        arrival_file = args.arrival_file,
        save_name    = args.save_name,
    )
    write_catalog_to_obs(params)


if __name__ == '__main__':
    main()
