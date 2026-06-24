"""
RESIF.py
============================
Fetch the RESIF FDSN earthquake catalog and convert it to the .obs bulletin format.

Events are queried year by year from the FDSN client. Only manual P/S picks
are written. The QuakeML file is also saved alongside the .obs output.

Usage
-----
    python fetch_obs/RESIF.py \\
        --file-name  org_catalogs/RESIF_20-25.xml \\
        --save-name  obs/RESIF_20-25.obs \\
        --client     RESIF \\
        --t1         2020-01-01T00:00:00 \\
        --t2         2026-01-01T00:00:00 \\
        --lat-min    41 --lat-max 44 \\
        --lon-min    -3 --lon-max 4 \\
        --mag-min    0  --event-type earthquake \\
        --mag-type   MLv
"""

import argparse
import logging
import os
from dataclasses import dataclass
from datetime import datetime as dt

from obspy import UTCDateTime, read_events
from obspy.clients.fdsn import Client


logger = logging.getLogger('fetch_obs.RESIF')

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
class RESIFParams:
    """
    Configuration for fetching and converting the RESIF catalog.

    Attributes
    ----------
    client_name : str        — FDSN client identifier (e.g. 'RESIF')
    t1          : UTCDateTime — query start time
    t2          : UTCDateTime — query end time
    lat_min     : float      — minimum latitude
    lat_max     : float      — maximum latitude
    lon_min     : float      — minimum longitude
    lon_max     : float      — maximum longitude
    mag_min     : float      — minimum magnitude
    event_type  : str        — event type filter (e.g. 'earthquake')
    file_name   : str        — path for the QuakeML output file
    save_name   : str        — path for the .obs output file
    mag_type    : str        — magnitude type to retain (e.g. 'MLv')
    """
    client_name: str
    t1:          UTCDateTime
    t2:          UTCDateTime
    lat_min:     float
    lat_max:     float
    lon_min:     float
    lon_max:     float
    mag_min:     float
    event_type:  str
    file_name:   str
    save_name:   str
    mag_type:    str


# ---------------------------------------------------------------------------
# FDSN helpers
# ---------------------------------------------------------------------------

def _fetch_year_slice(client, parameters, year_start, year_end):
    """Query the FDSN client for events in a given time slice."""
    return client.get_events(
        starttime=year_start, endtime=year_end,
        minlatitude=parameters.lat_min,   maxlatitude=parameters.lat_max,
        minlongitude=parameters.lon_min,  maxlongitude=parameters.lon_max,
        minmagnitude=parameters.mag_min,  eventtype=parameters.event_type,
        includeallorigins=False, includeallmagnitudes=False,
        includearrivals=True, orderby='time-asc',
    )


def _find_best_magnitude(event, mag_type):
    """
    Find the best available manual magnitude of a given type for an event.

    Returns
    -------
    (is_preferred, index)
        is_preferred=True, index=None  → use event.preferred_magnitude()
        is_preferred=False, index=int  → use event.magnitudes[index]
        is_preferred=False, index=None → no suitable magnitude found
    """
    try:
        pm = event.preferred_magnitude()
        if pm.magnitude_type == mag_type and 'auto' not in pm.creation_info.author:
            return True, None
        raise ValueError()
    except Exception:
        try:
            best = ()
            for i, mag in enumerate(event.magnitudes):
                if mag.magnitude_type != mag_type:
                    continue
                if 'auto' in getattr(getattr(mag, 'creation_info', None), 'author', ''):
                    continue
                unc = mag.mag_errors.uncertainty
                if not best or unc < best[1]:
                    best = (i, unc)
            if not best:
                for i, mag in enumerate(event.magnitudes):
                    if mag.magnitude_type == mag_type and getattr(mag, 'evaluation_status', None) == 'confirmed':
                        best = (i, None)
            return False, best[0] if best else None
        except Exception:
            return False, None


def _fetch_catalog(parameters):
    """Load an existing QuakeML catalog file and return it as an obspy Catalog."""
    return read_events(parameters.file_name, format='QUAKEML')


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_catalog(parameters, log_dir=None):
    """
    Query FDSN year by year, write the QuakeML file, and return the catalog.

    Parameters
    ----------
    parameters : RESIFParams
    log_dir    : str, optional — log directory; default: fetch_obs/console_output/

    Returns
    -------
    obspy.core.event.Catalog
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Log file    : {log_path}")
    logger.info(f"FDSN client : {parameters.client_name}")
    logger.info(f"Time range  : {parameters.t1} → {parameters.t2}")
    logger.info(f"Area        : lat [{parameters.lat_min}, {parameters.lat_max}]  "
                f"lon [{parameters.lon_min}, {parameters.lon_max}]")
    logger.info(f"mag_min     : {parameters.mag_min}  mag_type: {parameters.mag_type}")

    client   = Client(parameters.client_name)
    catalog  = None
    t1, t2   = parameters.t1, parameters.t2

    current_year = t1.year
    end_year     = t2.year

    if t1.month != 1 or t1.day != 1 or t1.hour != 0 or t1.minute != 0 or t1.second != 0:
        yr_start = t1
        yr_end   = min(UTCDateTime(f'{current_year + 1}-01-01T00:00:00'), t2)
        year_cat = _fetch_year_slice(client, parameters, yr_start, yr_end)
        catalog  = year_cat if catalog is None else catalog + year_cat
        logger.info(f"Fetched {len(year_cat)} event(s) from {yr_start} to {yr_end}")
        current_year += 1

    while current_year < end_year:
        yr_start = UTCDateTime(f'{current_year}-01-01T00:00:00')
        yr_end   = UTCDateTime(f'{current_year + 1}-01-01T00:00:00')
        year_cat = _fetch_year_slice(client, parameters, yr_start, yr_end)
        catalog  = year_cat if catalog is None else catalog + year_cat
        logger.info(f"Fetched {len(year_cat)} event(s) for year {current_year}")
        current_year += 1

    if t2.month != 1 or t2.day != 1 or t2.hour != 0 or t2.minute != 0 or t2.second != 0:
        yr_start = UTCDateTime(f'{end_year}-01-01T00:00:00')
        year_cat = _fetch_year_slice(client, parameters, yr_start, t2)
        catalog  = year_cat if catalog is None else catalog + year_cat
        logger.info(f"Fetched {len(year_cat)} event(s) from {yr_start} to {t2}")

    logger.info(f"Total events fetched: {len(catalog)}")
    catalog.write(parameters.file_name, format='QUAKEML')
    logger.info(f"QuakeML written: {parameters.file_name}")
    return catalog


def write_catalog_to_obs(parameters, log_dir=None):
    """
    Convert a QuakeML catalog to the .obs bulletin format.

    Only manual P/S picks are written. The QuakeML file is also re-saved.

    Parameters
    ----------
    parameters : RESIFParams
    log_dir    : str, optional — log directory; default: fetch_obs/console_output/

    Returns
    -------
    dict with key: output
    """
    if not logger.handlers:
        log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
        logger.info(f"Log file   : {log_path}")

    catalog = _fetch_catalog(parameters)
    logger.info(f"Catalog loaded: {parameters.file_name} ({len(catalog)} events)")

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

        for event in catalog:
            origin = event.origins[0]
            year   = origin.time.year
            month  = origin.time.month
            day    = origin.time.day
            hour   = origin.time.hour
            minute = origin.time.minute
            second = origin.time.second + origin.time.microsecond / 1e6

            try:
                i_val, i_mag = _find_best_magnitude(event, parameters.mag_type)
                if i_val:
                    magnitude      = event.preferred_magnitude().mag
                    magnitude_type = event.preferred_magnitude().magnitude_type
                    mag_author     = event.preferred_magnitude().creation_info.agency_id
                elif i_mag is not None:
                    magnitude      = event.magnitudes[i_mag].mag
                    magnitude_type = event.magnitudes[i_mag].magnitude_type
                    mag_author     = event.magnitudes[i_mag].creation_info.agency_id
                else:
                    raise ValueError()
            except Exception:
                n_skipped_mag += 1
                continue

            qual          = getattr(origin, 'quality', None)
            phases_count  = getattr(qual, 'associated_phase_count', None)
            h_uncertainty = getattr(qual, 'horizontal_uncertainty', None)
            v_uncertainty = getattr(qual, 'vertical_uncertainty', None)
            az_gap        = getattr(qual, 'azimuthal_gap', None)
            rms           = getattr(qual, 'standard_error', None)
            depth         = origin.depth / 1000.0

            f.write(
                f'# {year} {month} {day} {hour} {minute} {second} '
                f'{origin.latitude} {origin.longitude} {depth} {magnitude} '
                f'{magnitude_type} {mag_author} {phases_count} '
                f'{h_uncertainty} {v_uncertainty} {az_gap} {rms}\n'
            )
            n_events += 1

            for pick in event.picks:
                if (not pick.phase_hint.lower().startswith('p') and
                        not pick.phase_hint.lower().startswith('s')):
                    n_picks_skip_phase += 1
                    continue
                if getattr(pick, 'evaluation_mode', None) != 'manual':
                    n_picks_skip_manual += 1
                    continue

                network   = str(pick.waveform_id.network_code)
                station   = str(pick.waveform_id.station_code)
                phase     = str(pick.phase_hint)
                yr        = str(pick.time.year).zfill(4)
                mo        = str(pick.time.month).zfill(2)
                dy        = str(pick.time.day).zfill(2)
                hr        = str(pick.time.hour).zfill(2)
                mn        = str(pick.time.minute).zfill(2)
                sc        = str(pick.time.second).zfill(2)
                ms        = str(int(pick.time.microsecond / 1000)).zfill(3)
                error_mag = '0.05' if phase.lower().startswith('p') else '0.15'

                code       = (network + '.' + station).ljust(9)
                phase_type = phase[0].ljust(6)
                date       = yr + mo + dy
                hours      = hr + mn
                seconds    = sc + '.' + ms
                real_phase = phase.ljust(6)
                channel    = str(pick.waveform_id.channel_code).ljust(4)

                f.write(
                    f"{code} {'?'.ljust(4)} {'?'.ljust(4)} {'?'.ljust(1)} {phase_type} {'?'.ljust(1)} "
                    f"{date} {hours} {seconds} {'GAU'.ljust(3)} {error_mag.ljust(9)} "
                    f"{'-1.00e+00'.ljust(9)} {'-1.00e+00'.ljust(9)} {'-1.00e+00'.ljust(9)}"
                    f" # {real_phase} {channel} {'FDSN'.ljust(9)} {'None'.ljust(4)}\n"
                )
                n_picks += 1

            f.write('\n')

    logger.info(f"Events written                        : {n_events}")
    logger.info(f"Events skipped (no matching magnitude): {n_skipped_mag}")
    logger.info(f"Picks written                         : {n_picks}")
    logger.info(f"Picks skipped (wrong phase)           : {n_picks_skip_phase}")
    logger.info(f"Picks skipped (non-manual)            : {n_picks_skip_manual}")
    logger.info(f"Output: {parameters.save_name}")
    catalog.write(parameters.file_name, format='QUAKEML')
    logger.info(f"QuakeML re-saved: {parameters.file_name}")

    return {'output': parameters.save_name}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Fetch the RESIF FDSN catalog and write a .obs bulletin.'
    )
    parser.add_argument('--file-name',   required=True)
    parser.add_argument('--save-name',   required=True)
    parser.add_argument('--client',      default='RESIF')
    parser.add_argument('--t1',          required=True, help='Start time (ISO format)')
    parser.add_argument('--t2',          required=True, help='End time (ISO format)')
    parser.add_argument('--lat-min',     type=float, required=True)
    parser.add_argument('--lat-max',     type=float, required=True)
    parser.add_argument('--lon-min',     type=float, required=True)
    parser.add_argument('--lon-max',     type=float, required=True)
    parser.add_argument('--mag-min',     type=float, default=0)
    parser.add_argument('--event-type',  default='earthquake')
    parser.add_argument('--mag-type',    default='MLv')
    args = parser.parse_args()

    params = RESIFParams(
        client_name = args.client,
        t1          = UTCDateTime(args.t1),
        t2          = UTCDateTime(args.t2),
        lat_min     = args.lat_min,
        lat_max     = args.lat_max,
        lon_min     = args.lon_min,
        lon_max     = args.lon_max,
        mag_min     = args.mag_min,
        event_type  = args.event_type,
        file_name   = args.file_name,
        save_name   = args.save_name,
        mag_type    = args.mag_type,
    )
    generate_catalog(params)
    write_catalog_to_obs(params)


if __name__ == '__main__':
    main()
