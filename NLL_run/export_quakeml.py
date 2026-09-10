"""
export_quakeml.py
============================
Export the final relocated catalog as a QuakeML 1.2 bulletin.

Merges obs/SSST_result.obs (magnitudes + phase picks) with RESULT/SSST_result.csv
(full-precision hypocentres, confidence ellipsoid, location-PDF quality metrics)
on publicId, and writes RESULT/FINAL.xml.

Every event is exported. Each carries a boolean `pyr:usable` flag derived from
the location-PDF quality metrics, plus the metrics themselves so the cut can be
re-derived or re-tuned downstream without recomputing anything.

Two origins per event
---------------------
The preferred origin is the location-PDF **expectation** (`pyr:locationEstimator
= "expectation"`); the maximum-likelihood point rides along as a second origin
(`.../maxlike`, `pyr:locationEstimator = "maximum_likelihood"`).

The expectation is what NLLoc reports under `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`,
and it is the estimate the published error actually describes: NLLoc's covariance,
its confidence ellipsoid and this catalog's true_erh/true_erz are all second
moments *about the expectation*, so attaching them to the mode quoted an ellipsoid
that was not centred on the point it belonged to. The same mismatch ran through
the QC — pdf_metrics computes C68 on samples centred on their own mean, so the
metric gating `pyr:usable` was testing an ellipsoid around a point the catalog did
not publish. It now tests the published one.

`standard_error`, `azimuthal_gap` and `minimum_distance` on the preferred origin
are NLLoc's own re-evaluation at the expectation — it re-reads every arrival's
travel time from the time grids there and re-solves the origin time, RMS and
residuals — not the maximum-likelihood values. The secondary origin therefore
carries no quality, no uncertainty and no arrivals: those belong to the
expectation, and repeating them under the mode would restore the very mismatch
this arrangement removes.

Usability rule
--------------
An event is unusable when its location PDF cannot be trusted:

  c68_overconfident   (C68 - 0.6827) / C68_sigma_n < -2, i.e. the posterior mass
                      inside the nominal one-sigma ellipsoid falls short by more than
                      two null sigmas -> the reported ERH/ERZ are over-confident
  dip_bimodal         Hartigan's dip test rejects unimodality of the depth
                      marginal AND the modal interval is wider than true_erz
                      -> two competing depths, too far apart for the quoted
                      error to cover, so no defensible scalar depth
  no_pdf_metrics      the event has no .scat cloud, so nothing can be said

Psi (and its underlying negentropy J) is exported but deliberately never used to
reject: it is directionless. A tight, curved, perfectly informative PDF scores as
badly as a diffuse one, and 81% of this catalog exceeds its own Gaussian null on
J -- gating on it would discard 85% of the events.

Why dip_bimodal needs the separation term as well as the p-value: the dip is a
vertical sup-distance on the ECDF and is invariant under affine rescaling of the
depth axis, so it registers that a second mode exists but never how far away it
sits. Rejecting on p alone therefore discards events whose bimodality is real but
irrelevant -- both modes inside the quoted error, where depth +/- true_erz already
covers them honestly. Over this catalog that is the common case: the median modal
interval among p < 0.05 events is 0.59 x true_erz. Requiring the interval to
exceed true_erz moves the count from 5329 to 661.

Why C68 is compared to a simulated null rather than to 0.6827 directly: the
one-sigma ellipsoid is fitted to the very samples whose coverage is being measured,
so even a perfectly Gaussian cloud scatters around 0.6827, and that scatter is not
binomial -- sqrt(.68*.32/n) overestimates it by ~2x (0.021 vs ~0.012 at n=479).
pdf_metrics simulates the real spread per sample size into C68_sigma_n. With the
median sigma of ~0.0124 the -2 sigma cut sits at C68 < 0.658, which flags on the
order of 150 events; a raw C68 < 0.6827 cut would flag thousands, nearly all of
them within sampling noise of perfect coverage.

Usage
-----
    python NLL_run/export_quakeml.py \\
        --obs       obs/SSST_result.obs \\
        --csv       RESULT/SSST_result.csv \\
        --inventory stations/GLOBAL_inventory.xml \\
        --output    RESULT/FINAL.xml
"""

import argparse
import io
import logging
import os
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd
from obspy import UTCDateTime, read_inventory
from obspy.core.event import (Arrival, Catalog, Comment, ConfidenceEllipsoid,
                              CreationInfo, Event, Magnitude, Origin,
                              OriginQuality, OriginUncertainty, Pick,
                              QuantityError, ResourceIdentifier,
                              WaveformStreamID)

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR      = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT    = os.path.dirname(_MODULE_DIR)
_DEFAULT_LOG_DIR = os.path.join(_MODULE_DIR, 'console_output')

sys.path.insert(0, _PROJECT_ROOT)

from NLL_run.merge_regional_results import _ellipsoid_axis_to_xyz
from temp_picks.match_picks         import load_bulletin
from fetch_inventory.merge_station_inventories import station_priority

logger = logging.getLogger('export_quakeml')

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Custom namespace for everything QuakeML has no element for.
_NS    = 'http://shallow-depth-dl-catalog/quakeml/1.0'
_NSMAP = {'pyr': _NS}

_C68_Z_MIN        = -2.0    # null sigmas below Gaussian coverage before rejection
_NOMINAL_COVERAGE = 0.6827  # erf(1/sqrt(2)): one-sigma coverage, matching pdf_metrics.py
_CHUNK_SIZE       = 5000    # events serialized per pass (bounds peak memory)
_SPLIT_YEARS      = 5       # calendar length of each split part

_KM_PER_DEGREE    = 111.195 # NLL reports Dist in km, QuakeML wants degrees

# The .obs event header fields, counted after stripping the leading '# '.
_H_MAG        = 9
_H_MAG_TYPE   = 10
_H_MAG_AUTHOR = 11

# The .obs pick line always holds exactly 19 whitespace-separated tokens.
_P_CODE        = 0
_P_INS         = 1     # '*' marks a pick NLLoc used via S-P relative timing
_P_PHASE       = 4     # NLLoc phase code, always P or S
_P_DATE        = 6     # YYYYMMDD
_P_HHMM        = 7
_P_SECONDS     = 8
_P_ERR_MAG     = 10    # pick uncertainty, seconds (1 sigma)
_P_REAL_PHASE  = 15    # original label before reduction: Pg, Sn, ...
_P_CHANNEL     = 16
_P_PICK_ORIGIN = 17


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ExportQuakeMLParams:
    file_obs:       str
    file_csv:       str   # path to RESULT/SSST_result.csv
    file_inventory: str
    save_file:      str
    write_catalog:  bool = True   # also emit <save_file stem>_catalog.xml, picks stripped
    write_parts:    bool = True   # also emit year-range parts sized to _SPLIT_BUDGET_GB


def _companion_path(save_file, suffix):
    """Path of a companion file next to save_file, e.g. FINAL.xml -> FINAL_catalog.xml."""
    stem, extension = os.path.splitext(save_file)
    return f'{stem}_{suffix}{extension}'


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
# Custom namespace helpers
# ---------------------------------------------------------------------------

def _x(value):
    """Wrap a value as a custom-namespace QuakeML extra."""
    return {'value': str(value), 'namespace': _NS}


def _x_num(value):
    """Same, but None-safe for possibly-missing numeric columns."""
    return None if value is None or pd.isna(value) else _x(value)


def _extras(pairs):
    """Build an `extra` dict from (name, value) pairs, dropping the empty ones."""
    return {name: wrapped for name, wrapped in pairs if wrapped is not None}


# ---------------------------------------------------------------------------
# Maximum-likelihood origin
# ---------------------------------------------------------------------------

def _maxlike_time(expectation_time, ot_sec):
    """
    Rebuild the full ML origin time from its seconds-within-the-minute.

    NLLoc writes only `OT <sec>` on the .hyp MAXIMUM_LIKELIHOOD line: it calls
    hypotime2hrminsec before overwriting the origin time with the expectation's,
    so the ML hour and minute never reach the file (GridLib.c:3513). The minute is
    therefore taken from the expectation origin, and the result is snapped to the
    nearest minute — the two solutions can straddle a minute boundary, and over a
    500-event test run the raw seconds differed by up to 59 s for exactly that
    reason.
    """
    minute_start = UTCDateTime(expectation_time.year, expectation_time.month,
                               expectation_time.day, expectation_time.hour,
                               expectation_time.minute)
    candidate = minute_start + ot_sec
    offset    = candidate - expectation_time
    if offset > 30.0:
        candidate -= 60.0
    elif offset < -30.0:
        candidate += 60.0
    return candidate


def _maxlike_origin(row, public_id, expectation_origin):
    """
    The maximum-likelihood hypocentre as a secondary Origin, or None if absent.

    Deliberately carries no OriginQuality, OriginUncertainty or arrivals. Under
    `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION` every one of those is re-solved at the
    expectation — RMS, Gap, Dist and the residuals belong to that point, and the
    confidence ellipsoid is a second moment about it. Attaching them here would
    recreate exactly the mode-with-someone-else's-error mismatch this stage exists
    to remove.
    """
    if pd.isna(row.get('maxlike_latitude')):
        return None

    origin = Origin(
        resource_id  = ResourceIdentifier(f'smi:pyrenees/origin/{public_id}/maxlike'),
        time         = _maxlike_time(expectation_origin.time, float(row['maxlike_ot_sec'])),
        latitude     = float(row['maxlike_latitude']),
        longitude    = float(row['maxlike_longitude']),
        # Metres positive down; negative values are events above sea level.
        depth        = float(row['maxlike_depth']) * 1e3,
        depth_type   = 'from location',
        creation_info = CreationInfo(author='NonLinLoc SSST'),
    )
    origin.extra = _extras([('locationEstimator', _x('maximum_likelihood'))])
    return origin


# ---------------------------------------------------------------------------
# Station code resolution
# ---------------------------------------------------------------------------

def load_station_epochs(path):
    """
    Map each internal station code to the real stations it stands for.

    The .obs picks carry the project's unified code (e.g. 'FR.0041'), which is
    stored as `alternate_code` in GLOBAL_inventory.xml. 80 of ~2000 codes cover
    more than one real station (near-duplicates merged within 20 m), so each maps
    to a list of epochs rather than a single station.

    Candidates are ordered by the merge priority -- permanent network first, then
    RESIF/RENASS, then a station carrying a real elevation, the XX placeholder for
    uncalled or unknown networks last -- with the oldest start_date breaking ties.
    That is the same rule that decides which station lends the code its label and
    which one NonLinLoc is given the coordinates of, so all three agree. It only
    decides cases where several candidates match the same date -- disjoint epochs
    still resolve on the date alone.

    Returns
    -------
    dict[str, list[(network_code, station_code, start_date, end_date, elevation)]]
    """
    inventory = read_inventory(path)
    epochs    = defaultdict(list)
    for network in inventory:
        for station in network:
            if station.alternate_code:
                epochs[station.alternate_code].append(
                    (network.code, station.code, station.start_date, station.end_date,
                     station.elevation)
                )
    return {code: sorted(entries,
                         key=lambda entry: (station_priority(entry[0], entry[4], code),
                                            entry[2] or UTCDateTime(datetime.max)))
            for code, entries in epochs.items()}


def resolve_station(epochs, unified_code, when):
    """
    Resolve a unified station code to a real (network, station) pair.

    Picks the first candidate whose epoch contains `when`, and falls back to the
    first candidate when none matches (open-ended entries have start_date and
    end_date None). Candidates arrive in merge-priority order from
    load_station_epochs(), so the permanent network -- and never the XX placeholder
    -- wins whenever several cover the pick.

    The stations that were *not* picked are returned alongside, including those
    whose epoch does not cover the pick. They are exported with the pick so that
    waveform retrieval can fall back on them: the stations behind one unified
    code sit within 20 m of each other, so if the resolved one turns out to be
    the wrong label, the others are the places to look.

    Returns
    -------
    ((str, str), list[(str, str)]) or None
        The resolved (network, station) and the other candidates, in inventory
        order. None when the code is absent from the inventory.
    """
    candidates = epochs.get(unified_code)
    if not candidates:
        return None

    names = [(network_code, station_code) for network_code, station_code, _, _, _ in candidates]

    for index, (_, _, start, end, _) in enumerate(candidates):
        if (start is None or when >= start) and (end is None or when <= end):
            return names[index], names[:index] + names[index + 1:]

    return names[0], names[1:]


def _ambiguous_codes(epochs):
    """
    Unified codes whose real station cannot be recovered from the dates.

    Most codes covering several stations have disjoint epochs and resolve
    cleanly; a few list entries that all span all time, so the ordering decides
    (the merge priority, then inventory order). The stations behind one code sit
    within 20 m of each other, so this only picks between near-identical names.
    """
    def spans_all(epoch):
        _, _, start, end, _ = epoch
        return (start is None or start.year <= 1900) and (end is None or end.year >= 2500)

    return [code for code, entries in epochs.items()
            if len(entries) > 1 and sum(spans_all(e) for e in entries) > 1]


# ---------------------------------------------------------------------------
# Confidence ellipsoid
# ---------------------------------------------------------------------------

def _confidence_ellipsoid(row):
    """
    Convert the NLLoc confidence ellipsoid to its QuakeML representation.

    NLLoc reports azimuth/dip for axes 1 and 2 only; the third is their cross
    product (same construction as merge_regional_results._build_covariance).
    QuakeML instead requires the orientation of the *major* axis plus a rotation
    about it, so the three axes are sorted by length and re-expressed.

    All lengths carry NLLoc's own 3-DOF 68% chi-square scaling, hence the
    confidence_level of 68 reported alongside.
    """
    axes = [
        _ellipsoid_axis_to_xyz(row['EllipsoidAz1'], row['EllipsoidDip1'], row['EllipsoidLen1']),
        _ellipsoid_axis_to_xyz(row['EllipsoidAz2'], row['EllipsoidDip2'], row['EllipsoidLen2']),
    ]
    unit1 = axes[0] / row['EllipsoidLen1']
    unit2 = axes[1] / row['EllipsoidLen2']
    third = np.cross(unit1, unit2)
    axes.append(third / np.linalg.norm(third) * row['EllipsoidLen3'])

    # Sort by semi-axis length: minor, intermediate, major.
    lengths = [float(np.linalg.norm(axis)) for axis in axes]
    order   = np.argsort(lengths)
    minor, intermediate, major = (axes[i] for i in order)

    major_unit = major / np.linalg.norm(major)
    east, north, down = major_unit
    azimuth = np.degrees(np.arctan2(east, north)) % 360.0
    plunge  = np.degrees(np.arcsin(np.clip(down, -1.0, 1.0)))

    # Rotation of the minor axes about the major one, measured from the phi = 0
    # reference: the horizontal direction perpendicular to the major azimuth.
    azimuth_rad = np.radians(azimuth)
    reference   = np.array([np.cos(azimuth_rad), -np.sin(azimuth_rad), 0.0])
    orthogonal  = np.cross(major_unit, reference)
    intermediate_unit = intermediate / np.linalg.norm(intermediate)
    rotation = np.degrees(np.arctan2(float(np.dot(intermediate_unit, orthogonal)),
                                     float(np.dot(intermediate_unit, reference)))) % 360.0

    return ConfidenceEllipsoid(
        semi_major_axis_length        = float(np.linalg.norm(major)) * 1e3,
        semi_intermediate_axis_length = float(np.linalg.norm(intermediate)) * 1e3,
        semi_minor_axis_length        = float(np.linalg.norm(minor)) * 1e3,
        major_axis_azimuth            = float(azimuth),
        major_axis_plunge             = float(plunge),
        major_axis_rotation           = float(rotation),
    )


# ---------------------------------------------------------------------------
# Usability verdict
# ---------------------------------------------------------------------------

def classify(row):
    """
    Decide whether an event's location PDF is trustworthy enough to use.

    See the module docstring for why C68 is compared against its simulated null
    and why Psi never contributes.

    Returns
    -------
    (bool, str) — usable flag, and a comma-joined reason string ('' when usable)
    """
    reasons = []

    if pd.isna(row['n_scat']):
        reasons.append('no_pdf_metrics')
    else:
        c68_z = (row['C68'] - _NOMINAL_COVERAGE) / row['C68_sigma_n']
        if c68_z < _C68_Z_MIN:
            reasons.append('c68_overconfident')
        if bool(row['dip_reject']):
            reasons.append('dip_bimodal')

    return not reasons, ','.join(reasons)


# ---------------------------------------------------------------------------
# Event construction
# ---------------------------------------------------------------------------

def _build_pick(line, index, public_id, epochs, unresolved):
    """
    Convert one .obs pick line into a (Pick, Arrival) pair.

    Returns None when the line is malformed, which never happens on a bulletin
    written by this pipeline (all 1e6 lines hold exactly 19 tokens).
    """
    parts = line.split()
    if len(parts) <= _P_PICK_ORIGIN:
        logger.warning(f"{public_id}: malformed pick line, skipped: {line!r}")
        return None

    # The seconds field is added rather than passed to UTCDateTime so a value
    # of 60 or more (legal in some source bulletins) cannot raise.
    date = parts[_P_DATE]
    hhmm = parts[_P_HHMM]
    time = (UTCDateTime(f"{date[0:4]}-{date[4:6]}-{date[6:8]}T{hhmm[0:2]}:{hhmm[2:4]}:00")
            + float(parts[_P_SECONDS]))

    unified_code = parts[_P_CODE]
    resolved     = resolve_station(epochs, unified_code, time)
    if resolved is None:
        unresolved[unified_code] += 1
        network_code, station_code = unified_code.split('.', 1) if '.' in unified_code \
                                     else ('', unified_code)
        alternates = []
    else:
        (network_code, station_code), alternates = resolved

    channel     = parts[_P_CHANNEL]
    pick_origin = parts[_P_PICK_ORIGIN]

    pick = Pick(
        resource_id  = ResourceIdentifier(f'smi:pyrenees/pick/{public_id}/{index}'),
        time         = time,
        time_errors  = QuantityError(uncertainty=float(parts[_P_ERR_MAG])),
        waveform_id  = WaveformStreamID(
            network_code = network_code,
            station_code = station_code,
            channel_code = None if channel == 'None' else channel,
        ),
        phase_hint   = parts[_P_REAL_PHASE],
        # TEMP_* picks come from PhaseNet, everything else from analyst bulletins.
        evaluation_mode = 'automatic' if pick_origin.startswith('TEMP_') else 'manual',
    )
    pick.extra = _extras([
        ('unifiedCode',    _x(unified_code)),
        ('pickOrigin',     _x(pick_origin)),
        ('relativeTiming', _x('true' if parts[_P_INS] == '*' else 'false')),
        # Fallbacks for waveform retrieval, omitted when the code is unique.
        ('alternateStations',
         _x(','.join(f'{net}.{sta}' for net, sta in alternates)) if alternates else None),
    ])

    arrival = Arrival(
        resource_id = ResourceIdentifier(f'smi:pyrenees/arrival/{public_id}/{index}'),
        pick_id     = pick.resource_id,
        phase       = parts[_P_PHASE],
    )
    return pick, arrival


def build_event(obs_event, row, verdict, epochs, unresolved):
    """
    Build one QuakeML event from an .obs block and its CSV row.

    Hypocentre, uncertainties and quality come from the CSV at full precision;
    the .obs supplies only the magnitude and the picks.

    Parameters
    ----------
    obs_event  : temp_picks.match_picks.Event
    row        : pd.Series — the SSST_result.csv row for the same publicId
    verdict    : (bool, str) — the classify() result for that row
    epochs     : dict — from load_station_epochs()
    unresolved : Counter — accumulates unified codes absent from the inventory

    Returns
    -------
    obspy.core.event.Event
    """
    public_id = obs_event.public_id
    usable, reason = verdict

    event = Event(
        resource_id   = ResourceIdentifier(f'smi:pyrenees/event/{public_id}'),
        event_type    = 'earthquake',
        creation_info = CreationInfo(author='Shallow_Depth_DL_Catalog'),
    )
    event.comments = [Comment(text=f"usable={str(usable).lower()} reason={reason or 'none'}")]
    event.extra = _extras([
        ('publicId',     _x(public_id)),
        ('usable',       _x(str(usable).lower())),
        ('rejectReason', _x(reason)),
        ('sourceZone',   _x(row['source'])),
        ('nScat',        _x_num(row['n_scat'])),
        ('J',            _x_num(row['J'])),
        ('Psi',          _x_num(row['Psi'])),
        ('C68',          _x_num(row['C68'])),
        ('C68Z',         _x_num((row['C68'] - _NOMINAL_COVERAGE) / row['C68_sigma_n']
                                if not pd.isna(row['n_scat']) else None)),
        ('JNullP95',     _x_num(row['J_null_p95'])),
        ('C68SigmaN',    _x_num(row['C68_sigma_n'])),
        ('dipStat',      _x_num(row['dip_stat'])),
        ('dipPval',      _x_num(row['dip_pval'])),
        ('dipSeparationKm', _x_num(row['dip_sep_km'])),
        ('dipReject',    _x(str(bool(row['dip_reject'])).lower())),
    ])

    picks    = []
    arrivals = []
    stations = set()
    for index, line in enumerate(obs_event.picks):
        built = _build_pick(line, index, public_id, epochs, unresolved)
        if built is None:
            continue
        pick, arrival = built
        picks.append(pick)
        arrivals.append(arrival)
        stations.add((pick.waveform_id.network_code, pick.waveform_id.station_code))

    origin = Origin(
        resource_id  = ResourceIdentifier(f'smi:pyrenees/origin/{public_id}'),
        time         = UTCDateTime(row['date-time']),
        latitude     = float(row['latitude']),
        longitude    = float(row['longitude']),
        # Metres positive down; negative values are events above sea level.
        depth        = float(row['depth']) * 1e3,
        depth_errors = QuantityError(uncertainty=float(row['true_erz']) * 1e3),
        depth_type   = 'from location',
        quality      = OriginQuality(
            used_phase_count       = int(row['Nphs']),
            associated_phase_count = len(picks),
            used_station_count     = len(stations),
            standard_error         = float(row['RMS']),
            azimuthal_gap          = float(row['Gap']),
            minimum_distance       = float(row['Dist']) / _KM_PER_DEGREE,
        ),
        # The DOF scaling is deliberately mixed here: the ellipsoid axes keep
        # NLLoc's 3-DOF scaling, while horizontal_uncertainty and depth_errors
        # carry the catalog's own 2-DOF true_erh and 1-DOF true_erz (see
        # merge_regional_results._compute_true_erh). Neither is a projection of
        # the other, but those are the values the rest of the catalog uses.
        origin_uncertainty = OriginUncertainty(
            horizontal_uncertainty             = float(row['true_erh']) * 1e3,
            min_horizontal_uncertainty         = float(row['minHorUnc']) * 1e3,
            max_horizontal_uncertainty         = float(row['maxHorUnc']) * 1e3,
            azimuth_max_horizontal_uncertainty = float(row['azMaxHorUnc']),
            confidence_level                   = _NOMINAL_COVERAGE * 100,
            preferred_description              = 'horizontal uncertainty',
            confidence_ellipsoid               = _confidence_ellipsoid(row),
        ),
        creation_info = CreationInfo(author='NonLinLoc SSST'),
    )
    origin.arrivals = arrivals
    origin.extra = _extras([
        ('locationEstimator', _x('expectation')),
        ('pdfVolume',       _x_num(row['pdfVolume'])),
        ('ellipsoidVolume', _x_num(row['ellipsoidVolume'])),
        ('ellipsoidAz1',    _x_num(row['EllipsoidAz1'])),
        ('ellipsoidDip1',   _x_num(row['EllipsoidDip1'])),
        ('ellipsoidLen1',   _x_num(row['EllipsoidLen1'])),
        ('ellipsoidAz2',    _x_num(row['EllipsoidAz2'])),
        ('ellipsoidDip2',   _x_num(row['EllipsoidDip2'])),
        ('ellipsoidLen2',   _x_num(row['EllipsoidLen2'])),
        ('ellipsoidLen3',   _x_num(row['EllipsoidLen3'])),
    ])

    # The expectation origin is preferred; the maximum-likelihood one rides along
    # so a consumer that wants the mode is not forced back to the NLLoc output.
    maxlike = _maxlike_origin(row, public_id, origin)

    event.origins             = [origin] if maxlike is None else [origin, maxlike]
    event.picks               = picks
    event.preferred_origin_id = origin.resource_id

    # Magnitude is carried through verbatim from the .obs header, including the
    # 5010 events reading exactly 0.00. Those are OMP placeholders, not
    # measurements: the raw OMP .mag files hold the literal string '0.0' for
    # ~19% of events, the neighbouring 0.1 bin holds only 148, and the zero
    # fraction drops from ~40% to 0.3% in 2013 when OMP started computing ML
    # systematically. They are kept as-is because magnitudes are recomputed in a
    # later stage -- but they must not be fitted as data.
    header    = obs_event.header_line[2:].split()
    magnitude = Magnitude(
        resource_id    = ResourceIdentifier(f'smi:pyrenees/magnitude/{public_id}'),
        mag            = float(header[_H_MAG]),
        magnitude_type = header[_H_MAG_TYPE],
        origin_id      = origin.resource_id,
        creation_info  = CreationInfo(author=header[_H_MAG_AUTHOR]),
    )
    event.magnitudes             = [magnitude]
    event.preferred_magnitude_id = magnitude.resource_id

    return event


# ---------------------------------------------------------------------------
# Chunked serialization
# ---------------------------------------------------------------------------

def _split_envelope(text):
    """
    Split a serialized QuakeML document into (head, body, tail).

    head holds the XML declaration, the root element with its namespace
    declarations, and the opening <eventParameters> tag; body holds the events;
    tail closes both. Concatenating one head, many bodies and one tail yields a
    valid document -- obspy has no streaming writer, so this is how the catalog
    is written without holding all events in memory at once.
    """
    opening = text.index('>', text.index('<eventParameters')) + 1
    closing = text.rindex('</eventParameters>')
    return text[:opening], text[opening:closing], text[closing:]


def _serialize(events):
    """Serialize a list of events to a QuakeML string."""
    buffer = io.BytesIO()
    Catalog(events).write(buffer, format='QUAKEML', nsmap=_NSMAP)
    return buffer.getvalue().decode('utf-8')


class _QuakeMLWriter:
    """
    Assemble one QuakeML file from independently serialized chunks.

    The envelope of the first chunk is kept, every chunk's events are appended,
    and the closing tags are written by close().
    """

    def __init__(self, path):
        self.path     = path
        self.n_events = 0
        self._handle  = open(path, 'w', encoding='utf-8')
        self._tail    = None

    def append(self, serialized, n_events):
        head, body, tail = _split_envelope(serialized)
        if self._tail is None:
            self._handle.write(head)
            self._tail = tail
        self._handle.write(body)
        self.n_events += n_events

    def close(self):
        self._handle.write(self._tail or '')
        self._handle.close()


def year_period(year, span=_SPLIT_YEARS):
    """
    Calendar period a year belongs to, e.g. 2023 -> (2020, 2024).

    A read_events() call on the whole catalog peaks around 15 GB, because every
    event and pick becomes a Python object. The parts exist so a reader never
    has to pay that: each covers a fixed calendar span, so which file holds a
    given event is obvious from its date alone, and stays obvious as the catalog
    grows. Parts are not equal in size — pick density rises from ~7 picks/event
    in the 1980s to ~28 in the 2020s, so the recent ones are the heavy ones.
    """
    start = (year // span) * span
    return start, start + span - 1


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def write_quakeml(parameters, log_dir=None):
    """
    Merge the final .obs bulletin and result CSV into a QuakeML 1.2 file.

    Parameters
    ----------
    parameters : ExportQuakeMLParams
    log_dir    : str, optional — log directory (default: NLL_run/console_output/)

    Returns
    -------
    dict with keys: output, log, n_events, n_picks, n_usable, reasons, unresolved
    """
    log_path = _setup_logger(log_dir or _DEFAULT_LOG_DIR)
    logger.info(f"Obs file   : {parameters.file_obs}")
    logger.info(f"Result CSV : {parameters.file_csv}")
    logger.info(f"Inventory  : {parameters.file_inventory}")
    logger.info(f"Output     : {parameters.save_file}")

    epochs = load_station_epochs(parameters.file_inventory)
    logger.info(f"Station codes in inventory : {len(epochs)}")
    for code in _ambiguous_codes(epochs):
        chosen, *rest = (f'{net}.{sta}' for net, sta, _, _ in epochs[code])
        logger.warning(f"Ambiguous station code {code}: {chosen} used "
                       f"over {', '.join(rest)} (dates do not separate them)")

    _, obs_events = load_bulletin(parameters.file_obs)
    logger.info(f"Events in obs bulletin     : {len(obs_events)}")

    results = pd.read_csv(parameters.file_csv).set_index('publicId')
    logger.info(f"Events in result CSV       : {len(results)}")

    unresolved  = Counter()
    reasons     = Counter()
    n_events    = 0
    n_picks     = 0
    n_usable    = 0
    n_orphan    = 0

    parent = os.path.dirname(parameters.save_file)
    if parent:
        os.makedirs(parent, exist_ok=True)

    # The complete bulletin, plus the two companions that keep a reader from
    # ever having to load all of it at once.
    full_writer    = _QuakeMLWriter(parameters.save_file)
    catalog_writer = (_QuakeMLWriter(_companion_path(parameters.save_file, 'catalog'))
                      if parameters.write_catalog else None)
    part_writers   = {}

    chunk       = []
    chunk_group = None

    def flush():
        """Serialize the pending chunk into every output it belongs to."""
        nonlocal chunk
        if not chunk:
            return
        serialized = _serialize(chunk)
        full_writer.append(serialized, len(chunk))

        if chunk_group is not None:
            writer = part_writers.get(chunk_group)
            if writer is None:
                writer = _QuakeMLWriter(_companion_path(
                    parameters.save_file, f'{chunk_group[0]}_{chunk_group[1]}'))
                part_writers[chunk_group] = writer
            writer.append(serialized, len(chunk))

        if catalog_writer is not None:
            # Picks and arrivals are 83% of the load cost, so the catalog-only
            # companion drops them. The chunk is discarded right after, so the
            # events are stripped in place rather than copied.
            for event in chunk:
                event.picks = []
                for origin in event.origins:
                    origin.arrivals = []
            catalog_writer.append(_serialize(chunk), len(chunk))

        chunk = []

    for obs_event in obs_events:
        public_id = obs_event.public_id
        if public_id is None or public_id not in results.index:
            logger.warning(f"publicId {public_id!r} absent from result CSV — skipping")
            n_orphan += 1
            continue

        row     = results.loc[public_id]
        verdict = classify(row)
        event   = build_event(obs_event, row, verdict, epochs, unresolved)

        usable, reason = verdict
        reasons[reason or 'usable'] += 1
        n_usable += usable
        n_picks  += len(event.picks)
        n_events += 1

        # The bulletin is chronological, so a part is a contiguous run of
        # events; flushing on a part change keeps every chunk single-part.
        group = (year_period(event.origins[0].time.year)
                 if parameters.write_parts else None)
        if chunk and (group != chunk_group or len(chunk) >= _CHUNK_SIZE):
            flush()
        chunk_group = group
        chunk.append(event)

    flush()

    outputs = [full_writer] + ([catalog_writer] if catalog_writer else []) \
              + [part_writers[key] for key in sorted(part_writers)]
    for writer in outputs:
        writer.close()

    logger.info(f"Events written             : {n_events}")
    logger.info(f"Picks written              : {n_picks}")
    logger.info(f"Usable                     : {n_usable}")
    logger.info(f"Unusable                   : {n_events - n_usable}")
    for reason, count in reasons.most_common():
        logger.info(f"  {reason:<20s} : {count}")
    logger.info(f"Obs events with no CSV row : {n_orphan}")
    logger.info(f"CSV rows with no obs event : {len(results) - n_events - n_orphan}")
    if unresolved:
        logger.warning(f"Unresolved station codes   : {len(unresolved)}")
        for code, count in unresolved.most_common():
            logger.warning(f"  {code} ({count} picks)")
    else:
        logger.info("Unresolved station codes   : 0")

    logger.info("Files written:")
    for writer in outputs:
        size_mb = os.path.getsize(writer.path) / 1e6
        logger.info(f"  {os.path.basename(writer.path):<28s} "
                    f"{writer.n_events:6d} events  {size_mb:8.1f} MB")

    return {
        'files':      [{'path': w.path, 'n_events': w.n_events,
                        'size': os.path.getsize(w.path)} for w in outputs],
        'output':     parameters.save_file,
        'log':        log_path,
        'n_events':   n_events,
        'n_picks':    n_picks,
        'n_usable':   n_usable,
        'reasons':    dict(reasons),
        'unresolved': dict(unresolved),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Export the final relocated catalog as a QuakeML 1.2 bulletin.'
    )
    parser.add_argument('--obs',       required=True, help='Final .obs bulletin (obs/SSST_result.obs)')
    parser.add_argument('--csv',       required=True, help='Merged result CSV (RESULT/SSST_result.csv)')
    parser.add_argument('--inventory', required=True, help='Station inventory (stations/GLOBAL_inventory.xml)')
    parser.add_argument('--output',    required=True, help='Output QuakeML file (RESULT/FINAL.xml)')
    parser.add_argument('--no-catalog', action='store_true',
                        help='Skip the picks-free <output>_catalog.xml companion')
    parser.add_argument('--no-parts',   action='store_true',
                        help='Skip the year-range parts')
    parser.add_argument('--log-dir',   default=None,
                        help='Log directory (default: NLL_run/console_output/)')
    args = parser.parse_args()

    write_quakeml(ExportQuakeMLParams(
        file_obs       = args.obs,
        file_csv       = args.csv,
        file_inventory = args.inventory,
        save_file      = args.output,
        write_catalog  = not args.no_catalog,
        write_parts    = not args.no_parts,
    ), log_dir=args.log_dir)


if __name__ == '__main__':
    main()
