"""
build_theoretical_tables.py
============================
Compute theoretical P- and S-wave travel-time tables using Pyrocko's `cake` CLI.

For each combination of velocity model, source depth, and epicentral distance,
the script queries `cake arrivals` and extracts the first-arriving P and S times.
It then reduces the results across the "min" and "plu" models to obtain envelopes
(lower/upper bounds) representing the travel-time uncertainty due to velocity
variation.

Outputs
-------
- CSV file  : columns [distance, tp_low, tp_high, ts_low, ts_high]
- PNG figure : two-panel plot of P-wave and S-wave arrival bands (seaborn style)

Requirements
------------
- Pyrocko must be installed and `cake` must be available on the system PATH.
- Velocity model files (.nd format) must exist under temp_picks/models/.

Usage
-----
    python temp_picks/build_theoretical_tables.py

    # Override output paths
    python temp_picks/build_theoretical_tables.py \\
        --output temp_picks/tables_Pyr.csv \\
        --figure-output temp_picks/figures/tables_Pyr.png
"""

import argparse
import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

_MODULE_DIR     = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_OUTPUT = os.path.join(_MODULE_DIR, 'tables_Pyr.csv')
_DEFAULT_FIGURE = os.path.join(_MODULE_DIR, 'figures', 'tables_Pyr.png')

_MODELS = {
    # "ref" is the nominal model (100% velocities); stored for reference, not used in envelope.
    'ref': os.path.join(_MODULE_DIR, 'models', 'model_Pyr_100.nd'),
    'min': os.path.join(_MODULE_DIR, 'models', 'model_Pyr_95.nd'),
    'plu': os.path.join(_MODULE_DIR, 'models', 'model_Pyr_105.nd'),
}
_DEPTHS    = [0.0, 30.0]
_DISTANCES = list(np.arange(0, 100.5, 0.5))


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_arrivals(model, depth, distance):
    """
    Query `cake arrivals` for the first P and S arrival times.

    Runs the Pyrocko `cake` command-line tool as a subprocess and parses its
    fixed-width text output to extract arrival times.

    Parameters
    ----------
    model    : str   — path to the velocity model file (.nd format)
    depth    : float — source depth in km
    distance : float — epicentral distance in km

    Returns
    -------
    arrP : float — first P-wave arrival time in seconds, or np.nan if not found.
    arrS : float — first S-wave arrival time in seconds, or np.nan if not found.

    Notes
    -----
    The `cake arrivals` output is a fixed-width table. The relevant columns are:
      - Characters 14–20 : arrival time in seconds
      - Character  41    : phase code (P/p or S/s)
    These positions were determined from the cake output format and must match
    the installed version of Pyrocko.
    """
    cmd = [
        "cake", "arrivals",
        "--sdepth", str(depth),
        "--distances", str(distance),
        "--model", model,
        "--phase", "p,s,P,S"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Warning: cake exited with code {result.returncode} "
              f"(model={model}, depth={depth}, dist={distance})")

    arrP = np.nan
    arrS = np.nan

    for line in result.stdout.splitlines():
        if len(line) < 50:
            continue
        phase_char = line[41].lower()
        try:
            arrival_time = float(line[14:21])
        except (ValueError, IndexError):
            continue
        if phase_char == 'p':
            if np.isnan(arrP) or arrival_time < arrP:
                arrP = arrival_time
        elif phase_char == 's':
            if np.isnan(arrS) or arrival_time < arrS:
                arrS = arrival_time

    return arrP, arrS


def _get_times(models, depths, distances):
    """
    Compute P- and S-wave arrival times for all model/depth/distance combinations.

    Parameters
    ----------
    models    : dict         — mapping of model key to .nd file path
    depths    : list[float]  — source depths in km
    distances : list[float]  — epicentral distances in km

    Returns
    -------
    TP : dict — nested dict of P-wave times: TP[model_key][depth] = np.array
    TS : dict — nested dict of S-wave times: TS[model_key][depth] = np.array
    """
    TP = {m: {} for m in models}
    TS = {m: {} for m in models}

    for m in models:
        for depth in depths:
            tp = []
            ts = []
            print(f"  Computing arrivals — model: {m:>4s}  depth: {depth:>5.1f} km")
            for dist in distances:
                p, s = _get_arrivals(models[m], depth, dist)
                tp.append(p)
                ts.append(s)
            TP[m][depth] = np.array(tp)
            TS[m][depth] = np.array(ts)

    return TP, TS


def _compute_envelopes(distances):
    """
    Build a DataFrame of P/S arrival-time envelopes across all models and depths.

    The envelope (low/high bounds) is the element-wise minimum and maximum of
    the "min" and "plu" model arrivals, taken across all provided depths.

    Parameters
    ----------
    distances : list[float] — epicentral distances in km

    Returns
    -------
    pd.DataFrame
        Columns: distance, tp_low, tp_high, ts_low, ts_high.
    """
    TP, TS = _get_times(_MODELS, _DEPTHS, distances)

    # The model keys are velocity perturbations, not time bounds: "min" (95%
    # velocity) yields the *latest* arrivals and "plu" (105%) the earliest, so
    # neither key can be tied to one edge of the band. Take the envelope over
    # every (model, depth) combination.
    stack_p = np.array([TP[m][depth] for m in ('min', 'plu') for depth in _DEPTHS])
    stack_s = np.array([TS[m][depth] for m in ('min', 'plu') for depth in _DEPTHS])

    tp_low,  tp_high = stack_p.min(axis=0), stack_p.max(axis=0)
    ts_low,  ts_high = stack_s.min(axis=0), stack_s.max(axis=0)

    return pd.DataFrame({
        'distance': distances,
        'tp_low':   tp_low,
        'tp_high':  tp_high,
        'ts_low':   ts_low,
        'ts_high':  ts_high,
    })


def _build_figure(time_tables, output):
    """
    Save a two-panel figure showing P-wave and S-wave travel-time bands.

    Parameters
    ----------
    time_tables : pd.DataFrame — output of _compute_envelopes()
    output      : str          — file path for the saved PNG figure
    """
    sns.set_theme(style='whitegrid')
    palette = sns.color_palette('tab10')

    fig, (ax_p, ax_s) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax_p.fill_between(
        time_tables.distance, time_tables.tp_low, time_tables.tp_high,
        color=palette[0], alpha=0.4, label='P-wave band'
    )
    ax_p.set_ylim(0, 30)
    ax_p.set_ylabel('Arrival Time (s)')
    ax_p.set_title('P-wave arrivals — ±5% velocity variation')

    ax_s.fill_between(
        time_tables.distance, time_tables.ts_low, time_tables.ts_high,
        color=palette[1], alpha=0.4, label='S-wave band'
    )
    ax_s.set_ylim(0, 30)
    ax_s.set_xlabel('Distance (km)')
    ax_s.set_ylabel('Arrival Time (s)')
    ax_s.set_title('S-wave arrivals — ±5% velocity variation')

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def save_theoretical_tables(output=None, figure_output=None):
    """
    Compute P/S travel-time envelopes, save CSV and figure.

    Uses the Pyrenees velocity models (±5% variation) and source depths
    0–30 km, at 0.5 km distance steps up to 100 km.

    Parameters
    ----------
    output : str, optional
        Path for the output CSV. Defaults to temp_picks/tables_Pyr.csv.
    figure_output : str, optional
        Path for the output PNG. Defaults to temp_picks/figures/tables_Pyr.png.

    Returns
    -------
    dict
        Summary with keys: 'output', 'figure'.
    """
    output        = output        or _DEFAULT_OUTPUT
    figure_output = figure_output or _DEFAULT_FIGURE

    os.makedirs(os.path.dirname(output),        exist_ok=True)
    os.makedirs(os.path.dirname(figure_output), exist_ok=True)

    time_tables = _compute_envelopes(_DISTANCES)
    _build_figure(time_tables, figure_output)
    time_tables.to_csv(output, index=False)

    print(f"Tables saved to : {output}")
    print(f"Figure saved to : {figure_output}")

    return {'output': output, 'figure': figure_output}


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Compute theoretical P/S travel-time tables using Pyrocko cake.'
    )
    parser.add_argument(
        '--output', default=None,
        help='Output CSV path. Default: temp_picks/tables_Pyr.csv.'
    )
    parser.add_argument(
        '--figure-output', default=None,
        help='Output PNG path. Default: temp_picks/figures/tables_Pyr.png.'
    )
    args = parser.parse_args()
    save_theoretical_tables(args.output, args.figure_output)


if __name__ == '__main__':
    main()
