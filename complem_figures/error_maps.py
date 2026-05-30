"""
error_maps.py
============================
Generate per-period ERH/ERV error maps from an NLL result file.

For each time window, computes a windowed-median error grid (horizontal and
vertical) and overlays the raw event scatter. Figures are produced in
parallel using all available CPU cores and saved as PDFs in the output folder.

Usage
-----
    python complem_figures/error_maps.py \\
        --file        RESULT/FINAL.txt \\
        --map-folder  complem_figures/error_maps/ \\
        --time-range  5
"""

import argparse
import os
from dataclasses import dataclass
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from obspy import UTCDateTime

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ErrorMapsParams:
    file:       str
    mapFolder:  str
    time_range: int  # in years


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_file(file):
    """
    Parse an NLL result file into a DataFrame.

    Parameters
    ----------
    file : str — path to the NLL result bulletin

    Returns
    -------
    pd.DataFrame with columns: date, latitude, longitude, erh, erv, time
    """
    with open(file, 'r') as f:
        lines = f.readlines()

    events = []
    for line in lines:
        infos = line.split()
        year  = '19' + infos[0].rjust(2, '0') if float(infos[0]) > 75 else '20' + infos[0].rjust(2, '0')
        date  = UTCDateTime(f'{year}-{infos[1]}-{infos[2]}T00:00:00.00Z')
        events.append([date, float(infos[6]), float(infos[7]),
                        float(infos[12]), float(infos[13])])

    df         = pd.DataFrame(events, columns=['date', 'latitude', 'longitude', 'erh', 'erv'])
    df['time'] = df['date'].apply(lambda x: pd.Timestamp(x.datetime))
    return df


def _filter_dates(events, time_range):
    """
    Slice the events DataFrame into fixed-length time windows.

    Parameters
    ----------
    events     : pd.DataFrame — events with a 'time' column
    time_range : int          — window length in years

    Returns
    -------
    dict[str, pd.DataFrame] — keyed by period label e.g. '1976-1980'
    """
    filtered = {}
    for period_start in range(1976, 2026, time_range):
        period_end = period_start + 4
        mask = (
            (events['time'] >= pd.Timestamp(f'{period_start}-01-01')) &
            (events['time'] <= pd.Timestamp(f'{period_end}-12-31'))
        )
        filtered[f'{period_start}-{period_end}'] = events[mask]
    return filtered


def _add_subplot(events, ax, type):
    """
    Render a windowed-median error grid and event scatter onto a matplotlib axis.

    Parameters
    ----------
    events : pd.DataFrame — events with latitude, longitude, erh/erv columns
    ax     : matplotlib Axes
    type   : str          — error column to plot ('ERH' or 'ERV')

    Returns
    -------
    matplotlib QuadMesh
    """
    vmin, vmax       = 0, 5
    lat_min, lat_max = 42.0, 44.0
    lon_min, lon_max = -2.25, 3.5
    bins_lat, bins_lon = 400, 860

    lat_edges   = np.linspace(lat_min, lat_max, bins_lat + 1)
    lon_edges   = np.linspace(lon_min, lon_max, bins_lon + 1)
    median      = np.zeros((bins_lat, bins_lon))
    count       = np.zeros((bins_lat, bins_lon), dtype=int)
    window_size = 4

    for i in range(bins_lat):
        for j in range(bins_lon):
            lat_low  = max(lat_edges[i]   - window_size * (lat_edges[1] - lat_edges[0]), lat_min)
            lat_high = min(lat_edges[i+1] + window_size * (lat_edges[1] - lat_edges[0]), lat_max)
            lon_low  = max(lon_edges[j]   - window_size * (lon_edges[1] - lon_edges[0]), lon_min)
            lon_high = min(lon_edges[j+1] + window_size * (lon_edges[1] - lon_edges[0]), lon_max)
            mask = (
                (events['latitude']  >= lat_low)  & (events['latitude']  <= lat_high) &
                (events['longitude'] >= lon_low)  & (events['longitude'] <= lon_high)
            )
            window = events[mask]
            if len(window) > 0:
                median[i, j] = np.median(window[type.lower()])
                count[i, j]  = len(window)
            else:
                median[i, j] = np.nan

    median_masked = np.ma.masked_where(count < 10, median)
    mesh = ax.pcolormesh(lon_edges, lat_edges, median_masked,
                         vmax=vmax, vmin=vmin, cmap='rocket_r',
                         shading='auto', alpha=0.9)

    sns.scatterplot(x=events['longitude'], y=events['latitude'],
                    s=0.6, color='black', linewidth=0, ax=ax)

    col = type.lower()
    valid_medians = median[count >= 10]
    ax.text(0.99, 0.98,
            f"Mean Err.: {np.nanmean(valid_medians):.1f}\n"
            f"Max Err.: {np.nanmax(valid_medians):.1f}\n"
            f"Std Err.: {np.nanstd(valid_medians):.1f}\n"
            f"Median (all events): {np.nanmedian(events[col]):.1f}",
            transform=ax.transAxes,
            fontweight='bold', color='black', fontsize=8, ha='right', va='top')
    ax.text(0.01, 0.98, type, transform=ax.transAxes,
            fontweight='bold', color='black', ha='left', va='top')

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    return mesh


def _generate_plot(args, map_folder):
    """Produce and save one ERH/ERV error-map PDF for a single time period."""
    sns.set_theme()
    date, df    = args
    fig, axes   = plt.subplots(2, 1, figsize=(12, 10), layout='constrained')
    _add_subplot(df, axes[0], type='ERH')
    mesh = _add_subplot(df, axes[1], type='ERV')
    fig.colorbar(mesh, ax=[axes[0], axes[1]],
                 label='Median error (km) - 9x9 grid',
                 shrink=0.7, pad=0.025, aspect=50)
    plt.suptitle(f'Horizontal and vertical error maps\n{date}', fontweight='bold')
    plt.savefig(f'{map_folder}{date}.pdf')
    plt.close(fig)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(params):
    """
    Generate and save per-period ERH/ERV error map PDFs in parallel.

    Parameters
    ----------
    params : ErrorMapsParams

    Returns
    -------
    dict with keys: output, n_figures
    """
    events          = _read_file(params.file)
    events_filtered = _filter_dates(events, params.time_range)

    args_list = [(item, params.mapFolder) for item in events_filtered.items()]
    n_cores   = int(cpu_count() / 1.5)
    print(f'Using {n_cores}/{cpu_count()} CPU cores')

    with Pool(processes=n_cores) as pool:
        pool.starmap(_generate_plot, args_list)

    n_figures = len(args_list)
    print(f'Saved {n_figures} figures @ {params.mapFolder}')
    return {'output': params.mapFolder, 'n_figures': n_figures}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate per-period ERH/ERV error maps from an NLL result file.'
    )
    parser.add_argument('--file',       required=True,
                        help='NLL result file (e.g. RESULT/FINAL.txt)')
    parser.add_argument('--map-folder', required=True,
                        help='Output folder for PDF figures')
    parser.add_argument('--time-range', type=int, default=5,
                        help='Time-window length in years (default: 5)')
    args = parser.parse_args()

    generate_figure(ErrorMapsParams(
        file       = args.file,
        mapFolder  = args.map_folder,
        time_range = args.time_range,
    ))


if __name__ == '__main__':
    main()
