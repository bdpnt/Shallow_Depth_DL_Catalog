"""
depth_histogram.py
============================
Plot a histogram of earthquake depths from a .obs bulletin.

Depth axis is inverted (0 km at top, seismological convention) with
horizontal bars. Depths recorded as 'None' in the bulletin are ignored.

Usage
-----
    python complem_figures/depth_histogram.py \\
        --bulletin  obs/FINAL.obs \\
        --output    complem_figures/depth_histogram/FINAL.png
"""

import argparse
import os
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class DepthHistogramParams:
    file_bulletin: str
    fig_save:      str
    bin_size:      float = 2.0   # km per bin
    max_depth:     float = 60.0  # y-axis lower limit (km)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    """
    Plot and save the depth histogram for the bulletin.

    Parameters
    ----------
    parameters : DepthHistogramParams

    Returns
    -------
    dict with keys: output
    """
    sns.set_theme()

    with open(parameters.file_bulletin, 'r') as f:
        lines = f.readlines()

    bull = [line.lstrip('# ').rstrip('\n').split()
            for line in lines if line.startswith('# ')]
    bull_df = pd.DataFrame(bull, columns=[
        'Year', 'Month', 'Day', 'Hour', 'Min', 'Sec',
        'Lat', 'Lon', 'Dep', 'Mag', 'MagType', 'MagAuthor',
        'PhaseCount', 'HorUncer', 'VerUncer', 'AzGap', 'RMS',
    ])
    bull_df['Dep'] = pd.to_numeric(bull_df['Dep'], errors='coerce')
    depths = bull_df['Dep'].dropna()

    bins = np.arange(0, parameters.max_depth + parameters.bin_size, parameters.bin_size)

    os.makedirs(os.path.dirname(parameters.fig_save), exist_ok=True)

    fig, ax = plt.subplots(figsize=(5, 8), layout='constrained')
    ax.hist(depths, bins=bins, orientation='horizontal', color=sns.color_palette()[0])
    ax.invert_yaxis()
    ax.set_ylim(parameters.max_depth, 0)
    ax.set_ylabel('Depth (km)')
    ax.set_xlabel('Count')
    ax.set_title('Depth distribution')

    ax.text(0.98, 0.02,
            f'N = {len(depths)}\nMedian = {depths.median():.1f} km',
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9)

    plt.savefig(parameters.fig_save)
    plt.close(fig)

    print(f'Figure saved @ {parameters.fig_save}')
    return {'output': parameters.fig_save}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Plot a histogram of earthquake depths from a .obs bulletin.'
    )
    parser.add_argument('--bulletin', default=os.path.join(_PROJECT_ROOT, 'obs', 'GLOBAL.obs'),
                        help='Input .obs bulletin')
    parser.add_argument('--output',   default=os.path.join(_MODULE_DIR, 'depth_histogram', 'GLOBAL.png'),
                        help='Output figure path')
    parser.add_argument('--bin-size', type=float, default=2.0,
                        help='Bin size in km (default: 2.0)')
    parser.add_argument('--max-depth', type=float, default=60.0,
                        help='Y-axis lower limit in km (default: 60.0)')
    args = parser.parse_args()

    generate_figure(DepthHistogramParams(
        file_bulletin = args.bulletin,
        fig_save      = args.output,
        bin_size      = args.bin_size,
        max_depth     = args.max_depth,
    ))


if __name__ == '__main__':
    main()
