"""
station_map.py
============================
Generate a PyGMT map of seismic stations in the Pyrénées.

Reads station locations from a StationXML inventory (FDSN format) and plots
each station as an inverted triangle on a Pyrénées basemap. No station labels
are shown.

Usage
-----
    python complem_figures/station_map.py \\
        --inventory stations/GLOBAL_inventory.xml \\
        --output    complem_figures/station_map/station_map.png
"""

import argparse
import os
import tempfile
import xml.etree.ElementTree as ET
from dataclasses import dataclass

import pygmt as pg


# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR   = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_MODULE_DIR)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class StationMapParams:
    fileInventory: str
    figSave:       str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_inventory(path):
    """Parse a StationXML file and return list of {lat, lon, network} dicts."""
    tree = ET.parse(path)
    root = tree.getroot()
    ns = {'fdsn': 'http://www.fdsn.org/xml/station/1'}

    stations = []
    for network in root.findall('fdsn:Network', ns):
        net_code = network.get('code')
        for station in network.findall('fdsn:Station', ns):
            lat = float(station.find('fdsn:Latitude',  ns).text)
            lon = float(station.find('fdsn:Longitude', ns).text)
            stations.append({'lat': lat, 'lon': lon, 'network': net_code})
    return stations


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_figure(parameters):
    """
    Parse a StationXML inventory and save a PyGMT station map.

    Parameters
    ----------
    parameters : StationMapParams

    Returns
    -------
    dict with key: output
    """
    stations = _read_inventory(parameters.fileInventory)
    print(f"Inventory read @ {parameters.fileInventory} ({len(stations)} stations)")

    _PRIMARY = {'FR', 'RD', 'RA', 'ES', 'CA', 'LC'}

    primary = [s for s in stations if s['network'] in _PRIMARY]
    am      = [s for s in stations if s['network'] == 'AM']
    other   = [s for s in stations if s['network'] not in _PRIMARY and s['network'] != 'AM']

    region = [-4.0, 4.0, 41.0, 45.0]

    fig = pg.Figure()
    with pg.config(MAP_FRAME_TYPE='fancy+'):
        fig.basemap(region=region, projection='M6i', frame='af')
    fig.coast(water='skyblue', land='#777777', resolution='i',
              area_thresh='0/0/1', borders='1/0.75p,black')

    if other:
        fig.plot(x=[s['lon'] for s in other], y=[s['lat'] for s in other],
                 style='i0.12c', fill="#72a04e", pen='0.08p, black', transparency=5)
    if am:
        fig.plot(x=[s['lon'] for s in am], y=[s['lat'] for s in am],
                 style='i0.14c', fill='#d9876f', pen='0.1p, black', transparency=10)
    if primary:
        fig.plot(x=[s['lon'] for s in primary], y=[s['lat'] for s in primary],
                 style='i0.16c', fill="#79b2d4", pen='0.15p, black', transparency=20)

    legend_lines = (
        "S 0.15c i 0.12c #72a04e 0.08p,black 0.3c Temporary networks\n"
        "S 0.15c i 0.14c #d9876f 0.1p,black  0.3c RaspberryShake network\n"
        "S 0.15c i 0.16c #79b2d4 0.15p,black 0.3c Permanent networks\n"
    )
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp:
        tmp.write(legend_lines)
        tmp_path = tmp.name
    fig.legend(spec=tmp_path, position="JBL+jBL+o0.2c", box="+gwhite+p0.5p")
    os.unlink(tmp_path)

    os.makedirs(os.path.dirname(parameters.figSave), exist_ok=True)
    fig.savefig(parameters.figSave, dpi=300)

    print(f"Figure saved @ {parameters.figSave}")
    return {'output': parameters.figSave}


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate a PyGMT map of Pyrénées seismic stations.'
    )
    parser.add_argument('--inventory', required=True,
                        help='Input StationXML inventory file')
    parser.add_argument('--output', required=True,
                        help='Output figure path (PNG or PDF)')
    args = parser.parse_args()

    generate_figure(StationMapParams(
        fileInventory = args.inventory,
        figSave       = args.output,
    ))


if __name__ == '__main__':
    main()
