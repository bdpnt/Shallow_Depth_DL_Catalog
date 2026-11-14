"""
crop_quakeml.py
============================
Crop the final QuakeML bulletin to a latitude/longitude box.

Streams an already exported RESULT/FINAL.xml (see export_quakeml.py) and copies
every event whose preferred origin -- the location-PDF expectation, i.e. the
published hypocentre -- falls inside the box into a new, independently valid
QuakeML file. Bounds are inclusive. Events are kept whatever their `pyr:usable`
flag, which travels with them, as do their picks and maximum-likelihood origin.

The input is walked with lxml iterparse, one event at a time, so memory stays
flat whatever the size of the input -- read_events() on FINAL.xml peaks ~15 GB.

Usage
-----
    python NLL_run/crop_quakeml.py --output RESULT/FINAL_crop.xml \\
        --lat-min 42.5 --lat-max 43.0 --lon-min 0.0 --lon-max 1.0
"""

import argparse
import os
import sys

from lxml import etree

# ---------------------------------------------------------------------------
# Module paths
# ---------------------------------------------------------------------------

_MODULE_DIR    = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT  = os.path.dirname(_MODULE_DIR)
_DEFAULT_INPUT = os.path.join(_PROJECT_ROOT, 'RESULT', 'FINAL.xml')

sys.path.insert(0, _PROJECT_ROOT)

from NLL_run.export_quakeml import _split_envelope

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

_BED = '{http://quakeml.org/xmlns/bed/1.2}'

_HEAD_BYTES = 65536   # enough to hold the XML declaration and the root element
_TAIL_BYTES = 4096    # enough to hold the closing tags


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _envelope(path):
    """Head and tail of a QuakeML file, read from its two ends."""
    with open(path, 'rb') as handle:
        head = handle.read(_HEAD_BYTES)
        handle.seek(max(0, os.path.getsize(path) - _TAIL_BYTES))
        tail = handle.read()
    # Undecodable bytes can only sit at the cut points, inside the event body
    # that _split_envelope discards.
    head, _, tail = _split_envelope(head.decode('utf-8', errors='ignore')
                                    + tail.decode('utf-8', errors='ignore'))
    return head, tail


def _preferred_point(event):
    """(latitude, longitude) of an event element's preferred origin."""
    preferred = event.findtext(f'{_BED}preferredOriginID')
    for origin in event.iterchildren(f'{_BED}origin'):
        if origin.get('publicID') == preferred:
            return (float(origin.findtext(f'{_BED}latitude/{_BED}value')),
                    float(origin.findtext(f'{_BED}longitude/{_BED}value')))
    raise ValueError(f"Event {event.get('publicID')} has no preferred origin")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def crop_quakeml(input_path, output_path, lat_min, lat_max, lon_min, lon_max):
    """
    Write the events of input_path whose preferred origin lies in the box.

    Returns
    -------
    dict with keys: output, n_scanned, n_kept, size
    """
    if lat_min >= lat_max or lon_min >= lon_max:
        raise ValueError(f"Empty box: lat {lat_min}..{lat_max}, lon {lon_min}..{lon_max}")

    head, tail = _envelope(input_path)
    n_scanned  = 0
    n_kept     = 0

    parent = os.path.dirname(output_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    with open(output_path, 'w', encoding='utf-8') as out:
        out.write(head)
        for _, event in etree.iterparse(input_path, events=('end',), tag=f'{_BED}event'):
            n_scanned += 1
            lat, lon = _preferred_point(event)
            if lat_min <= lat <= lat_max and lon_min <= lon <= lon_max:
                out.write(etree.tostring(event, encoding='unicode'))
                n_kept += 1
            # Free the event and everything already parsed before it.
            event.clear()
            while event.getprevious() is not None:
                del event.getparent()[0]
        out.write(tail)

    return {
        'output':    output_path,
        'n_scanned': n_scanned,
        'n_kept':    n_kept,
        'size':      os.path.getsize(output_path),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Crop the final QuakeML bulletin to a latitude/longitude box.'
    )
    parser.add_argument('--input',   default=_DEFAULT_INPUT,
                        help='Exported QuakeML bulletin (default: RESULT/FINAL.xml)')
    parser.add_argument('--output',  required=True, help='Cropped QuakeML file')
    parser.add_argument('--lat-min', required=True, type=float)
    parser.add_argument('--lat-max', required=True, type=float)
    parser.add_argument('--lon-min', required=True, type=float)
    parser.add_argument('--lon-max', required=True, type=float)
    args = parser.parse_args()

    summary = crop_quakeml(args.input, args.output,
                           args.lat_min, args.lat_max, args.lon_min, args.lon_max)
    print(f"Events scanned : {summary['n_scanned']}")
    print(f"Events kept    : {summary['n_kept']}")
    print(f"Output         : {summary['output']} ({summary['size']/1e6:.1f} MB)")


if __name__ == '__main__':
    main()
