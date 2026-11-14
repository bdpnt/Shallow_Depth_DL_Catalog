# Weekly Code Audit — 2026-11-01

## Files reviewed
- NLL_run/crop_quakeml.py (new, 757321e)
- complem_figures/basemap_lines.py (821ad6a)
- complem_figures/error_maps.py (821ad6a)
- fetch_inventory/_fill_missing_elevations.py (e9ab793)
- fetch_inventory/merge_station_inventories.py (e9ab793)

## Findings

### fetch_inventory/merge_station_inventories.py
- **[High]** Lines 250 / 308 / 538–543: `_arbitrate_group_elevation` calls `_get_elevation` (OpenTopoData) from `_combine_close_stations`, which runs unconditionally — `parameters.fill_elevations` is checked only around `_fill_station_elevations`. So `--no-fill-elevations` / `fill_elevations=False` no longer gives an offline run: every disagreeing group still triggers an HTTP request (10 s timeout each when offline, plus a 1 s sleep). The module docstring and CLAUDE.md both still state this flag skips "the only step of the pipeline that reaches the network". Fix: pass `parameters.fill_elevations` into `_combine_close_stations` and skip arbitration (log the spread as a warning, keep declared values) when it is False — or give arbitration its own flag.
- **[Low]** Line 244 / 262: when the representative is chosen on a higher-ranked term (permanent network) but carries no elevation, while two *other* members disagree by > 100 m, arbitration assigns it `chosen` and counts it as "changed by EU-DEM arbitration" rather than as an in-group fill. Behaviour is reasonable; only the logged counts are mislabelled. Consider counting `best[1].elevation == 0` separately.

### fetch_inventory/_fill_missing_elevations.py
- **[Low]** Line 40: `int(response['results'][0]['elevation'])` truncates toward zero (EU-DEM returns floats, e.g. 659.9 → 659). Use `round()` — sub-metre, but the truncation is systematic (always biases down for land stations).
- **[Low]** Line 38: a quota/rate-limit response (HTTP 429, possibly non-JSON) is caught by the broad `except` and returns 0, which every caller treats as "lookup failed". That is the same silent-zero failure mode that motivated leaving Open-Elevation; with the new 1 000 calls/day cap, a re-run on the same day (156 fills + arbitrations each time) can exhaust it. Suggest logging `response.status_code` explicitly and/or aborting after N consecutive failures instead of degrading every station silently.

### NLL_run/crop_quakeml.py
- **[Low]** Lines 91–101: if `--output` equals `--input`, the input is truncated by `open(output_path, 'w')` after the envelope is read, and `iterparse` then reads an empty/partial file — the bulletin is destroyed. Add `if os.path.abspath(input_path) == os.path.abspath(output_path): raise ValueError(...)`.
- **[Low]** Line 88 vs CLI: the empty-box check rejects `lon_min >= lon_max`, so a box crossing the antimeridian is impossible — irrelevant for the Pyrenees, noted only because the error message calls it "empty".

### complem_figures/error_maps.py
- **[Low]** Line 244: `draw_relief` is called once per panel (ERH and ERV), so the `.npz` is re-read, decompressed and re-hillshaded twice per page, and once more per page in `--time-range` batch mode. Cheap at 15'' but avoidable: load + shade once (e.g. `functools.lru_cache` on a helper returning `rgb, extent`) and only `imshow` per axis.

### complem_figures/basemap_lines.py
- No genuine issues. (`draw_relief` assumes ascending `lat`/`lon`, which `load_earth_relief` guarantees for the cache it writes.)

## Summary
5 files reviewed, 7 issues found (1 High, 6 Low).
