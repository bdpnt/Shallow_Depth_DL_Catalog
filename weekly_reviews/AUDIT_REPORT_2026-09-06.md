# Weekly Code Audit — 2026-09-06

## Files reviewed
- NLL_run/export_quakeml.py
- NLL_run/generate_regional_runfiles.py
- complem_figures/cross_section.py
- complem_figures/station_colocation.py
- fetch_inventory/merge_station_inventories.py

## Findings

### NLL_run/export_quakeml.py
- **[High]** Line 731: `chosen, *rest = (f'{net}.{sta}' for net, sta, _, _ in epochs[code])`
  unpacks each `epochs[code]` entry into **4** targets, but as of commit `ed62e9c`
  every entry is a **5**-tuple `(network, station, start_date, end_date, elevation)`
  (see `load_station_epochs`, lines 291–294). This raises
  `ValueError: too many values to unpack (expected 4)` as soon as
  `_ambiguous_codes(epochs)` returns any code — i.e. whenever the inventory holds a
  unified code covering several stations that all span all time, which the docstring
  and the warning path treat as a normal, expected case. The export aborts before
  writing any QuakeML.
  Fix: `for net, sta, _, _, _ in epochs[code]` (5 targets). The neighbouring
  unpackings (`resolve_station`, `_ambiguous_codes.spans_all`) were already updated
  to 5; this one was missed.

### NLL_run/generate_regional_runfiles.py
- **[Medium]** Lines 157 & 161 (`_find_station_info`, rewritten in `ed62e9c`):
  `alternate_code_map.get(alternate_code)` is called with no default and then
  iterated. If a station code present in `GLOBAL.obs` picks is absent from
  `GLOBAL_code_map.txt`, `.get` returns `None` and the `for` loop raises
  `TypeError: 'NoneType' object is not iterable`. Similarly, if the code is in the
  map but `inventory.select(...)` yields nothing, `candidates` stays empty and
  `min(candidates, ...)` (line 161) raises `ValueError: min() arg is an empty
  sequence`. Both surface as opaque errors mid–GTSRCE generation rather than a clear
  "station code X not found" message. Suggested fix: guard for a missing/empty
  lookup and log the offending `alternate_code` (skip it or raise a descriptive
  error).

### complem_figures/station_colocation.py
- **[Low]** Lines 261 & 267 (`_panel_sensitivity`): `np.sort(pair_distances)` is
  computed twice on the full upper-triangle pair array (~n²/2 ≈ 2.9M elements for
  ~2400 stations). Sort once into a local (e.g. `ordered = np.sort(pair_distances)`)
  and reuse it for both the `searchsorted` over `radii` and the single
  `at_threshold` lookup.

## Summary
5 files reviewed, 3 issues found (1 High, 1 Medium, 1 Low).
