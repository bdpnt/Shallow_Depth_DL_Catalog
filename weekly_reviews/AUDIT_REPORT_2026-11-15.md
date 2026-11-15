# Weekly Code Audit — 2026-11-15

Commits covered: `0c8f841`..`971c13a` (PRs #51 and #52): `--mag_version` in `network_history.py`, QAFI faults and profile insets in `event_maps.py`, the new `workflow_diagram.py` and `pick_cross_section.py`, the cross-section location map and grey "no information" fill, `--split` in `error_maps.py`, `--background` in `zone_map.py`, the Okabe-Ito colours in `station_map.py`, and the `failles/` → `faults/` rename.

## Files reviewed
- complem_figures/cross_section.py
- complem_figures/cross_section_panels.py
- complem_figures/error_maps.py
- complem_figures/event_maps.py
- complem_figures/network_history.py
- complem_figures/pick_cross_section.py
- complem_figures/station_map.py
- complem_figures/workflow_diagram.py
- complem_figures/zone_map.py
- zone_Arette/plot_cross_section.py

## Findings

### complem_figures/pick_cross_section.py
- **[Medium]** Lines 300–304: `print_command()` passes only the geometry (`--lon0/--lat0/--azimut/--length/--width/--depth-max`). It leaves out `--uncert-h`, `--uncert-v` and `--solution`, but the picker applies all three to the preview. A profile picked with `--uncert-h 1.5 --uncert-v 1.5` (the docstring's own example) or `--solution maxlike` prints a command that `cross_section_panels.py` runs with the 1.0 km / expectation defaults. The figure then shows a different event set from the one previewed, and nothing warns about it. Fix: append `--uncert-h {p.uncert_h:g} --uncert-v {p.uncert_v:g} --solution {p.solution}` to the printed command. The picker's `--bulletin` has no panels equivalent, so at least note in the output that the panels use their own default `--panel` files.

### complem_figures/cross_section_panels.py
- **[Low]** Line 209: `_location_map` hard-codes `faults=True`. The `.gmt` files it then reads (`faults/QAFI_v4_traces.gmt`, `faults/QFIv4_Debated_traces.gmt`) are in `faults/`, which is gitignored (`.gitignore:5`), and they have to be made by hand with `ogr2ogr`. On a fresh checkout, or before that conversion has run, every `cross_section_panels.py` run fails in `event_maps.generate_figure` *after* the panel PNG/PDF have been written. The traceback blames the location map, not the panels, and there is no way to turn it off. `event_maps.py` makes `--faults` opt-in, but here it is always on. Fix: pass `faults=os.path.exists(event_maps._QAFI_FILE)` (print a note when it is false), or add a `--no-faults` flag.

### complem_figures/cross_section.py
- **[Low]** Lines 493–509: the new NaN split removes events without the metric. If **every** projected event lacks it (for example `--use-err psi` on a format-6 CSV not yet annotated by `pdf_metrics.py`), `X`, `Z` and `cval` become empty. The code still calls `pygmt.makecpt`, `fig.plot` with empty arrays and `fig.colorbar` for a metric that was never drawn, which can make GMT fail or produce a misleading colour bar. Fix: after the split, skip the coloured plot and the colour bar when `cval.size == 0` (the grey layer already shows the events).

### complem_figures/network_history.py
- **[Low]** Lines 213–216 and 281–284 (`--mag_version`):
  - `float(header[10])` raises `ValueError` on any event line whose magnitude field is the literal `None`. `.obs` headers use that placeholder for missing values; `zone_map._read_bulletin` replaces it for exactly this reason.
  - If no event reaches `MAG_MIN`, `erz_big.min()` raises on an empty array before the figure is drawn.
  
  Fix: parse with `float(v) if v != 'None' else np.nan`, and guard the summary print with `if len(erz_big)`.

### complem_figures/workflow_diagram.py
- **[Low]** Line 35: `_PROJECT_ROOT` is defined but never used. The script reads no project file and builds its outputs from `_MODULE_DIR`. Fix: remove the constant.

## Summary
10 files reviewed, 5 issues found (1 Medium, 4 Low). The following checked out:
- The `length='limit'` fix to `cross_section_panels._project`.
- The visible-count print.
- `compass_point` rounding, including 348.75°+ → N.
- `label_profile_ends` left/right orientation.
- The inset CPT reuse in `event_maps._draw_profile_inset`.
- The B/F end-triangle colours of the new `JBC+e` colour bars.
- The `error_maps.py` split/two-panel refactor, which behaves the same as before.
- The `zone_map._read_bulletin` extraction.
- `cross_section.py` importing `event_maps`: no circular import, and `event_maps` keeps PyGMT out of its top level, so `pick_cross_section.py` can import it in `seisbench_env`.
- The `failles/` → `faults/` path rename. It is consistent in `cross_section.py` and `zone_Arette/plot_cross_section.py`, but the local directory has to be renamed to match, since it is untracked.
