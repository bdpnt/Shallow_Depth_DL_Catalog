# Weekly Code Audit — 2026-09-27

## Files reviewed
- complem_figures/aoi_map_report.py
- complem_figures/basemap_lines.py
- complem_figures/cross_section_report.py
- complem_figures/depth_histogram_compare_report.py
- complem_figures/error_maps_report.py
- complem_figures/error_vs_quality_report.py
- complem_figures/network_history_report.py
- complem_figures/profile_location_map.py
- complem_figures/source_depth_maps_report.py
- complem_figures/ssst_corrections_report.py
- complem_figures/temporary_network_impact_report.py
- global_obs/magnitude_regression_report.py

All twelve are new report-figure siblings added in commits `995ca80` and
`930b7d8` (merged as PR #48). None modifies an existing pipeline module; each
reuses the parsing/geometry of its non-report counterpart.

## Findings

### complem_figures/depth_histogram_compare_report.py
- **[High]** Lines 62 and 153 — the default maximum-likelihood column is wrong,
  so panel (c) silently duplicates panel (b). `generate_figure(..., column='expect_z',
  column_ml='depth', ...)` and `--column-ml` default `'depth'` both point at the
  CSV's `depth` column. But in `RESULT/*.csv` the `depth` column **is** the PDF
  expectation and equals `expect_z` exactly — confirmed in
  `NLL_run/merge_regional_results.py:133-137` ("the CSV's latitude/longitude/depth
  are identical to its expect_lat/expect_lon/expect_z") and CLAUDE.md §NLL output.
  The maximum-likelihood depth lives in `maxlike_depth` (added by
  `_MAXLIKE_COLUMNS`, `merge_regional_results.py:181-182`, and read as
  `maxlike_depth` by `export_quakeml.py:253`). With the current default, panel
  (b) uses `expect_z` and panel (c) uses `depth == expect_z`, so the green
  "maximum-likelihood" panel and its printed "median max-likelihood" statistic
  are identical to the orange expectation panel — the exact comparison the figure
  exists to show is lost, and a reader would wrongly conclude the two estimators
  give the same depth distribution. The docstring bug is baked in too: line 11
  describes panel (c) as "maximum-likelihood hypocentres (default column: depth)".
  Suggested fix: default `column_ml='maxlike_depth'` in both the function
  signature (line 62) and the argparse default (line 153).

### complem_figures/error_maps_report.py
- **[Low]** Lines 252-254 — crashes on an all-empty median field. When no
  9×9 window reaches `_MIN_N = 10` events (a sparse `--start`/`--end` period, or a
  narrow region), `windowed_median` returns an all-NaN grid, so `valid =
  median[np.isfinite(median)]` is empty and `np.max(valid)` raises `ValueError:
  zero-size array to reduction operation maximum` (with `np.mean(valid)` also
  emitting a RuntimeWarning). The figure would otherwise draw fine (an empty map).
  Suggested fix: guard the summary print, e.g. only compute mean/max when
  `valid.size` and print "0 cells" otherwise.

### complem_figures/network_history_report.py
- **[Low]** Line 49 vs. lines 278-279 — the `PERMANENT` set and the panel (b)
  legend disagree on what the blue curve counts. `PERMANENT = {"FR", "RA", "RD",
  "ES", "CA", "LC", "G", "GE", "IU", "MN", "WM"}` drives `contrib_perm`, but the
  curve is labelled "of which permanent networks (FR, RA, RD, ES, CA, LC)" — six
  networks, not eleven. If any contributing unified code carries a `G`/`GE`/`IU`/
  `MN`/`WM` prefix (plausible for ORFEUS/GFZ-sourced codes), the plotted
  "permanent" curve exceeds what the legend claims to count, with no cue to the
  reader. (Note also that `AM`, which CLAUDE.md §1 counts as permanent, is
  deliberately absent here and instead shown as a semi-permanent deployment in
  panel (c) — that exclusion is self-consistent with the legend; the extra five
  global-network codes are the mismatch.) Suggested fix: either trim `PERMANENT`
  to the six networks named in the legend, or extend the legend text to match the
  set actually counted.

## Notes (reviewed, no action)
- `aoi_map_report.py`: the half-plane sign convention (`_cross`/`_kept`) and the
  Sutherland-Hodgman clip (`_clip_halfplane`) are internally consistent, and the
  shaded band is correctly the intersection of the two universal lines with the
  RESIF and IGN/ICGC source lines.
- `magnitude_regression_report.py`: the `mtype not in line` substring test in
  `extract_headers` (line 89) can over-match (e.g. `"ML"` inside `"MLv"`), but it
  is reproduced verbatim from `generate_magnitude_models.py` and the source
  bulletins are single-magtype-dominant, so no divergence from the pipeline.
- `basemap_lines.py`, `cross_section_report.py`, `error_vs_quality_report.py`,
  `profile_location_map.py`, `source_depth_maps_report.py`,
  `ssst_corrections_report.py`, `temporary_network_impact_report.py`: parsing
  indices, colour-scale rules and reused helpers match their non-report
  counterparts; no correctness issues found.

## Summary
12 files reviewed, 3 issues found (1 High, 0 Medium, 2 Low).
