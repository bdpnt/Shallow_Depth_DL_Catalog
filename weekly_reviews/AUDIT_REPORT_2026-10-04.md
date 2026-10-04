# Weekly Code Audit — 2026-10-04

Scope: Python files modified on the `claude` branch in the last 7 days.
Two commits account for all changes:

- `df5990b` — *refactor: fold the report figure scripts into complem_figures modules*
  (the `*_report.py` scripts become suffix-less modules; most of the diff is renames
  plus shared-helper imports in place of copied code).
- `4a7dbba` — *feat: add --before-year filter to the depth histogram figures*.

All changed files are figure-generation scripts under `complem_figures/`; none touch
the relocation pipeline, the catalog, or any published value.

## Files reviewed
- complem_figures/aoi_map.py
- complem_figures/basemap_lines.py
- complem_figures/cross_section_panels.py
- complem_figures/depth_histogram.py
- complem_figures/depth_histogram_compare.py
- complem_figures/depth_histogram_report.py *(deleted — folded into depth_histogram.py)*
- complem_figures/error_maps.py
- complem_figures/error_maps_report.py *(deleted — folded into error_maps.py)*
- complem_figures/error_vs_quality.py
- complem_figures/event_maps.py
- complem_figures/event_maps_report.py *(deleted — folded into event_maps.py)*
- complem_figures/magnitude_regression.py
- complem_figures/network_history.py
- complem_figures/plot_pdf_cloud_2d.py
- complem_figures/profile_location_map.py
- complem_figures/source_depth_maps.py
- complem_figures/ssst_corrections.py
- complem_figures/ssst_corrections_report.py *(deleted — folded into ssst_corrections.py)*
- complem_figures/temporary_network_impact.py
- complem_figures/temporary_network_impact_report.py *(deleted — folded into temporary_network_impact.py)*

## Findings

### complem_figures/error_maps.py
- **[Medium]** Lines 263–264: `np.max(valid)` / `np.mean(valid)` are computed on
  `valid = median[np.isfinite(median)]` with no empty-array guard. `windowed_median`
  only fills a cell when `keep.sum() >= _MIN_N` (10 events in the window, line 181),
  so a period that has events but never reaches 10 in any ~10 km window yields an
  all-NaN `median` and an empty `valid`. `np.max([])` then raises
  `ValueError: zero-size array to reduction operation maximum which has no identity`.
  The `if period.empty: continue` guard (line 313) does **not** cover this case — the
  period is non-empty, just sparse.
  Failure scenario: the pipeline driver `generate_complem_figures.py` runs this in
  batch over 1976–2025 in 5-year periods on `RESULT/NLL_result.csv`. An early sparse
  period (e.g. 1976–1980, a handful of events) raises, and because the exception
  propagates out of `plot_period` → `generate_figure`, the **entire** error-maps batch
  aborts before the dense, important later periods are rendered. (`np.mean([])` only
  warns; the hard crash is `np.max`.) The pre-refactor code was also vulnerable, but
  this commit rewrote the line and the gap remains.
  Suggested fix: guard the stats, e.g. `if valid.size:` print the mean/max, else print
  that 0 cells reached the `_MIN_N` minimum and the map is empty.

### complem_figures/ssst_corrections.py
- **[Low]** Lines 697–698: `_draw` now pins the frame exactly to the grid
  (`set_xlim(lons[0] - 0.5*(lons[1]-lons[0]), lons[-1] + 0.5*(lons[1]-lons[0]))`, same
  for y). Two minor consequences:
  (a) *Cosmetic regression* — the old `_draw_map` autoscaled to include the station
  triangle plotted at `sta_lon/sta_lat`. With the default per-station framing
  (`_map_grid` pads only 0.15° around the station's event cloud), an edge-of-network
  station whose events all lie to one side can sit just outside the grid, and its
  reference triangle is now silently clipped off-frame. Affects only the per-station
  atlas pages; shared-extent / `_PAGE_EXTENT` pages still cover the whole region.
  (b) *Latent `IndexError`* — `lons[1]` / `lats[1]` assume ≥2 nodes. A degenerate
  `--extent` with equal bounds (e.g. `--extent=0,1,43,43`) makes `np.arange` return a
  single-element axis and `lats[1]` raises. Not reachable through normal use, only via
  invalid `--extent` input.
  Suggested fix: derive the half-cell from the node spacing with a `len(...) > 1` guard,
  and if the station marker must stay visible, expand the limits to include
  `sta_lon/sta_lat` (or draw it with `clip_on=False`).

### complem_figures/depth_histogram.py
- **[Low]** `generate_figure` (title at ~line 144): when `--before-year` excludes every
  event, `depths` is empty, `median = depths.median()` is `NaN`, and the title renders
  `median = nan km`. No crash (an empty histogram is drawn), so impact is cosmetic, but
  it is an unhandled empty-selection boundary. Suggested fix: if `n_total == 0`, note
  "no events before <year>" in the title instead of formatting `NaN`.

## Verified correct (not reported)
- `--before-year` parsing in both depth-histogram files: `depth_histogram.py` uses
  `lstrip('# ')` + split (Year at field 0, Dep at field 8) and `depth_histogram_compare.py`
  uses raw `split()` (Year at field 1, Dep at field 9) — each internally consistent with
  its own depth index; the CSV branch filters on `date-time` year correctly.
- `aoi_map.py` `_cross` / `_kept` sign convention and all four boundary-line definitions
  match `global_obs/filter_events_by_aoi._is_in_aoi` exactly.
- Cross-module imports after the refactor all resolve: `event_maps` exports
  (`EventMapsParams`, `REPORT_REGION`, `_read_catalog`, `profile_geometry`,
  `hypocentre_columns`), `pdf_metrics` symbols pulled by `error_vs_quality.py`
  (`_ERROR_ROWS`, `_ERROR_PREDICTORS`, `_ERROR_CONTROL`, `_error_ylim`, `_partial_spearman`,
  `_quantile_bins`, `_CLIP_PCT`) and `read_scat`, `cross_section._dest_point`, and the
  `plot_pdf_cloud` helpers imported by `plot_pdf_cloud_2d.py`.
- `error_maps.py` period-end fix (`range(1976, 2026, time_range)`, inclusive
  `{end}-12-31 23:59:59.999999`) and the vectorized `windowed_median` rewrite match the
  original edge convention; the CSV is read once per `generate_figure`, not per period.
- `ssst_corrections.py` numerical core (`correction_field`, `event_correction_spread`,
  `event_displacement`, `_scan_bulletin`) is unchanged by the commit; the
  `_plot_timeline` return-type change to a list is matched at its only caller.

## Summary
20 files reviewed (14 current modules + 6 folded/deleted `*_report.py` siblings),
3 issues found: 1 Medium, 2 Low. No High-severity defects. All issues are confined to
figure-generation scripts; none affect the catalog or any published hypocentre/metric.
