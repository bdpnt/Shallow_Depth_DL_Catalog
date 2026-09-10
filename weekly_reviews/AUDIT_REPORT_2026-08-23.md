# Weekly Code Audit — 2026-08-23

## Files reviewed
- NLL_run/export_quakeml.py
- NLL_run/pdf_metrics.py
- NLL_run/merge_regional_results.py
- complem_figures/event_ranking.py
- complem_figures/plot_pdf_cloud.py
- final_steps.py
- run_SSST.py

Window: 2026-08-16 → 2026-08-23 (`claude` branch). Commits reviewed:
`eb9d348`, `8670ae8`, `3948603`, `8321bfa`, `1a7879b`, `39bc77e`, `397c13d`.
This week added the QuakeML export stage (`export_quakeml.py`, `final_steps.py`),
the location-PDF quality metrics (`pdf_metrics.py`) and wired them into
`run_SSST.py`, plus a shared gridding helper reused by `event_ranking.py`.

## Findings

### NLL_run/pdf_metrics.py
- **[Low]** Line 83: `_GAUSS_ENTROPY_3D = 0.5 * 3 * np.log(2 * np.pi * np.e)` is
  defined but never referenced. `negentropy()` recomputes the same quantity
  inline at line 239 (`0.5 * d * np.log(2 * np.pi * np.e)`). Dead constant
  introduced with the new file. *Suggested fix:* delete the constant, or have
  `negentropy()` use it (the inline form is `d`-generic, so deletion is the
  cleaner option).
- **[Low]** Lines 783 / 788 (`_generate_voxel_html`): when no voxel reaches
  `_VOXEL_MIN_COUNT`, `values`/`counts` are empty and both
  `float(np.nanmax(values))` (line 783) and `counts.max()` (line 788) raise on
  the empty selection, so the whole HTML panel is lost rather than rendered
  blank. Only reachable on degenerate/very sparse input, and it is swallowed by
  the `generate_figures` try/except in `compute_metrics` (line 954), so the CSV
  is safe — but the figure fails silently instead of drawing an empty scene.
  *Suggested fix:* early-return (skip the panel) when `lat_i.size == 0`.

### NLL_run/export_quakeml.py
- **[Low]** Line 126: the `write_parts` comment says the year-range parts are
  "sized to `_SPLIT_BUDGET_GB`", but no such constant exists. Parts are actually
  fixed-length calendar spans set by `_SPLIT_YEARS = 5` (see `year_period`,
  line 588). Misleading comment left over from an earlier byte-budget design.
  *Suggested fix:* reword to reference the fixed 5-year `_SPLIT_YEARS` span.
- **[Low]** Lines 521–527 (`build_event`): the magnitude header is parsed
  positionally (`header[_H_MAG]` … `header[_H_MAG_AUTHOR]`, indices 9–11) with no
  length guard, whereas `_build_pick` (line 352) explicitly guards pick-line
  token count before indexing. A short/malformed event header raises IndexError
  instead of being logged and skipped like a bad pick line. The pipeline
  guarantees the 12-token header format, so this is defensive-only, but the
  asymmetry with the pick-line handling is worth closing. *Suggested fix:* guard
  `len(header) > _H_MAG_AUTHOR` and log-and-skip (or emit a magnitude-less
  event) on failure.

### NLL_run/merge_regional_results.py
- No issues. The only change adds `ellipsoidVolume = 4/3·π·L1·L2·L3`; correct and
  consistent with the identical formula in `pdf_metrics.compute_metrics`.

### complem_figures/event_ranking.py
- No issues. The refactor to the shared `windowed_stat_grid` preserves the map
  semantics; dropping non-finite samples up front (new helper) is a mild
  improvement over the previous per-cell `np.median`, not a regression. The PDF
  metrics are merged defensively (`how='left`, only columns actually present).

### complem_figures/plot_pdf_cloud.py
- No issues. Swaps obspy's `read_nlloc_scatter` for `pdf_metrics.read_scat`,
  which correctly recovers the 3 samples per cloud that obspy's `[4:]` header
  slice silently drops.

### run_SSST.py
- No issues. The `compute_metrics` step is correctly wrapped as non-fatal and
  idempotent, matching the file's own docstring and the surrounding stages.

### final_steps.py
- No issues. Validates required inputs before running and reports clearly.

## Summary
7 files reviewed, 4 issues found (0 High, 0 Medium, 4 Low). No correctness or
performance defects that affect catalog output; all findings are dead code,
misleading comments, or defensive-guard asymmetries.
