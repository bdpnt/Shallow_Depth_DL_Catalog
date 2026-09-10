# Weekly Code Audit — 2026-08-30

## Files reviewed
- NLL_run/export_quakeml.py
- NLL_run/generate_regional_runfiles.py
- NLL_run/generate_ssst_runfiles.py
- NLL_run/merge_regional_results.py
- NLL_run/pdf_metrics.py
- NLL_run/run_ssst.py
- complem_figures/cross_section.py
- complem_figures/plot_pdf_cloud.py
- complem_figures/ssst_corrections.py

Window: 2026-08-23 → 2026-08-30 (`claude` branch). Commits reviewed:
`235438d`, `977c268`, `528c807`, `9381bb1`, `4d758ac`, `041015f`, `e7662a3`,
`dc2d91c`, `e04ea22`, `8bdf5b7`.

This week's theme is the **maximum-likelihood / expectation split**: NLLoc now
runs with `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION` (`generate_regional_runfiles.py`,
`generate_ssst_runfiles.py`, `run_ssst.py`), so the published hypocenter is the
PDF expectation — the point the confidence ellipsoid is actually a second moment
about. `merge_regional_results.py` recovers the ML point into new `maxlike_*`
columns; `export_quakeml.py` publishes it as a secondary origin; `plot_pdf_cloud.py`
draws it as the mode marker; `pdf_metrics.py` retunes the `dip_reject` rule to
require a physical mode separation. Plus one large new diagnostic,
`complem_figures/ssst_corrections.py` (reconstructs the deleted Loc2ssst
correction grids), and a new SSST cross-section format in `cross_section.py`.

## Findings

### complem_figures/ssst_corrections.py
- **[Low]** Line 228 (`parse_hyp_file`): the phase-line branch indexes
  `f_[0]`, `f_[4]`, `float(f_[17])` with no token-count guard. A truncated or
  malformed phase line (fewer than 18 whitespace fields) raises `IndexError`,
  which aborts the parse of the entire `.hyp` file and, because `extract()`
  wraps no per-file try/except, kills the whole cached `(zone, iteration)`
  extraction rather than skipping the bad line. NLLoc's phase format is fixed
  (the `>` separator makes the residual reliably field 17), so this is
  defensive-only — but it mirrors the same guard asymmetry flagged in
  `export_quakeml.build_event` last week. *Suggested fix:* guard
  `len(f_) > 17` (skip the arrival otherwise), matching the guarded pick
  parsing elsewhere in the project.

### NLL_run/pdf_metrics.py
- **[Low]** Lines ~1210-1213 (`compute_metrics`): the new reject rule
  `(df['dip_pval'] < _DIP_ALPHA) & (df['dip_sep_km'] > _DIP_SEP_ERZ_FACTOR * df['true_erz'])`
  silently evaluates to `False` (not rejected) whenever `true_erz` is NaN, since
  any comparison against NaN is False. An event with a genuinely bimodal,
  widely-separated depth PDF but a missing/degenerate `true_erz` therefore flips
  from *rejected* (old p-value-only rule) to *kept/usable*. In the real pipeline
  `true_erz` is written by `merge_regional_results.py` for every located event
  immediately upstream, so this edge is effectively unreachable — reported only
  because the failure mode is a silent *loosening* of the usability gate, the
  direction that matters least visibly. *Suggested fix:* none required for the
  pipeline; if hardening is wanted, treat a NaN `true_erz` on an evaluated cloud
  (`dip_pval` finite) as reject rather than keep.

### NLL_run/merge_regional_results.py
- No issues. `_read_maxlike` parses `MAXIMUM_LIKELIHOOD` fields 2/4/6/8
  (lat/lon/depth/OT) correctly, keys on `PUBLIC_ID`, and raises a clear,
  actionable error (not a silent skip) when the summary `.hyp` is missing or was
  produced without expectation mode — the right call for a multi-day campaign.
  Missing per-event blocks degrade to NaN with a warning.

### NLL_run/export_quakeml.py
- No issues. The `_maxlike_time` reconstruction is sound: NLLoc writes only
  seconds-within-the-minute for the ML origin, and snapping to the nearest
  minute against the expectation origin (±30 s window) correctly handles the
  minute-boundary straddle the docstring describes. The secondary ML origin
  deliberately carries no quality/uncertainty/arrivals, which is consistent —
  those quantities are re-solved at the expectation and would be wrong on the
  mode.

### NLL_run/generate_regional_runfiles.py · generate_ssst_runfiles.py · run_ssst.py
- No issues. The `SAVE_NLLOC_EXPECTATION` additions are consistent across all
  three run-file generators, and the flag ordering respects NLLoc's "last
  LOCHYPOUT wins" rule (`run_ssst.py` still appends the authoritative line with
  the octree/fmamp flags on the final iteration only).

### complem_figures/plot_pdf_cloud.py
- No issues. Ellipsoid centres correctly switch to the expectation
  (`longitude`/`latitude`/`depth`) and the diamond markers to the new
  `maxlike_*` columns; the `_load_nll_reference` `KeyError` guard for CSVs that
  predate the columns is a genuine migration safeguard.

### complem_figures/cross_section.py
- No issues. Format 6 wiring is correct: the `pygmt.project` input dropped from
  5 to 4 columns (`lon, lat, depth, cvals`) and the `pz` convention output is
  read back consistently as `X=col0, Z=col1, cval=col2`. The `--usable` mask
  mirrors `export_quakeml`'s `pyr:usable` rule (C68 z-score ≥ −2, not
  `dip_reject`, metrics present), and the out-of-range `np.clip` before
  `makecpt` correctly avoids the reversed-CPT-turns-black artifact the comment
  cites.

### Cross-cutting note (not a defect)
The switch to `SAVE_NLLOC_EXPECTATION` changes what the result-CSV
`latitude`/`longitude`/`depth` mean for *every* consumer, including unmodified
maps (`event_maps.py`, `depth_maps.py`, `error_maps.py`). This is intended and
correct — those plots now show the published expectation, and the only scripts
that needed the mode were updated to read `maxlike_*`. Recorded here only so the
semantic shift is on the record; no change is warranted.

## Summary
9 files reviewed, 2 issues found (0 High, 0 Medium, 2 Low). Both are narrow,
defensive edge cases at parse/boundary points that the real pipeline never
reaches; no correctness or performance defect affects catalog output. The
maximum-likelihood/expectation split is implemented coherently across the run-file
generators, the merge/rematch stage, the QuakeML export and the diagnostics.
