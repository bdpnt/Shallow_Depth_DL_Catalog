# Weekly Code Audit — 2026-09-20

## Files reviewed
- complem_figures/event_maps_report.py
- NLL_run/pdf_metrics.py

## Findings

### complem_figures/event_maps_report.py
New file (commit `f89b9cc`): report-ready PyGMT depth map with cross-section
overlays. The swath geometry (`_dest_point` / `profile_geometry`) is a faithful
copy of `cross_section.py`'s selection math, and the `.obs` column drop/rename
indices match `event_maps.py` exactly, so the drawn strip matches the projected
one and there is no off-by-one in the parser. Two minor issues:

- **[Low]** Lines 163-171 — asymmetric guard on `.obs` input. `--column expect_z`
  raises a clear `SystemExit` for a `.obs` bulletin (only the maximum-likelihood
  depth exists there), but the parallel option `--position expect` has no such
  guard: it is only read in the CSV branch (lines 170-171) and is silently
  ignored for a `.obs` file. A user asking for expectation epicentres on a `.obs`
  input gets maximum-likelihood positions with no warning. Suggested fix: add a
  symmetric check in the `.obs` branch, e.g. `if parameters.position != 'maxlike':
  raise SystemExit('--position is meaningful only for a RESULT/*.csv ...')`.

- **[Low]** Lines 192-197 — the "Quality filter applied — N events kept" count is
  printed before `dropna(subset=['Latitude', 'Longitude', 'Depth'])` on line 197,
  so the reported count overstates the number of events actually plotted whenever
  coordinates are missing. Cosmetic, but the printed figure should match what is
  drawn. Suggested fix: drop the NaN coordinate rows before the count is printed,
  or print the final `len(events_df)` after the `dropna`.

### NLL_run/pdf_metrics.py
Change in commit `75f073f`: adds a `('Dist', 'Nearest station distance (km)',
'linear')` entry to `_ERROR_PREDICTORS` and updates two docstrings. No issues.
The subplot grid derives its width from `n_cols = len(_ERROR_PREDICTORS) + 1`
(line 910), so the extra panel is placed automatically with no hardcoded column
count to fall out of sync. `_plot_error_panel` drops NaN on the predictor
(`dropna(subset=[predictor, error_col])`, line 816) and guards the empty case, so
the new `Dist` column is handled robustly; `Dist` is present in the SSST/NLL
result CSVs alongside the existing `RMS`/`Gap` predictors.

## Summary
2 files reviewed, 2 issues found (0 High, 0 Medium, 2 Low).
