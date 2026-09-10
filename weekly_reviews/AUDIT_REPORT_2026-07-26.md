# Weekly Code Audit — 2026-07-26

## Files reviewed
- add_temp_picks.py
- complem_figures/event_ranking.py
- complem_figures/plot_pdf_cloud.py
- temp_picks/convert_picks.py

Scope: commits on `claude` from the last 7 days —
`17d59f7` (TEMP_STB parquet pick format), `191bd6e` (PDF-cloud volume table),
`eec800e` (event_ranking gridmap), `7d3fe3b` (event_ranking correlation matrix).
Review focused on the diffs introduced by these commits plus their surrounding
context.

## Findings

### complem_figures/event_ranking.py
- **[Medium]** Lines 239–251 (`_add_gridmap_subplot`): the windowed-median grid is
  filled with a nested Python loop over every cell (`bins_lat × bins_lon` ≈
  100 × 288 ≈ 28,800 cells at the default `bin_size=0.02`), and each iteration
  recomputes a fresh pandas boolean mask over the full event `Series`
  (`(lat >= lat_low) & (lat <= lat_high) & (lon >= lon_low) & (lon <= lon_high)`)
  followed by `values[mask]` index-aligned slicing. On the full Pyrenees SSST
  catalog this is O(cells × events) with per-iteration pandas overhead and can
  run into minutes. The cell widths `lat_edges[1]-lat_edges[0]` /
  `lon_edges[1]-lon_edges[0]` are also constants recomputed inside the loop.
  Suggested fix: convert `lat`, `lon`, and `values` to NumPy arrays once before
  the loop (drop pandas alignment), and hoist the two constant cell-width terms
  out of the loop. This is an opt-in diagnostic (`--figures true`), so it is not
  on the critical path, but the slowdown scales poorly with catalog size.

- **[Low]** Lines 295 and 339 (`_generate_gridmap_figure` /
  `_generate_corr_matrix_figure`): `os.makedirs(os.path.dirname(output_path), exist_ok=True)`
  raises `FileNotFoundError: ''` when the user passes a bare filename (no
  directory component) via `--figure-output` / `--corr-matrix-output`, because
  `os.path.dirname('foo.pdf') == ''` and `os.makedirs('', exist_ok=True)` fails.
  The default paths always include a directory, so this only bites on a
  user-supplied bare filename. Suggested fix: guard with
  `if os.path.dirname(output_path):` or use `os.path.dirname(output_path) or '.'`.
  (The same pattern appears at line 380 in `generate_ranking`.)

- **[Low]** Lines 309–320 (`_generate_corr_matrix_figure`): when the
  `--max-erh` / `--max-erz` cutoffs exclude every event (or `ranking_df` is
  empty), `filtered[cols].corr()` returns an all-NaN matrix and the heatmap is
  rendered with `nan` annotations instead of failing or warning. A short guard
  that skips the figure and prints a warning when `filtered` is empty would make
  the degenerate case obvious rather than producing a misleading PDF.

### temp_picks/convert_picks.py
- No blocking issues. The new TEMP_STB row-based path (`convert_temp_stb`,
  `_convert_parquet_dir`) reads each leaf parquet individually to bound memory,
  filters by `phase_score`, and accesses row fields by attribute name (correctly
  guarding against `pd.read_parquet(columns=...)` column-order assumptions). Note
  (informational, not a defect): for the parquet path `n_input` counts *all* rows
  including those later dropped for low `phase_score`, whereas for the line-based
  path `n_input` counts only real pick lines — the two "Input pick lines" log
  values therefore mean slightly different things, but the parquet path logs
  `Rows dropped (low phase_score)` separately so the totals still reconcile.

### add_temp_picks.py
- No issues. The `_converted_output_path` helper correctly redirects converted
  output under `_PICK_FILES` regardless of the (gitignored) source location, and
  the TEMP_STB entries are ordered after TEMP_OMP so `match_picks.py`'s per-event
  dedup absorbs any OMP duplicates. If the `all_picks/PICKS_MARC/*.pq` datasets
  are absent, `_convert_parquet_dir`'s glob yields nothing and the pipeline
  degrades gracefully (empty output, no crash).

### complem_figures/plot_pdf_cloud.py
- No issues. `pdf_volume` / `ellipsoid_volume` are extracted after `row.iloc[0]`,
  so they are scalars by the time `_fmt` applies the `:.3g` format spec, and the
  `None` fallback (`'n/a'`) covers the missing-CSV case.

## Summary
4 files reviewed, 3 issues found (1 Medium, 2 Low) — all confined to the opt-in
figure-generation paths of `event_ranking.py`; no correctness defects in the
pick-ingestion pipeline changes.
