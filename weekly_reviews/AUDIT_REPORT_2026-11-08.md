# Weekly Code Audit — 2026-11-08

Commits covered: `9c81eaa` (shaded relief on every map), `160adb4` and `0d29ae4` (new `d90_maps.py`).

## Files reviewed
- NLL_run/pdf_metrics.py
- complem_figures/aoi_map.py
- complem_figures/basemap_lines.py
- complem_figures/cross_section.py
- complem_figures/d90_maps.py
- complem_figures/event_maps.py
- complem_figures/event_ranking.py
- complem_figures/profile_location_map.py
- complem_figures/source_depth_maps.py
- complem_figures/ssst_corrections.py
- complem_figures/station_map.py
- complem_figures/zone_map.py
- zone_Arette/plot_arette_nll_map.py
- zone_Arette/plot_cross_section.py

## Findings

### NLL_run/pdf_metrics.py
- **[Low]** Line 71: `from complem_figures.basemap_lines import draw_relief` now runs at module import time. `run_SSST.py:65` imports `pdf_metrics`, so a broken or missing `complem_figures/basemap_lines.py` would stop `run_SSST.py` at startup, before any relocation. CLAUDE.md says figure generation here is non-fatal, and this import goes around that guarantee. Fix: import `draw_relief` lazily inside `_draw_map_panel` (as `aoi_map._render_pygmt` does), so a failure lands in the existing `try/except` around `generate_figures`.

### complem_figures/profile_location_map.py
- **[Low]** Lines 27 and 162: the docstring and the comment say the above-sea-level (magenta) triangle is at the **top** of the colour bar and the deeper-than-`max_depth` one at the bottom. The bar is vertical (`JMR`) with no `+r`, and `makecpt(reverse=True)` reverses the colours, not z. So z increases upward: the background triangle (z < 0, magenta) sits at the **bottom** and the foreground (black) at the top. `source_depth_maps.py`'s "left triangle" is correct for its horizontal bar. Fix: swap top/bottom in both texts (check against a rendered figure first).

### complem_figures/d90_maps.py
- **[Low]** Lines 109–112: the guard only checks for `C68`, but `classify()` also reads `n_scat`, `C68_sigma_n` and `dip_reject`. A partially annotated CSV fails with a bare `KeyError` partway through the row loop, not with the explanatory `ValueError`. Fix: check the full column set `{'n_scat', 'C68', 'C68_sigma_n', 'dip_reject'}`.
- **[Low]** Line 112: `classify(row)` over `df.iterrows()` builds a pandas Series for each of ~47 k rows, which is the slowest part of a ~10 s script. The same pattern is in `seismicity_3d.py:120`. Fix: vectorise as `df.n_scat.notna() & ((df.C68 - 0.6827) / df.C68_sigma_n >= -2) & ~df.dip_reject.astype(bool)`, or add a vectorised `classify_frame` next to `classify` in `export_quakeml.py` so both callers share one definition.

## Summary
14 files reviewed, 4 issues found (all Low). The relief change itself checks out:
- The `_relief_grey` flip is applied and undone symmetrically, so the grid keeps its orientation.
- The crossed `COLOR_BACKGROUND`/`COLOR_FOREGROUND` matches GMT's `-I`, which swaps B/F.
- Every `--no-dem` path keeps the former behaviour.
