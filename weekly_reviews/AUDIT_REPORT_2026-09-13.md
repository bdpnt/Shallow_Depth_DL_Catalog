# Weekly Code Audit — 2026-09-13

Scope: Python files modified in the last 7 days on the `claude` branch.

## Files reviewed
- NLL_run/export_quakeml.py
- NLL_run/merge_regional_results.py
- NLL_run/pdf_metrics.py
- complem_figures/cross_section.py
- complem_figures/depth_histogram_report.py
- complem_figures/plot_pdf_cloud.py
- complem_figures/plot_pdf_cloud_report.py
- complem_figures/temporary_network_impact.py
- temp_picks/build_theoretical_tables.py

## Findings

### complem_figures/temporary_network_impact.py
- **[Low]** Lines 317–330 (`_plot_distance`): the panel-title statistics are computed with
  `np.median(actual["before"])` / `np.median(actual["after"])` and `(actual[...] < _NEAR_KM).mean()`
  on the raw `min_dist_km` arrays. `_scan_bulletin` (line 166) assigns `np.nan` to an event whose
  picks all land on station codes absent from the GTSRCE files. `np.median` propagates NaN (unlike
  pandas `.median()`, which skips it), so a single all-unresolved event would print `nan` in the
  title, while the value returned in `generate_figures`' summary dict (lines 405–406, pandas
  `.median()`) would be a real number — an internal inconsistency. Reachability is very low
  (an event would need *every* station unresolved, and events carry ~20 picks each), so this is a
  latent-robustness note rather than an observed failure. Suggested fix: use `np.nanmedian(...)` and
  compute the "within N km" fractions over the non-NaN subset (mask NaN before the comparisons) so
  the annotations match the skipna summary values.

### temp_picks/build_theoretical_tables.py
- **[Low]** Lines 169–176 (`_compute_envelopes` / `_get_times`): the envelope now correctly
  stacks only the `'min'` and `'plu'` models (a good fix — the previous code wrongly tied the
  `'min'` key to the low time bound). However `_get_times(_MODELS, ...)` still runs the Pyrocko
  `cake` CLI for the `'ref'` (100 %) model across every depth × distance, and those arrivals are
  then discarded (the comprehension hardcodes `('min', 'plu')`). This is wasted external-process
  work on each table build. Pre-existing (not introduced this week), but the recent change made
  `'ref'`'s exclusion explicit. Suggested fix: pass only the envelope models to `_get_times`
  (e.g. `{k: _MODELS[k] for k in ('min', 'plu')}`), keeping the `'ref'` path in `_MODELS` for
  documentation only.

### complem_figures/plot_pdf_cloud_report.py
- **[Low–Medium]** Lines 494 / 498 (`_build_figure`) and `_set_limits`/`_span`: an empty-but-valid
  `.scat` crashes opaquely. `generate_figure` guards only for a *missing* `.scat` (line 665), but
  `read_scat` returns an empty array for a valid file whose header declares 0 samples. With `n == 0`,
  `pdf.max()` (`--color log`/`linear`) and the `v.min()/v.max()` in the limit helpers raise a bare
  `ValueError: zero-size array ...` instead of a clear message. Suggested fix: after reading the
  cloud, `if len(cloud) == 0: raise ValueError(f'{scat_path}: empty scatter cloud')`.
- **[Low]** Lines 597–598: the title mixes two confidence scalings. The horizontal semi-axes are
  NLLoc's fixed `min_hor_unc`/`max_hor_unc` (~68% QML scaling), but the vertical semi-axis is
  `sqrt(k3 * cov[2,2])` with `k3 = chi2.ppf(--confidence, 3)`, as is the drawn ellipse. When
  `--confidence` ≠ 0.6827 the printed horizontal pair no longer describes the plotted ellipse.
  Suggested fix: derive all three semi-axes from `cov` at the chosen `k3`, or state that the
  horizontal pair is NLLoc's fixed value independent of `--confidence`.
- **[Low]** Line 413 (`_build_figure`): `to_local, to_geo = info['to_local'], info['to_geo']`
  unpacks `to_local`, which is never used in the function (the real use is at lines 670–671).
  Suggested fix: `to_geo = info['to_geo']`.
- **[Low]** Line 762: stale `--cmap` help text reads "Colormap for --color pdf", but `--color`
  choices are `log|linear|none` — there is no `pdf` mode. Suggested fix: reword to
  "Colormap for `--color log/linear`".

## Notes (reviewed, no action needed)
- The `0.68 → 0.6827` (`erf(1/sqrt(2))`, one-sigma) coverage change is consistent across
  `export_quakeml.py`, `pdf_metrics.py`, `cross_section.py` and `plot_pdf_cloud.py`; the recomputed
  `c68_z` in `export_quakeml._verdict` (line 427) and the `C68_z` from `pdf_metrics._add_null_columns`
  now use the same constant, so `--usable`/`pyr:usable` selections stay aligned.
- `merge_regional_results.py` deliberately keeps `_S3_3DOF = 3.53` (NLLoc's literal
  `DELTA_CHI_SQR_68_3`, `matrix_statistics.h`) rather than `chi2.ppf(0.6827, 3)`. This is correct:
  the factor exists only to divide back out the scaling the binary actually applied to the ellipsoid
  axes, so it must be the value NLLoc used, not the theoretical one. The 1-DOF/2-DOF factors that
  are applied fresh correctly use `0.6827`.
- `temporary_network_impact.py`'s event-header parse (`line[2:].split()`, lat=`tokens[6]`,
  lon=`tokens[7]`) matches the authoritative `match_picks._parse_event_header` for `NLL_result.obs`
  exactly; empty-month handling, the shared colour scale, and the last-bin clipping are all sound.
- `cross_section.py`: the swath edges are now built with `_dest_point(..., azimut ± 90, largeur_coupe)`
  (perpendicular to the profile, correct for any azimuth — the old due-east/west offset also lacked a
  proper degree conversion), the region is the bounding box of the four corners (valid for a
  southward profile), and the 1:1 panel aspect (`panel_height ∝ depth range / length`) is sound.
- `depth_histogram_report.py`: column indices (Dep=8, PhaseCount=12) match the rest of the pipeline;
  imports all used. Minor code smell only: `line.lstrip('# ')` (line 61) is a character-set strip,
  not a prefix strip, but produces the correct result because every event header begins with a
  numeric year — no action needed.
- `plot_pdf_cloud.py`'s `--confidence` help text was corrected from the stale "(default: 0.9)" to
  match the actual `0.6827` default.

## Summary
9 files reviewed. The week's changes are predominantly a principled, internally consistent
one-sigma-coverage (0.6827) harmonization, a genuine travel-time envelope bug fix, correct
cross-section swath/aspect fixes, and new report-figure variants. 6 issues found — 5 Low and
1 Low–Medium (an opaque crash on an empty `.scat`); no High or Medium correctness bugs.
