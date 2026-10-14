# Consolidated Review — 2026-10-01

## Scope
Every finding of the 15 weekly audits (`weekly_reviews/AUDIT_REPORT_2026-06-21.md`
→ `2026-09-27.md`) and of `review/full-audit-2026-07-12.md` was re-checked
against HEAD (`4a7dbba`).

The next campaign restarts the whole workflow from stage 1 (station inventory),
so this review keeps only what **breaks, crashes or corrupts results on that
fresh run**:
- Anything wrong only because an on-disk product predates a fix is left out.
- Latent defects (code wrong, current inputs don't trigger it) are in an appendix.
- Figure-script defects have their own section.
- Performance and cleanup items are dropped (they remain in the original reports).

Reach was measured with read-only scripts on the **raw inputs**
(`org_catalogs/`, `temp_picks/pick_files/`), which a fresh run consumes
unchanged, or on outputs whose format the code reproduces. Every count below
comes from those checks, not from the original reports.

Status tags: **still open** (unchanged since the report), **escalated**
(measured impact higher than reported), **new** (not in any report).

## 1. Pipeline — breaks or corrupts on a fresh run (stage order)

### fetch_obs/OMP.py
- **[High] escalated** (from 2026-07-19, Low). Lines 203-206: a phase line whose
  P quality is `4` is skipped whole (`continue`), so the S reading on the same
  line is lost too. On the raw archives: **38 825 S picks (17 489 events) in
  `OMP_78-19.mag` + 3 219 (1 267 events) in `OMP_2016.mag`**, every one with S
  quality < 4 and a time present. The S time is absolute (pick hh:mm from
  cols 25:29, seconds from 105:110), not an offset from the placeholder P. It is
  physically sound: e.g. ATE at 20.8 km, S 6.25 s after a 1979-06-02 origin.
  This affects the whole 1978–2019 record, which is all OMP.
  *Fix:* on P quality 4, skip only the P write and still parse S.
- **[Medium] new.** Lines 225 and 261: a negative pick second is rolled back
  from the **event's** minute (`minute_p = minute - 1`) instead of the pick's
  own (`minute_p - 1`). In `OMP_78-19.mag` 21 P picks have a negative second,
  and all 21 have a pick minute ≠ event minute, so each is written 60 s early.
  Example, 1989-07-06 04:41:52.92: ESCF `42 -5.74` → written 04:40:54.26
  instead of 04:41:54.26. In the last run 2 such picks survived to
  `SSST_result.obs` / `FINAL.xml`, ~55 s before the origin. The S branch
  (line 261) has the same bug but 0 occurrences in the raw files.
  *Fix:* `minute_p -= 1` / `minute_s -= 1`.

### fetch_obs/LDG.py
- **[Medium] new.** Lines 146-156 copy `arrtime` verbatim, and LDG's own export
  dates some picks **one day ahead** while hour, minute and second are right.
  The raw line `ATE;11/02/2021 22:51:38.40;Sg;…;timeres 0.32` belongs to an
  event on 10/02/2021 22:51:16; LDG's 0.32 s residual shows only the date is
  wrong. That is **93 picks in 77 events** in `LDG_20-25.obs`, all ~86 400 s off.
  EDT locates around them, but they inflate the published RMS
  (`PYRENEES_052660`: RMS 4 912 s; `PYRENEES_050434`: 1 071 s, both with tight
  ellipsoids) and are exported to `FINAL.xml`.
  *Fix:* when a pick lies ≈ 1 day from its origin, rebuild its date from the
  origin date (or drop picks outside a sanity window).

### build_global_bulletin.py · global_obs/fuse_bulletins.py
- **[High] new.** The stage-3 file selection picks up the previous run's
  outputs. `remap` (`build_global_bulletin.py:48`) and fusion (`:150`) both take
  `obs/*.obs`. Fusion's only exclusion is `'GLOBAL' not in p`
  (`fuse_bulletins.py:961`). So on a restart with `obs/` untouched:
  1. `remap` rewrites all 22 matching files in place. On already-remapped files
     it deletes the picks (the non-idempotence noted in CLAUDE.md §3).
  2. Fusion then ingests **`NLL_result.obs`, `NLL_result_augmented.obs` and
     `SSST_result.obs`** as three extra source bulletins.

  The result is relocated copies of the same events matched against their
  originals (or duplicated where they moved beyond the thresholds), plus
  external temp picks leaking into `GLOBAL.obs` before §4.
  *Fix:* use the explicit six-source list already written out at
  `build_global_bulletin.py:131-136` for both remap and fusion.

### NLL_run/reformate_obs.py · run_SSST.py
- **[High] new.** `reformate_obs.py:140` creates `obs/nlloc_obs/GLOBAL_<N>/` with
  `exist_ok=True` and never empties it (zone 1 currently holds 11 518 files).
  `run_SSST.py:193` then locates every `*.nlloc_obs` in it. The new fusion
  reassigns `PYRENEES_%06d` in chronological order, so last campaign's files
  whose ID isn't rewritten in that zone are relocated alongside the new ones.
  They enter the zone CSVs, can win the cross-zone dedup (lowest `pdfVolume`),
  and are then re-matched by `publicId` to a **different** event's magnitude and
  picks in `SSST_result.obs`.
  *Fix:* empty `outputDir` before writing, on a fresh start only
  (`iteration_start == 0`, `run_SSST.py:157`).

### NLL_run/export_quakeml.py
- **[High] still open** (from 2026-09-06). Line 731:
  `for net, sta, _, _ in epochs[code]` unpacks 4 targets from the 5-tuples that
  `load_station_epochs` returns. Running `_ambiguous_codes` on the current
  inventory returns `['FR.0013']`. The fresh inventory is built from the same
  StationXMLs, so `final_steps.py` will raise `ValueError: too many values to
  unpack` before writing any QuakeML. (Today's `FINAL.xml`, dated 08-24, was
  written before the change that broke this line.)
  *Fix:* `for net, sta, *_ in epochs[code]`.

### temp_picks/convert_picks.py
- **[Low] still open** (from 2026-07-05). Lines 459-469 (`convert_temp_oth`):
  the S block is not gated on the P quality, contrary to the docstring. In raw
  `pick_files/merged_other.txt`, **~148 S picks** come from quality-4 P lines,
  103 of them anchored on the `0.00` placeholder. The `match_picks` residual
  band probably rejects most of them, but they should not be generated.
  *Fix:* gate the S block on `quality_p < 4 or quality_p == 9`.

## 2. Figure scripts — wrong on a fresh run

### complem_figures/zone_map.py
- **[Medium] escalated** (from 2026-06-21, latent). Line 102 loads zones from
  `<root>/loc/`, which no longer exists (renamed `run/nll_loc/`). `_ZONES` is
  `{}`, so the zone overview map silently draws no zones.
  *Fix:* point it at `run/nll_loc`, and raise if nothing is found.

### complem_figures/plot_pdf_cloud_2d.py
- **[Low] still open** (from 2026-09-13, then `plot_pdf_cloud_report.py`).
  Line 660: the `--cmap` help text says "Colormap for --color pdf", but the
  `--color` choices are `log|linear|none`.

## 3. Appendix — latent defects, not triggered by current inputs
The code is wrong, but measured on today's inputs it never fires (evidence in
brackets). Fix them opportunistically when touching the file.

- `fuse_bulletins.py:274,511`: unbounded loops in `find_pick_lines` / `_add_phases_to_lines` [every `.obs` ends in `\n\n`]. `:687`: `mean([])` when every magnitude is NaN [0 NaN]. `:1312,1314`: CLI loose defaults 30 km / 10 s instead of 50 / 30 [the orchestrator passes 50 / 30].
- `merge_regional_results.py:235-248`: duplicate log uses `sort_values().iloc[0]`, selection uses `idxmin` [0 ties in 2 664 duplicate groups]. `:91`: no zero-length-axis guard [0]. `:221-222`: no guard for a missing or empty zone CSV.
- `pdf_metrics.py:1215-1217`: `dip_reject` keeps an event whose `true_erz` is NaN [0 NaN].
- `fetch_obs/OMP.py:154-163`: `23:59:60` cascades to hour 24 after the clamp [0 in raw; 1 hour ≥ 24 is clamped correctly].
- `generate_regional_runfiles.py:157,161`: `_find_station_info` on a `None` / empty lookup [remap removes unknown codes upstream].
- `apply_magnitude_models.py:187`: `StopIteration` is reported as "model missing".
- `filter_events_by_aoi.py:241`: non-atomic in-place rewrite.
- `remap_picks_to_unified_codes.py:138-153`: a lone candidate is accepted outside its epoch [not measured].
- `reformate_obs.py:97`: `publicId` collisions overwrite silently [IDs unique by construction].
- `generate_ssst_runfiles.py:208`: `nz_ls` hard-coded [covers the 761-node LOCGRID].
- `generate_regional_runfiles.py:183`: fixed 4-line header copy [matches the format].
- Token-count guards: `filter_distant_picks.py:148-149`, `match_picks.py:147`, `export_quakeml.py:626-629`, `ssst_corrections.py:244` [all formats fixed].
- `run_zone.py:29,89`: hard-coded `_NLL_BIN`, `sys.exit` inside a library function [works on this machine].
- `run_SSST.py:295`: `--zones` values not stripped [only bites with spaces].
- `depth_maps.py:84`: `period_start + 4` ignores `time_range` [only 5 is used].
- `error_maps.py:262`: crash on an all-NaN field; `pdf_metrics.py:1062,1067`: empty voxel panel [default periods populated].
- `ssst_evolution.py:49,167`: one shared `_YLIM` [plotted median/IQR lie inside it]; `:136`: log of non-positive values.
- `plot_pdf_cloud.py:194-201`: `.hdr` not checked [all 322 980 present].
- `plot_pdf_cloud_2d.py`: empty cloud crashes; title semi-axes mix scalings (`:520`) [only with a non-default `--confidence`]; unused `to_local` (`:337`).
- `event_maps.py:237-239`: count printed before `dropna`; `:228`: silent `source_filter` when the column is missing; malformed `.obs`/CSV columns [no NaN coordinates].
- `event_ranking.py:163,219`: log/ratio of `pdfVolume = 0` [4 in `NLL_result.csv`, none survive to SSST]; `:306,358,399`: `makedirs('')` on a bare filename; all-NaN `corr()`.
- `temporary_network_impact.py:352,370`: `np.median` on NaN [0 all-unresolved events].
- `network_history.py:49` vs `:279`: `PERMANENT` set wider than its legend [G/GE/IU/MN/WM codes exist but carry 0 picks].

Dropped as performance or cleanup only: `run_ssst` re-copying chunks and
`copytree`, `depth_maps` cell loop, `station_colocation` double sort, the unused
`ref` model in `build_theoretical_tables`, dead imports and constants, stale
comments and docstrings, missing-file messages, `parse_nll_output` chdir.

## 4. Closed since the reports
- `fuse_bulletins` `found_possible` and match-frame performance
- code-map block leak in `_build_alternate_code_map`
- `remap` concat and O(P×S) lookup
- `filter_distant_picks` per-pick lookup
- `append_station_delays` +3 offset (the old offset was also checked correct for the existing `stat_totcorr` layout)
- OMP non-digit qualities and second-overflow cascade
- AOI block-wise rebuild
- magnitude models selected by key
- elevation API timeout
- `error_maps` window length and performance
- `event_ranking` gridmap loop
- `plot_pdf_cloud` help text, `KeyError` guard and cross-zone label
- `depth_histogram_compare` maxlike default (2026-09-27 High)
- `event_maps` `--solution` guard on `.obs`
- `finalize_nll_catalog.py` and `nll_phase_1.py` deleted

## Suggested fix order
1. Stage-3 file selection and `nlloc_obs` cleanup — they decide whether the restart is clean.
2. OMP and LDG parsers — before the fetch.
3. `export_quakeml.py:731` — one line.
4. TEMP_OTH S gating, then the figure scripts.

## Summary
16 reports, ~85 findings re-checked:
- **7 pipeline defects affect a fresh run:** 4 High, 2 Medium, 1 Low. Of those, 4 are new and 1 is escalated.
- **2 figure defects**.
- **~25 latent defects**, in the appendix.
- **15 closed**.
