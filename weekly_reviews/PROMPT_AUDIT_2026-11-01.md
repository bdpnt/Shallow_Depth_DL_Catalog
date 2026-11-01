# Instruction-File Audit — 2026-11-01

Audited at `claude` @ 77df86a.

## Files audited
- CLAUDE.md (80 lines / 7 809 characters)
- .claude/rules/complem_figures.md (44 / 9 992)
- .claude/rules/export.md (31 / 5 280)
- .claude/rules/fetch_harmonize.md (37 / 3 683)
- .claude/rules/inventory.md (21 / 4 119)
- .claude/rules/nll.md (52 / 6 909)
- .claude/rules/python_patterns.md (189 / 5 370)
- .claude/rules/results_format.md (28 / 4 131)
- .claude/rules/ssst.md (39 / 6 534)
- .claude/rules/ssst_corrections.md (29 / 5 968)
- .claude/rules/temp_picks.md (30 / 5 930)

No other CLAUDE.md / CLAUDE.local.md / AGENTS.md is tracked, and nothing else under `.claude/` is tracked (no skills, agents or commands). Every rules file is under 40 000 characters, and the root CLAUDE.md is under 200 lines.

Checks that came back clean:
- Every backticked `.py` path exists.
- Every function, constant and CLI flag named in the files exists (`station_priority`, `_arbitrate_group_elevation`, `_DIP_SEP_ERZ_FACTOR`, `_SPLIT_YEARS`, `timeGridSpacing`, `windowed_median_map`, `draw_relief`, `--keep-grid`, `--no-dem`, `--extract-from`, crop_quakeml's `--lat-min`…).
- The anchors `fetch_all_bulletins.py:46-49` and `fuse_bulletins.py:979` still point at what the text describes.
- Only three `.py` files were added in the last 35 days: `NLL_run/crop_quakeml.py`, `NLL_run/keep_awake.py` and `complem_figures/seismicity_3d.py`. All three are documented.
- No instruction file hard-codes a Claude model ID.

## Findings

### CLAUDE.md
- **[High]** Line 48: "Stage 1's elevation fill is the only step that reaches the network (~3 min); `--no-fill-elevations` skips it." This is wrong in three ways:
  1. Since e9ab793 (2026-10-29), `_combine_close_stations` calls `_arbitrate_group_elevation` (`fetch_inventory/merge_station_inventories.py:308`), which queries OpenTopoData through `_get_elevation` (`:250`). `merge_inventory` calls `_combine_close_stations` unconditionally (`:538`), and only `_fill_station_elevations` is gated by `fill_elevations` (`:540`). So `--no-fill-elevations` no longer gives an offline run.
  2. Stage 2 also downloads. `fetch_obs/RESIF.py:171` builds an FDSN `Client`, and `fetch_obs/ICGC.py:113,209` calls `requests.get` on `sismocat.icgc.cat`. Both are run by `fetch_all_bulletins.py:55,71-72`.
  3. `--no-fill-elevations` exists only on the `merge_station_inventories.py` CLI (`:573`). The stage script `build_global_inventory.py` has no argparse at all: `main()` calls `run_pipeline()`, and `MergeInventoryParams` is built at `:34-39` without `fill_elevations`. Passing the flag to the stage script has no effect.

  Proposed edit:
  - old: `- ⚠ Stage 1's elevation fill is the only step that reaches the network (~3 min); `--no-fill-elevations` skips it.`
  - new: `- ⚠ Network access: stage 2 downloads RESIF (FDSN) and ICGC (sismocat web) bulletins, and stage 1 queries OpenTopoData (EU-DEM) twice — for the elevation fill (~3 min, skipped by `MergeInventoryParams.fill_elevations=False` / `merge_station_inventories.py --no-fill-elevations`; `build_global_inventory.py` takes no flags) and for the per-group elevation arbitration (`_arbitrate_group_elevation`), which no flag disables.`
- **[Medium]** Line 37: "Every `complem_figures` reader defaults to the expectation; `--solution maxlike` switches." This contradicts `.claude/rules/complem_figures.md:11`, which makes `cross_section.py` an exception (it reads `df['latitude']`, `cross_section.py:276`). There are other undocumented readers without `--solution`: `depth_maps.py` reads `latitude/longitude/depth` from `RESULT/NLL_result.csv` (`:63-66`), and `error_maps.py` and `event_ranking.py` have no `--solution` either. On a CSV written under `SAVE_NLLOC_EXPECTATION`, `latitude` equals `expect_lat`, so these readers plot the expectation. On an older CSV they plot the maximum-likelihood point. The pre-split text said "today's `RESULT/*.csv` are such files"; a813d13 dropped that clause. Whether the current local CSVs predate expectation mode is UNVERIFIABLE here, because `RESULT/` is gitignored.

  Proposed edit:
  - old: `Every `complem_figures` reader defaults to the expectation; `--solution maxlike` switches.`
  - new: `Readers with `--solution` (`event_maps`, `cross_section_panels`, `depth_histogram`, `profile_location_map`, `seismicity_3d`) default to the expectation and switch with `--solution maxlike`; `cross_section.py`, `depth_maps.py`, `error_maps.py` and `event_ranking.py` read `latitude`/`longitude`/`depth`, which is the expectation only on a CSV written in expectation mode.`
- **[Medium]** Line 38: "**0.05 s P / 0.15 s S** everywhere (fetch, `LOCQUAL2ERR`, added picks)". This contradicts `.claude/rules/fetch_harmonize.md:21` ("OMP is the exception: 0.05 s normally / 0.10 s when flagged, for P and S alike"). The code confirms the exception: `fetch_obs/OMP.py:239,275` write `0.05`/`0.10` for both phases. OMP supplies 38 514 of the 57 856 `GLOBAL.obs` events.

  Proposed edit:
  - old: `**0.05 s P / 0.15 s S** everywhere (fetch, `LOCQUAL2ERR`, added picks).`
  - new: `**0.05 s P / 0.15 s S** (fetch, `LOCQUAL2ERR`, added picks) — except OMP picks at fetch, 0.05 s / 0.10 s when flagged, for P and S alike.`

### .claude/rules/inventory.md
- **[High]** Line 19: "⚠ **This is the only step of the pipeline that reaches the network** (~1 s/station, so ~3 min); `--no-fill-elevations` / `MergeInventoryParams.fill_elevations=False` skips it". The evidence is the same as for CLAUDE.md:48. Line 18 of this same file introduces the EU-DEM arbitration, which also calls the network and runs whatever `fill_elevations` says (`merge_station_inventories.py:250,308,538-543`). Stage 2 reaches RESIF and ICGC.

  Proposed edit:
  - old: `⚠ **This is the only step of the pipeline that reaches the network** (~1 s/station, so ~3 min); `--no-fill-elevations` / `MergeInventoryParams.fill_elevations=False` skips it`
  - new: `⚠ **Reaches the network** (~1 s/station, so ~3 min); `--no-fill-elevations` (`merge_station_inventories.py` CLI only — `build_global_inventory.py` takes no flags) / `MergeInventoryParams.fill_elevations=False` skips this fill but **not** `_arbitrate_group_elevation`, which still queries EU-DEM for every group spreading over `ELEVATION_TOLERANCE`. Stage 2 (RESIF FDSN, ICGC web) also downloads`
- **[Medium]** Line 5: the glob `stations/**` matches no tracked file, because `stations/` is in `.gitignore`. It probably fires as intended on the user's local data. Proposed edit: none; keep the glob if matching local untracked files is intended.

### .claude/rules/nll.md
- **[Medium]** Line 19: the anchor `NLL_run/generate_regional_runfiles.py:303-353` has drifted (the file last changed in 7cad9e8, 2026-09-05). Lines 303-353 now hold the `generate_run` docstring and the start of the control-file code. The control file is now written at lines 326-381: `lines = []` at 326, `CONTROL` at 328, `LOCQUAL2ERR` at 377, `LOCPHSTAT` at 378, and the write at 380-381.

  Proposed edit:
  - old: `` `NLL_run/generate_regional_runfiles.py:303-353` ``
  - new: `` `NLL_run/generate_regional_runfiles.py:326-381` (`generate_run`) ``
- **[Medium]** Line 38: the anchor `run_NLL.py:62-69` has drifted. `_ZONES = {` is now at line 65, the six boxes at 66-71, and the closing brace at 72 (`run_NLL.py` last changed in 2e5de46, 2026-10-14).

  Proposed edit:
  - old: `` `run_NLL.py:62-69` ``
  - new: `` `run_NLL.py:65-72` (`_ZONES`) ``
- **[Medium]** Lines 51-52 contradict each other:
  - Line 51 says `obs/NLL_result.obs` has "50 918 events / 711 721 picks".
  - Line 52 (c0ea4e0) says "(same events, **1 059 360 picks** on the 2026-10-23 run, 51 568 events from 727 742 picks …)".
  - Elsewhere the old figures remain: `complem_figures.md:33,37` (50 918 events, 1 020 452 picks), `ssst.md:20` (50 918 → 46 224) and the CLAUDE.md:27 table (50 918).

  Line 52 describes a newer NLL run, but calls its events "same events" while giving a different count. A reader cannot tell which run any of the other figures belong to.

  Proposed edit:
  - old: `(same events, **1 059 360 picks** on the 2026-10-23 run, 51 568 events from 727 742 picks: +331 618, +45.6 %; SANIMS alone +13 466)`
  - new: `(same events as its input. 2026-10-23 run: 51 568 events / 727 742 → **1 059 360 picks**, +331 618, +45.6 %, SANIMS alone +13 466. The 50 918-event figures above and in ssst.md / complem_figures.md are from the earlier NLL run)`
- **[Medium]** Line 5: the glob `run/nll/**` matches no tracked file (`run/` is gitignored). Proposed edit: none if intentional.
- **[Low]** Line 24: "`VGGRID` / `LOCGRID` | … `nz = 761` → 35 km". Only `LOCGRID` uses `nz = 761` (`_compute_grid_corners`, `:285`). `VGGRID` is written as `2 {VGGRID[0]} {VGGRID[1]}` with the default `[9000, 800]` (`:71`, `:335`), so the velocity grid has 800 depth nodes, which reach about 36.95 km. This does not change any decision.

### .claude/rules/results_format.md
- **[Medium]** Line 3: the glob `RESULT/**` matches no tracked file (`RESULT/` is gitignored). Proposed edit: none if intentional.

### .claude/rules/ssst.md
- **[Medium]** Line 9: the glob `run/ssst/**` matches no tracked file (`run/` is gitignored). Proposed edit: none if intentional.
- **[Low]** Line 19:
  - The CLI list omits `--rebuild-grids` (`run_SSST.py:312`, "Rebuild the P/S grids even when they exist").
  - "`--decaff` (allow sleep, §4)" points to §4, but since the split the sleep paragraph lives in CLAUDE.md:52, not in `nll.md`.
  - Line 17 calls SSST the "Final stage", although stage 7 (the export) now follows it.

### .claude/rules/complem_figures.md
- **[Low]** These cross-references dangle after the split:
  - Line 40 says `NLL_run/pdf_metrics.py` is described "(see above)", but that description now lives in `ssst.md:25-39`.
  - Line 42 refers to "the two driver scripts above". `generate_complem_figures.py` and `generate_complem_maps.py` appear only in this file's frontmatter and in CLAUDE.md:55-56.

### .claude/rules/ssst_corrections.md
- **[Low]** Lines 10 and 18 give different numbers for what reads as the same quantity, the iteration-0 → final ("pure SSST") displacement. Line 10 says median 1.09 km, p90 4.74 km. Line 18 says median 1.33 km, p90 6.85 km over 53 754 events. Separately, `complem_figures.md:33` uses "~1.33 km median" for the NLL_result.csv → SSST_result.csv shift, which is a different comparison. The figures may come from different event subsets (line 10's spread map is P-only). This is UNVERIFIABLE without the run outputs.

### Coverage
- **[Low]** `run_gamma_detection.py` (top level, last touched be47d88 2026-04-04) matches no stage rules file; only the generic `**/*.py` of `python_patterns.md` applies. No instruction file mentions it. It is a standalone PhaseNet + GaMMA detection script that reaches FDSN, and it sits outside the 7-stage pipeline. A one-line mention in CLAUDE.md would stop it being mistaken for a pipeline step.

## Summary
11 files audited, 16 findings (2 high, 9 medium, 5 low).

Both High findings describe one hazard that is no longer accurate: the network-access warning. Since e9ab793, elevation arbitration queries OpenTopoData whatever `--no-fill-elevations` says. Stage 2 has always downloaded RESIF and ICGC. And the flag cannot be passed through `build_global_inventory.py`.

Note: skills, agents and settings in ~/.claude/ are local to the user's machine and not covered here; run `/doctor prompt-audit` locally for those.
