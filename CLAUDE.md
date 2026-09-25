# CLAUDE.md — Shallow_Depth_DL_Catalog

## Project Overview

Shallow_Depth_DL_Catalog is an earthquake catalog processing and relocation pipeline focused on the **Pyrenees region**, covering seismicity from **1978 to 2025**.

The goal is to produce a unified, publication-quality earthquake catalog with harmonized magnitudes and improved hypocenter locations, by integrating data from 5 independent seismic networks.

**Area of interest** — three different extents are in play, and they are not the same:
- the RESIF FDSN query box, lat **41–44** / lon −3 to 4 (`fetch_all_bulletins.py:46-49`) — the only geographic filter applied at fetch time, and only to RESIF;
- the AOI filter proper (`global_obs/filter_events_by_aoi.py`), which is **not a rectangle**: two oblique half-planes applied to every source — south of (44.00 N, −0.25 E)→(43.25 N, 3.50 E) and north of (42.50 N, −2.25 E)→(42.00 N, 0.25 E) — plus a per-source line that keeps **RESIF only north** of (43.00 N, −2.25 E)→(42.00 N, 2.25 E) and **IGN/ICGC only south** of (43.75 N, −2.25 E)→(42.00 N, 6.25 E). LDG and OMP get no per-source line;
- the 6 NLL zone boxes (§4), whose union is what actually gets relocated.

`obs/GLOBAL.obs` spans lat 41.30–44.16, lon −2.99 to 4.06.

---

## Pipeline Summary

The workflow follows 7 main stages:

### 1. Station Inventory Fusion
- Source: FDSN XML files + OMP CSV files (in `stations/`)
- Script: `build_global_inventory.py` → calls modules in `fetch_inventory/`
- Output: `stations/GLOBAL_inventory.xml` + `stations/GLOBAL_code_map.txt`
- Each station gets a unique code (`NET.NNNN`); duplicates are removed by distance threshold (20m). Co-located stations inherit the code of the group's **best-ranked** member under `station_priority()` — lowest wins: **permanent network** (`FR RA RD ES CA LC AM`, the `perm` rows of `stations/all_networks.xlsx`, AM/RaspberryShake counted as permanent) before temporary; **RESIF/RENASS** (`FR RA RD`) among the permanent; a **real elevation** before one left at 0; `XX` last; oldest `start_date` breaking ties. `all_networks.xlsx` is *not* read at runtime — `openpyxl` is absent from `seisbench_env`, so the two sets are module constants and the sheet stays the human-maintained reference
- **That one rule is used in all three places a group must yield one representative**: the code it lends the group (`_combine_close_stations`), the lat/lon/elevation written into the NLL `GTSRCE` lines (`generate_regional_runfiles._find_station_info`, §4) and the station name in the QuakeML export (`export_quakeml.load_station_epochs`, §7). Before it, each of the three picked a different member. **They share the ranking, not the outcome**: the export applies it only after the pick's date, so for 37 codes the name it prints is a different member than GTSRCE's donor (§7). `_find_station_info` resolves the code-map block's members **against the inventory**, not by their order in the file, so it stays correct against a `GLOBAL_code_map.txt` written before the rule existed
- Stations left without an elevation inherit one from their co-located group, then from **Open-Elevation** (`_fill_station_elevations`, reusing `_fill_missing_elevations._get_elevation`). 156 of 2 151 carried none — `XI` 99, `24` 31, `X7` 15, `ZU` 8, `XX` 2, `6E` 1, whole nodal deployments whose StationXML never populated the field. ⚠ **This is the only step of the pipeline that reaches the network** (~1 s/station, so ~3 min); `--no-fill-elevations` / `MergeInventoryParams.fill_elevations=False` skips it
- Sources: 4 FDSN inventories (RESIF 346, ICGC 262, ORFEUS 103, GFZ 12) + 12 OMP/temporary deployments (1 691 stations). 2 414 in → 2 151 stations / 1 997 unique codes
- ⚠ **This stage is interactive and cannot run unattended.** `merge_station_inventories.check_inventory()` prompts at the terminal for every station code appearing in more than one network, and the operator's answers are **not logged anywhere** — `GLOBAL_inventory.xml` is not reproducible from code + inputs alone
- `_remove_fdsn_duplicates.py` and `_convert_csv_to_stationxml.py` are **standalone preprocessing**, run by hand per deployment; `build_global_inventory.py` imports neither. They assign network code `XX` to OMP stations (see the `XX`-last term above). The third `_`-prefixed module, `_fill_missing_elevations.py`, is still hand-run on the CSVs, but its `_get_elevation()` **is** now imported by `merge_station_inventories`

### 2. Catalog Fetching & Conversion
- Sources: RESIF (FDSN), ICGC, IGN, LDG, OMP (in `org_catalogs/`)
- Script: `fetch_all_bulletins.py` → calls modules in `fetch_obs/`
- Output: **six** `.obs` files in `obs/` — OMP contributes two, `OMP_78-19.obs` and `OMP_2016.obs` (a re-picked 2016: 2 178 events against 1 568 for the same year in the archive file; both are fused, the overlap resolved by the fusion's loose+pick-validated matching)
- 75 926 events / 1 784 332 picks total. RESIF/ICGC/IGN/LDG cover **2020–2025 only**; the whole 1978–2019 record comes from OMP, so the catalog is not homogeneous in source composition across time
- Only **manual** P/S picks are kept (RESIF `evaluation_mode`, IGN/ICGC GSE2 `m__`). An event with no magnitude of the expected type is dropped entirely, picks included
- Pick uncertainty is **assigned, not measured** — no source publishes it: **0.05 s for P, 0.15 s for S**, reused as `LOCQUAL2ERR` in §4. OMP is the exception: 0.05 s normally / 0.10 s when flagged, **for P and S alike**

### 3. Catalog Harmonization
- Script: `build_global_bulletin.py` → calls modules in `global_obs/`
- Steps:
  1. **remap_picks_to_unified_codes.py** — associates picks with unified station codes
  2. **generate_magnitude_models.py** — builds regression models to convert all magnitude types to ML
  3. **apply_magnitude_models.py** — applies the models to all `.obs` files
  4. **filter_events_by_aoi.py** — filters each source bulletin to the area of interest (run in `pygmt_env`)
  5. **fuse_bulletins.find_and_merge_doubles** — de-duplicates each source catalog individually before fusion (1 s / 50 km)
  6. **fuse_bulletins.py** — spatially/temporally merges all catalogs into `obs/GLOBAL.obs`
  7. **plot_global_catalog_map.py** — generates a map of the merged catalog
- Matching thresholds: strict 15 km / 2 s / 1.5 mag units; loose 50 km / 30 s confirmed by ≥1 shared P-phase pick
- Result: 75 926 → **57 856 events / 1 285 838 picks** in `obs/GLOBAL.obs`. A merged event keeps the **first** location (so RESIF's, whenever it contributes — it is the main bulletin) and the **mean** magnitude; duplicated picks collapse to the **first** occurrence (392 610 removed). `PUBLIC_ID PYRENEES_%06d` is assigned here in chronological order and is the join key of every later stage
- **Only three sources need a conversion model.** `apply_magnitude_models` runs on `obs/*_20-25.obs` only, looking up a model as `"<MagType> <FileAuthor>"`. In `GLOBAL.obs`: 38 514 events `ML OMP` — **already equivalent to ML LDG**, so no model is applied and none is needed — 18 382 `ML LDG`, and 960 `MD LDG` (duration magnitude, kept as `MD`, the one unharmonized group). Note that no OMP→LDG model *could* be fitted in any case: OMP ends in 2019 and LDG starts in 2020, so no event pair exists to regress on
- ⚠ **Step 5 is interactive** — `find_and_merge_doubles` prompts per group, once per source catalog (6 times). Auto-merge only for pairs within 0.15 s / 10 km. Unlike §1 these choices *are* logged (kept/dropped bulletin IDs)
- ⚠ **`build_global_bulletin.py` must not be re-run on already-processed bulletins.** `remap_picks_to_unified_codes` rewrites `obs/*.obs` **in place** and is not idempotent: on a second run it parses the already-substituted code (`FR.0041` → station name `0041`), finds no match, and **deletes the pick**. The glob now matches 22 files, including `GLOBAL.obs`, every `GLOBAL_<N>*.obs`, `NLL_result*.obs` and `SSST_result.obs`
- `_remove_magnitudes_under_1` is **commented out** (`fuse_bulletins.py:982`) although `fuse_bulletins`'s docstring still advertises it; 31.3 % of `GLOBAL.obs` is below ML 1.0
- `global_obs/add_temporary_picks.py` and `list_magnitude_types.py` are imported by nothing — out of the pipeline

### 4. Earthquake Relocation (NonLinLoc)
The study area is too large for a single NLL run, so it is split into **6 geographic zones**, processed with up to **3 zones running concurrently** (zones are independent; NLL runs are the bottleneck).

- **`run_NLL.py`** — for each zone: generates one `.obs` file and one `.in` run file `run/nll/run_<N>_DELAYS.in` (plus GTSRCE station file), runs Vel2Grid → Grid2Time → NLLoc via `NLL_run/run_zone.py`, cleans up `.hdr` files, then generates a second-pass run file `run/nll/run_<N>_NLL.in` by appending per-station delay corrections derived from first-run arrival-time residuals and reruns NLLoc via `run_zone.py` with `--corrections-pass` (grids already built), cleaning up `.hdr` files again right after — this per-zone `.hdr` cleanup keeps at most a few zones' worth on disk at once instead of waiting for all 6 zones to finish. Once all zones are done, exports the locdelay summary via `export_locdelay_info` to `run/nll/locdelays/`. NLL working folders live inside `run/`: `nll_model/` (velocity grids), `nll_time/` (travel-time grids), `nll_loc/` (per-zone NLLoc output).
- **`NLL_run/run_zone.py`** — runs Vel2Grid → Grid2Time → NLLoc in sequence for a given `.in` file; `--corrections-pass` skips Vel2Grid/Grid2Time; accepts a `zone_label` to prefix log output when zones run concurrently

> ⚠️ `run_NLL.py` has **no CLI** — `python run_NLL.py` (with any argument, `--help` included) starts the full pipeline immediately: it rewrites `obs/GLOBAL.obs` in place, regenerates the zone `.obs`/`.in` files, and begins building travel-time grids. Read the file instead of invoking it.

#### The NLL control file — what actually defines the catalog

`NLL_run/generate_regional_runfiles.py:303-353` writes the whole control file as hardcoded literals. These values condition every hypocenter, so they are listed here rather than left buried:

| Statement | Value | Note |
|---|---|---|
| `TRANS` | `LAMBERT WGS-84 <lat_sw> <lon_sw> 42 44 0.0` | standard parallels bracket the range; origin is per-zone |
| `VGGRID` / `LOCGRID` | **0.05 km spacing**, top at **−3 km**, `nz = 761` → **35 km** | −3 km encloses the summit stations; 35 km sits just below the Moho |
| `LAYER` (×5) | 0.0 → 5.50/3.20 · 1 → 5.60/3.26 · 4 → 6.10/3.55 · 11 → 6.40/3.72 · 34 → 8.00/4.50 (Vp/Vs ≈ 1.72, no gradients) | **`Pyrenees_1D`, from Souriau & Pauchet (1998).** 1-D and identical in all 6 zones. The same model is duplicated as Pyrocko `.nd` files in `temp_picks/models/` for the pick-association bands, with nothing checking the two copies stay in sync |
| `LOCSEARCH` | `OCT 50 50 5 0.001 50000 500 1 0` | oct-tree; the 50 000 samples are what populate the `.scat` clouds `pdf_metrics.py` reads |
| `LOCMETH` | `EDT_OT_WT 9999 4 -1 -1 1.72 145 -1.0 0` | Equal Differential Time (robust to outlier picks — right for a 5-agency merge); **minimum 4 phases**, i.e. zero redundancy; Vp/Vs 1.72 because `GTMODE GRID2D` builds P grids only |
| `LOCGAU` / `LOCGAU2` | `0.05 0.0` / `0.01 0.01 2.0` | model error comparable to the pick error |
| `LOCQUAL2ERR` | `0.05 0.15 0.05 0.15 99999.9` | the same 0.05/0.15 s written at fetch time (§2); the alternation is the P/S distinction, not a quality ladder |
| `LOCPHASEID` | `P → P p G PN PG`, `S → S s G SN SG` | collapses `Pn` onto the direct-`P` grid; the 80 km cut removes most distances where that matters. `G` appears in **both** lines, which cannot be right |
| `CONTROL` | `1 54321` | fixed seed → the relocation is deterministic |

Grids are extended **100 km beyond the zone box** on each side, to cover stations out to the 80 km pick limit.

#### Data loss at this stage — not accounted for anywhere else

- `filter_distant_picks` drops picks from stations more than **80 km** from the event: **512 361 of 1 285 838 picks, 39.8 %**. It rewrites `obs/GLOBAL.obs` in place, but is idempotent (a rerun removes 0).
- The 6 zone boxes are at `run_NLL.py:62-69`. Their union does **not** cover the catalog: **5 510 events (9.5 %) fall outside every zone and are silently never relocated** (nothing west of −2.00 or east of 3.50; in the west only lat 42.50–43.50 is covered). A further 8 508 fall in two zones and are relocated twice (dedup by lowest `pdfVolume`); 1 428 more produce no solution. **Total 6 938 events, 12.0 % of `GLOBAL.obs`, absent from `NLL_result.csv`.**
- Second pass: a station delay is applied only with **≥ 100 phases** and `StdDev ≥ 0`.

**Both stages report the location-PDF expectation, not the maximum-likelihood point.** Every NLLoc control file carries `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`, which sets the hypocenter to the expectation and **re-solves origin time, RMS, per-pick residuals, Gap and Dist there** — it re-reads each arrival's travel time from the time grids at the new point, so the reported quality belongs to the reported location. The reason is a consistency break, not a preference: NLLoc's `CovXX..ZZ`, its confidence ellipsoid and this project's `true_erh`/`true_erz` are second moments **about the expectation**, so pairing them with the mode quoted an ellipsoid that was not centred on the point it described. `pdf_metrics.py` had the same break — it whitens samples about their own mean, so `C68`, the metric gating `pyr:usable`, was testing an ellipsoid around a point the catalog did not publish.

Measured over the ML-based `ssst_run1`: median mode→mean offset 0.12 km horizontal / 0.26 km depth, p90 1.72 / 2.72 km, and `dh/true_erh` vs `Psi` gives Spearman ρ = −0.59 — the two estimates separate exactly where the PDF is least Gaussian. The clearest case is the search-grid ceiling at −3 km: **789 events sat pinned on it** (3 116 in the −4..−2 km bin) because an argmax on a truncated grid piles up on the boundary. No expectation is pinned; that bin falls to 816 and events above sea level go from 14.1 % to 8.9 %.

The maximum-likelihood solution is not lost — NLLoc writes it to the `.hyp` as `MAXIMUM_LIKELIHOOD  MaxLikeLat .. Long .. Depth .. OT ..`, verified bit-identical to what a baseline run reports. Costs ~25 % wall clock (measured, 500 events); no file-descriptor growth. The `.scat` clouds are byte-identical to a maximum-likelihood run, so `pdf_metrics.py` is unaffected.

### 5. Post-relocation Processing
- **`run_NLL.py`** (after all 6 zones complete):
  1. Reads the 6 per-zone NLL CSV summaries, deduplicates zone-overlap events (kept: lowest `pdfVolume`), writes → `RESULT/NLL_result.csv`
  2. Rematches relocated events back to `obs/GLOBAL.obs` via `publicId` to recover metadata not present in NLL output (e.g. magnitude)
  3. Saves matched events to `obs/NLL_result.obs` (50 918 events / 711 721 picks)
- **`add_temp_picks.py`** (**required**, run after — `run_SSST.py` reads its output as its only input): augments `obs/NLL_result.obs` with picks from external sources → `obs/NLL_result_augmented.obs` (same 50 918 events, **1 020 452 picks: +308 731, +43.4 %**)

### 6. SSST Relocation
Final stage, run after `add_temp_picks.py`: relocates `obs/NLL_result_augmented.obs` with iterative Source-Specific Station Terms (NonLinLoc `Loc2ssst`). Zones run **strictly sequentially** (memory-bound), everything inside a zone in parallel (`NLLOC_CORES` NLLoc chunks / `LOC2SSST_CORES` Loc2ssst instances).

- **`run_SSST.py`** — orchestrator; campaign configuration (RUN_NAME, CHAR_DISTS, VPVS, core counts, LSPHSTAT) at the top of the file; CLI `--zones`, `--iteration-start`, `--iteration-stop` (partial campaigns/resume). Per zone: cuts `obs/GLOBAL_<N>_SSST.obs` + `stations/GTSRCE_SSST_<N>.txt`, derives `run/ssst/run_<N>_NLL.in` (NLLoc) and `run/ssst/run_<N>_SSST.in` (Loc2ssst) from `run/nll/run_<N>_DELAYS.in` (`NLL_run/generate_ssst_runfiles.py`), splits the bulletin into per-event files in `obs/nlloc_obs/GLOBAL_<N>/` (`NLL_run/reformate_obs.py`), builds P+S grids in `run/ssst_model|ssst_time` (VpVs −9.99, real S grids), and runs the iteration loop (`NLL_run/run_ssst.py`: len(CHAR_DISTS) SSST iterations + final NLLoc-only relocation, outputs under `run/ssst_loc/<RUN_NAME>/`).
- After all zones: merges the final-iteration CSVs → `RESULT/SSST_result.csv` (dedup by lowest `pdfVolume`), rematches against `obs/NLL_result_augmented.obs` via `publicId` → `obs/SSST_result.obs` (same modules as the NLL stage). 50 918 → **46 224 events**; the 4 694 lost are mostly the `LSPHSTAT NRdgsMin = 6` phase minimum, against NLL's 4.
- `LSPHSTAT [0.15, 6, 200.0, 0.3, 0.5, 10.0]` (RMS ≤ 0.15 s, ≥ 6 phases, gap ≤ 200°, P/S residual ≤ 0.3/0.5 s, Len3 ≤ 10 km) selects which events *teach* the corrections — deliberately stricter than what *receives* them. `NRdgsMin` doubles as the NLLoc min-phases threshold, read back from the generated control file so the two cannot diverge.
- The SSST control files are **derived from `run/nll/run_<N>_DELAYS.in`** (the first-pass file, so the NLL static `LOCDELAY`s are deliberately not inherited — iteration 0 at L = 9999 recomputes that same quantity). `LSGRID`/`LSOUTGRID` are re-gridded to a coarse **1 km** spacing: Loc2ssst holds two full buffers in RAM per instance.
- **`NLL_run/pdf_metrics.py`** — last step of `run_SSST.py` (after the merge and the rematch): reads the per-event `.scat` clouds of each zone's final SSST iteration and rewrites `RESULT/SSST_result.csv` in place with the location-PDF quality columns, then saves the diagnostic figures below. Non-fatal there (the results are already written) and idempotent, so it can be rerun standalone on any campaign; ~65 s for ~46 k events (~50 s metrics, ~15 s figures).
  ```bash
  python NLL_run/pdf_metrics.py --run-name ssst_run1
  python NLL_run/pdf_metrics.py --run-name ssst_run1 --figures false   # metrics only
  ```
  Uses its own `.scat` reader — **not** `obspy.io.nlloc.util.read_nlloc_scatter`, which does `np.fromfile(...)[4:]` and silently drops the first 3 samples of every file (4 records dropped where the header is 1).

  Figures land in `complem_figures/pdf_metrics/` (plotting is wrapped: a failure there never touches the already-written CSV):
  - `<run>_metrics_vs_quality.pdf` — Ψ / `C68` / `dip_stat` against depth, RMS, Nphs, Gap, Dist (hexbin density + median/IQR over equal-count bins + the per-bin Gaussian null)
  - `<run>_error_vs_quality.pdf` — the inverse view: the published error `true_erh` / `true_erz` (log axis) against RMS, Gap, Nphs, nearest-station distance, Ψ, `C68`, and a box split on `dip_reject`. Each panel reports Spearman ρ **and** ρ given Nphs, because the marginal ρ of RMS against ERH is +0.00 while the Nphs-controlled one is +0.62 — RMS rises with the phase count and ERH falls with it, and the two cancel. Ψ is the strongest single correlate of the published error (ρ = −0.87 ERH, −0.85 ERZ); Ψ and `C68` are computed on whitened samples, so neither is definitionally tied to the size of the ellipsoid
  - `<run>_gridmap.pdf` — the same three metrics as lon/lat windowed-median maps
  - `<run>_gridmap_depth.pdf` — those maps split into depth quartiles, colour scale shared per row
  - `<run>_voxels_<metric>.html` — interactive 3-D lat/lon/depth cells (Plotly)

  The maps colour by `C68_z = (C68 − 0.6827) / C68_sigma_n` rather than raw `C68`, and the 1-D nulls are per bin: both are read from each event's own `J_null_p95` / `C68_sigma_n`, so nothing breaks when `n_scat` is not near-constant as it happens to be in `ssst_run1`.

### 7. Final Bulletin Export
Last stage, run after `run_SSST.py`. Merges `obs/SSST_result.obs` (magnitudes + picks) and `RESULT/SSST_result.csv` (full-precision hypocentres + PDF metrics) on `publicId` into `RESULT/FINAL.xml`, a QuakeML 1.2 bulletin.

- **`final_steps.py`** — orchestrator; the post-SSST finalization stage. Currently one step, further steps get appended as `[n/N]`.
- **`NLL_run/export_quakeml.py`** — the export. All 46 224 events, all 1 001 095 picks. ~5 min, ~3.4 GB peak RAM, ~0.86 GB output. Serializes in chunks of 5 000 events and splices the bodies into one file (obspy has no streaming QuakeML writer).
  ```bash
  python final_steps.py
  python NLL_run/export_quakeml.py --obs obs/SSST_result.obs --csv RESULT/SSST_result.csv \
      --inventory stations/GLOBAL_inventory.xml --output RESULT/FINAL.xml
  ```
- **Two origins per event.** The preferred origin is the PDF expectation (`pyr:locationEstimator = "expectation"`), carrying the time, the `OriginQuality`, the `OriginUncertainty` + ellipsoid and the arrivals — all of which NLLoc evaluated at that point. The maximum-likelihood hypocentre rides along as a second origin (`smi:pyrenees/origin/<pid>/maxlike`, `pyr:locationEstimator = "maximum_likelihood"`) with lat/lon/depth and a reconstructed time, and **deliberately no quality, no uncertainty and no arrivals** — those belong to the expectation, and repeating them under the mode would restore the very mismatch this arrangement removes (§4). Emitted only when the CSV carries `maxlike_*`.
- Each event carries `pyr:usable` = `(C68 − 0.6827)/C68_sigma_n ≥ −2` **and** `not dip_reject` **and** metrics present → 45 415 usable / 809 unusable, plus a `pyr:rejectReason` string. **Ψ is exported but never rejects** (directionless; 81 % of the catalog exceeds its own `J` null). The `C68` cut is against the simulated null, not raw 0.6827 — that difference is ~160 events flagged instead of thousands.
- Station codes are resolved to real `NET.STA` via inventory `alternate_code`, epoch-aware; when several candidates cover the same date the winner is `station_priority()` from §1 — the merge's own rule, so the permanent network wins and `XX` loses last (currently decides only `FR.0013` → `FR.MTHF`). But the ranking applies **after** the date, not instead of it — 79 of 80 multi-station codes have disjoint epochs, so the pick's date usually decides, and should. For **37 codes** the exported name is therefore a different member of the group than the one GTSRCE took coordinates from: same site to within 20 m, and to within 61 m in elevation now that groups share elevations. Everything that re-derives the representative *after* the merge passes the group's code to `station_priority()`, which adds a last term preferring the member whose network matches the code's prefix. That is needed because the in-group elevation fill erases criterion 3 — once a member has inherited the winner's elevation the term goes neutral and the date could hand the group elsewhere — and the prefix is the durable record of who actually won. The merge itself must not pass it: it assigns the codes, so the term would be circular. Resolution is never destructive: every pick keeps `pyr:unifiedCode`, and codes covering several stations also carry `pyr:alternateStations` (the non-chosen ones, **including date-mismatched entries**) as waveform-retrieval fallbacks — ~15% of picks.
- The `pyr:` prefix is bound to `http://shallow-depth-dl-catalog/quakeml/1.0`. QuakeML 1.2's schema is closed, so every non-standard field lives in that namespace; validators skip foreign namespaces, which is why the file still validates.
- `read_events('RESULT/FINAL.xml')` works (46 224 events in 255 s, 0 warnings) but peaks at **~15 GB RAM** — ~2M objects, of which picks+arrivals are 83% (12.5 KB per pick, 55.8 KB per event). So the stage writes the catalog three ways: `FINAL.xml` (everything, ~15 GB), `FINAL_catalog.xml` (all events, picks stripped, measured 2.5 GB / 33 s — the FDSN `includearrivals=false` convention), and `FINAL_<from>_<to>.xml` (5-year calendar periods set by `_SPLIT_YEARS`, each independently valid; together an exact partition of `FINAL.xml`). Parts are very uneven — `FINAL_2020_2024.xml` holds 40% of all picks and measures 6.0 GB, vs 0.05 GB for 1975-1979. Disable either companion with `--no-catalog` / `--no-parts`.
- Prefer `lxml.etree.iterparse` for anything that does not need the whole catalog resident; it walks the full file in seconds at flat memory.
- `Mag 0.00` (5 010 events) is an OMP placeholder, written through unchanged because magnitudes are recomputed later — do not fit it as data.

---

## Complementary Analysis

Scripts in `complem_figures/` for visualization and statistics:
- `event_maps.py` — geographic maps of seismicity
- `gutenberg_richter.py` — magnitude-frequency distribution
- `depth_maps.py` — depth distribution
- `depth_histogram.py` — histogram of event depths
- `error_maps.py` — location uncertainty maps
- `cross_section.py` — vertical cross-sections; `--format 6` reads `RESULT/SSST_result.csv` directly (`true_erh`/`true_erz` as the error filter), colours the section by `--use-err erh|erv|psi|c68z`, and `--usable` applies the same acceptance test as `pyr:usable` in the QuakeML export
- `station_map.py` — map of seismic stations
- `zone_map.py` — overview map of the 6 NLL zones
- `event_ranking.py` — ranks events by pdfVolume/ellipsoidVolume change (NLL → SSST); flags multi-zone/zone-changed events; carries the PDF-quality columns (post-SSST only) when `pdf_metrics.py` has annotated the CSV, and runs without them otherwise; optional gridmap PDF
- `plot_pdf_cloud.py` — interactive 3D PDF scatter-cloud of one event across SSST iterations (Plotly). Each ellipsoid is centred on the PDF expectation, which is also the published hypocentre; the diamond markers and the path joining them are the **maximum-likelihood** point, read from the `.hyp` `MAXIMUM_LIKELIHOOD` line. The gap between a diamond and its ellipsoid centre is that iteration's mode-vs-mean distance. Its NLL reference panel needs the `maxlike_*` columns in `RESULT/NLL_result.csv`
- `ssst_evolution.py` — per-zone pdfVolume/EllipsoidLen3/RMS evolution across SSST iterations (convergence QC)
- `station_colocation.py` — the one diagnostic aimed at **stage 1**: is 20 m the right co-location radius? Reads the 2 414 pre-merge stations (`stations/*/*.xml`) against the groups in `GLOBAL_inventory.xml`. Separations are bimodal with a real trough between 20 m and 200 m, so the threshold does sit in a gap, but it is not free: 651 pairs within 20 m, 811 within 30 m, the extras inside the dense `AM`/`8M`/`XI` nodal arrays. "Merged" is read from the **outcome**, not from the distance — which is what exposes the **423 pairs that merged from beyond the radius** through the single-linkage chaining of `_combine_close_stations`, and the 38 inside it that did not because `check_inventory` removed one by hand. Three panels (nearest-neighbour split merged/not, threshold sensitivity, per-network), `seisbench_env`, ~20 s, → `complem_figures/station_colocation/`
- `temporary_network_impact.py` — what the external picks of §5's `add_temp_picks.py` bought, from the one place the question is cleanly answerable: `obs/NLL_result.obs` and `obs/NLL_result_augmented.obs` hold the **same 50 918 events with byte-identical hypocentres** — their event lines differ in `PhaseCount` alone — so anything that moves between them is the added picks, with no relocation mixed in. That is why it does **not** compare `NLL_result.csv` against `SSST_result.csv`, which would confound the two (different event sets, and hypocentres moved ~1.33 km median). Two figures in `complem_figures/temporary_network_impact/`, `seisbench_env`, ~5 s:
  - `phase_count_timeline.pdf` — two chronological strips 1978→2025 coloured by the median phases per event per month, before above and after below on a **shared** colour scale (that sharing is the figure). Empty months (8 of 565) stay grey rather than closing the gap. The default `vmax` is the p99 of the pooled monthly medians because exactly one bin in 1 114 reaches 52 and would otherwise eat half the ramp; the colorbar declares the clip with `extend='max'`. Median 12 → 15 phases overall, but the gain is **entirely post-2010** — only 585 of the 17 077 pre-2010 events gain anything
  - `nearest_station_distance.pdf` — overlapping before/after histograms of the epicentral distance to the closest recording station, all events beside the 26 150 that gained picks. Median 10.80 → 9.89 km, and events with a station within 5 km go 14.7 % → 19.3 % (14.1 % → 23.0 % on the gained-picks panel). Distances beyond the axis limit fold into the last bin, but every printed statistic is computed on the **unclipped** values

  Station positions come from `stations/GTSRCE_*.txt`, i.e. the coordinates NLLoc itself located with, so the result is checkable against the binary: the computed median reproduces NLLoc's own `Dist` column to **0.02 km**. One pick in 1 020 452 sits on a code absent from those files (`FR.0002`) and is reported, never silently dropped.
- `ssst_corrections.py` — reconstructs and maps the SSST travel-time corrections themselves. The `.ssst` grids Loc2ssst wrote are deleted by `run_SSST.py`'s own cleanup, so the field is recomputed from the surviving per-event `.hyp` locations with Loc2ssst's own formula (`corr = Σ res·w / Σ w`, `w = exp(−d²/L²) + floor`, **d measured event-to-grid-node, not to the station**). `res = obs − pred` and the correction is *added* to the predicted time, so red = arrivals later than the 1-D model predicts. Each iteration is an increment; the final grids carry their sum. Outputs to `complem_figures/ssst_corrections/`:
  - `<run>_station_atlas.pdf` — one page per station×phase field (300 with ≥100 usable arrivals): five increments (9999→50→15→5→1 km) + cumulative total, depth-slice maps on a shared diverging scale, unsupported nodes greyed. A station in two zones is drawn in its better-sampled one — the zones ran separate Loc2ssst passes, so their fields are independent.
  - `<run>_spread_map.pdf` — catalog-wide map of the **across-station** spread of the total correction per event (P only). A correction common to all of an event's stations is absorbed by the origin time and cannot move the hypocentre, so the dispersion — not the mean — is the part that *can* relocate. **Not an impact map**: against the pure SSST displacement (iteration-0 vs final location, same picks) ρ = +0.10, or +0.20 given Nphs; median displacement rises only 0.95 → 1.27 km across the whole spread range. What moves a hypocentre is the *azimuthal pattern* of the differential correction, which a std discards. Nphs suppresses the marginal ρ (spread ρ = +0.42 with Nphs, displacement ρ = −0.18) — the same confound as RMS vs ERH in `pdf_metrics.py`. For reference the pure SSST displacement has median 1.09 km, p90 4.74 km.

  **Nothing is interpolated between map nodes** — one flat `pcolormesh` cell per independently evaluated node; the smoothness on the page is the Gaussian kernel itself. So node spacing is a real constraint: the 0.02° default is ~2.2 km lat / ~1.6 km lon at 43°N, which **undersamples the L = 1 km panel**. The 300-page atlas keeps it coarse on purpose (213 s); use `--stations FR.0041:P --map-spacing 0.005` (~55 s/page, writes `*_station_atlas_selection.pdf`) for pages that go on a slide.

  `--extent=lon0,lon1,lat0,lat1` (the `=` form is required when lon0 is negative) draws every page on one shared frame so pages can be compared, clipped to each zone's `LSGRID` — read from `run/ssst/run_<N>_SSST.in`, since off that box Loc2ssst computed nothing. Default frames each page on its own station's events, which is tighter but not comparable page to page.

  **`--depth` defaults to each station's own median event depth**, not a fixed value. The kernel is 3-D, so a slice `dz` off the event mass damps every event by `exp(−dz²/L²)` — negligible at L = 15 km, fatal at L = 1 km. A fixed 10 km slice against a catalog median of 6.6–7.3 km cut `FR.0041`'s median summed weight in the finest panel from 5.7 to 2.0. Pass `--depth` only to put several stations on a common plane.

  - `<run>_displacement_map.pdf` — **the impact map**: each event's iteration-0 hypocentre vs its final one (same picks, same parameters, only the corrected grids differ; matched by `publicId`), so the difference is the corrections alone. Median distance moved + arrows of the median horizontal displacement per 0.25° cell, and a signed median depth change (red = deeper). 53 754 events, median move 1.33 km, p90 6.85 km — ~1 km through the dense central network but **3–4 km at the western and eastern margins**.

  Both catalog-wide maps share `windowed_median_map`: a **median** over a circular window (longitude scaled by `cos(lat)`), because displacement is heavy-tailed (max 99.5 km) and a mean lets single events paint neighbourhoods.

  ```bash
  python complem_figures/ssst_corrections.py                 # all three, ~7 min
  python complem_figures/ssst_corrections.py --product displacement
  python complem_figures/ssst_corrections.py --extract-only  # fill the cache only
  ```
  Parsing the ~276 k per-event `.hyp` files takes ~100 s, cached as `.npz` per (zone, iteration) in `run/ssst_corrections_cache/` (33 MB) — **delete it after a new campaign**, the cache keys on (zone, iteration) and not on the run that produced them. The parser prefix-matches `HYPOCENTER` / `GEOGRAPHIC` and ignores everything else, so the `MAXIMUM_LIKELIHOOD` line added under `SAVE_NLLOC_EXPECTATION` passes straight through; the source positions it reads become the expectations, which is right — Loc2ssst reads the same line and used the same points. **Validated against the binary**: Loc2ssst re-run on zone 6 (`FR.0047`, P and S) agrees node-by-node over all 303×313×40 nodes to 0.0000 ms at L=9999 and 0.0736 ms at L=1, against ±0.36 s amplitudes; its own "163 accepted" matches the module's LSPHSTAT selection exactly. Omitting `LOCFILES` from the control file is what keeps that check cheap — `ihave_time_input_grids = flag_out_grid * flag_nlloc_outfile` (`Loc2ssst.c:590`), so Loc2ssst writes the correction grid and skips the travel-time grids, needing no Grid2Time rebuild.

  Note when reading the figures: increments are **not** monotonically decreasing — they peak at L = 15 km, the scale where the residual field has real spatial structure. And LSPHSTAT admits only a subset of the catalog: 23 727 events at iteration 0, rising to 34 674 by the final relocation.

`complem_figures/pdf_metrics/` holds figures too, but they are written by `NLL_run/pdf_metrics.py` (see above), not by a module living here.

> **Environments**: `seisbench_env` is the project default — the whole pipeline (`build_global_inventory.py` → `run_SSST.py`) runs in it unprefixed. `pygmt_env` is the only exception, for the modules importing PyGMT or `xarray`.
> - `seisbench_env` → `generate_complem_figures.py` (Gutenberg-Richter, depth maps, error maps)
> - `pygmt_env`     → `generate_complem_maps.py` (event maps for each zone and final catalog), `cross_section.py`, and the two `build_global_bulletin.py` steps it launches itself via `conda run` (`filter_events_by_aoi.py`, `plot_global_catalog_map.py`)

`event_ranking.py`, `plot_pdf_cloud.py`, and `ssst_evolution.py` are standalone diagnostics, not wired into the two driver scripts above; they read `RESULT/*.csv` and `run/nll_loc/` / `run/ssst_loc/<run-name>/` directly. `station_colocation.py` is standalone too, but reads `stations/` rather than any relocation output, and `temporary_network_impact.py` likewise, reading only `obs/NLL_result*.obs` and `stations/GTSRCE_*.txt`.

`zone_Arette/` — focused analysis of the Arette seismic zone.

---

## External Pick Ingestion (temp_picks/)

A self-contained sub-pipeline for ingesting picks from external sources into `obs/NLL_result.obs`, producing `obs/NLL_result_augmented.obs`. All scripts live in `temp_picks/` and are importable as a package (`from temp_picks.<module> import <function>`). Log files are written to `temp_picks/console_output/`.

The root-level script **`add_temp_picks.py`** orchestrates the full pipeline in sequence (6 steps). Only steps 1–2 (theoretical tables, QC figure) are skipped when their output exists; the merge, convert and match steps always run.

Eight external datasets are ingested, **in a fixed order that encodes precedence** — `TEMP_STB` after `TEMP_OMP` so `match_picks`'s per-event `(station, phase)` dedup silently drops the overlap. Each source is matched against the *previous* output, so the stage is **not re-runnable** against an existing `NLL_result_augmented.obs`. Measured over the eight: **308 731 picks added**, 31.9 M skipped for no event within 60 s (these are continuous PhaseNet detection streams, not bulletins), 534 312 as duplicates, 201 200 outside the travel-time band, 1 661 ambiguous, **0** for an unknown station.

Added picks get the project's standard **0.05 s P / 0.15 s S** — i.e. machine picks are given the same declared precision as manual readings, and the available `phase_score` is not used to modulate it.

The travel-time bands come from the **same velocity model as NonLinLoc**, as Pyrocko `.nd` files, perturbed **±5 %** and enveloped over source depths 0 and 30 km. The two copies of the model are maintained separately and nothing checks that they still agree.

| Script | Role |
|--------|------|
| `build_theoretical_tables.py` | Runs Pyrocko's `cake` CLI to compute P/S travel-time envelopes (±5% velocity, 0–100 km) → `temp_picks/tables_Pyr.csv` |
| `merge_omp_picks.py` | Merges all yearly OMP/PhaseNet CSV files from `picks_OMP/` subdirectories → `pick_files/merged_omp.csv`; station `SMC` and year `2026` excluded by default |
| `merge_pyrenees_picks.py` | Concatenates RaspberryShake/PhaseNet `.txt` files from `picks_station_pyrenees/` and `picks_station_pyrenees2/` → `pick_files/merged_pyrenees.txt` and `pick_files/merged_pyrenees2.txt` |
| `convert_picks.py` | Converts external pick files to `.obs` pick line format; maps station names to internal codes via `GLOBAL_code_map.txt`. Formats `TEMP_OBS`, `TEMP_RSB`, `TEMP_OMP`, `TEMP_OTH`, and `TEMP_STB` are supported; new formats are registered in `FORMAT_HANDLERS`. `TEMP_STB` (Strasbourg/RENASS-OMP parquet picks) reads directories of parquet files instead of a single text file and drops picks below a `phase_score` threshold. Unresolved stations are logged as an end-of-run summary. |
| `match_picks.py` | Matches converted picks to bulletin events: 60 s time window + residual filter (±0.1 s P, ±0.3 s S, plus ±2.5 s t0-error margin); appends new picks and updates `PhaseCount`; chains against `obs/NLL_result.obs` → `obs/NLL_result_augmented.obs`; auto-sorts output via `sort_picks`. |
| `sort_picks.py` | Sorts pick lines within each event block by ascending arrival time. |
| `plot_travel_times.py` | QC figure: scatter of observed (distance, travel time) picks over theoretical P/S bands. |

---

## Key Data Formats

### `.obs` (custom seismic bulletin)
- One block per event, separated by blank lines
- Event line starts with `# `: location, magnitude, quality parameters (azimuth gap, RMS, horizontal/vertical uncertainty)
- Following lines: one pick per station (station code, phase P/S, arrival time, uncertainties)

### NLL output
- Per-zone CSV summary: `run/nll_loc/GLOBAL_<N>/GLOBAL_<N>.obs.sum.grid0.loc.csv` — relocated hypocenter parameters including `publicId` (links back to the input `.obs` event), `pdfVolume` (location PDF volume; smaller = tighter), the confidence-ellipsoid axes, etc.
- Merged result: `RESULT/NLL_result.csv` — deduplicated across all 6 zones; horizontal/vertical uncertainties `true_erh` / `true_erz` are derived from the 3-D confidence ellipsoid, rescaled to DOF-appropriate one-sigma (p = 0.6827) confidence factors (not axis-aligned `errH`/`errZ` projections, and not NLL's raw 3-DOF ellipsoid scaling): `true_erz` is a 1-DOF marginal standard deviation, `true_erh` is a 2-DOF horizontal error ellipse reduced via the geometric mean of its semi-axes
- Does **not** contain magnitude or full pick metadata → rematching to `obs/GLOBAL.obs` via `publicId` is necessary
- `ellipsoidVolume` = 4/3·π·Len1·Len2·Len3, written alongside `true_erh`/`true_erz` for direct comparison against `pdfVolume`
- `latitude` / `longitude` / `depth` are the **PDF expectation** and equal `expect_lat` / `expect_lon` / `expect_z` exactly — NLLoc runs with `SAVE_NLLOC_EXPECTATION`, so `RMS`, `Gap`, `Dist` and the origin time in `date-time` are its own re-evaluation at that point (see §4). This is what puts the published location and the published error on the same point
- `maxlike_latitude` / `maxlike_longitude` / `maxlike_depth` / `maxlike_ot_sec` — the maximum-likelihood hypocenter, added by `merge_regional_results.py` from the `MAXIMUM_LIKELIHOOD` line of the summary `.hyp` sitting beside each zone CSV (the CSV itself never carries it). `maxlike_ot_sec` is **seconds within the minute only** — NLLoc discards the ML hour/minute before writing — so the full timestamp has to be rebuilt against the expectation origin and snapped to the nearest minute; the two solutions do straddle minute boundaries. The merge **raises** if that line is absent, since its absence means the run was not made in expectation mode

### Location-PDF quality columns (`RESULT/SSST_result.csv` only)
Added by `NLL_run/pdf_metrics.py` from the `.scat` scatter clouds.
- `J` / `Psi` — negentropy `KL(p ‖ N(μ,Σ))` in nats and `Psi = exp(-J)`, the effective-volume ratio. `Psi = 1` is exactly Gaussian. **Directionless — never a rejection criterion on its own.**
- `C68` — posterior mass inside the nominal one-sigma ellipsoid (`chi2.ppf(0.6827, 3) = 3.5267`, not `chi2.ppf(0.68, 3)`). The only directional metric, so it drives keep/reject: `> 0.6827` conservative (safe), `< 0.6827` over-confident.
- `dip_stat` / `dip_pval` / `dip_sep_km` / `dip_reject` — Hartigan's dip test on the depth marginal, computed at a common subsample of 400 so p-values are comparable across events. `dip_sep_km` is the width of Hartigan's modal interval, i.e. how far apart the competing depths are. **`dip_reject` needs both conditions**: `dip_pval < 0.05` **and** `dip_sep_km > true_erz`. The dip statistic is a vertical distance on the ECDF and is invariant under rescaling of the depth axis, so it registers that a second mode exists but never how far away — on its own it discards events whose bimodality is real but sits entirely inside the quoted error, where `depth ± true_erz` already covers both modes. Requiring the separation term moves the count 5 329 → 661. The threshold is `_DIP_SEP_ERZ_FACTOR` in `pdf_metrics.py`; at 2.0 it would be vacuous (the `dip_sep_km`/`true_erz` ratio tops out at 2.24 over this catalog).
- `J_null_p95` / `C68_sigma_n` — per-event n-matched Gaussian nulls. **`J` and `C68` are uninterpretable without them**: the k-NN estimator bias depends on n, and the `C68` null spread is not binomial (μ/Σ are fitted on the samples being tested).
- `n_scat` — sample count backing every metric above.
- These measure **consistency, not accuracy** — velocity-model bias is invisible to all of them.

---

## Git Workflow

- All AI-generated commits go to the **`claude` branch**, never to `main`
- Commit messages must be clean and descriptive so changes are understandable without reading the diff
- Use the format `type: description` — e.g. `fix: ...`, `feat: ...`, `docs: ...`. Never use scoped form `fix(module): ...`
- **Never push automatically to main** — commit and push to the **`claude` branch** without asking, but notice the user
- The user reviews changes locally and decides when to merge or push to `main`

### Pull Request Format
- PR body: **Summary section only** (bullet points of what changed and why). No "Test plan" section.
