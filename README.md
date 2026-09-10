# Shallow_Depth_DL_Catalog

A pipeline for building a unified, publication-quality earthquake catalog for the **Pyrenees region**, integrating data from five independent seismic networks over the period **1978–2025**. The merged catalog spans lat 41.30–44.16, lon −2.99 to 4.06; the area of interest is an oblique polygon rather than a box, and differs by source — see [Area of interest](#area-of-interest).

The workflow covers catalog fetching, station inventory fusion, magnitude harmonization, event merging, two-stage probabilistic earthquake relocation with NonLinLoc (per-station delay corrections, then iterative SSST), and result visualization.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Pipeline](#pipeline)
  - [1. Station Inventory Fusion](#1-station-inventory-fusion)
  - [2. Catalog Fetching & Conversion](#2-catalog-fetching--conversion)
  - [3. Catalog Harmonization](#3-catalog-harmonization)
  - [4. Earthquake Relocation (NonLinLoc)](#4-earthquake-relocation-nonlinloc)
  - [5. Post-relocation Processing](#5-post-relocation-processing)
  - [6. External Pick Ingestion (temp_picks)](#6-external-pick-ingestion-temp_picks)
  - [7. SSST Relocation](#7-ssst-relocation)
  - [8. Final Bulletin Export](#8-final-bulletin-export)
- [Complementary Analysis](#complementary-analysis)
- [Dependencies](#dependencies)
- [A note on AI assistance](#a-note-on-ai-assistance)

---

## Overview

Five seismic catalogs are integrated:

| Source | Network | Format | Period |
|--------|---------|--------|--------|
| RESIF | French national network | FDSN / QuakeML | 2020–2025 |
| ICGC | Catalan network | Web / text | 2020–2025 |
| IGN | Spanish national network | Text | 2020–2025 |
| LDG | French seismological bulletin | Text | 2020–2025 |
| OMP | Pyrenean Observatory | Text / .mag | 1978–2019 |

All catalogs are converted to a common `.obs` format, magnitudes are harmonized towards a common **ML (LDG)** scale, and events are merged into a single `GLOBAL.obs` bulletin.

> **Only three of the five sources need a conversion model.** `MLv` (RESIF), `mb_Lg` (IGN) and ICGC's `ML` are regressed onto ML LDG; **the OMP ML is already equivalent to the LDG ML** for the pre-2020 events it covers, so no model is applied to it and none is needed. In `GLOBAL.obs`: 38 514 events `ML OMP`, 18 382 `ML LDG`, and 960 `MD LDG` — the last being LDG's own duration magnitude, kept as `MD` and identifiable by its `MagType`. Note that no OMP→LDG regression could have been fitted in any case, since OMP ends in 2019 and LDG is fetched only from 2020, leaving no overlapping event pair.

### Area of interest

Three different extents are in play:

- **The fetch box**, lat 41–44 / lon −3 to 4, applied as an FDSN query parameter to **RESIF only** (`fetch_all_bulletins.py:46-49`). The other four sources are regional by construction.
- **The AOI filter** (`global_obs/filter_events_by_aoi.py`), which is **not a rectangle**. Two oblique half-planes apply to every source — events must lie south of the line (44.00 N, −0.25 E) → (43.25 N, 3.50 E) and north of (42.50 N, −2.25 E) → (42.00 N, 0.25 E), so the bounds follow the trend of the range rather than parallels. A third, **per-source** line implements a division of labour: **RESIF is kept only north** of (43.00 N, −2.25 E) → (42.00 N, 2.25 E), **IGN and ICGC only south** of (43.75 N, −2.25 E) → (42.00 N, 6.25 E), with a deliberate overlap band where all three contribute and the fusion arbitrates. LDG and OMP get no per-source line.
- **The 6 NLL zone boxes**, whose union is what actually gets relocated — and which do not cover the whole AOI (see [§4](#4-earthquake-relocation-nonlinloc)). Earthquakes are then relocated using **NonLinLoc** across 6 geographic sub-zones, and final results are compiled into `RESULT/NLL_result.csv` and `obs/NLL_result.obs`. After external picks are ingested (`obs/NLL_result_augmented.obs`), a final **SSST** relocation stage (`run_SSST.py`, iterative NLLoc + Loc2ssst) produces `RESULT/SSST_result.csv` and `obs/SSST_result.obs`.

---

## Project Structure

```
Shallow_Depth_DL_Catalog/
│
├── fetch_all_bulletins.py        # Entry point: fetch & convert all catalogs
├── build_global_inventory.py     # Entry point: fuse all station inventories
├── build_global_bulletin.py      # Entry point: harmonize & merge catalogs
├── run_NLL.py                # Entry point: run NLL relocation (6 zones) and finalize the catalog
├── add_temp_picks.py             # Entry point: augment NLL_result.obs with external picks
├── run_SSST.py               # Entry point: SSST relocation of the augmented catalog
├── final_steps.py                # Entry point: export the final QuakeML bulletin
├── generate_complem_figures.py   # Entry point: matplotlib figures
├── generate_complem_maps.py      # Entry point: PyGMT event maps
│
├── fetch_obs/                # Catalog fetching & .obs conversion modules
│   ├── RESIF.py
│   ├── ICGC.py
│   ├── IGN.py
│   ├── LDG.py
│   └── OMP.py
│
├── fetch_inventory/          # Station inventory fusion modules
│   ├── merge_station_inventories.py
│   ├── _remove_fdsn_duplicates.py
│   ├── _fill_missing_elevations.py
│   └── _convert_csv_to_stationxml.py
│
├── global_obs/               # Catalog harmonization modules
│   ├── remap_picks_to_unified_codes.py
│   ├── list_magnitude_types.py       (standalone; not in the pipeline)
│   ├── generate_magnitude_models.py
│   ├── apply_magnitude_models.py
│   ├── add_temporary_picks.py        (superseded by temp_picks/; not in the pipeline)
│   ├── filter_events_by_aoi.py
│   ├── fuse_bulletins.py
│   └── plot_global_catalog_map.py
│
├── NLL_run/                  # NonLinLoc workflow modules
│   ├── run_zone.py                  # Run Vel2Grid → Grid2Time → NLLoc for one zone
│   ├── generate_regional_runfiles.py
│   ├── append_station_delays.py
│   ├── export_locdelay_info.py
│   ├── parse_nll_output.py          (deprecated)
│   ├── filter_distant_picks.py
│   ├── match_pre_post_relocation.py
│   ├── merge_regional_results.py
│   ├── generate_ssst_runfiles.py    # SSST: derive control files from the NLL run files
│   ├── reformate_obs.py             # SSST: split a bulletin into per-event .nlloc_obs
│   ├── run_ssst.py                  # SSST: iterative NLLoc + Loc2ssst workflow (one zone)
│   ├── pdf_metrics.py               # SSST: location-PDF quality metrics from .scat clouds + their figures
│   └── export_quakeml.py            # Final: merge .obs + result CSV into a QuakeML bulletin
│
├── complem_figures/          # Visualization & statistical analysis
│   ├── event_maps.py
│   ├── depth_maps.py
│   ├── error_maps.py
│   ├── depth_histogram.py
│   ├── gutenberg_richter.py
│   ├── cross_section.py
│   ├── station_map.py
│   ├── zone_map.py
│   ├── event_ranking.py
│   ├── plot_pdf_cloud.py
│   ├── ssst_evolution.py
│   ├── ssst_corrections.py
│   └── station_colocation.py
│
├── zone_Arette/              # Focused analysis of the Arette seismic zone
│
├── temp_picks/               # External pick ingestion & QC sub-pipeline
│   ├── build_theoretical_tables.py  # Compute P/S travel-time bands (Pyrocko/cake)
│   ├── merge_omp_picks.py           # Merge yearly OMP/PhaseNet CSV files
│   ├── merge_pyrenees_picks.py      # Merge RaspberryShake/PhaseNet text files
│   ├── convert_picks.py             # Convert external pick files to .obs format
│   ├── match_picks.py               # Match converted picks to bulletin events
│   ├── sort_picks.py                # Sort picks by arrival time within each event
│   ├── plot_travel_times.py         # Plot theoretical bands vs observed picks
│   ├── models/                      # Velocity model files (.nd)
│   ├── pick_files/                  # Input pick files (raw) + merged outputs
│   ├── tables_Pyr.csv               # Computed travel-time table
│   ├── figures/                     # Output figures
│   └── console_output/              # Log files
│
├── org_catalogs/              # Raw input catalogs (not modified)
├── obs/                      # .obs bulletin files (source + merged)
├── stations/                 # Station inventories (XML + unified)
├── run/                      # NLL working directory
│   ├── nll/                  # NLL run configuration files (.in) + locdelays/
│   ├── nll_model/            # Velocity model grids (NLL)
│   ├── nll_time/             # Travel time grids (NLL)
│   ├── nll_loc/              # NLL output files per zone
│   ├── ssst/                 # SSST control files (.in), per-zone tmp/, log/ (run journals)
│   ├── ssst_model/           # Velocity model grids (SSST, P+S)
│   ├── ssst_time/            # Travel time grids (SSST, P+S)
│   └── ssst_loc/             # SSST output per campaign (iterations + final locations)
├── RESULT/                   # Merged relocation results (.csv)
└── mag_model/                # Saved magnitude conversion models
```
(`obs/nlloc_obs/GLOBAL_<N>/` holds the per-event `.nlloc_obs` files of the SSST stage.)

---

## Conventions

All Python scripts in this project share the same interface contract:
- **CLI**: every script accepts `--help` and can be run directly from the command line
- **Public API**: every module is importable as a Python package (e.g. `from temp_picks.match_picks import match_picks`)
- **Logging**: timestamped log files are written to a `console_output/` directory local to each sub-pipeline

---

## Pipeline

Entry points, in execution order:

```
build_global_inventory.py → fetch_all_bulletins.py → build_global_bulletin.py
→ run_NLL.py → add_temp_picks.py → run_SSST.py → final_steps.py
```

### 1. Station Inventory Fusion

**Script:** `build_global_inventory.py`  
**Module:** `fetch_inventory/merge_station_inventories.py` → `merge_inventory()`

Merges all station XML inventories into a single unified inventory: 4 FDSN sources (RESIF 346 stations, ICGC 262, ORFEUS 103, GFZ 12) and 12 OMP/temporary deployments (1 691 stations). Each station receives a unique code of the form `NET.NNNN`; stations within 20 m of each other share a code, inheriting that of the group's **best-ranked** member. 2 414 stations in → 2 151 retained across 46 networks, grouped into 1 997 unique codes.

The ranking is `station_priority()` in `merge_station_inventories.py`, and it is the **single rule** used everywhere a group of co-located stations has to yield one representative — the code it lends the group, the coordinates and elevation written into the NonLinLoc `GTSRCE` lines ([§4](#4-earthquake-relocation-nonlinloc)), and the station name resolved in the QuakeML export ([§8](#8-final-bulletin-export)). Lowest wins: **permanent network** (`FR RA RD ES CA LC AM`, the `perm` rows of `stations/all_networks.xlsx`) before temporary; **RESIF/RENASS** (`FR RA RD`) among the permanent ones; a station carrying a **real elevation** before one left at 0; the `XX` placeholder last; oldest `start_date` breaking ties. Before this rule the three stages each picked a different member of the same group.

One wrinkle the rule has to survive: the in-group elevation fill below **erases its own third criterion**. Once a member has inherited the winner's elevation, that term goes neutral, and a later reader re-deriving the ranking could hand the group to a different station than the merge did. The record of who won is not lost, though — it is the code's own prefix, `AM.0043` being an `AM` code because an `AM` station won it. So every consumer that re-derives the representative *after* the merge passes the code to `station_priority()`, which adds a last term preferring the member whose network matches that prefix. The merge itself must not pass it: it is what assigns the codes, so the term would be circular there. Without it exactly one group of 80 ends up labelled after one station and described by another.

Stations left without an elevation inherit one from their co-located group, then from the Open-Elevation API — 156 of 2 151 carried none, whole nodal deployments (`XI` 99, `24` 31, `X7` 15) whose StationXML never populated the field, and a station's elevation enters the travel-time computation directly. This is the **only step of the pipeline that reaches the network**; `--no-fill-elevations` skips it.

OMP CSV data is pre-processed with `_remove_fdsn_duplicates.py` and `_convert_csv_to_stationxml.py`. These two are **standalone tools run by hand** per deployment — `build_global_inventory.py` imports neither. They assign the placeholder network code `XX` to OMP stations, which is what the `XX`-last term of the ranking resolves. The third, `_fill_missing_elevations.py`, is still run by hand on the CSVs, but its `_get_elevation()` is now imported by `merge_station_inventories` for the fill above.

> ⚠️ **This stage is interactive and cannot run unattended.** `merge_inventory` calls `check_inventory()` unconditionally, which prompts at the terminal for every station code appearing in more than one network and asks which network(s) to remove. There is no flag and no non-interactive path. The operator's answers are **not logged, saved, or committed**, so `stations/GLOBAL_inventory.xml` is not reproducible from the code and the inputs alone — re-running the stage means re-making every one of those judgements. It is the only step of the pipeline that is not reproducible; the relocation itself is deterministic (fixed NonLinLoc seed).

**Outputs:**
- `stations/GLOBAL_inventory.xml` — unified QuakeML inventory
- `stations/GLOBAL_code_map.txt` — mapping between original and unified station codes

---

### 2. Catalog Fetching & Conversion

**Script:** `fetch_all_bulletins.py`  
**Modules:** `fetch_obs/` (one module per source)

Downloads or reads each catalog and converts it to the `.obs` format. RESIF and ICGC are fetched dynamically; IGN, LDG, and OMP are read from local files in `org_catalogs/`.

**Outputs:** **six** `.obs` files in `obs/` — five sources, but OMP contributes two:

| Source | File | Events | Picks |
|---|---|---|---|
| OMP | `OMP_78-19.obs` | 37 859 | 673 716 |
| RESIF | `RESIF_20-25.obs` | 11 406 | 227 432 |
| ICGC | `ICGC_20-25.obs` | 9 511 | 289 203 |
| LDG | `LDG_20-25.obs` | 7 968 | 280 638 |
| IGN | `IGN_20-25.obs` | 7 004 | 159 333 |
| OMP | `OMP_2016.obs` | 2 178 | 48 126 |
| | **Total** | **75 926** | **1 784 332** |

`OMP_2016.obs` is a **re-picked version of 2016** — 2 178 events against 1 568 for the same year inside `OMP_78-19.obs`. Both files are fused; the fusion resolves 1 523 of the overlapping pairs by loose matching confirmed on shared P picks.

Three selection decisions are made here and are not revisited later:

- **Only manual P/S picks are kept** — by `evaluation_mode` for RESIF, by the GSE2 `m__` flag for IGN and ICGC. LDG and OMP are manual by construction.
- **An event with no magnitude of the expected type is dropped entirely**, picks included.
- **RESIF, ICGC, IGN and LDG cover 2020–2025 only**; the entire 1978–2019 record comes from OMP. The catalog is therefore not homogeneous in source composition across time, which any completeness or rate analysis has to account for.

**Pick uncertainty is assigned, not measured.** No source bulletin publishes a per-pick timing error, so one is supplied: **0.05 s for P, 0.15 s for S**. The same two values reappear as `LOCQUAL2ERR` in the NonLinLoc control file and are reused for externally-ingested picks, so the error budget is consistent end to end — with one exception: **OMP assigns 0.05 s (or 0.10 s when flagged) to S picks as well as P**, using its own analyst quality digits instead. Since OMP supplies 40 % of the catalog's picks and effectively all of them before 2020, pre-2020 S arrivals are declared three times more precise than everyone else's.

#### .obs format

Each event occupies one block separated by a blank line:

```
# YYYY MM DD HH MM SS.ss  Lat  Lon  Dep  Mag  MagType  Author  Nph  ErrH  ErrV  Gap  RMS
STA.CODE  INS  CMP  ONSET  PHASE  DIR  YYYYMMDD  HHMM  S.MS  Err  ErrMag  Coda  Amp  Period  # Phase  Chan  Origin  PGV
...
```

---

### 3. Catalog Harmonization

**Script:** `build_global_bulletin.py`  
**Module:** `global_obs/`

Runs the following steps in sequence:

| Step | Module | Function | Description |
|------|--------|----------|-------------|
| 1 | `remap_picks_to_unified_codes.py` | `remap_picks_to_unified_codes()` | Associates picks with unified station codes from global inventory |
| 2 | `generate_magnitude_models.py` | `convert_magnitudes()` | Builds piecewise ODR regression models (breakpoint at M=2): MLv RESIF, mb_Lg IGN, ML ICGC → ML LDG |
| 3 | `apply_magnitude_models.py` | `apply_magnitude_models()` | Applies the models to convert all `.obs` magnitudes to ML LDG |
| 4 | `filter_events_by_aoi.py` | `filter_events_by_aoi()` | Removes events outside the area of interest |
| 5 | `fuse_bulletins.py` | `find_and_merge_doubles()` | Deduplicates each source catalog individually before fusion |
| 6 | `fuse_bulletins.py` | `fuse_bulletins()` | Matches and merges all cleaned catalogs into `GLOBAL.obs` |
| 7 | `plot_global_catalog_map.py` | `plot_global_catalog_map()` | Generates a map of the merged catalog |

**Matching thresholds (step 6 fusion):** strict ≤15 km / ≤2 s / ≤1.5 mag units (ML–ML pairs); loose ≤50 km / ≤30 s, confirmed by ≥1 shared P-phase pick (same station, Δt ≤ 1 s).

**How a merged event is built.** When events from several sources match, latitude, longitude and depth are taken from the **first** source that has them — which is RESIF whenever it contributes, since it is the main bulletin — while the magnitude becomes the **mean** across contributing sources. The picks of all sources are pooled, then de-duplicated by (unified station code, phase), keeping the **first** occurrence; 392 610 pick lines are removed this way, so the rule decides the arrival time of a large fraction of the catalog. `PUBLIC_ID PYRENEES_%06d` is assigned here, in chronological order, and is the join key of every later stage — note it is a *position in a list*, not a content hash, so re-running the fusion with one event inserted renumbers everything after it.

> ⚠️ **`build_global_bulletin.py` must not be re-run on already-processed bulletins.** Step 1 rewrites `obs/*.obs` **in place** and is not idempotent: on a second run it parses the already-substituted code (`FR.0041` → station name `0041`), finds no inventory match, and **deletes the pick**. The glob now matches 22 files, including `GLOBAL.obs`, every `GLOBAL_<N>*.obs`, `NLL_result.obs`, `NLL_result_augmented.obs` and `SSST_result.obs`, so a rerun would strip the picks from every later stage's output too.

> ⚠️ **Step 5 is interactive.** `find_and_merge_doubles` presents every group of near-identical events for a decision, once per source catalog (six times). Only pairs within 0.15 s and 10 km are merged automatically. Unlike the inventory stage, these choices *are* logged with the kept and dropped bulletin IDs.

Note also that `_remove_magnitudes_under_1` is **commented out** (`fuse_bulletins.py:982`) even though `fuse_bulletins`'s own docstring still advertises the ML 1.0 cut; 31.3 % of `GLOBAL.obs` lies below ML 1.0.

**Outputs:**
- `obs/GLOBAL.obs` — unified catalog: **57 856 events / 1 285 838 picks** (18 070 cross-source duplicates merged)
- `obs/MAPS/` — statistics figures
- `mag_model/` — serialized magnitude models. The fitted relations, and how well they fit:

| Model | M ≥ 2 | R² | M < 2 | R² | pairs |
|---|---|---|---|---|---|
| `MLv RESIF` → ML LDG | `y = 1.118x + 0.089` | 0.624 | `y = 0.762x + 0.800` | 0.511 | 3 200 |
| `mb_Lg IGN` → ML LDG | `y = 1.147x − 0.197` | 0.762 | `y = 0.802x + 0.494` | **0.343** | 3 101 |
| `ML ICGC` → ML LDG | `y = 0.810x + 0.982` | 0.765 | `y = 0.689x + 1.225` | 0.732 | 2 068 |

Events are paired for the fit at a deliberately strict 10 km / 2 s — much tighter than the fusion's own thresholds, since a regression is corrupted by a single bad pair far more than a catalog merge is. Two properties are worth knowing: every low-magnitude branch fits worse than its high-magnitude counterpart, so the magnitudes of the small events the project targets are the least well constrained; and the two branches use **different estimators** — orthogonal distance regression above M = 2, but ordinary least squares under a continuity constraint below it, because SciPy's ODR has no constrained form. The estimator that ignores error in x is therefore applied exactly where that error is largest. The models are applied as point estimates with no uncertainty propagated.

---

### 4. Earthquake Relocation (NonLinLoc)

The study area is too large for a single NLL run, so it is divided into **6 geographic zones**. Each zone is processed independently — up to **3 zones run concurrently** — then results are merged.

**Script:** `run_NLL.py`

#### The NonLinLoc control file

`NLL_run/generate_regional_runfiles.py:303-353` writes the entire control file as hardcoded literals. These values condition every hypocentre in the catalog, so they are recorded here. The velocity model is **`Pyrenees_1D`, after Souriau & Pauchet (1998)**:

| Statement | Value | |
|---|---|---|
| `TRANS` | `LAMBERT WGS-84 <lat_sw> <lon_sw> 42 44 0.0` | Lambert conformal conic, standard parallels bracketing the range; origin per zone |
| `VGGRID` / `LOCGRID` | **0.05 km spacing**, top at **−3 km**, `nz = 761` → bottom at **35 km** | −3 km encloses the summit stations; 35 km sits just below the Moho of the model |
| `LAYER` ×5 | 0.0 → 5.50/3.20 · 1 → 5.60/3.26 · 4 → 6.10/3.55 · 11 → 6.40/3.72 · 34 → 8.00/4.50 km/s, Vp/Vs ≈ 1.72, no gradients | **`Pyrenees_1D`, from Souriau & Pauchet (1998).** 1-D and identical in all six zones. The same model is duplicated as Pyrocko `.nd` files under `temp_picks/models/` for the pick-association bands, and nothing checks that the two copies stay in sync |
| `LOCSEARCH` | `OCT 50 50 5 0.001 50000 500 1 0` | oct-tree; the 50 000 samples populate the `.scat` clouds the PDF metrics are computed from |
| `LOCMETH` | `EDT_OT_WT 9999 4 -1 -1 1.72 145 -1.0 0` | Equal Differential Time — robust to outlier picks, which suits a five-agency merge. **Minimum 4 phases**, i.e. a formally determined solution with no redundancy; such events are located and judged later by the PDF metrics rather than rejected here |
| `LOCGAU` / `LOCGAU2` | `0.05 0.0` / `0.01 0.01 2.0` | model error comparable in size to the pick error |
| `LOCQUAL2ERR` | `0.05 0.15 0.05 0.15 99999.9` | the same 0.05/0.15 s written at fetch time; the alternation is the P/S distinction re-expressed, not a quality ladder |
| `LOCPHASEID` | `P → P p G PN PG` · `S → S s G SN SG` | collapses `Pn` onto the direct-P grid — acceptable only because the 80 km cut removes most distances where `Pn` arrives first. `G` appears in **both** lines, which cannot be correct |
| `CONTROL` | `1 54321` | fixed seed, so the relocation is reproducible bit for bit |

Grids are extended **100 km beyond the zone box** on each side, so that stations out to the 80 km pick limit fall inside them.

#### What this stage removes

Two large reductions happen here, neither of them previously recorded:

- **Picks from stations more than 80 km from the event are discarded** — **512 361 of 1 285 838, or 39.8 %**. The cut keeps the inversion inside the range where the 1-D model's direct Pg/Sg arrivals are the first arrivals, and sharpens depth resolution, since distant arrivals constrain depth barely at all while still pulling the epicentre. `filter_distant_picks` rewrites `obs/GLOBAL.obs` in place but is idempotent (a second run removes 0).
- **The six zone boxes do not cover the catalog.** Their union admits 52 346 of 57 856 events, so **5 510 events — 9.5 % — fall outside every zone and are silently never relocated**: nothing west of −2.00° or east of 3.50°, and in the west only latitudes 42.50–43.50. A further 8 508 events fall inside two zones and are located twice (the duplicate resolved by lowest `pdfVolume`), and 1 428 more enter a zone but yield no solution. **In total 6 938 events, 12.0 % of `GLOBAL.obs`, are absent from `RESULT/NLL_result.csv`.**

| Zone | Latitude | Longitude | Events |
|---|---|---|---|
| 1 | 42.50 – 43.50 | −2.00 – −0.75 | 11 822 |
| 2 | 42.50 – 43.25 | −1.00 – 0.50 | 23 333 |
| 3 | 42.00 – 43.25 | 0.25 – 1.00 | 5 976 |
| 4 | 42.00 – 43.00 | 0.75 – 2.25 | 14 079 |
| 5 | 42.00 – 43.00 | 2.00 – 3.50 | 4 704 |
| 6 | 42.75 – 43.75 | 2.25 – 3.50 | 940 |

Neighbouring zones overlap by 0.25° of longitude on purpose: an event near a boundary is located twice, in two grids with two different station sets, and the tighter of the two PDFs wins.

For each zone:

1. **First pass.** Calls `NLL_run/generate_regional_runfiles.py` → `generate_run()`:
   - Generates `obs/GLOBAL_<N>.obs` (regional subset, with far picks removed — done once globally before any zone starts)
   - Generates `stations/GTSRCE_<N>.txt` (station list)
   - Generates `run/nll/run_<N>_DELAYS.in` (NLL configuration file)

   NLL is then launched automatically via `NLL_run/run_zone.py`, which runs Vel2Grid → Grid2Time → NLLoc in sequence. Before each run it checks free disk space on the output filesystem, and after each step it scans the subprocess output for fatal disk/memory errors (which NLL programs can otherwise report as a silent exit code 0), aborting the pipeline if one is detected. Console log lines are prefixed with the zone label (e.g. `[zone 3]`) since multiple zones may log concurrently.
   - Cleans up `.hdr` files left by NLL in the `run/nll_loc/GLOBAL_<N>/` folder.

2. **Second (corrections) pass.** Calls `NLL_run/append_station_delays.py` → `append_station_delays()` to read per-station average residuals (LOCDELAY entries) from the first pass and append qualifying station delay corrections to a second-pass run file `run/nll/run_<N>_NLL.in`. NLLoc is then launched automatically via `NLL_run/run_zone.py` with `--corrections-pass` (Vel2Grid and Grid2Time are skipped since the grids are already built).
   - Cleans up `.hdr` files left by NLL in the `run/nll_loc/GLOBAL_<N>/` folder again, right after this zone's second pass, and deletes the zone's travel-time grids (`run/nll_time/Pyrenees_<N>/`) — nothing reads them after the second pass and they dominate disk usage. This per-zone cleanup keeps at most a few zones' worth on disk at a time instead of waiting for all 6 zones to finish.

`generate_run()` and `append_station_delays()` reconfigure a shared logger on every call, which isn't safe if two zones call them at the same instant, so `run_NLL.py` serializes just those two (fast, non-NLL) calls behind a lock; the slow NLL subprocess runs themselves are not serialized and run fully in parallel across zones.

#### Diagnostic — `NLL_run/export_locdelay_info.py`

Called automatically by `run_NLL.py` once all zones have completed both passes. Reads the LOCDELAY station corrections from all second-pass run files and exports them to `run/nll/locdelays/locdelay_summary.txt`, keeping only entries with |residual| > 0.3 s. Useful for identifying stations with systematically biased travel-time residuals.

---

### 5. Post-relocation Processing

**Script:** `run_NLL.py` (runs after all 6 zones complete)
**Modules:** `NLL_run/merge_regional_results.py`, `NLL_run/match_pre_post_relocation.py`

1. Reads the 6 per-zone NLL CSV summaries (`run/nll_loc/GLOBAL_<N>/GLOBAL_<N>.obs.sum.grid0.loc.csv`), deduplicates events that appear in multiple overlapping zones (kept: lowest `pdfVolume`), and writes → `RESULT/NLL_result.csv`. True horizontal/vertical errors (`true_erh` / `true_erz`) are derived from the 3-D confidence ellipsoid, rescaled to DOF-appropriate one-sigma (p = 0.6827) confidence factors: `true_erz` is a 1-DOF marginal standard deviation, `true_erh` is a 2-DOF horizontal error ellipse reduced to the geometric mean of its semi-axes. `ellipsoidVolume` (4/3·π·Len1·Len2·Len3) is also written, so the Gaussian ellipsoid volume can be compared directly against `pdfVolume`, NLLoc's OCT-tree-integrated PDF volume — the two diverge for non-Gaussian events.
2. Rematches relocated events back to `obs/GLOBAL.obs` via the `publicId` field to recover metadata absent from NLL output (magnitude, pick details, etc.).
3. Saves matched events → `obs/NLL_result.obs`

#### Which hypocentre the catalog reports

A NonLinLoc location is a probability density, not a point, and NLLoc offers two ways to collapse it: the **maximum-likelihood** point (the mode of the PDF) and the **expectation** (its mean). This catalog reports the expectation, at both the NLL and the SSST stage, via `LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`.

The reason is a consistency requirement rather than a preference. The uncertainty the catalog publishes — NLLoc's covariance, its confidence ellipsoid, and the `true_erh` / `true_erz` derived from them above — is a **second moment about the expectation**. Quoting it next to the mode therefore attached an error ellipsoid to a point that ellipsoid was not centred on. The same break ran through the quality control: `pdf_metrics.py` measures `C68` on samples whitened about their own mean, so the one metric that decides whether an event is usable was testing an ellipsoid around a location the catalog did not publish. It now tests the published one.

`SAVE_NLLOC_EXPECTATION` is not a relabelling of coordinates. NLLoc re-solves the event at the expectation: it re-reads every arrival's travel time from the time grids at the new point and recomputes the **origin time, the RMS, the per-pick residuals**, then redoes the azimuthal gap and the station distances. That is what makes the fix possible at all — the origin time is analytically marginalised out of the spatial PDF, so no post-processing of the output files could have produced a matching one.

How much it matters, measured over the previous maximum-likelihood catalog (46 224 events):

| | median | p90 | max |
|---|---|---|---|
| horizontal mode→mean offset | 0.12 km | 1.72 km | 88.9 km |
| depth offset | 0.26 km | 2.72 km | 31.0 km |

Small for most events, and that is the point: where the PDF is close to Gaussian the two estimates agree and the choice is immaterial. They separate where it is not — `dh / true_erh` against `Ψ` gives Spearman ρ = **−0.59**, so the mode drifts furthest from the mean, relative to the quoted error, precisely on the events whose PDF is least Gaussian and whose mode therefore means least.

The sharpest case is the top of the search grid, at −3 km. **789 events sat pinned exactly on that boundary**, and 3 116 fell in the −4..−2 km bin: an argmax over a truncated grid piles up against the edge, because the single highest cell can be the last one before the grid stops. No expectation is pinned there — the mean integrates over the whole density instead of picking its peak — and the bin drops to 816, with events placed above sea level falling from 14.1 % to 8.9 %. For a catalog built to study shallow seismicity, that spike was an artefact of the estimator, not a result.

Nothing is discarded. The maximum-likelihood solution is written to each `.hyp` as a `MAXIMUM_LIKELIHOOD` line — verified bit-identical to what a maximum-likelihood run reports — and recovered into the merged CSV as `maxlike_latitude` / `maxlike_longitude` / `maxlike_depth` / `maxlike_ot_sec`, then published as a second origin in `FINAL.xml`. The scatter clouds are byte-identical either way, so every PDF metric is unaffected. The cost is ~25 % of NLLoc wall clock.

---

### 6. External Pick Ingestion (temp_picks)

**Script:** `add_temp_picks.py`
**Modules:** `temp_picks/`

Sits between the NLL and SSST relocation stages: a self-contained sub-pipeline for ingesting picks from external sources into `obs/NLL_result.obs`, producing `obs/NLL_result_augmented.obs` (the input of `run_SSST.py`). **This stage is required, not optional** — `run_SSST.py` has no other input. The root-level script **`add_temp_picks.py`** orchestrates steps 1–6 automatically; only steps 1 and 2 are skipped when their output already exists.

It runs *between* the two relocations because the association test needs an event location to predict a travel time against: the NLL pass produces locations good enough to associate new picks, and the SSST pass then exploits them.

**Eight external datasets** are ingested in a fixed order, and the order encodes a precedence — `TEMP_STB` sources come after `TEMP_OMP` so that the per-event `(station, phase)` de-duplication silently drops whatever they duplicate. Each source is matched against the *previous* output, which is what makes that work, and which also means the stage **cannot be re-run** against an existing `NLL_result_augmented.obs`.

Measured across the eight sources:

| Outcome | Picks |
|---|---|
| **Added** | **308 731** (711 721 → 1 020 452, +43.4 %) |
| Skipped — no bulletin event within 60 s | 31 881 596 |
| Skipped — duplicate `(station, phase)` on the event | 534 312 |
| Skipped — travel time outside the theoretical band | 201 200 |
| Skipped — ambiguous, several events matched | 1 661 |
| Skipped — station not in inventory | 0 |

The 99 % rejection rate is expected and is the point: these are **continuous PhaseNet detection streams**, not curated bulletins, so the stage is an association problem rather than a merge. Ambiguous picks are **discarded rather than guessed** — a wrong association would inject a false arrival into the SSST relocation and corrupt both the hypocentre and the station correction derived from it, whereas a discarded pick costs only information.

Added picks receive the project's standard **0.05 s P / 0.15 s S**, i.e. machine picks are given the same declared precision as an analyst's manual reading. The `phase_score` is available and is *not* used to modulate that.

| Step | Script | Description |
|------|--------|-------------|
| 1 | `build_theoretical_tables.py` | Uses Pyrocko's `cake` CLI to compute P/S travel-time envelopes across ±5% velocity models and source depths of 0 and 30 km, for epicentral distances 0–100 km → `tables_Pyr.csv`. The model is the **same one NonLinLoc uses**, as `.nd` files; expressing the tolerance as a velocity envelope rather than a fixed time window makes it grow with distance, as travel-time model error actually does |
| 2 | `plot_travel_times.py` | QC figure: overlays all observed (distance, travel time) picks from a bulletin on top of the theoretical P/S bands. Skipped if the figure already exists. Also usable as a standalone script. |
| 3 | `merge_omp_picks.py` | Merges all yearly OMP/PhaseNet CSV files from `picks_OMP/` subdirectories → `pick_files/merged_omp.csv`. Station `SMC` and year `2026` are excluded by default; configurable via `--drop-years`. Rows with a PhaseNet `phase_score` below 0.5 (default) are dropped. |
| 4 | `merge_pyrenees_picks.py` | Concatenates RaspberryShake/PhaseNet `.txt` files from `picks_station_pyrenees/` and `picks_station_pyrenees2/` → `pick_files/merged_pyrenees.txt` and `pick_files/merged_pyrenees2.txt`. Lines with a `prob=` value below 0.5 (default) are dropped. |
| 5 | `convert_picks.py` | Converts external pick files to the project's `.obs` pick line format; maps short station names to internal codes via `GLOBAL_code_map.txt`. Supports formats `TEMP_OBS`, `TEMP_RSB`, `TEMP_OMP`, `TEMP_OTH`, and `TEMP_STB`; new formats are added as handler functions. `TEMP_STB` (Strasbourg/RENASS-OMP picks, `temp_picks/all_picks/PICKS_MARC/`) reads directories of parquet files instead of a single text file, and drops picks below a `phase_score` threshold (default 0.5). Unresolved stations are reported as an end-of-run summary. |
| 6 | `match_picks.py` | For each converted pick, finds candidate events within a 60 s origin-time window, filters by theoretical travel-time residual (±0.1 s P, ±0.3 s S, plus ±2.5 s t0-error margin), and appends matched picks to the bulletin. Chains against `obs/NLL_result.obs` → `obs/NLL_result_augmented.obs`. Runs `sort_picks` automatically on the output. |
| — | `sort_picks.py` | Sorts all pick lines within each event block by ascending arrival time. Invoked automatically by step 6; also usable as a standalone script on any bulletin. |

---

### 7. SSST Relocation

**Script:** `run_SSST.py`
**Modules:** `NLL_run/generate_ssst_runfiles.py`, `NLL_run/reformate_obs.py`, `NLL_run/run_ssst.py`

Final relocation stage, run **after** `run_NLL.py` and `add_temp_picks.py`. It relocates the augmented catalog (`obs/NLL_result_augmented.obs`) with **Source-Specific Station Terms** (NonLinLoc's `Loc2ssst`): instead of the single static `LOCDELAY` correction per station of the NLL stage, each station/phase gets a 3-D correction grid, smoothed with a characteristic distance that shrinks over iterations (`CHAR_DISTS`, default `[9999, 50, 15, 5, 1]` km — the first, huge value is equivalent to static terms). Each iteration relocates the catalog with the previous iteration's corrected travel-time grids, then recomputes the corrections from the new residuals; a final NLLoc-only pass relocates the catalog with the last grids (saving oct-tree and fmamp output).

Zones are processed **strictly sequentially** (Loc2ssst holds two full correction-grid buffers in RAM per instance), but everything inside a zone runs in parallel: the catalog is split round-robin across `NLLOC_CORES` concurrent NLLoc processes, and the station list across `LOC2SSST_CORES` concurrent Loc2ssst instances (a per-station split is exact — each station's correction depends only on its own residuals).

Per zone:

1. **`NLL_run/generate_ssst_runfiles.py`** — cuts the zone bulletin `obs/GLOBAL_<N>_SSST.obs` from the augmented catalog, writes `stations/GTSRCE_SSST_<N>.txt` (only the stations picked in the zone), and derives both control files from the NLL-stage `run/nll/run_<N>_DELAYS.in` (TRANS, grids, and localization parameters are reused verbatim, so the zone information entered in `run_NLL.py` stays the single source of truth): `run/ssst/run_<N>_NLL.in` (NLLoc) and `run/ssst/run_<N>_SSST.in` (Loc2ssst, with LSGRID/LSOUTGRID derived from the LOCGRID extent at a coarse 1 km spacing).
2. **`NLL_run/reformate_obs.py`** — splits the zone bulletin into one `.nlloc_obs` file per event in `obs/nlloc_obs/GLOBAL_<N>/` (named and headed by `publicId`, which NLLoc propagates into its output — the rematch key).
3. **`NLL_run/run_ssst.py`** — builds the initial **P and S** travel-time grids (`VpVs = -9.99`: real S grids, so S station terms are gridded independently; skipped when already present), then runs the iteration loop. Fails fast on any non-zero NLLoc/Loc2ssst exit, on fatal disk/memory errors in the subprocess logs, and on empty obs/grid globs. Supports partial campaigns and resuming (`--iteration-start` / `--iteration-stop`).
4. **Grid cleanup** — the zone's travel-time grids (the big disk consumers) are deleted automatically once the zone is done: every `ssst_corr<i>/` grid set except the last, which is kept with its symlinks materialized (the resume input of a partial campaign, the reusable final SSST model of a finished one). The initial grids (`run/ssst_time/Pyrenees_<N>/`) go too as soon as iteration 0 has run.

Once all zones are complete, the final-iteration CSVs feed the same chain as the NLL stage: `merge_regional_results.merge_bulletins` → `RESULT/SSST_result.csv` (zone-overlap duplicates resolved by lowest `pdfVolume`), then `match_pre_post_relocation.save_bulletin` rematching against `obs/NLL_result_augmented.obs` → `obs/SSST_result.obs`. **50 918 → 46 224 events**; the 4 694 lost are mostly the `LSPHSTAT` phase minimum of 6, against NLL's 4.

`LSPHSTAT = [0.15, 6, 200.0, 0.3, 0.5, 10.0]` — RMS ≤ 0.15 s, ≥ 6 phases, gap ≤ 200°, P/S residual ≤ 0.3/0.5 s, ellipsoid `Len3` ≤ 10 km — selects which events *teach* the corrections, and is deliberately stricter than what *receives* them: a residual from a badly-located event measures the location error, not the velocity-model error. Its `NRdgsMin` doubles as the NLLoc min-phases threshold, read back out of the generated Loc2ssst control file so the two selections cannot drift apart. The selection widens as locations improve — 23 727 events admitted at iteration 0, 34 674 by the final relocation.

The SSST control files are derived from `run/nll/run_<N>_DELAYS.in`, the **first-pass** NLL file, so the NLL stage's static `LOCDELAY` corrections are deliberately not inherited: iteration 0 at L = 9999 km recomputes that same quantity from scratch, and keeping both would double-count. `LSGRID`/`LSOUTGRID` are re-gridded from the `LOCGRID` extent to a coarse **1 km** spacing, since Loc2ssst holds two full buffers in RAM per instance — and since the correction field is intrinsically smooth at the kernel scale, sampling it finer than the smallest `L` would add nothing.

Every SSST iteration reports the expectation, not only the final relocation — see [Which hypocentre the catalog reports](#which-hypocentre-the-catalog-reports). Loc2ssst builds its corrections from the residuals in the `.hyp` files, so this keeps the residuals it consumes consistent with the hypocentre the catalog goes on to publish.

#### Location-PDF quality metrics — `NLL_run/pdf_metrics.py`

Run automatically by `run_SSST.py` as its last step, once the merge and the rematch are done. It reads the per-event `.scat` scatter clouds of each zone's final SSST iteration, joins them to the merged catalog on `(source, publicId)`, rewrites `RESULT/SSST_result.csv` in place with eleven extra columns, and saves the diagnostic figures described below. Runtime is ~65 s for ~46 000 events (~50 s for the metrics, ~15 s for the figures).

The step is deliberately **non-fatal**: at that point the campaign's results are already written, so a metrics failure prints a warning instead of making a multi-hour run look like it failed. It is also idempotent — re-running replaces the columns rather than duplicating them — so it can be repeated standalone at any time, on any campaign.

The metrics exist to decide **which events are trustworthy enough to keep**:

- **Ψ** (`Psi`, with `J = −ln Ψ`) — how much of the PDF's shape the confidence ellipsoid actually captures. `Ψ = 1` is an exactly Gaussian PDF; smaller means the true PDF is curved, clustered or heavy-tailed and the ellipsoid is a poor stand-in for it. Published as a quality indicator only: it says the ellipsoid is *wrong*, not whether it is too big or too small, so it never rejects an event on its own.
- **C_68** (`C68`) — the fraction of the PDF actually inside the nominal one-sigma ellipsoid. This is the one metric that carries a direction, so it drives the keep/reject decision: `C_68 > 0.6827` means the quoted `ERH`/`ERZ` are **conservative** (safe — the true region is tighter than stated), while `C_68 < 0.6827` means they are **over-confident**, i.e. the location is less well constrained than the error bars claim.
- **`dip_stat` / `dip_pval` / `dip_sep_km`** — Hartigan's dip test on the depth marginal, plus the width of its modal interval in km. A small p-value means the depth PDF has two competing solutions; `dip_sep_km` says how far apart they are. **Rejection requires both** — `dip_pval < 0.05` *and* `dip_sep_km > true_erz` — because the dip statistic is a vertical distance on the ECDF and is invariant under rescaling of the depth axis: modes 0.5 km apart and modes 50 km apart give the same statistic and the same p-value. Rejecting on the p-value alone therefore discards events whose bimodality is real but physically irrelevant, both modes sitting inside the quoted error where `depth ± true_erz` already covers them honestly. That is the common case here: among `p < 0.05` events the median modal interval is 0.59 × `true_erz`, and adding the separation term moves the reject count from 5 329 to 661. The threshold lives in `_DIP_SEP_ERZ_FACTOR`.

`J` and `C_68` cannot be read on their own: the k-NN entropy estimator is biased in a sample-size-dependent way, and the `C_68` null spread is *not* binomial (μ and Σ are fitted on the very samples being tested, which suppresses the scatter — 0.013 simulated vs 0.021 binomial at n≈479). Each event therefore also carries `J_null_p95` and `C68_sigma_n`, simulated from Gaussian clouds of that event's own sample count, so a cut is a one-line comparison: `C68 < 0.6827 − 3·C68_sigma_n` (over-confident) or `J > J_null_p95` (non-Gaussian).

> All three measure whether the reported uncertainty is **self-consistent**, not whether the location is **accurate**. Velocity-model error is a bias none of them can see: an event can score perfectly and still be systematically mislocated. Thresholds are not established values.

Observed over the current catalog: `C_68` is over-confident for only 0.1 % of events and conservative for 57 %, while the two-condition dip rule rejects 1.43 % (661 events; the p-value alone would reject 11.5 %). Both degrade monotonically with azimuthal gap (median Ψ 0.90 at gap < 90° falling to 0.17 above 270°) and with falling phase count — and so does the dip reject rate, which is what one expects if the events with genuinely far-apart depth modes are the poorly-constrained ones.

##### Diagnostic figures

The same run writes four outputs to `complem_figures/pdf_metrics/`, which is what makes the statements above checkable rather than asserted. Plotting is wrapped in its own guard, so a failure there leaves the annotated CSV untouched; `--figures false` skips it entirely.

| File | Answers |
|------|---------|
| `<run>_metrics_vs_quality.pdf` | How do Ψ, `C_68` and the dip statistic relate to the classical quality indicators? 3 metrics × depth, RMS, Nphs, Gap, Dist — hexbin density, median + IQR over equal-count bins of the indicator, and the Gaussian null of that same bin |
| `<run>_gridmap.pdf` | Where in the Pyrenees is the reported uncertainty least trustworthy? Windowed-median lon/lat maps, 0.02° cells |
| `<run>_gridmap_depth.pdf` | Does that pattern change with depth? The same maps split into depth quartiles, colour scale shared across the slices of a row |
| `<run>_voxels_<metric>.html` | Interactive lat/lon/depth cells for each metric (Plotly, self-contained) |

Two details keep the figures honest when the catalog changes. The maps colour by `C68_z = (C_68 − 0.6827) / C68_sigma_n` instead of raw `C_68` — a cell median mixes events with different sample counts, and only the z-score puts them on a common scale — and the 1-D null overlays are computed per bin from the `J_null_p95` / `C68_sigma_n` of the events *in that bin*, so they bend by themselves if `n_scat` correlates with the indicator. Nothing is hardcoded from the near-constant `n_scat` of `ssst_run1`.

Latitude and longitude are deliberately absent from the 1-D grid: a 1-D latitude panel marginalizes over longitude and hides exactly the structure the maps are there to show.

The campaign configuration (`RUN_NAME`, `CHAR_DISTS`, `VPVS`, core counts, `LSPHSTAT` — whose `NRdgsMin` doubles as the NLLoc min-phases threshold, read back from the generated Loc2ssst control so the two selections cannot diverge) lives at the top of `run_SSST.py`. Per-zone run journals are written to `run/ssst/log/`; intermediate location outputs (`loc_ssst_corr<i>/`) are deletable after validation (each journal ends with the ready-to-paste commands), and the chunk directories of each iteration are deleted automatically after each merge.

---

### 8. Final Bulletin Export

`final_steps.py` — the last entry point, run once `run_SSST.py` has finished. Everything it needs already exists but is split across two files, each holding half the catalog:

| Input | Holds | Missing |
|-------|-------|---------|
| `obs/SSST_result.obs` | magnitudes, all 1 001 095 phase picks | numbers rounded to the header's display precision |
| `RESULT/SSST_result.csv` | full-precision hypocentres, confidence ellipsoid, PDF quality metrics | magnitude, picks |

The stage merges them on `publicId` (1:1 across 46 224 events) into **`RESULT/FINAL.xml`**, a QuakeML 1.2 bulletin — the first output of this project in a format readable outside it. Hypocentres, uncertainties and quality come from the CSV at full precision; the `.obs` supplies only the magnitude and the picks. Runtime ~5 min, peak ~3.4 GB RAM, ~0.86 GB output.

Every event is exported, each carrying a boolean **`pyr:usable`** flag plus the metrics behind it, so a consumer can filter without recomputing anything and can re-tune the cut without rerunning the stage.

```
usable = (C68 − 0.6827) / C68_sigma_n ≥ −2   and   not dip_reject   and   metrics present
```

Measured on the current catalog: **45 415 usable / 809 unusable** — 661 `dip_bimodal`, 155 `c68_overconfident`, 4 `no_pdf_metrics` (11 events fail on two counts). The reason string is exported alongside the flag.

Three things that rule deliberately does *not* do. **Ψ never rejects**: it is directionless, and 81 % of this catalog exceeds its own Gaussian null on `J`, so gating on it would discard 85 % of the events — it is exported as an indicator only. And the `C_68` cut is taken **against the simulated null, not against 0.6827**: with the median `C68_sigma_n ≈ 0.0124` the −2 σ cut sits at `C_68 < 0.658`, flagging ~160 events where a raw `C_68 < 0.6827` would flag thousands — the events in between are within sampling noise of perfect coverage. And **the dip p-value never rejects on its own**: it cannot see how far apart the modes are, so it is paired with `dip_sep_km > true_erz` (above), which is the difference between 661 events flagged and 5 329.

Station codes are resolved back to real network/station names: the picks carry the project's unified code (`FR.0041`), which is matched against `alternate_code` in `GLOBAL_inventory.xml` and resolved by the pick's own date. 80 codes cover more than one station (near-duplicates merged within 20 m); 79 have disjoint epochs and resolve on the date alone. When several candidates cover the same date, the winner is the one the merge itself would keep — `station_priority()` from [§1](#1-station-inventory-fusion), so the permanent network wins and the **`XX` placeholder loses last** (this currently decides one code, `FR.0013` → `FR.MTHF`, 1 331 picks; it is logged). The export applies that ranking **after** the date, though, not instead of it: 79 of the 80 multi-station codes have disjoint epochs, so the pick's own date usually decides, and it should — a 2015 pick must not be named after a station installed in 2019. For 37 codes this means the name printed here and the member GTSRCE took coordinates from are different stations of the same group. They sit within 20 m of each other and, since the merge now shares elevations within a group, within 61 m in elevation.

Resolution is a judgement, so it is never destructive. Every pick keeps its unified code in `pyr:unifiedCode`, and whenever the code covers more than one station the **stations that were not chosen are exported too**, in `pyr:alternateStations` — including those whose epoch does not cover the pick, since an epoch that excludes a real pick is itself a sign the metadata may be wrong. The stations behind one unified code sit within 20 m of each other, so if the resolved label turns out to be wrong, those are exactly the places to fetch waveforms from instead. This affects ~15 % of picks (150 395 on multi-station codes).

#### QuakeML output

Standard QuakeML carries the hypocentre, the `OriginQuality`, the `OriginUncertainty` with its full `ConfidenceEllipsoid`, the magnitude, and one `Pick` + `Arrival` per phase.

Each event carries **two origins**, which QuakeML supports directly. The preferred one — the one `preferredOriginID` points at, and the only one most consumers will ever look at — is the PDF expectation, and it holds everything: the origin time, the quality, the uncertainty and the arrivals. The maximum-likelihood hypocentre is the second origin, at `smi:pyrenees/origin/<publicId>/maxlike`, carrying only its coordinates and its own origin time. The two are told apart by `pyr:locationEstimator`, `"expectation"` or `"maximum_likelihood"`.

The second origin is deliberately bare. Under `SAVE_NLLOC_EXPECTATION` the RMS, the azimuthal gap, the station distances and the confidence ellipsoid are all evaluated at the expectation, and the ellipsoid is a moment about it; copying them onto the maximum-likelihood origin would recreate exactly the mismatch described in [Which hypocentre the catalog reports](#which-hypocentre-the-catalog-reports). A consumer who wants the mode gets the mode, and is not handed an error bar that belongs to a different point.

Everything QuakeML has no element for goes into a custom namespace, `http://shallow-depth-dl-catalog/quakeml/1.0`, carried by the prefix **`pyr`** (for Pyrenees). This is not decoration: QuakeML 1.2's schema is closed, so a bare `<usable>` element inside the standard namespace would make the file invalid. Foreign namespaces are skipped under lax processing, which is how `FINAL.xml` validates against the QuakeML 1.2 schema while carrying fields the standard has never heard of. The prefix itself is cosmetic — consumers match on the URI — and the URI is an identifier, not a resolvable address.

| Level | Field | Meaning |
|-------|-------|---------|
| event | `usable` | `true` / `false` — the rule above |
| event | `rejectReason` | `''`, or `dip_bimodal` / `c68_overconfident` / `no_pdf_metrics`, comma-joined |
| event | `Psi`, `J`, `C68`, `C68Z` | the quality metrics; `C68Z` is the null-normalized z-score the cut uses |
| event | `JNullP95`, `C68SigmaN`, `nScat` | the per-event Gaussian null and the sample count backing it |
| event | `dipStat`, `dipPval`, `dipSeparationKm`, `dipReject` | Hartigan's dip test on the depth marginal, the width of its modal interval in km, and the two-condition flag built from them |
| event | `publicId`, `sourceZone` | the catalog's join key, and which NLL zone won the dedup |
| origin | `locationEstimator` | `expectation` (preferred origin) or `maximum_likelihood` (secondary) |
| origin | `pdfVolume`, `ellipsoidVolume` | the two volume measures, for direct comparison |
| origin | `ellipsoidAz1` … `ellipsoidLen3` | the raw NLLoc ellipsoid columns, so the QuakeML conversion stays auditable |
| pick | `unifiedCode` | the project's internal station code, before resolution |
| pick | `alternateStations` | the other real stations sharing that unified code, comma-joined — fallbacks for waveform retrieval; absent when the code maps to one station |
| pick | `pickOrigin` | provenance: `OMP`, `FDSN`, `LDG`, `IGN`, `ICGC`, `TEMP_*` |
| pick | `relativeTiming` | `true` when NLLoc used the pick via S−P relative timing (the `*` flag), not its absolute time |

#### Which file to read

`read_events('RESULT/FINAL.xml')` loads the whole catalog in one call — 46 224 events in 255 s, zero warnings, every field surviving the round trip — **but it peaks at ~15 GB of RAM**, because 46 224 events and 1 001 095 picks become roughly two million Python objects, each with a registered `ResourceIdentifier`. Measured per object: 12.5 KB per pick+arrival and 55.8 KB per event, so **the picks are 83 % of the cost**.

The stage therefore writes the same catalog three ways. They are redundant on purpose: pick whichever matches what you are doing.

| File | Holds | Size | `read_events` cost |
|------|-------|------|--------------------|
| `FINAL.xml` | everything, one file | 834 MB | 15.1 GB / 255 s |
| `FINAL_catalog.xml` | all 46 224 events, **no picks or arrivals** | 178 MB | 2.5 GB / 33 s |
| `FINAL_<from>_<to>.xml` | the full bulletin cut into 5-year calendar periods | 0.1–322 MB | 0.05–6.0 GB |

- **Want a catalog?** `FINAL_catalog.xml` — origins, magnitudes, quality, uncertainties and the `pyr:` usability fields, minus the phases. This is the FDSN convention (`includearrivals=false`), and it is what most analyses need.
- **Want phases for a period?** the matching `FINAL_<from>_<to>.xml`. Which file holds an event follows from its date alone. The parts partition `FINAL.xml` exactly — no event appears twice, none is missing — and each is independently schema-valid, so `read_events` accepts a glob to recombine any subset.
- **Want everything at once?** `FINAL.xml`, on a machine with the memory for it.

Parts are far from equal, and five-year bins do not make them so: pick density rises from ~7 picks/event in the 1980s to ~28 in the 2020s, and the network densified on top of that. `FINAL_2020_2024.xml` alone holds **14 537 events and 401 295 picks — 40 % of every pick in the catalog** — and measures 6.0 GB to load, against 88 KB and 0.05 GB for `FINAL_1975_1979.xml`. So the split removes the 15 GB requirement but leaves one genuinely heavy file. If a harder ceiling is ever needed, `_SPLIT_YEARS` in `NLL_run/export_quakeml.py` is the single knob: setting it to 2 or 1 re-cuts the dense years without touching anything else.

```python
from obspy import read_events
cat = read_events('RESULT/FINAL_catalog.xml')   # ~2.6 GB, seconds
cat[0].extra['usable'].value                    # 'true' / 'false'

phases = read_events('RESULT/FINAL_2020_2024.xml')
phases[0].picks[0].extra['alternateStations'].value
```

For anything that does not need objects in memory, **stream instead**: `lxml.etree.iterparse` walks all 46 224 events of the full file in seconds at flat memory, and the `pyr:` fields are plain child elements of `<event>` and `<pick>` (namespace `{http://shallow-depth-dl-catalog/quakeml/1.0}`). That is how every verification of this file was run.

Three caveats a reader of `FINAL.xml` needs:

- **The hypocentre is the PDF expectation, and the mode is right beside it.** The preferred origin is the mean of the location PDF; the second origin is its maximum-likelihood point, and for a non-Gaussian event the two can sit further apart than the quoted error. That is not an inconsistency in the file — it is the honest shape of the answer, and `Ψ`, `C_68` and the dip columns say which events it applies to. A reader who wants one number should take the preferred origin, which is the one the uncertainty and the quality actually describe. The reasoning, and what changed when the catalog switched, is in [Which hypocentre the catalog reports](#which-hypocentre-the-catalog-reports).
- **The DOF scaling is mixed inside `OriginUncertainty`, by design.** The `ConfidenceEllipsoid` axes keep NLLoc's 3-DOF 68 % scaling, while `horizontalUncertainty` and the depth uncertainty carry the catalog's own `true_erh` (2-DOF) and `true_erz` (1-DOF). Those are the values the rest of the catalog uses, but neither is a projection of the other. NLLoc reports azimuth/dip for axes 1 and 2 only, so the major-axis orientation QuakeML requires is reconstructed from their cross product.
- **`Mag 0.00` on 5 010 events is a placeholder, not a measurement**, and is written through unchanged because magnitudes are recomputed later. The raw OMP `.mag` files hold the literal string `0.0` for ~19 % of events; the neighbouring `0.1` bin holds 148, and the zero fraction falls from ~40 % to 0.3 % in 2013 when OMP began computing ML systematically. All 5 010 are `(ML, OMP)` — `LDG.py` and `RESIF.py` skip magnitude-less events instead of filling them. Do not fit that bin as data.

---

## Complementary Analysis

Two driver scripts run the `complem_figures/` modules:
- `generate_complem_figures.py` — matplotlib figures: depth histograms, Gutenberg-Richter distributions, and per-period depth and error maps
- `generate_complem_maps.py` — PyGMT event maps for each of the 6 NLL zones and the final catalog

Each module can also be run standalone:

| Script | Description |
|--------|-------------|
| `event_maps.py` | Geographic maps of seismicity (from `.obs`, `.txt`, or `.csv` NLL summary) |
| `depth_maps.py` | Per-period windowed-median depth maps |
| `error_maps.py` | Per-period spatial distribution of location uncertainties (ERH, ERV) |
| `depth_histogram.py` | Histogram of event depths |
| `gutenberg_richter.py` | Magnitude-frequency distribution (Gutenberg-Richter law) |
| `cross_section.py` | Vertical cross-sections of seismicity |
| `station_map.py` | Map of seismic stations |
| `zone_map.py` | Overview map of the 6 NLL relocation zones |
| `event_ranking.py` | Ranks events by pdfVolume/ellipsoidVolume change (NLL → SSST); flags multi-zone and zone-changed events; carries the PDF-quality columns (Ψ, C_68, dip test) when present; optional whole-region gridmap PDF |
| `plot_pdf_cloud.py` | Interactive 3D Plotly visualization of one event's NLLoc PDF scatter-cloud across SSST iterations |
| `ssst_evolution.py` | Per-zone plot of pdfVolume/EllipsoidLen3/RMS evolution across SSST iterations (convergence QC) |
| `ssst_corrections.py` | Reconstructs and maps the SSST travel-time corrections themselves — a per-station/phase atlas, the across-station spread, and the displacement SSST actually produced (see [below](#mapping-the-ssst-corrections--complem_figuresssst_correctionspy)) |
| `station_colocation.py` | Diagnostic for the 20 m co-location radius of [§1](#1-station-inventory-fusion): the observed distribution of station separations against what the merge actually merged |

`event_ranking.py`, `plot_pdf_cloud.py`, and `ssst_evolution.py` are standalone diagnostics for the SSST stage, run directly rather than wired into `generate_complem_figures.py` / `generate_complem_maps.py`; they read `RESULT/NLL_result.csv`, `RESULT/SSST_result.csv`, and the per-zone `run/nll_loc/` / `run/ssst_loc/<run-name>/` outputs directly. `event_ranking.py` picks up the PDF-quality columns (Ψ, C_68, dip test) automatically once `NLL_run/pdf_metrics.py` has annotated `SSST_result.csv`, and runs without them otherwise — they are post-SSST only, since `run/nll_loc/` holds no `.scat` clouds to compare against. It also shares that module's windowed-grid helper (`windowed_stat_grid`), so its gridmaps and the PDF-metric maps bin space identically.

`station_colocation.py` is the one diagnostic that looks at stage 1 rather than at the relocation. It reads the 2 414 pre-merge stations from `stations/*/*.xml` and the groups recorded in `GLOBAL_inventory.xml`, and asks whether 20 m is the right radius. The separations are bimodal with a real trough between 20 m and 200 m, so the threshold does sit in a gap — but it is not a free choice: 651 pairs fall within 20 m, 811 within 30 m, and the extra ones are inside the dense nodal arrays (`AM`, `8M`, `XI`), where genuinely distinct nodes would start collapsing into one code. Crucially it reads "merged" from the **outcome** rather than from the distance, which is what exposes the 423 pairs that merged from *beyond* the radius through the single-linkage chaining of `_combine_close_stations` and the 38 inside it that did not, because `check_inventory` removed one of them by hand.

```bash
python complem_figures/station_colocation.py             # → complem_figures/station_colocation/
python complem_figures/station_colocation.py --threshold 30
```

One folder here has no script of its own: `complem_figures/pdf_metrics/` receives the figures written by `NLL_run/pdf_metrics.py` (see [above](#location-pdf-quality-metrics--nll_runpdf_metricspy)).

### Mapping the SSST corrections — `complem_figures/ssst_corrections.py`

Standalone diagnostic that answers "what did SSST actually correct, and where". `Loc2ssst` writes an explicit correction grid per station and phase (`ssst_corr<i>/*.ssst.*`), but `run_SSST.py` deletes those directories to reclaim disk, so nothing of the correction field survives a finished campaign. The module recomputes it from the per-event `.hyp` locations, which do survive, using the formula `Loc2ssst` itself uses (`Loc2ssst.c:1091-1102`):

```
corr_i(sta, pha, x,y,z) = SUM_e res_e * w_e / SUM_e w_e
w_e                     = exp(-|x_e - x|² / L_i²) + weight_floor
```

The distance is **event-to-grid-node** — the station's own position never enters the weight, which is what makes the correction *source*-specific. `res = obs − pred` and the correction is **added** to the predicted travel time, so **red = arrivals later than the 1-D model predicts** (the real path is slower), blue = earlier. Each iteration's field is an **increment** on the previous grids; the correction the final grids carry is the sum of the five.

Two products, from one reconstruction:

- `<run>_station_atlas.pdf` — one page per station×phase field (all 300 with ≥ 100 usable arrivals, `--min-picks`): the five increments (`9999 → 50 → 15 → 5 → 1 km`) then their cumulative total, as depth-slice maps on a shared diverging scale. Nodes with no event inside the smoothing kernel are greyed — the field there is only the station static term. A station falling in two zones is drawn in the zone where it has the most arrivals, since the two zones ran separate `Loc2ssst` passes over different event sets and their fields are independent.
- `<run>_spread_map.pdf` — one catalog-wide map: at each event, the spread **across its recording stations** of the total correction applied to its P picks. A correction common to every station of an event is absorbed exactly by the origin time and cannot move the hypocentre, so the across-station dispersion, not the mean, is the part of the field that *can* relocate. It is **not** a map of where SSST actually did relocate — see the caveat below.

```bash
python complem_figures/ssst_corrections.py                   # all three, ~7 min
python complem_figures/ssst_corrections.py --product displacement   # the impact map, <1 s
python complem_figures/ssst_corrections.py --product atlas --min-picks 300
python complem_figures/ssst_corrections.py --extract-only    # fill the parse cache
# presentation-quality pages for two chosen fields (~55 s each)
python complem_figures/ssst_corrections.py --product atlas \
    --stations FR.0041:P,RD.0038:P --map-spacing 0.005
```

Parsing the ~276 k per-event `.hyp` files takes ~100 s and is cached as `.npz` per (zone, iteration) in `run/ssst_corrections_cache/` (33 MB), so only the first run pays for it.

**How the maps are drawn.** Nothing is interpolated between map nodes: `pcolormesh` draws one flat cell per node, and every node is an independent evaluation of the formula above. The smoothness on the page *is* the Gaussian kernel — the field is intrinsically smooth at scale `L`. That makes the node spacing a real constraint: the `--map-spacing` default of 0.02° is ~2.2 km in latitude and ~1.6 km in longitude at 43°N, so the `L = 1 km` panel is **undersampled** and part of its speckle is aliasing. The full 300-page atlas keeps the coarse spacing deliberately (212 s); drop to `--map-spacing 0.005` with `--stations` for the pages that need to be publication-grade, which resolves the finest panel properly at ~55 s per page.

**A shared map frame is optional.** By default each page is framed on its own station's events — tight, but not comparable from page to page. `--extent=lon0,lon1,lat0,lat1` (the `=` form is required when `lon0` is negative) puts every page on one frame, clipped to each zone's `LSGRID` as read from `run/ssst/run_<N>_SSST.in`. That clipping matters: zone 2's correction grid spans only lon −2.22 to 1.57, so a full-Pyrenees frame is largely territory where `Loc2ssst` computed nothing, and it is greyed rather than filled with the station's static term. On a `-2.25,3.5,41.75,43.75` frame each zone's grid covers ~68 % of the page.

**The slice depth is data-driven, not fixed.** `--depth` defaults to each station's own median event depth (~6.9 km for `FR.0041`) rather than a round number. This matters because the kernel is 3-D: a slice `dz` off the event mass damps every contribution by `exp(−dz²/L²)`. That is harmless at `L = 15 km` but devastating at `L = 1 km` — an earlier fixed 10 km slice, against a catalog median depth of 6.6–7.3 km, cut `FR.0041`'s median summed weight in the finest panel from 5.7 to 2.0 and left it looking far emptier than the data warrants. Pass `--depth` only to compare stations on a common plane, accepting that cost.

**Validated against the binary.** `Loc2ssst` was re-run for real on zone 6 (station `FR.0047`, P and S) at both ends of the schedule and compared node by node over all 303 × 313 × 40 nodes: max |diff| = 0.0000 ms at `L = 9999 km` (exact to float32) and 0.0736 ms at `L = 1 km`, against correction amplitudes of ±0.36 s. `Loc2ssst` independently reported "652 location files read, 163 accepted", matching the module's `LSPHSTAT` selection exactly. The check is cheap because `ihave_time_input_grids = flag_out_grid * flag_nlloc_outfile` (`Loc2ssst.c:590`): omitting `LOCFILES` from the control file makes `Loc2ssst` write the correction grid and skip the travel-time grids, so no `Grid2Time` rebuild is needed.

- `<run>_displacement_map.pdf` — **the impact map**: where SSST actually moved events, and which way. Each event's iteration-0 hypocentre against its final one — same picks, same parameters, only the corrected grids differing — so the difference is the corrections and nothing else, matched by `publicId`. Top panel: median distance moved, with arrows for the median horizontal displacement per 0.25° cell. Bottom panel: median change in depth, signed (red = pushed deeper). Over 53 754 events the median move is **1.33 km**, p90 **6.85 km** — but it is ~1 km through the dense central network and **3–4 km at the western and eastern margins**, which is the headline: SSST moves poorly-constrained edge events most.

**The spread map is not an impact map.** Measured against the *pure* SSST displacement — each event's iteration-0 location versus its final location, same picks, only the grids differing — over 39 080 events:

| | |
|---|---|
| `ρ(spread, displacement)` | **+0.10** |
| `ρ(spread, displacement \| Nphs)` | **+0.20** |
| median displacement, bottom → top spread quintile | 0.95 km → 1.27 km |

Disagreement between stations is *necessary* for an event to move but nowhere near sufficient: what displaces a hypocentre is the **azimuthal pattern** of the differential correction, and a standard deviation discards all of that geometry — spread distributed evenly around the azimuths largely cancels. The marginal ρ understates even this, because `Nphs` suppresses it (well-recorded events carry more spread, ρ = +0.42, yet resist moving, ρ = −0.18), the same confound documented for RMS against ERH in `pdf_metrics.py`. Use `--product displacement` for impact. (Those ρ figures are computed over the 39 080 events carrying a spread — ≥ 5 P picks — whose median displacement is 1.09 km; the displacement map itself uses all 53 754 matched events, median 1.33 km.)

Both catalog-wide maps share `windowed_median_map`, so they bin space identically. It is a **median** over a circular window — longitude scaled by `cos(lat)` so the window is round in km — because displacement is heavy-tailed (one event moves 99.5 km) and a mean would let single events paint whole neighbourhoods.

Two further properties worth knowing before reading the figures: the increments are **not** monotonically decreasing — they peak at `L = 15 km`, the scale at which the residual field carries genuine spatial structure — and the events feeding the corrections are a strict subset of the catalog, since `LSPHSTAT` (RMS ≤ 0.15 s, Nphs ≥ 6, gap ≤ 200°, Len3 ≤ 10 km) admits 23 727 events at iteration 0, rising to 34 674 by the final relocation as the locations improve.

Map modules apply a quality filter (erh ≤ 3 km, erv ≤ 3 km, gap ≤ 300°, rms ≤ 0.5 s) by default; use `--no-filter` for pre-relocation catalogs where errors are unavailable.

`zone_Arette/` contains a focused analysis of the Arette seismic zone, including gap/RMS statistics across different station distance cutoffs and yearly temporal analysis.

---

## Dependencies

| Package | Use |
|---------|-----|
| `obspy` | Seismic data I/O, FDSN client, inventory management |
| `pandas`, `numpy` | Data manipulation |
| `pyarrow` | Parquet I/O (`temp_picks/convert_picks.py`'s `TEMP_STB` format) |
| `scipy` | ODR regression, spatial queries (KDTree), statistics |
| `scikit-learn` | Regression diagnostics (R²) for magnitude models |
| `matplotlib`, `seaborn` | Plotting |
| `plotly` | Interactive 3D PDF-cloud visualization — `complem_figures/plot_pdf_cloud.py` |
| `diptest` | Hartigan's dip test for depth multimodality — `NLL_run/pdf_metrics.py` |
| `xarray` | Grid handling for cross-sections |
| `pygmt` | Geographic maps |
| `joblib` | Magnitude model serialization |
| `requests` | ICGC catalog fetching |
| `pyproj` | Coordinate transformations |
| **NonLinLoc** | Probabilistic earthquake location — Vel2Grid, Grid2Time, NLLoc, Loc2ssst (external tool, invoked automatically by `run_NLL.py` / `run_SSST.py`) |
| **Pyrocko** / **cake** | Theoretical travel-time computation (`temp_picks/build_theoretical_tables.py`) |

---

## A note on AI assistance

Parts of this codebase were written or modified with the help of **[Claude Code](https://claude.ai/code)** (Anthropic). As a researcher, I believe in being transparent about the use of AI tools in scientific work. All AI-generated code in this project has been reviewed before being committed to the main branch.
