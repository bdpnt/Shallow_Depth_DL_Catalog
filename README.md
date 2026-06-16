# Shallow_Depth_DL_Catalog

A pipeline for building a unified, publication-quality earthquake catalog for the **Pyrenees region** (lat 41–45°N, lon -3 to 4°E), integrating data from five independent seismic networks over the period **1978–2025**.

The workflow covers catalog fetching, station inventory fusion, magnitude harmonization, event merging, probabilistic earthquake relocation with NonLinLoc, and result visualization.

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

All catalogs are converted to a common `.obs` format, magnitudes are harmonized to a common **ML (LDG)** scale, and events are merged into a single `GLOBAL.obs` bulletin. Earthquakes are then relocated using **NonLinLoc** across 6 geographic sub-zones, and final results are compiled into `RESULT/FINAL.csv` and `obs/FINAL.obs`.

---

## Project Structure

```
Shallow_Depth_DL_Catalog/
│
├── fetch_all_bulletins.py        # Entry point: fetch & convert all catalogs
├── build_global_inventory.py     # Entry point: fuse all station inventories
├── build_global_bulletin.py      # Entry point: harmonize & merge catalogs
├── prepare_nll_inputs.py         # Entry point: prepare NLL run files (6 zones)
├── generate_nll_corrections.py   # Entry point: prepare second-pass NLL run files
├── finalize_nll_catalog.py       # Entry point: compile and match final catalog
├── add_temp_picks.py             # Entry point: augment FINAL.obs with external picks
├── generate_complem_figures.py   # Entry point: matplotlib figures (seisbench_env)
├── generate_complem_maps.py      # Entry point: PyGMT event maps (pygmt_env)
├── run_gamma_detection.py        # Entry point: standalone PhaseNet/GaMMA detection
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
│   ├── list_magnitude_types.py
│   ├── generate_magnitude_models.py
│   ├── apply_magnitude_models.py
│   ├── add_temporary_picks.py
│   ├── filter_events_by_aoi.py
│   ├── fuse_bulletins.py
│   └── plot_global_catalog_map.py
│
├── NLL_run/                  # NonLinLoc workflow modules
│   ├── generate_regional_runfiles.py
│   ├── append_station_delays.py
│   ├── export_locdelay_info.py
│   ├── parse_nll_output.py          (deprecated)
│   ├── filter_distant_picks.py
│   ├── match_pre_post_relocation.py
│   └── merge_regional_results.py
│
├── complem_figures/          # Visualization & statistical analysis
│   ├── event_maps.py
│   ├── depth_maps.py
│   ├── error_maps.py
│   ├── depth_histogram.py
│   ├── gutenberg_richter.py
│   ├── cross_section.py
│   ├── station_map.py
│   └── zone_map.py
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
├── ORGCATALOGS/              # Raw input catalogs (not modified)
├── obs/                      # .obs bulletin files (source + merged)
├── stations/                 # Station inventories (XML + unified)
├── run/                      # NLL run configuration files (.in)
├── loc/                      # NLL output files per zone
├── RESULT/                   # Merged relocation results (.csv)
├── model/                    # Velocity model grids (NLL)
├── time/                     # Travel time grids (NLL)
└── MAGMODELS/                # Saved magnitude conversion models
```

---

## Conventions

All Python scripts in this project share the same interface contract:
- **CLI**: every script accepts `--help` and can be run directly from the command line
- **Public API**: every module is importable as a Python package (e.g. `from temp_picks.match_picks import match_picks`)
- **Logging**: timestamped log files are written to a `console_output/` directory local to each sub-pipeline

---

## Pipeline

### 1. Station Inventory Fusion

**Script:** `build_global_inventory.py`  
**Module:** `fetch_inventory/merge_station_inventories.py` → `merge_inventory()`

Merges all station XML inventories (FDSN networks + OMP) into a single unified inventory. Each station receives a unique code; duplicates within 20 m are removed. OMP CSV data is pre-processed with `_remove_fdsn_duplicates.py`, `_fill_missing_elevations.py`, and `_convert_csv_to_stationxml.py` before fusion.

**Outputs:**
- `stations/GLOBAL_inventory.xml` — unified QuakeML inventory
- `stations/GLOBAL_code_map.txt` — mapping between original and unified station codes

---

### 2. Catalog Fetching & Conversion

**Script:** `fetch_all_bulletins.py`  
**Modules:** `fetch_obs/` (one module per source)

Downloads or reads each catalog and converts it to the `.obs` format. RESIF and ICGC are fetched dynamically; IGN, LDG, and OMP are read from local files in `ORGCATALOGS/`.

**Outputs:** individual `.obs` files in `obs/`  
(e.g. `RESIF_20-25.obs`, `IGN_20-25.obs`, `OMP_78-19.obs`, …)

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

**Outputs:**
- `obs/GLOBAL.obs` — unified catalog
- `obs/MAPS/` — statistics figures
- `MAGMODELS/` — serialized magnitude models

---

### 4. Earthquake Relocation (NonLinLoc)

The study area is too large for a single NLL run, so it is divided into **6 geographic zones**. Each zone is processed independently, then results are merged.

#### Pre-run — `prepare_nll_inputs.py`

Calls `NLL_run/generate_regional_runfiles.py` → `generate_run()` for each zone:
- Generates `obs/GLOBAL_1.obs` … `obs/GLOBAL_6.obs` (regional subsets, with far picks removed)
- Generates `stations/GTSRCE_1.txt` … `stations/GTSRCE_6.txt` (station lists)
- Generates `run/run_1.in` … `run/run_6.in` (NLL configuration files)

Then run NLL externally for each zone:

```bash
Vel2Grid run/run_<N>.in
Grid2Time run/run_<N>.in
NLLoc run/run_<N>.in
```

#### Second pass — `generate_nll_corrections.py`

For each zone, this script does two things:
1. Cleans up `.hdr` files left by NLL in the `loc/GLOBAL_<N>/` folder.
2. Reads per-station average residuals (LOCDELAY entries) from the first run and appends qualifying station delay corrections to generate second-pass run files `run/run_<N>_PR.in`.

Run NLL again:

```bash
Vel2Grid run/run_<N>_PR.in
Grid2Time run/run_<N>_PR.in
NLLoc run/run_<N>_PR.in
```

#### Diagnostic — `NLL_run/export_locdelay_info.py`

Called automatically by `generate_nll_corrections.py` at the end of the second-pass run. Reads the LOCDELAY station corrections from all second-pass run files and exports them to `run/locdelays/locdelay_summary.txt`, keeping only entries with |residual| > 0.3 s. Useful for identifying stations with systematically biased travel-time residuals.

---

### 5. Post-relocation Processing

**Script:** `finalize_nll_catalog.py`  
**Modules:** `NLL_run/merge_regional_results.py`, `NLL_run/match_pre_post_relocation.py`

1. Cleans up `.hdr` files left by NLL in each `loc/GLOBAL_<N>/` folder.
2. Reads the 6 per-zone NLL CSV summaries (`loc/GLOBAL_<N>/GLOBAL_<N>.obs.sum.grid0.loc.csv`), deduplicates events that appear in multiple overlapping zones (kept: lowest `pdfVolume`), and writes → `RESULT/FINAL.csv`. True horizontal/vertical errors (`true_erh` / `true_erz`) are derived from the 3-D confidence ellipsoid.
3. Rematches relocated events back to `obs/GLOBAL.obs` via the `publicId` field to recover metadata absent from NLL output (magnitude, pick details, etc.).
4. Saves matched events → `obs/FINAL.obs`

---

## Complementary Analysis

Two driver scripts run the `complem_figures/` modules in **different conda environments**:
- `generate_complem_figures.py` (`seisbench_env`) — matplotlib figures: depth histograms, Gutenberg-Richter distributions, and per-period depth and error maps
- `generate_complem_maps.py` (`pygmt_env`) — PyGMT event maps for each of the 6 NLL zones and the final catalog

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

Map modules apply a quality filter (erh ≤ 3 km, erv ≤ 3 km, gap ≤ 300°, rms ≤ 0.5 s) by default; use `--no-filter` for pre-relocation catalogs where errors are unavailable.

`zone_Arette/` contains a focused analysis of the Arette seismic zone, including gap/RMS statistics across different station distance cutoffs and yearly temporal analysis.

---

## External Pick Ingestion (temp_picks)

Scripts in `temp_picks/` implement a self-contained sub-pipeline for ingesting picks from external sources into `obs/FINAL.obs`, producing `obs/FINAL_augmented.obs`. The root-level script **`add_temp_picks.py`** orchestrates steps 1–6 automatically.

| Step | Script | Description |
|------|--------|-------------|
| 1 | `build_theoretical_tables.py` | Uses Pyrocko's `cake` CLI to compute P/S travel-time envelopes across ±5% velocity models and source depths of 0–30 km, for epicentral distances 0–100 km → `tables_Pyr.csv` |
| 2 | `plot_travel_times.py` | QC figure: overlays all observed (distance, travel time) picks from a bulletin on top of the theoretical P/S bands. Run automatically as step 2; skipped if the figure already exists. Also usable as a standalone script. |
| 3 | `merge_omp_picks.py` | Merges all yearly OMP/PhaseNet CSV files from `picks_OMP/` subdirectories → `pick_files/merged_omp.csv`. Station `SMC` and year `2026` are excluded by default; configurable via `--drop-years`. |
| 4 | `merge_pyrenees_picks.py` | Concatenates RaspberryShake/PhaseNet `.txt` files from `picks_station_pyrenees/` and `picks_station_pyrenees2/` → `pick_files/merged_pyrenees.txt` and `pick_files/merged_pyrenees2.txt`. |
| 5 | `convert_picks.py` | Converts external pick files to the project's `.obs` pick line format; maps short station names to internal codes via `GLOBAL_code_map.txt`. Supports formats `TEMP_OBS`, `TEMP_RSB`, and `TEMP_OMP`; new formats are added as handler functions. Unresolved stations are reported as an end-of-run summary. |
| 6 | `match_picks.py` | For each converted pick, finds candidate events within a 60 s origin-time window, filters by theoretical travel-time residual (±0.1 s P, ±0.3 s S, plus ±2.5 s t0-error margin), and appends matched picks to the bulletin. Chains against `obs/FINAL.obs` → `obs/FINAL_augmented.obs`. Runs `sort_picks` automatically on the output. |
| 7 | `sort_picks.py` | Sorts all pick lines within each event block by ascending arrival time. Also usable as a standalone script on any bulletin. |

---

## Dependencies

| Package | Use |
|---------|-----|
| `obspy` | Seismic data I/O, FDSN client, inventory management |
| `pandas`, `numpy` | Data manipulation |
| `scipy` | ODR regression, spatial queries (KDTree), statistics |
| `scikit-learn` | Regression diagnostics (R²) for magnitude models |
| `matplotlib`, `seaborn` | Plotting |
| `xarray` | Grid handling for cross-sections |
| `pygmt` | Geographic maps (requires separate `pygmt_env` conda environment) |
| `joblib` | Magnitude model serialization |
| `requests` | ICGC catalog fetching |
| `seisbench`, `torch` | PhaseNet phase detection — `run_gamma_detection.py` |
| `gamma` | GaMMA event association — `run_gamma_detection.py` |
| `pyproj` | Coordinate transformations |
| **NonLinLoc** | Probabilistic earthquake location (external tool, run manually) |
| **Pyrocko** / **cake** | Theoretical travel-time computation (`temp_picks/build_theoretical_tables.py`) |

---

## A note on AI assistance

Parts of this codebase were written or modified with the help of **[Claude Code](https://claude.ai/code)** (Anthropic). As a researcher, I believe in being transparent about the use of AI tools in scientific work. All AI-generated code in this project has been reviewed and verified line-by-line before being committed to the main branch.
