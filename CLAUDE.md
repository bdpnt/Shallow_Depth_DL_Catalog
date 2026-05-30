# CLAUDE.md — Seisbench2025

## Project Overview

Seisbench2025 is an earthquake catalog processing and relocation pipeline focused on the **Pyrenees region** (latitude 41–45°N, longitude -3 to 4°E), covering seismicity from **1978 to 2025**.

The goal is to produce a unified, publication-quality earthquake catalog with harmonized magnitudes and improved hypocenter locations, by integrating data from 5 independent seismic networks.

---

## Pipeline Summary

The workflow follows 5 main stages:

### 1. Station Inventory Fusion
- Source: FDSN XML files + OMP CSV files (in `stations/`)
- Script: `build_global_inventory.py` → calls modules in `fetch_inventory/`
- Output: `stations/GLOBAL_inventory.xml` + `stations/GLOBAL_code_map.txt`
- Each station gets a unique code; duplicates are removed by distance threshold (20m)

### 2. Catalog Fetching & Conversion
- Sources: RESIF (FDSN), ICGC, IGN, LDG, OMP (in `ORGCATALOGS/`)
- Script: `fetch_all_bulletins.py` → calls modules in `fetch_obs/`
- Output: individual `.obs` files per source in `obs/`

### 3. Catalog Harmonization
- Script: `build_global_bulletin.py` → calls modules in `global_obs/`
- Steps:
  1. **remap_picks_to_unified_codes.py** — associates picks with unified station codes
  2. **generate_magnitude_models.py** — builds regression models to convert all magnitude types to ML
  3. **apply_magnitude_models.py** — applies the models to all `.obs` files
  4. **fuse_bulletins.py** — spatially/temporally merges all catalogs into `obs/GLOBAL.obs`
  5. **plot_global_catalog_map.py** — generates a map of the merged catalog
- Matching thresholds: 15 km distance, 2 s time, 1.5 magnitude units, ≥2 picks

### 4. Earthquake Relocation (NonLinLoc)
The study area is too large for a single NLL run, so it is split into **6 geographic zones**.

- **`prepare_nll_inputs.py`** — generates one `.obs` file and one `.in` run file per zone, plus GTSRCE station files
- External (run manually in terminal):
  ```
  Vel2Grid run/<runfile.in>
  Grid2Time run/<runfile.in>
  NLLoc run/<runfile.in>
  ```
- **`generate_nll_corrections.py`** — generates second-pass run files by appending per-station delay corrections derived from first-run arrival-time residuals; also exports the locdelay summary via `export_locdelay_info`
- Second pass: same external commands repeated

### 5. Post-relocation Processing
- **`finalize_nll_catalog.py`**:
  1. Cleans NLL output files
  2. Merges the 6 regional results into `RESULT/FINAL.txt`
  3. Rematches relocated events back to `obs/GLOBAL.obs` to recover metadata not present in NLL output (e.g. magnitude)
  4. Saves matched events to `obs/FINAL.obs`
- **`add_temp_picks.py`** (optional, run after): augments `obs/FINAL.obs` with picks from external sources → `obs/FINAL_augmented.obs`

---

## Complementary Analysis

Scripts in `complem_figures/` for visualization and statistics:
- `event_maps.py` — geographic maps of seismicity
- `gutenberg_richter.py` — magnitude-frequency distribution
- `depth_maps.py` — depth distribution
- `error_maps.py` — location uncertainty maps
- `cross_section.py` — vertical cross-sections

> **Environments**:
> - `seisbench_env` → `generate_complem_figures.py` (Gutenberg-Richter, depth maps, error maps)
> - `pygmt_env`     → `generate_complem_maps.py` (event maps for each zone and final catalog)

`zone_Arette/` — focused analysis of the Arette seismic zone.

---

## External Pick Ingestion (temp_picks/)

A self-contained sub-pipeline for ingesting picks from external sources into `obs/FINAL.obs`, producing `obs/FINAL_augmented.obs`. All scripts live in `temp_picks/` and are importable as a package (`from temp_picks.<module> import <function>`). Log files are written to `temp_picks/console_output/`.

The root-level script **`add_temp_picks.py`** orchestrates the full pipeline (steps 1–5 below) in sequence.

| Script | Role |
|--------|------|
| `build_theoretical_tables.py` | Runs Pyrocko's `cake` CLI to compute P/S travel-time envelopes (±5% velocity, 0–100 km) → `temp_picks/tables_Pyr.csv` |
| `merge_omp_picks.py` | Merges all yearly OMP/PhaseNet CSV files from `picks_OMP/` subdirectories → `pick_files/merged_omp.csv`; station `SMC` and year `2026` excluded by default |
| `merge_pyrenees_picks.py` | Concatenates RaspberryShake/PhaseNet `.txt` files from `picks_station_pyrenees/` and `picks_station_pyrenees2/` → `pick_files/merged_pyrenees.txt` and `pick_files/merged_pyrenees2.txt` |
| `convert_picks.py` | Converts external pick files to `.obs` pick line format; maps station names to internal codes via `GLOBAL_code_map.txt`. Formats `TEMP_OBS`, `TEMP_RSB`, and `TEMP_OMP` are supported; new formats are registered in `FORMAT_HANDLERS`. Unresolved stations are logged as an end-of-run summary. |
| `match_picks.py` | Matches converted picks to bulletin events: 60 s time window + residual filter (±0.1 s P, ±0.3 s S, plus ±2.5 s t0-error margin); appends new picks and updates `PhaseCount`; chains against `obs/FINAL.obs` → `obs/FINAL_augmented.obs`; auto-sorts output via `sort_picks`. |
| `sort_picks.py` | Sorts pick lines within each event block by ascending arrival time. |
| `plot_travel_times.py` | QC figure: scatter of observed (distance, travel time) picks over theoretical P/S bands. |

---

## Key Data Formats

### `.obs` (custom seismic bulletin)
- One block per event, separated by blank lines
- Event line starts with `# `: location, magnitude, quality parameters (azimuth gap, RMS, horizontal/vertical uncertainty)
- Following lines: one pick per station (station code, phase P/S, arrival time, uncertainties)

### NLL output
- `.hyp` / `FINAL.txt` — relocated hypocenter parameters
- Does **not** contain magnitude or full pick metadata → rematching to `.obs` is necessary

---

## Git Workflow

- All AI-generated commits go to the **`claude` branch**, never to `main`
- Commit messages must be clean and descriptive so changes are understandable without reading the diff
- Use the format `type: description` — e.g. `fix: ...`, `feat: ...`, `docs: ...`. Never use scoped form `fix(module): ...`
- **Never push automatically to main** — commit and push to the **`claude` branch** without asking, but notice the user
- The user reviews changes locally and decides when to merge or push to `main`
