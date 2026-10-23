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

Seven stages, each with an orchestrator script. Per-stage detail (parameters, measured counts, rationale) lives in path-scoped rules under `.claude/rules/`, loaded automatically when a matching file is read or edited — **read the relevant rules file before answering a question about a stage without opening its code**.

| # | Stage | Script | Output | Rules file |
|---|---|---|---|---|
| 1 | Station inventory fusion | `build_global_inventory.py` (`fetch_inventory/`) | `stations/GLOBAL_inventory.xml`, `GLOBAL_code_map.txt` | `inventory.md` |
| 2 | Catalog fetching & conversion | `fetch_all_bulletins.py` (`fetch_obs/`) | six `obs/*.obs` (75 926 events) | `fetch_harmonize.md` |
| 3 | Catalog harmonization | `build_global_bulletin.py` (`global_obs/`) | `obs/GLOBAL.obs` (57 856 events) | `fetch_harmonize.md` |
| 4–5 | NLL relocation (6 zones) + rematch | `run_NLL.py` (`NLL_run/`) | `RESULT/NLL_result.csv`, `obs/NLL_result.obs` (50 918) | `nll.md`, `results_format.md` |
| 5b | External pick ingestion (**required**) | `add_temp_picks.py` (`temp_picks/`) | `obs/NLL_result_augmented.obs` | `temp_picks.md` |
| 6 | SSST relocation + PDF metrics | `run_SSST.py` | `RESULT/SSST_result.csv`, `obs/SSST_result.obs` (46 224) | `ssst.md`, `results_format.md` |
| 7 | QuakeML export | `final_steps.py` | `RESULT/FINAL*.xml` | `export.md` |

Figures and diagnostics: `complem_figures/` → `complem_figures.md` (and `ssst_corrections.md`).

### Cross-cutting conventions
- `PUBLIC_ID PYRENEES_%06d`, assigned in stage 3 in chronological order, is the join key of every later stage.
- **Published hypocentres are the location-PDF expectation, not the maximum-likelihood point** (`LOCHYPOUT ... SAVE_NLLOC_EXPECTATION`): `latitude`/`longitude`/`depth` in `RESULT/*.csv` equal `expect_*`, and RMS/Gap/Dist/origin time are re-evaluated there. The ML point is kept as `maxlike_*`. Every `complem_figures` reader defaults to the expectation; `--solution maxlike` switches.
- Pick uncertainty is assigned, not measured: **0.05 s P / 0.15 s S** everywhere (fetch, `LOCQUAL2ERR`, added picks).
- `pyr:usable` = `(C68 − 0.6827)/C68_sigma_n ≥ −2` **and** `not dip_reject` **and** metrics present. Ψ never rejects.
- `station_priority()` (stage 1) is the single rule for picking a co-located group's representative everywhere.
- The 1-D velocity model `Pyrenees_1D` exists twice (NLL control file and `temp_picks/models/*.nd`); nothing checks they agree.

### ⚠ Hazards — read before running anything
- ⚠ Stage 1 — **This stage is interactive and cannot run unattended.** `merge_station_inventories.check_inventory()` prompts at the terminal for every station code appearing in more than one network, and the operator's answers are **not logged anywhere** — `GLOBAL_inventory.xml` is not reproducible from code + inputs alone
- ⚠ Stage 3 — **Step 5 is interactive** — `find_and_merge_doubles` prompts per group, once per source catalog (6 times). Auto-merge only for pairs within 0.15 s / 10 km. Unlike §1 these choices *are* logged (kept/dropped bulletin IDs)
- ⚠ **`build_global_bulletin.py` must not be re-run on already-processed bulletins.** `remap_picks_to_unified_codes` rewrites `obs/*.obs` **in place** and is not idempotent: on a second run it parses the already-substituted code (`FR.0041` → station name `0041`), finds no match, and **deletes the pick**. Remap and fusion take the six source bulletins as an explicit list (`_source_catalogs` in `build_global_bulletin.py`, whose order is the fusion order), so later-stage outputs left in `obs/` are never touched — but the six themselves must come fresh from stage 2
- ⚠ **`add_temp_picks.py` is not re-runnable** against an existing `obs/NLL_result_augmented.obs` (each source matches against the previous output), and needs `temp_picks/pick_files/Pamplona_PN_picks.txt` from `../PhaseNet_Pamplona` first.
- ⚠ Stage 1's elevation fill is the only step that reaches the network (~3 min); `--no-fill-elevations` skips it.

> ⚠️ `python run_NLL.py` starts the full pipeline immediately: it rewrites `obs/GLOBAL.obs` in place, regenerates the zone `.obs`/`.in` files, and begins building travel-time grids. Its only option is `--decaff`; `--help` is safe.

**Sleep.** `run_NLL.py` and `run_SSST.py` keep the Mac awake by default (`NLL_run/keep_awake.py`: `caffeinate -i -w <pid>`, released by itself when the process ends, crash or kill included); `--decaff` turns it off. `-i` blocks idle sleep only — closing the lid still sleeps a laptop, which freezes the run (no timeouts anywhere, so it resumes on wake) rather than killing it.

> **Environments**: `seisbench_env` is the project default — the whole pipeline (`build_global_inventory.py` → `run_SSST.py`) runs in it unprefixed. `pygmt_env` is the only exception, for the modules importing PyGMT or `xarray`.
> - `seisbench_env` → `generate_complem_figures.py` (Gutenberg-Richter, depth maps, error maps)
> - `pygmt_env`     → `generate_complem_maps.py` (event maps for each zone and final catalog), `cross_section.py`, `cross_section_panels.py`, `profile_location_map.py`, `aoi_map.py`, `source_depth_maps.py`, and the two `build_global_bulletin.py` steps it launches itself via `conda run` (`filter_events_by_aoi.py`, `plot_global_catalog_map.py`)

---

## Key Data Formats

### `.obs` (custom seismic bulletin)
- One block per event, separated by blank lines
- Event line starts with `# `: location, magnitude, quality parameters (azimuth gap, RMS, horizontal/vertical uncertainty)
- Following lines: one pick per station (station code, phase P/S, arrival time, uncertainties)

Result CSV columns: see `.claude/rules/results_format.md`.

---

## Git Workflow

- All AI-generated commits go to the **`claude` branch**, never to `main`
- Commit messages must be clean and descriptive so changes are understandable without reading the diff
- Use the format `type: description` — e.g. `fix: ...`, `feat: ...`, `docs: ...`. Never use scoped form `fix(module): ...`
- **Never push automatically to main** — commit and push to the **`claude` branch** without asking, but notice the user
- The user reviews changes locally and decides when to merge or push to `main`

### Pull Request Format
- PR body: **Summary section only** (bullet points of what changed and why). No "Test plan" section.
