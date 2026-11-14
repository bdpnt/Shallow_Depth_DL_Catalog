---
paths:
  - "fetch_all_bulletins.py"
  - "fetch_obs/**"
  - "build_global_bulletin.py"
  - "global_obs/**"
  - "complem_figures/magnitude_regression.py"
  - "complem_figures/aoi_map.py"
  - "complem_figures/source_depth_maps.py"
---

# Stages 2–3 — Catalog Fetching & Harmonization

### 2. Catalog Fetching & Conversion
- Sources: RESIF (FDSN), ICGC, IGN, LDG, OMP (in `org_catalogs/`)
- Script: `fetch_all_bulletins.py` → calls modules in `fetch_obs/`
- Output: **six** `.obs` files in `obs/` — OMP contributes two, `OMP_78-19.obs` and `OMP_2016.obs` (a re-picked 2016: 2 178 events against 1 568 for the same year in the archive file; both are fused, the overlap resolved by the fusion's loose+pick-validated matching)
- 75 926 events / 1 784 332 picks total. RESIF/ICGC/IGN/LDG cover **2020–2025 only**; the whole 1978–2019 record comes from OMP, so the catalog is not homogeneous in source composition across time
- Only **manual** P/S picks are kept (RESIF `evaluation_mode`, IGN/ICGC GSE2 `m__`). An event with no magnitude of the expected type is dropped entirely, picks included
- LDG's export dates **93 picks (77 events) one day after their event** with the right hour/minute/second (LDG's own residual stays small); `fetch_obs/LDG.py` moves any pick lying 1 day ± 300 s after its origin back one day and logs the count
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
- `_remove_magnitudes_under_1` is **commented out** (`fuse_bulletins.py:979`) although `fuse_bulletins`'s docstring still advertises it; 31.3 % of `GLOBAL.obs` is below ML 1.0
- `global_obs/add_temporary_picks.py` and `list_magnitude_types.py` are imported by nothing — out of the pipeline
