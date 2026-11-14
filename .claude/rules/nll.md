---
paths:
  - "run_NLL.py"
  - "NLL_run/**"
  - "run/nll/**"
---

# Stages 4–5 — NonLinLoc relocation and post-processing

### 4. Earthquake Relocation (NonLinLoc)
The study area is too large for a single NLL run, so it is split into **6 geographic zones**, processed with up to **3 zones running concurrently** (zones are independent; NLL runs are the bottleneck).

- **`run_NLL.py`** — for each zone: generates one `.obs` file and one `.in` run file `run/nll/run_<N>_DELAYS.in` (plus GTSRCE station file), runs Vel2Grid → Grid2Time → NLLoc via `NLL_run/run_zone.py`, cleans up `.hdr` files, then generates a second-pass run file `run/nll/run_<N>_NLL.in` by appending per-station delay corrections derived from first-run arrival-time residuals and reruns NLLoc via `run_zone.py` with `--corrections-pass` (grids already built), cleaning up `.hdr` files again right after — this per-zone `.hdr` cleanup keeps at most a few zones' worth on disk at once instead of waiting for all 6 zones to finish. Once all zones are done, exports the locdelay summary via `export_locdelay_info` to `run/nll/locdelays/`. NLL working folders live inside `run/`: `nll_model/` (velocity grids), `nll_time/` (travel-time grids), `nll_loc/` (per-zone NLLoc output).
- **`NLL_run/run_zone.py`** — runs Vel2Grid → Grid2Time → NLLoc in sequence for a given `.in` file; `--corrections-pass` skips Vel2Grid/Grid2Time; accepts a `zone_label` to prefix log output when zones run concurrently


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
- **`add_temp_picks.py`** (**required**, run after — `run_SSST.py` reads its output as its only input): augments `obs/NLL_result.obs` with picks from external sources → `obs/NLL_result_augmented.obs` (same events, **1 059 360 picks** on the 2026-10-23 run, 51 568 events from 727 742 picks: +331 618, +45.6 %; SANIMS alone +13 466)
