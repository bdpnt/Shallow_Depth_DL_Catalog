---
paths:
  - "run_SSST.py"
  - "NLL_run/run_ssst.py"
  - "NLL_run/generate_ssst_runfiles.py"
  - "NLL_run/reformate_obs.py"
  - "NLL_run/pdf_metrics.py"
  - "NLL_run/keep_awake.py"
  - "run/ssst/**"
  - "complem_figures/ssst_evolution.py"
  - "complem_figures/plot_pdf_cloud*.py"
---

# Stage 6 — SSST Relocation

### 6. SSST Relocation
Final stage, run after `add_temp_picks.py`: relocates `obs/NLL_result_augmented.obs` with iterative Source-Specific Station Terms (NonLinLoc `Loc2ssst`). Zones run **strictly sequentially** (memory-bound), everything inside a zone in parallel (`NLLOC_CORES` NLLoc chunks / `LOC2SSST_CORES` Loc2ssst instances).

- **`run_SSST.py`** — orchestrator; campaign configuration (RUN_NAME, CHAR_DISTS, VPVS, core counts, LSPHSTAT) at the top of the file; CLI `--zones`, `--iteration-start`, `--iteration-stop` (partial campaigns/resume), `--keep-grid` (keep each finished zone's final grid set), `--decaff` (allow sleep, §4). Per zone: cuts `obs/GLOBAL_<N>_SSST.obs` + `stations/GTSRCE_SSST_<N>.txt`, derives `run/ssst/run_<N>_NLL.in` (NLLoc) and `run/ssst/run_<N>_SSST.in` (Loc2ssst) from `run/nll/run_<N>_DELAYS.in` (`NLL_run/generate_ssst_runfiles.py`), splits the bulletin into per-event files in `obs/nlloc_obs/GLOBAL_<N>/` (`NLL_run/reformate_obs.py`, which first deletes any `.nlloc_obs` left by an earlier run — NLLoc locates every file in the folder), builds P+S grids in `run/ssst_model|ssst_time` (VpVs −9.99, real S grids), and runs the iteration loop (`NLL_run/run_ssst.py`: len(CHAR_DISTS) SSST iterations + final NLLoc-only relocation, outputs under `run/ssst_loc/<RUN_NAME>/`).
- After all zones: merges the final-iteration CSVs → `RESULT/SSST_result.csv` (dedup by lowest `pdfVolume`), rematches against `obs/NLL_result_augmented.obs` via `publicId` → `obs/SSST_result.obs` (same modules as the NLL stage). 50 918 → **46 224 events**; the 4 694 lost are mostly the `LSPHSTAT NRdgsMin = 6` phase minimum, against NLL's 4.
- `LSPHSTAT [0.15, 6, 200.0, 0.3, 0.5, 10.0]` (RMS ≤ 0.15 s, ≥ 6 phases, gap ≤ 200°, P/S residual ≤ 0.3/0.5 s, Len3 ≤ 10 km) selects which events *teach* the corrections — deliberately stricter than what *receives* them. `NRdgsMin` doubles as the NLLoc min-phases threshold, read back from the generated control file so the two cannot diverge.
- The SSST control files are **derived from `run/nll/run_<N>_DELAYS.in`** (the first-pass file, so the NLL static `LOCDELAY`s are deliberately not inherited — iteration 0 at L = 9999 recomputes that same quantity). `LSGRID`/`LSOUTGRID` are re-gridded to a coarse **1 km** spacing: Loc2ssst holds two full buffers in RAM per instance.
- **Disk, not RAM, is what the iterations accumulate** (every NLLoc/Loc2ssst is a subprocess, its RAM is returned on exit). Each iteration writes a full `ssst_corr<i>/` grid set, so `run_ssst.py` retires the set an iteration read (`run/ssst_time/Pyrenees_<N>/` at 0, `ssst_corr<i-1>/` after) as soon as its Loc2ssst succeeded, moving the uncorrected stations' grids forward instead of copying them: **at most two grid sets on disk**, against initial + 5 before. Once a zone has run its final relocation that last set is deleted too (30–70 GB per zone, read by nothing downstream) unless `--keep-grid` is given; a partial campaign keeps it as the resume input. `loc_ssst_corr<i>/` (per-iteration locations, read by `ssst_evolution.py`, `plot_pdf_cloud.py`, `ssst_corrections.py`) are kept.
- **NLLoc's RAM scales with an event's arrival count**: for 2-D grids it holds one full grid per arrival while locating an event. The SSST grids are therefore re-gridded at **0.1 km** (`GenSSSTRunParams.timeGridSpacing`, same 450 × 40 km extent) against the NLL stage's 0.05: ~29 → ~7.4 MB per arrival. At 0.05 the 2017 Maupasacq events (up to 726 picks, 205 events > 200) needed ~21 GB each, and round-robin chunks hit them together since they cluster in time — this is what exhausted 48 GB. Measured at 0.1: 5.4 GB / 89 s for the 726-arrival event, and a 101-arrival event moves 6/18/40 m (x/y/z), 6 ms. `build_grids()` rebuilds grids whose header does not match the control file's `VGGRID`, so grids left at another spacing are never silently reused.
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
