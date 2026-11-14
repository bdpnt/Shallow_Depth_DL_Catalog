---
paths:
  - "RESULT/**"
  - "NLL_run/**"
  - "run_NLL.py"
  - "run_SSST.py"
  - "final_steps.py"
  - "complem_figures/**"
---

# Result CSV formats (NLL output, location-PDF quality columns)

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
