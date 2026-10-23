---
paths:
  - "final_steps.py"
  - "NLL_run/export_quakeml.py"
  - "complem_figures/seismicity_3d.py"
  - "complem_figures/cross_section.py"
---

# Stage 7 — Final Bulletin Export

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
