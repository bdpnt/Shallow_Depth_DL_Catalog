# Weekly Code Audit — 2026-10-18

Scope: Python files touched on the `claude` branch in the last 7 days
(commit range `cc6b7be..HEAD`, 2026-10-09 → 2026-10-16). Reviewed the actual
diffs plus surrounding context; all 20 files compile.

## Files reviewed
- NLL_run/generate_ssst_runfiles.py
- NLL_run/keep_awake.py
- NLL_run/reformate_obs.py
- NLL_run/run_ssst.py
- NLL_run/run_zone.py
- add_temp_picks.py
- build_global_bulletin.py
- complem_figures/depth_histogram.py
- complem_figures/error_maps.py
- complem_figures/plot_pdf_cloud_2d.py
- complem_figures/ssst_corrections.py
- complem_figures/zone_map.py
- fetch_obs/LDG.py
- fetch_obs/OMP.py
- global_obs/fuse_bulletins.py
- global_obs/remap_picks_to_unified_codes.py
- run_NLL.py
- run_SSST.py
- temp_picks/convert_picks.py
- temp_picks/match_picks.py

## Findings

### NLL_run/run_zone.py / NLL_run/run_ssst.py
- **[Low]** `run_zone.py:29` and `run_ssst.py:55`: `_NLL_BIN` now points to a
  **temporary** locally-patched build
  (`/Users/bdupont/Desktop/Codes/NonLinLoc_fixed/src/bin`). It is clearly
  self-documented as TEMPORARY with a revert note, so this is a tracking item,
  not a defect: it must be reverted to `NonLinLoc/src/bin` before the branch
  merges to `main`, otherwise every machine without that side-build (CI, other
  operators) runs the wrong binary or fails outright. As with the previous
  value, the path is also machine-specific absolute — a longer-term `NLL_BIN`
  env-var / config override would remove both the revert risk and the
  hardcoded home directory. No change required this week; flagged so it is not
  forgotten at merge time.

### build_global_bulletin.py
- **[Low]** Lines 46-54: the new explicit `_source_catalogs` list (correctly
  replacing the non-deterministic `glob.glob("obs/*.obs")` for remap and
  fusion — a genuine improvement that also protects later-stage outputs in
  `obs/`) orders the sources `RESIF, ICGC, IGN, LDG, OMP_2016, OMP_78-19`. The
  inline comment justifies this as "ICGC precedes IGN as before," but the
  prior fusion order came from an **unsorted** `glob.glob`, so prior
  precedence between ICGC and IGN was filesystem-dependent, not a guaranteed
  "as before." Because a merged event keeps the *first* location, this order
  decides which hypocenter survives for ICGC/IGN duplicate events that RESIF
  does not cover. The order is reasonable (it matches the source order in
  CLAUDE.md §2) — suggestion is only to reword the comment so it does not
  imply byte-for-byte continuity with previous runs, and to confirm
  ICGC-over-IGN is the intended precedence.

### complem_figures/ssst_corrections.py
- **[Low]** Lines 1181-1185: the new `--extent` validation runs
  `[float(v) for v in args.extent.split(',')]` *before* checking arity/order,
  so a malformed non-numeric or wrong-length value (e.g. `--extent=a,b,c,d` or
  `--extent=1,2,3`) raises a bare `ValueError`/`IndexError` traceback instead
  of the intended clean `p.error(...)` message. Suggested fix: wrap the parse
  in `try/except (ValueError)` and route failures (including wrong count)
  through `p.error`, so all malformed-extent cases produce the same usage
  error.

## Notes (reviewed, no issue)
- `fetch_obs/LDG.py` day-shift fix: `yr/mo/dy` are set for every pick and only
  overridden on the shift branch (no undefined-variable path); the
  one-day-±300 s detection via datetime arithmetic is correct.
- `fetch_obs/OMP.py` `minute_p/minute_s -= 1`: both are defined before the
  decrement; the fix correctly rolls the *pick's own* minute back on negative
  seconds, and the result feeds a `timedelta` so negative minutes normalize.
- `temp_picks/match_picks.py`: the `n_added_near` increment (line 414) and the
  new decrement (line 434) use the same 20.0 km threshold — symmetric and
  correct; no KeyError risk in the `'*'`-pick drop loop.
- `temp_picks/convert_picks.py`: day-granularity epoch comparison, network
  preference, the HYPO71 regex (reduced `P X ` and full ` EPUX ` phase
  blocks), the quality-4 early return (now dropping the S pick too), and
  card-minute-based S seconds are all correct. TEMP_OTH correctly omits the
  `network=` argument (HYPO71 cards carry no network).
- `global_obs/{fuse_bulletins,remap_picks_to_unified_codes}.py`: the
  `folder_path`/`folder_bulletin` → `files` rename is consistent across the
  dataclass, CLI, and `build_global_bulletin.py` call sites; no stale
  references remain in these modules.
- `run_ssst.py`/`run_SSST.py` rolling grid retirement and end-of-zone cleanup,
  `keep_awake.py`, and the `depth_histogram.py` / `error_maps.py` /
  `zone_map.py` empty-input guards are correct.

## Summary
20 files reviewed, 3 issues found (all Low severity). No High or Medium issues.
This week's changes are high quality: the pick-date corrections (LDG, OMP,
HYPO71 TEMP_OTH), the glob→explicit-files conversion, and the new edge-case
guards are all sound.
