# Weekly Code Audit — 2026-10-11

Review window: the 12 commits on `claude` since the 2026-10-04 audit
(`f4cd784` … `cc6b7be`).

## Files reviewed
- NLL_run/export_quakeml.py
- NLL_run/merge_regional_results.py
- add_temp_picks.py
- build_global_bulletin.py
- complem_figures/seismicity_3d.py
- temp_picks/convert_picks.py
- temp_picks/match_picks.py

## Findings

### temp_picks/match_picks.py
- **[Low]** Line 433: the new unpaired-`*`-pick drop loop (step 6) decrements
  `n_added` for every dropped pick but never decrements `n_added_near`.
  `n_added_near` is incremented at line 415 for any added pick whose matched
  station is within 20 km of the event; a dropped pick that was within 20 km
  therefore leaves `n_added_near` overcounted. The stat is both logged
  ("Added - station within 20 km of event", line 452) and returned in the
  result dict (`n_added_near`, line 471), so both can exceed the true count.
  This affects reporting only — the output `.obs` picks are correct.
  Suggested fix: when dropping a pick in step 6, also decrement `n_added_near`
  if that pick's `dist_km_matched` was ≤ 20 km (track the distance alongside
  the line in `added_relative`), or recompute `n_added_near` from the surviving
  picks after the drop loop.

### complem_figures/seismicity_3d.py
- **[Low]** Line 120: `out['usable'] = [classify(row)[0] for _, row in df.iterrows()]`
  runs a per-row `df.iterrows()` over the full SSST result (~46 k rows).
  `iterrows()` materializes a Series per row and is the slow path; here it adds
  a few seconds every time the page is built. `classify()` reads only
  `n_scat`, `C68`, `C68_sigma_n`, and `dip_reject`, all of which are already
  present as columns, so the `usable` flag can be computed vectorized
  (`c68_z = (df['C68'] - _C68_NOMINAL) / df['C68_sigma_n]`; reject where
  `n_scat` is NaN, `c68_z < _C68_Z_MIN`, or `dip_reject`). Non-critical for a
  standalone figure script; noting it because the per-row reuse of a scalar
  classifier is an avoidable cost on the catalog-sized input.

## Notes (reviewed, no issue)
- `export_quakeml.py` (line 731): the 4→5-field epoch unpack fix is correct and
  complete. `load_station_epochs` now yields 5-tuples
  `(net, sta, start, end, elevation)`; every other consumer (`resolve_station`,
  the sort key at line 296, the alternates builder at line 491) indexes fields
  rather than destructuring, so no other unpack site was left stale.
- `merge_regional_results.py`: dropping no-hypocentre rows (NaN lat/lon/depth)
  and `pdfVolume <= 0` rows before the smallest-`pdfVolume` dedup is sound —
  both would otherwise win the `idxmin` dedup spuriously. `idxmin` skips NaN
  `pdfVolume`, and the pre-filter removes the degenerate zero case, so the
  dedup now selects a genuine solution.
- `convert_picks.py`: the new `instrument`/`relative_timing` plumbing is
  backward compatible — the added parameters default to `'?'` / `False`
  (prior behaviour), `ROW_BASED_FORMATS` contains only `TEMP_STB` whose handler
  accepts `instrument`, and the pick-line field indices used downstream in
  `match_picks` (station `parts[0]`, instrument `parts[1]`, phase `parts[4]`)
  match `_format_pick_line`'s output.
- `match_picks.py` step 6 pairing logic: verified it never keeps a genuinely
  unpaired `*` pick nor drops a valid S-P pair — `relative_phases` is built
  from all `*` picks in the event (pre-existing + added), the subset test
  `{'P','S'} <= relative_phases[station]` is monotonic, and no `KeyError` is
  possible since each added pick's own station/phase is always present.
- `add_temp_picks.py` / `build_global_bulletin.py`: wiring and `os.makedirs`
  guards (`MAPS`, `STATS`) are correct; `os` is already imported.

## Summary
7 files reviewed, 2 issues found (0 High, 0 Medium, 2 Low). Both are
non-blocking: one diagnostic-stat overcount in `match_picks.py` and one
avoidable per-row loop in the `seismicity_3d.py` figure builder. No
correctness defects in the catalog-producing paths this week.
