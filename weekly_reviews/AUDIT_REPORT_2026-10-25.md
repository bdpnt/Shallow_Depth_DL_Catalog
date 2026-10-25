# Weekly Code Audit — 2026-10-25

Scope: Python files changed on `claude` since 2026-10-18 — one commit, `c0ea4e0` (feat: ingest SANIMS QuakeML picks). `a813d13` changed only docs. The SANIMS source file (`temp_picks/all_picks/SANIMS/`, gitignored) is not in this checkout, so the findings below come from reading the code, not from running it on the data.

## Files reviewed
- add_temp_picks.py
- temp_picks/convert_picks.py

## Findings

### temp_picks/convert_picks.py
- **[Medium]** Line 675 (`_convert_quakeml`) / line 569 (`convert_temp_qml`): the QuakeML `phaseHint` is passed through unchanged. Every downstream consumer assumes exactly `'P'` or `'S'`. Line 573 gives any value other than `'P'` the S error (0.15 s). `match_picks.interpolate_travel_time` (line 216) checks a `Pg`/`Pn` pick against the **S** travel-time band with the S tolerance. The `(station, phase)` dedup key `'Pg'` never matches an existing `'P'`, so the slot is filled twice. SeisComP bulletins commonly write `Pg`/`Pn`/`Sg`/`Sn`, and NLLoc's `LOCPHASEID` only knows the upper-case `PG`/`PN`. Fix: normalise to the first letter (`hint[:1].upper()`) and drop anything that isn't P/S. Also log a `Counter` of the raw hints once, so the SANIMS phase vocabulary gets checked against the 13 466 picks already added.
- **[Medium]** Line 675: if a pick has no `phaseHint` (it is optional in QuakeML 1.2), `phase_type` is `None`. `_format_pick_line` then raises `AttributeError` on `phase.ljust(6)` (line 273), which stops the whole `add_temp_picks.py` run at step 5. The fix is the same as above: skip the pick and count it as dropped, or fall back to the `phase` of the preferred origin's arrival for that `pickID`, which is already read in the loop above.
- **[Low]** Line 672: `wid` is `None` for a pick without `waveformID`, so `wid.get(...)` crashes in the same way. Add a guard that counts the pick as dropped.
- **[Low]** Line 685: `event.clear()` without removing processed siblings from the parent still leaves one empty `<event>` element per event attached to the tree. So the "flat memory" claim in the docstring is not strictly true. This doesn't matter at 1 470 events, but it would for a large bulletin. The usual idiom is `while event.getprevious() is not None: del event.getparent()[0]`.

### add_temp_picks.py
No issues. The precedence order, with SANIMS first, matches the comment and `.claude/rules/temp_picks.md`.

## Summary
2 files reviewed, 4 issues found (2 Medium, 2 Low).
