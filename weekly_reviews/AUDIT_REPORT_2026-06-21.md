# Weekly Code Audit — 2026-06-21

Window: Python files modified on the `claude` branch in the last 7 days
(commits from 2026-06-15 through 2026-06-18).

## Files reviewed
- NLL_run/filter_distant_picks.py
- NLL_run/generate_regional_runfiles.py
- NLL_run/match_pre_post_relocation.py
- NLL_run/merge_regional_results.py
- NLL_run/parse_nll_output.py
- complem_figures/depth_maps.py
- complem_figures/error_maps.py
- complem_figures/event_maps.py
- complem_figures/zone_map.py
- finalize_nll_catalog.py
- generate_complem_figures.py
- generate_complem_maps.py
- generate_nll_corrections.py
- global_obs/fuse_bulletins.py
- global_obs/remap_picks_to_unified_codes.py
- temp_picks/convert_picks.py
- temp_picks/match_picks.py
- temp_picks/merge_omp_picks.py
- temp_picks/merge_pyrenees_picks.py
- temp_picks/sort_picks.py

## Findings

### global_obs/fuse_bulletins.py
- **[High]** Line 603: `found_possible.append(possible_row.index[0])` always records the **first** candidate's index, but the loop (`for _, row in possible_row.iterrows()`, line 577) can match a *later* candidate after `continue`-ing past earlier ones. `possible_row` is a boolean-filtered view of `possible_match` (line 550) and retains its original index, so when the match is not the first candidate, the wrong row is promoted into `validated` (line 629) and the wrong row is dropped from `possible_match` (line 632) — corrupting the loose-match accounting (the fused `new_lines` themselves use the correct `event_idx2`, so output lines are unaffected). Fix: iterate `for orig_idx, row in possible_row.iterrows():` and `found_possible.append(orig_idx)`.
- **[Medium]** Lines 274-281 (`find_pick_lines`): the `while True` loop advances `curr += 1` and indexes `all_lines[curr]`, breaking only on a line starting with `'\n'`. If the target event is the last block and the file has no trailing blank line, this runs off the end with `IndexError`. Add a `curr < len(all_lines)` bound (as `_remove_duplicate_picks` at line 715 already does).
- **[Medium]** Lines 513-516 (`_add_phases_to_lines`): same unbounded-loop hazard — `while not old_lines[curr].startswith('\n')` with `curr += 1` will `IndexError` on the final block when there is no trailing blank line. Add a `curr < len(old_lines)` guard.
- **[Medium]** Lines 686-687 (`_replace_mean_magnitudes`): `mean(mags)` is computed after filtering out `'Nan'` entries. If every source magnitude for an event is `'Nan'`, `mags` is empty and `numpy.mean([])` returns `nan` with only a RuntimeWarning, silently writing a `nan` magnitude into the catalog. Guard the empty case explicitly (skip / sentinel).
- **[Low]** Line 37: `from numpy import mean` is redundant with `import numpy as np` (line 40) and only used once (line 687). Remove it and use `np.mean`.
- **[Low]** Line 18: `import sys` is never used. Dead import.

### NLL_run/generate_regional_runfiles.py
- **[Medium]** Lines 105-118 (`_build_alternate_code_map`): `network_code`/`station_code` are assigned only inside the inner `while` when a `  Station` line is found. If an `Alternate` block contains no `  Station` line before the blank-line break, line 118 reads stale values from a previous iteration — or raises `NameError` on the first block. Initialize both to `None` before the inner loop and `continue` if still `None`.
- **[Medium]** Lines 137-138 (`_find_station_info`): `alternate_code_map.get(alternate_code)` returns `None` for an absent code, so `codes[0]` raises `TypeError`; `inventory.select(...).networks[0]` raises `IndexError` on no match. Since every bulletin station code is fed here, one unknown station aborts all GTSRCE generation. Guard the missing-code / empty-select case and log+skip.
- **[Low]** Lines 156-158 (`_gen_child_obs`): `f.writelines(lines[:4])` unconditionally copies a fixed 4-line header with no validation; if the source bulletin's header length differs (or the file has <4 lines) the child `.obs` header is silently wrong. Validate or derive the header span, and add an explanatory comment.

### complem_figures/depth_maps.py
- **[Medium]** Line 84 (`_filter_dates`): `period_end = period_start + 4` hardcodes a 5-year window while the loop steps by `time_range` (line 83). Correct only for the default `time_range=5`; any other `--time-range` produces overlapping windows (e.g. 3 → 1976-1980, 1979-1983) or gaps (e.g. 10 → 1976-1980, 1986-1990) silently. Fix: `period_end = period_start + time_range - 1`.
- **[Medium]** Lines 118-131: the windowed-median grid is built with a pure-Python double loop over `bins_lat * bins_lon` (≈344,000) cells, each masking the full events DataFrame — O(cells × events), the dominant cost, run once per period. Replace with `scipy.stats.binned_statistic_2d` or a binned/indexed approach.
- **[Medium]** Line 63 (`_read_file`): no guard for a missing or empty `FINAL.csv`. A missing file raises an unguarded `FileNotFoundError`; an empty DataFrame propagates to `np.nanmean` over an all-NaN array (RuntimeWarning + NaN annotations). Add an existence/non-empty check that fails clearly.

### complem_figures/error_maps.py
- **[Medium]** Lines 84-91 (`_filter_dates`): identical hardcoded-window bug as depth_maps.py — `period_end = period_start + 4` ignores `time_range`. Fix: `period_end = period_start + time_range - 1`.
- **[Medium]** Lines 119-132: identical O(cells × events) pure-Python double-loop performance hotspot as depth_maps.py.
- **[Medium]** Line 63 (`_read_file`): identical missing-file / empty-DataFrame boundary gap as depth_maps.py.

### complem_figures/event_maps.py
- **[Medium]** Lines 112-122 (`.obs` path): the drop/rename pipeline assumes every `# ` event line splits into ≥17 whitespace fields. A short/malformed line yields a DataFrame missing columns 13-16, so `.rename`/`.astype(float)`/`_remove_high_err` fail with an opaque `KeyError`. Validate field count at this boundary.
- **[Medium]** Lines 142-143 (`.csv` path): `_remove_high_err` (default on) accesses `df.erh`/`df.erv`/`df.gap`/`df.rms`, which exist only after a rename that assumes `FINAL.csv` contains `true_erh`/`true_erz`/`Gap`/`RMS`. A partial/pre-relocation CSV lacking any of these raises `AttributeError` rather than a clear message.
- **[Low]** Line 139: `if parameters.source_filter and 'source' in events_df.columns:` silently plots **all** events when `source_filter` is set but the `source` column is missing — a per-zone map would then silently render every event. Log a warning when the column is absent but a filter was requested.

### NLL_run/merge_regional_results.py
- **[Medium]** Line 121: `pd.read_csv(path)` has no guard for a missing zone CSV (`FileNotFoundError`) or an empty one (`EmptyDataError`, e.g. a zone with zero relocations), aborting the merge. Catch/skip and warn; also handle `frames` being empty before `pd.concat`.
- **[Medium]** Lines 130-147: the duplicate-resolution **log** sorts each group by `pdfVolume` and reports `iloc[0]` as "kept", while the row actually kept is chosen via `groupby(...).idxmin()`. On `pdfVolume` ties these two paths can disagree, so the log can name a different source than what was written. Drive both log and selection from the same `idxmin` result.
- **[Low]** Lines 79-80 (`_build_covariance`): `v1 / len1` and `v2 / len2` divide by ellipsoid axis lengths with no guard against a zero-length (degenerate) axis, which would silently yield NaN in `true_erh`/`true_erz`. Add a guard or a comment that NLL never emits zero-length axes.

### NLL_run/filter_distant_picks.py
- **[Medium]** Lines 134, 144-145: a `# ` header with fewer than 9 whitespace tokens makes `float(parts[7])`/`parts[8]` raise `IndexError`/`ValueError`, aborting mid-run during the in-place overwrite; a missing `parameters.fileBulletin` raises an unguarded `FileNotFoundError`. Validate `len(parts) >= 9` (or try/except + skip) and confirm the file exists before rewriting it in place.
- **[Low]** Lines 144-145: `event_lat`/`event_lon` are never reset between event blocks. Today a pre-header pick is correctly skipped, but if header parsing is later hardened to "skip malformed", subsequent picks would be silently measured against the *previous* event's coordinates. Reset to `None` on block boundary / failed parse, with a comment on the stateful carry-over.

### global_obs/remap_picks_to_unified_codes.py
- **[Medium]** Lines 86-97 (`find_unique_stations`): `pd.concat([unique_sta, pd.DataFrame([new_row])])` is called once per station inside the loop — quadratic (each concat copies the growing frame). Accumulate rows in a list and build one `pd.DataFrame(rows)` after the loop.
- **[Low]** Lines 138-153: when exactly one candidate matches, `alternate_code` is set (line 139) before the date-range check, and the date filter only narrows when `len(working) == 1` — it never rejects the lone candidate. An out-of-service-window pick is thus silently remapped. Add a comment if intended, or reject when the date is out of range.
- **[Low]** Line 18: `import sys` is never used. Dead import.

### complem_figures/zone_map.py
- **[Medium]** Line 102: `_ZONES = _load_zones(...)` runs at **import time**; if `loc/` is missing or has no parseable `GLOBAL_*/last.in`, `_ZONES` is silently `{}`, the module still imports, and `generate_figure` draws zero zones with no warning (failure surfaces far from cause). Defer the load into `generate_figure` or raise on empty.
- **[Low]** Line 168: `zip(_ZONES.items(), _ZONE_COLORS)` truncates to the shorter sequence. `_ZONE_COLORS` has exactly 6 entries, so a 7th discovered zone (`GLOBAL_7`) would be silently dropped. Document the 6-zone assumption or extend the palette.

### temp_picks/match_picks.py
- **[Medium]** Line 148 (`_update_phase_count`): `parts[12]` is hardcoded as the PhaseCount field with no length check; a header with <13 tokens raises `IndexError` and aborts the run. Validate length (log+skip or raise clearly) and add a comment tying the index to the documented header layout.
- **[Low]** Lines 108-117 / 419-425: only the *last* `PUBLIC_ID` line per block is retained on read, and any non-pick metadata line with <5 tokens is kept in `picks` but excluded from the `pick_keys` dedup set (line 117), so it can never be deduplicated. Low practical impact; add a guard/comment that only true pick lines are expected here.

### temp_picks/merge_pyrenees_picks.py
- **[Low]** Lines 25-26: `_PROJECT_ROOT` is computed but never referenced. Dead code — remove.
- **[Low]** Lines 66-78 (`_merge_directory`): files are concatenated without ensuring each source ends in a newline. An input `.txt` lacking a trailing newline joins its last line onto the next file's first line, producing a malformed pick line downstream. Normalize newlines between files.

### finalize_nll_catalog.py
- **[Low]** Lines 46-50: the six per-zone CSV paths are passed to `merge_bulletins` with no existence check; a zone that produced no `*.loc.csv` triggers `FileNotFoundError` and aborts the whole finalize. Filter `csv_files` to existing paths and log the missing zones.

### temp_picks/convert_picks.py
- **[Low]** Lines 73, 452: neither the code-map load nor the input file open guards against a missing/empty file — `FileNotFoundError` propagates uncaught. Since `convert_file` is also imported as a library API, a clearer boundary error would help.

### temp_picks/merge_omp_picks.py
- **[Low]** Lines 121-124: `os.listdir(input_dir)` raises a raw `FileNotFoundError` if `input_dir` does not exist, with no friendly message. Add an existence check that logs clearly. (Output is not at risk — listdir precedes the write-open.)

### NLL_run/parse_nll_output.py
- **[Low]** Lines 88-128 (`_parse_hypo71`): `os.chdir` into the loc folder is restored only on the success path; if `open(hypo_file)` raises, the `os.chdir(current_folder)` at line 128 is skipped, leaving the process in the wrong directory. Use `try/finally`. (Module is marked DEPRECATED at line 4 and no longer called — low priority, but it is retained dead code worth confirming as intentional.)

## Summary
20 files reviewed, 35 issues found (1 High, 18 Medium, 16 Low).
No issues found in: NLL_run/match_pre_post_relocation.py, generate_nll_corrections.py,
generate_complem_maps.py (minor import-time-empty-`loc/` caveat only),
generate_complem_figures.py, temp_picks/sort_picks.py.

Highest-priority item: the wrong-index bug in `fuse_bulletins.py:603`, which can
mis-promote/mis-drop loose P-phase matches in the fusion bookkeeping. The two
hardcoded `+ 4` date windows (`depth_maps.py:84`, `error_maps.py:84`) are latent
logic errors that only surface with a non-default `--time-range`.
