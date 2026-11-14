---
paths:
  - "add_temp_picks.py"
  - "temp_picks/**"
  - "complem_figures/temporary_network_impact.py"
---

# External Pick Ingestion (`temp_picks/`)

A self-contained sub-pipeline for ingesting picks from external sources into `obs/NLL_result.obs`, producing `obs/NLL_result_augmented.obs`. All scripts live in `temp_picks/` and are importable as a package (`from temp_picks.<module> import <function>`). Log files are written to `temp_picks/console_output/`.

The root-level script **`add_temp_picks.py`** orchestrates the full pipeline in sequence (6 steps). Only steps 1–2 (theoretical tables, QC figure) are skipped when their output exists; the merge, convert and match steps always run.

Ten external datasets are ingested, **in a fixed order that encodes precedence** — the **SANIMS** bulletin first, so its manual analyst picks win every `(station, phase)` slot over the machine-picked sources, and `TEMP_STB` after `TEMP_OMP` so `match_picks`'s per-event `(station, phase)` dedup silently drops the overlap. Each source is matched against the *previous* output, so the stage is **not re-runnable** against an existing `NLL_result_augmented.obs`. **SANIMS** (`TEMP_QML`) is read from `temp_picks/all_picks/SANIMS/Bulletin_Sanims.xml`, the QuakeML of the ICTJA/SeisComP SANIMS experiment (2019-09 → 2020-11, 1 470 events, 44 349 picks), copied from the bulletin folder delivered with it. Its `Cata_Sanims.txt` / `Sanims_bulletin.txt` / `Sanims_report.txt` exports are a stale 874-event subset and are not read. Only manual picks are kept, and picks whose arrival has `timeWeight 0` are dropped (36 in total). The analysts' own pick uncertainties (median 0.2 s) are **deliberately overwritten** by the project's 0.05/0.15 s. 10 stations are absent from the inventory and stay unresolved (1 181 picks): `H2.TO01/CA01/CA04`, `FB.FIJA`, `CA.VILA`, and the Barcelona stations `YY.CASN/TIAN/JPLA/BALM/ISAB`. The last, **`Pamplona_PN_picks.txt`** (`TEMP_RSB`, prob ≥ 0.5), is not produced here: it is written into `temp_picks/pick_files/` by the sibling project `../PhaseNet_Pamplona` (`python main.py export`), PhaseNet picks for the Pamplona region incl. the YK temporary stations (C015–C018, F001–F003, Oct–Dec 2020). It comes **last**, so it only fills `(station, phase)` slots no other source filled, and it must exist before `add_temp_picks.py` runs. Figures that follow were measured over the eight in between, before SANIMS and Pamplona were added, and before the 2026-10-14 `TEMP_OTH` fixes (S read from the card's HH:MM, the original maupasacq, quality-4 lines dropped whole), which change them. Measured over the eight: **308 731 picks added**, 31.9 M skipped for no event within 60 s (these are continuous PhaseNet detection streams, not bulletins), 534 312 as duplicates, 201 200 outside the travel-time band, 1 661 ambiguous, **0** for an unknown station.

`pick_files/merged_other.txt` (`TEMP_OTH`) is **built by hand, not by code**: the plain concatenation of the HYPO71 phase files in `all_picks/PICKS_MATTHIEU/` in the order argeles, data_bigorre_naturels, nat_pisco, stryde, spf, maupasacq_seismotech_reprocessed (CR stripped). Use that original maupasacq file, never a reformatted copy: the old `maupasacq.txt` had carried P seconds ≥ 60 into the next minute without adjusting the S, putting 5 950 S picks 60 s late (deleted 2026-10-14). In HYPO71, both P and S seconds count from the card's HH:MM and may exceed 60.

Added picks get the project's standard **0.05 s P / 0.15 s S** — i.e. machine picks are given the same declared precision as manual readings, and the available `phase_score` is not used to modulate it.

The travel-time bands come from the **same velocity model as NonLinLoc**, as Pyrocko `.nd` files, perturbed **±5 %** and enveloped over source depths 0 and 30 km. The two copies of the model are maintained separately and nothing checks that they still agree.

| Script | Role |
|--------|------|
| `build_theoretical_tables.py` | Runs Pyrocko's `cake` CLI to compute P/S travel-time envelopes (±5% velocity, 0–100 km) → `temp_picks/tables_Pyr.csv` |
| `merge_omp_picks.py` | Merges all yearly OMP/PhaseNet CSV files from `picks_OMP/` subdirectories → `pick_files/merged_omp.csv`; station `SMC` and year `2026` excluded by default |
| `merge_pyrenees_picks.py` | Concatenates RaspberryShake/PhaseNet `.txt` files from `picks_station_pyrenees/` and `picks_station_pyrenees2/` → `pick_files/merged_pyrenees.txt` and `pick_files/merged_pyrenees2.txt` |
| `convert_picks.py` | Converts external pick files to `.obs` pick line format; maps station names to internal codes via `GLOBAL_code_map.txt` — by short name, restricted to the pick's network when the format carries one (`TEMP_RSB`/`TEMP_OMP`/`TEMP_STB`/`TEMP_QML`; short names such as `C015` are reused by YK 2020 and YS 2021), then by epoch compared **at day granularity** (the pick passes only its date). Formats `TEMP_OBS`, `TEMP_RSB`, `TEMP_OMP`, `TEMP_OTH`, `TEMP_STB` and `TEMP_QML` are supported; new formats are registered in `FORMAT_HANDLERS`. `TEMP_STB` (Strasbourg/RENASS-OMP parquet picks) reads directories of parquet files instead of a single text file and drops picks below a `phase_score` threshold. `TEMP_QML` reads a QuakeML 1.2 bulletin event by event (`lxml` iterparse), keeping manual, non-zero-weight picks. Unresolved stations are logged as an end-of-run summary. |
| `match_picks.py` | Matches converted picks to bulletin events: 60 s time window + residual filter (±0.1 s P, ±0.3 s S, plus ±2.5 s t0-error margin); appends new picks and updates `PhaseCount`; chains against `obs/NLL_result.obs` → `obs/NLL_result_augmented.obs`; auto-sorts output via `sort_picks`. |
| `sort_picks.py` | Sorts pick lines within each event block by ascending arrival time. |
| `plot_travel_times.py` | QC figure: scatter of observed (distance, travel time) picks over theoretical P/S bands. |
