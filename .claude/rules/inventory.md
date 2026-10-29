---
paths:
  - "build_global_inventory.py"
  - "fetch_inventory/**"
  - "stations/**"
  - "complem_figures/station_colocation.py"
  - "NLL_run/generate_regional_runfiles.py"
  - "NLL_run/export_quakeml.py"
---

# Stage 1 — Station Inventory Fusion

- Source: FDSN XML files + OMP CSV files (in `stations/`)
- Script: `build_global_inventory.py` → calls modules in `fetch_inventory/`
- Output: `stations/GLOBAL_inventory.xml` + `stations/GLOBAL_code_map.txt`
- Each station gets a unique code (`NET.NNNN`); duplicates are removed by distance threshold (20m). Co-located stations inherit the code of the group's **best-ranked** member under `station_priority()` — lowest wins: **permanent network** (`FR RA RD ES CA LC AM`, the `perm` rows of `stations/all_networks.xlsx`, AM/RaspberryShake counted as permanent) before temporary; **RESIF/RENASS** (`FR RA RD`) among the permanent; a **real elevation** before one left at 0; `XX` last; oldest `start_date` breaking ties. `all_networks.xlsx` is *not* read at runtime — `openpyxl` is absent from `seisbench_env`, so the two sets are module constants and the sheet stays the human-maintained reference
- **That one rule is used in all three places a group must yield one representative**: the code it lends the group (`_combine_close_stations`), the lat/lon/elevation written into the NLL `GTSRCE` lines (`generate_regional_runfiles._find_station_info`, §4) and the station name in the QuakeML export (`export_quakeml.load_station_epochs`, §7). Before it, each of the three picked a different member. **They share the ranking, not the outcome**: the export applies it only after the pick's date, so for 37 codes the name it prints is a different member than GTSRCE's donor (§7). `_find_station_info` resolves the code-map block's members **against the inventory**, not by their order in the file, so it stays correct against a `GLOBAL_code_map.txt` written before the rule existed
- **Declared elevations are trusted**: a DEM fills gaps and arbitrates conflicts, it never overrides a lone declared value — only clear errors get corrected (by hand, in the source file). A co-located group whose declared elevations spread over **100 m** (`ELEVATION_TOLERANCE`) has its **representative** set to the declared value closest to EU-DEM at its position (`_arbitrate_group_elevation`, logged; a warning when even that value is > 100 m from the DEM). Other members keep their own value, and the arbitration runs after the code is chosen, so it cannot change a code. The only group concerned once `HO` is added is 24.SG07 (373 m, a typo — EU-DEM 660 m) / HO.V07 (662 m). The 20 m grouping is purely horizontal: nothing checks the vertical separation of a group (e.g. `LC.CANF`, in the Canfranc tunnel, and `IB.E149` carry identical coordinates and 1 205.5 m, ~283 m below the surface)
- Stations left without an elevation inherit one from their co-located group, then from **EU-DEM 25 m** via the OpenTopoData public API (`_fill_station_elevations`, reusing `_fill_missing_elevations._get_elevation`; 1 call/s, 1 000/day, no key). It replaced Open-Elevation, whose anonymous monthly quota ran out (HTTP 429) and made every lookup silently return 0. 156 of 2 151 carried none — `XI` 99, `24` 31, `X7` 15, `ZU` 8, `XX` 2, `6E` 1, whole nodal deployments whose StationXML never populated the field. ⚠ **This is the only step of the pipeline that reaches the network** (~1 s/station, so ~3 min); `--no-fill-elevations` / `MergeInventoryParams.fill_elevations=False` skips it
- Sources: 4 FDSN inventories (RESIF 346, ICGC 262, ORFEUS 103, GFZ 12) + 12 OMP/temporary deployments (1 691 stations). 2 414 in → 2 151 stations / 1 997 unique codes
- `_remove_fdsn_duplicates.py` and `_convert_csv_to_stationxml.py` are **standalone preprocessing**, run by hand per deployment; `build_global_inventory.py` imports neither. They assign network code `XX` to OMP stations (see the `XX`-last term above). The third `_`-prefixed module, `_fill_missing_elevations.py`, is still hand-run on the CSVs, but its `_get_elevation()` **is** now imported by `merge_station_inventories`
