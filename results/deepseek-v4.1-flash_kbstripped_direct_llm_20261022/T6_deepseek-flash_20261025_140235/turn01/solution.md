## Bottom line

The run reports `status: success`, but **the assessment itself did not succeed**. Every quantity the question asks for — global and basin-mean nitracline depths, model-minus-WOA23 differences, spatial correlation, RMSE, basin ordering — came back as `NaN` in the results object. Two independent bugs in the script are responsible, and both must be fixed before any of the requested numbers exist:

1. **The WOA23 mask files were never parsed.** `landsea_rows` and `basin_rows` are both empty lists, so `landsea` has `n_valid = 0`, `ocean_mask` is `False` everywhere, `basin_map` is entirely `NaN`, and `basin_group` is `0` everywhere. Because `area_weighted_mean()` requires a non-empty valid mask, it returns `NaN` for every basin of every source; likewise `corr_rmse` finds zero valid grid cells and returns `NaN`.
2. **Two of the three model grids were built from index coordinates, not geographic coordinates.** For CESM2, `model_lats`/`model_lons` came out as 1…384 and 1…320 (the `nlat`/`nlon` index variables), and for MPI-ESM1-2-LR as 0…219 and 0…255 (the `j`/`i` indices). Only GFDL-ESM4 has true coordinates (−89.5…89.5, 0.5…359.5). Consequently the nearest-neighbour regridding for CESM2 and MPI matched nonsense positions, and the equatorial-Pacific section is empty for CESM2 (`sections['CESM2'].n_valid = 0`) and has only 132 valid cells out of 6400 for MPI (a spurious strip rather than a section).

## What the numbers that *do* exist show — and why they cannot be used

The only finite fields are the raw, unmasked, unweighted nitracline arrays after regridding:

| Source | n valid (of 64 800) | min (m) | max (m) | mean (m) |
|---|---|---|---|---|
| WOA23 | 41 071 | 2.5 | 1 025 | 138.0 |
| CESM2 | 24 367 | 10.0 | 673.2 | 149.1 |
| GFDL-ESM4 | 49 464 | 6.3 | 950.0 | 104.9 |
| MPI-ESM1-2-LR | 16 570 | 11.5 | 522.5 | **28.3** |

These are *not* the requested statistics: they are not area-weighted, not restricted to the WOA ocean mask, and they are computed on a regridded grid whose coordinates are wrong for two models. The MPI value (global-mean nitracline at ~28 m, versus 138 m observed) is a clear artefact — a nitracline that shallow everywhere is physically implausible for the oligotrophic subtropics and reflects the broken index-coordinate regridding rather than a model bias. The WOA field itself looks sane (mean 138 m, maximum 1025 m, consistent with the deepest searched gradient being centred between the 1000 m and 1100 m levels, as specified).

Because the masks are empty, I cannot report:
- global means or any model-minus-WOA difference (all `NaN`);
- Pacific/Atlantic/Indian means or differences (all `NaN`);
- correlation to 3 d.p. or RMSE to 1 d.p. (all `NaN`);
- the observed basin ordering or any model's reproduction of it — the lists returned (`['Atlantic','Pacific','Indian']` for all four sources) are simply the alphabetical default produced by sorting `NaN` values, not a physical result.

There is therefore **no defensible answer to "which model reproduces the observed nitracline best"**: with no correlation or RMSE against observations, and with two of the three model fields geometrically corrupted, no ranking can be supported.

## The figures

Both figures were written to disk, but neither is a valid assessment figure.

- **Figure 1 (2×2 nitracline maps)** — the four panels use one shared `YlGnBu` scale (2nd–98th percentile of the pooled values), deeper = darker, which is the requested convention. However, the land mask applied is `ocean_mask`, which is `False` everywhere, so no land masking is drawn: the WOA panel shows only where the climatology has data, and the CESM2/MPI panels show fields placed at wrong geographic positions. Only the GFDL-ESM4 panel is on a trustworthy grid.
- **Figure 2 (equatorial Pacific sections, 2°S–2°N, 120°E–80°W, 0–1000 m, depth increasing downward)** — the WOA23 and GFDL-ESM4 panels are populated (WOA 14 245 valid cells; GFDL 4 941), share a 0–99th-percentile `viridis` scale in µmol kg⁻¹, and each carries its own red nitracline line (WOA mean along the section 123 m, GFDL 118 m — closely comparable). The CESM2 panel is **completely blank** (0 valid cells, no nitracline line), and the MPI panel contains only a small spurious fragment (132 cells, nitracline defined at only 11 of 160 longitudes, mean 45 m). Note also that the model sections were converted from mol m⁻³ to µmol kg⁻¹ by a fixed factor 10⁶/1025; that is a reasonable approximation but is not an in-situ density conversion, so a small density-dependent bias remains in the model panels relative to WOA.

## What needs fixing to actually answer the question

1. **Read the masks correctly.** Inspect `landsea_01.msk` and `basinmask_01.msk` directly (they may be binary, fixed-width, or differently headed than the `#`-comment ASCII assumed by `read_woa_msk`, which silently returned zero rows without raising). A robust fallback is to define the ocean mask from the WOA nitrate field itself (`np.isfinite(woa_no3[...,0])`), and to take the basin codes from the mask's NetCDF equivalent or a known WOA13 basin-code array. Add explicit assertions (`assert ocean_mask.sum() > 30000`, `assert set(np.unique(basin_group)) == {0,1,2,3}`) so a silent mask failure cannot recur. Marginal-sea folding should then be re-verified against the WOA code list (1 Atlantic, 2 Pacific, 3 Indian; 4–6, 9 → Atlantic; 7–8 → Indian; 10 Southern and 11 Arctic excluded).
2. **Use the real 2-D/1-D geographic coordinates** for CESM2 and MPI-ESM1-2-LR (the `lat`/`lon` or `nav_lat`/`nav_lon`/`latitude`/`longitude` variables), not `nlat`/`nlon`/`j`/`i`. Then re-do the native-grid nitracline computation and the regridding; the section extraction will follow automatically (a 2°S–2°N band should contain several model rows and yield a continuous longitude–depth field).
3. **Weight with each model's own `areacello`** (as the task specifies) rather than `cos(lat)` on the WOA grid, or at minimum document the approximation; this matters most for the strongly non-uniform CESM2 grid.
4. Optionally convert model nitrate with in-situ density (`rho = ` EOS-80/TEOS-10 from the models' theta and salinity) so the Figure 2 shared colour scale is strictly comparable between models and WOA.
5. Re-run and report only after `summary`, `corr_rmse` and `model_basin_orders` contain finite values; the basin-ordering test should compare ranks of finite means, not `NaN` sorts.

Until steps 1 and 2 are done, the honest answer is that the simulated nitracline depth of CESM2, GFDL-ESM4 and MPI-ESM1-2-LR has **not** been successfully assessed against WOA23 in this run.
