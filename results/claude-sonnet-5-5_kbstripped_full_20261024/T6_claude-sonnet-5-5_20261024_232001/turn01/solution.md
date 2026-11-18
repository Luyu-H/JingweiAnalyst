## Summary

MPI-ESM1-2-LR reproduces the observed nitracline best on the shared WOA23 1° grid. It has the lowest RMSE (151.8 m) and the highest correlation (0.617), and it is the only model that reproduces the WOA23 basin ordering. It also has the largest negative bias (about −33 m), and all three models are shallow-biased and only modestly skilful.

## Methods as run

- **Data.**
  - CMIP6 historical `no3` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, 1971-01 to 2000-12, each on its own native grid.
  - WOA23 `nitrate` from `load_climatology_data` with `time_slice="annual"`. The comparison-mask attributes record this as the all-period composite.
  - Each load used `depth_range=[0, first level >1000 m]`, with the limit taken from `load_depth_levels`: CESM2 1106 m, GFDL-ESM4 1100 m, MPI 1085 m, WOA23 1050 m.
  - Values ≥1e10 (the 1e20 fill) were set to NaN.
- **Grid labels (from the task notes; the code did not check them).** GFDL-ESM4 is `gr`, already on a 1° archive grid. CESM2 and MPI-ESM1-2-LR are `gn`. MPI's coarse grid limits its resolution of the equatorial Pacific.
- **Order of operations.**
  1. Each model's monthly `no3` was averaged over 1971-2000 with `temporal_aggregate(freq="mean")`.
  2. The nitracline depth was then computed from that mean state with `compute_nutricline_depth(min_depth=0, max_depth=1000)`.
  3. WOA23 is already a climatology, so step 1 was skipped for it.
  - The code did not call `compute_vertical_gradient` separately; the gradient is handled inside `compute_nutricline_depth`.
  - Columns with fewer than 3 valid levels were set to NaN, never 0.
- **Units.** Models are per volume (mol m⁻³) and WOA23 is per mass (µmol kg⁻¹). The depth maps were computed on the native fields with no conversion, and the nitracline depth is scale-invariant.
- **Basin and global means.**
  - `select_basin` was used with `include_marginal_seas=True` and `depth_dependent=False`. The labels are WOA-style, not the OMIP 10-region ones.
  - Pacific, Atlantic and Indian were kept. Southern and Arctic were excluded from the basin results.
  - Means come from `compute_basic_stat(op="area_weighted_mean")` on native grids.
  - The "global" mean is over all valid ocean cells, so it includes the Southern and Arctic Oceans.
- **Shared-grid comparison.**
  - Model maps were regridded to the WOA23 1° grid with `regrid_horizontal(method="conservative")`.
  - `build_comparison_mask(ocean_only=True)` gave a shared mask of **40,462 cells**. The smallest input had 41,019 valid cells, and the land-sea mask removed none.
  - RMSE, bias and `mean_data`/`mean_reference` come from `compute_difference_metrics`, and the correlation from `analyze_correlation`.
  - The native-grid basin and global means and the shared-mask RMSE and correlation are different quantities. The native-grid global difference is not the shared-mask bias.

## Global mean nitracline depth (native grids, m)

| Source | Mean | Difference from WOA23 |
|---|---|---|
| WOA23 | 166 | – |
| CESM2 | 147 | −18 |
| GFDL-ESM4 | 134 | −32 |
| MPI-ESM1-2-LR | 133 | −33 |

## Basin means (native grids, m; difference from WOA23 in brackets)

| Source | Pacific | Atlantic | Indian |
|---|---|---|---|
| WOA23 | 201 | 165 | 180 |
| CESM2 | 163 (−38) | 110 (−56) | 166 (−14) |
| GFDL-ESM4 | 160 (−41) | 91 (−74) | 170 (−10) |
| MPI-ESM1-2-LR | 184 (−17) | 112 (−53) | 121 (−58) |

**Basin ordering, deepest to shallowest**
- WOA23: Pacific > Indian > Atlantic.
- CESM2: Indian > Pacific > Atlantic. This does not match WOA23.
- GFDL-ESM4: Indian > Pacific > Atlantic. This does not match WOA23.
- MPI-ESM1-2-LR: Pacific > Indian > Atlantic. This matches WOA23, but its Indian mean (121 m) is only about 9 m deeper than its Atlantic mean (112 m).

## Shared 1° mask (40,462 cells)

| Model | RMSE (m) | Bias (m) | `mean_data` (m) | Pearson r | Area-weighted r |
|---|---|---|---|---|---|
| CESM2 | 163.1 | −16.7 | 150.6 | 0.396 | 0.403 |
| GFDL-ESM4 | 161.4 | −31.6 | 135.7 | 0.502 | 0.498 |
| MPI-ESM1-2-LR | 151.8 | −32.5 | 134.8 | 0.617 | 0.598 |

`mean_reference` (WOA23) is 167.3 m in all three comparisons. The unweighted mean of the masked WOA23 map is 140.6 m, so these tool-reported means appear to be area-weighted. That is my inference from the numbers, not something the code set explicitly.

The spatial standard-deviation ratios (model/WOA23) are 0.56 for CESM2, 0.76 for GFDL-ESM4 and 0.88 for MPI-ESM1-2-LR.

## Interpretation

- **Bias and RMSE.**
  - All three models have a shallow bias, concentrated in the Atlantic (−53 to −74 m).
  - CESM2 has the smallest global and shared-mask bias, but the weakest spatial correlation and the weakest spatial variability.
  - MPI has the largest mean bias but the best pattern agreement and RMSE.
  - RMSE values of 150–163 m are comparable to the mean depth itself, so the skill of all three is modest.
- **Best model.** MPI-ESM1-2-LR ranks first on both RMSE and correlation, so it is judged best overall, with the caveat of its large negative bias.
- **Basin ordering.** CESM2 and GFDL-ESM4 put the Indian Ocean deeper than the Pacific, mainly because their Pacific nitracline is too shallow (−38 and −41 m).
- **MPI regional errors.** MPI's good Pacific mean (−17 m) contrasts with a very shallow Indian mean (−58 m).

## Figures

- **Figure 1.** A 2×2 map panel of nitracline depth for WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, with the models regridded to the WOA23 1° grid.
  - All panels share one sequential `YlGnBu` scale, where deeper is darker. The scale runs from 0 m to the pooled 98th percentile (about 584 m), with an "extend max" arrow.
  - NaN cells (land and invalid columns) are left blank.
  - Panel titles give the native-grid global mean.
- **Figure 2.** A 2×2 longitude-depth section of nitrate through the equatorial Pacific, in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR.
  - It covers 120°E–80°W and 0–1000 m, with depth increasing downward.
  - All panels use one `viridis` scale in mmol m⁻³, from 0 to the pooled 99th percentile (about 40.9).
  - WOA23 was converted from µmol kg⁻¹ with `transform_unit` under a constant seawater density of 1.025 kg L⁻¹ (the recorded `density_assumption`, which is approximate). The models were converted from mol m⁻³ to mmol m⁻³ with no density needed.
  - The in-situ density route via `compute_potential_density_anomaly` was therefore not used; this is the "everything in mmol m⁻³" alternative.
  - Model sections were bilinearly regridded to the WOA23 grid.
  - The 2°S–2°N average was a cos(latitude)-weighted mean on the 1° grid, not `compute_basic_stat`.
  - The red line on each panel is the 2°S–2°N mean of that source's own nitracline map (WOA23 native; models regridded conservatively). It is not the nitracline of the averaged profile.
  - The equatorial-band mean of that line is 115 m for WOA23, 122 m for CESM2, 120 m for GFDL-ESM4 and 46 m for MPI-ESM1-2-LR. MPI's much shallower value is consistent with its coarse equatorial resolution.

## Caveats

- The discrete maximum-gradient diagnostic depends on each dataset's vertical grid, which differs between WOA23 (about 25 m spacing near the surface, coarser deeper) and the models.
- The diagnostic is less well defined in weakly stratified high-latitude columns.
