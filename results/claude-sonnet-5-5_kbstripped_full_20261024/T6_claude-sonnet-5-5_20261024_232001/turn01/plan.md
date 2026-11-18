# Plan

1. Load the data.
   - Use `load_regular_data` for cmip6_historical `no3` from CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, 1971-01 to 2000-12. Load each model on its own native grid and keep the profile down to the level just below 1000 m, so the gradient at the deepest level inside the top 1000 m can be computed.
   - Use `load_climatology_data` for woa23 `nitrate` (the all-period composite) on the 1° grid, also down to the level below 1000 m.
   - Load the WOA23 land-sea mask with `load_landsea_mask` (ocean is value > 1) and the basin mask with `load_basin_mask`.
   - Load `load_areas` for each model's areacello and for the WOA23 grid.
   - Load `load_depth_levels` for each dataset so the 1000 m search window and the layer boundaries are known.
   - Note [cmip6-grid-labels]: GFDL-ESM4 is `gr` (already on a 1° archive grid with a regridded land-sea mask), while CESM2 and MPI-ESM1-2-LR are `gn`. State this in the answer. MPI's coarser grid limits its resolution of the equatorial Pacific.
   - Note [missing-values-and-land-cells]: Build each model's valid mask from its own `no3`. Treat values near 1e20 as unmasked fill. Cells with no valid profile or too few levels must stay NaN in the nitracline map, never 0.

2. Average to the mean state first.
   - Use `temporal_aggregate` with `freq='mean'` on each model's monthly `no3` for 1971-2000, so month lengths are weighted. Do not use `compute_basic_stat` with a time mean. WOA23 is already a climatology.
   - Note [monthly-means-are-time-averages]: Take the 1971-2000 mean state before the gradient, because the nitracline is a nonlinear diagnostic. State this order in the answer.

3. Compute the nitracline depth.
   - Run `compute_vertical_gradient` on the native nitrate fields (models in mol m-3, WOA23 in µmol kg-1) and on the WOA23 field.
   - Use `compute_nutricline_depth` for the depth of the steepest increase of nitrate with depth. Search only the levels within the top 1000 m; the deepest one uses the level just below 1000 m.
   - Note [cmip6-bgc-units-per-volume]: The nitracline depth is scale-invariant, so take it from the native fields with no density conversion and no ×1000 rescaling.

4. Compute global and basin means.
   - Use `select_basin` with the WOA23 basin mask, `include_marginal_seas=True` and `depth_dependent=False`. Keep Pacific, Atlantic and Indian. Exclude Southern and Arctic as separate basins. Carry the same basin mask onto each model's grid.
   - Note [basin-and-region-definitions]: Use the surface labels because the nitracline depth is a 2-D quantity of the water column. Say that the mask follows WOA labelling and not the OMIP 10-region one.
   - Use `compute_basic_stat` with the area-weighted mean over each source's valid cells, using the areas from step 1. Do this globally and for each of the three basins, for the three models and WOA23. Report the values to the nearest whole metre, and each model's difference from WOA23 globally and by basin.
   - Note [omip-native-grid-vs-regridded]: Compute these means on the native grids. Regrid only for the maps and the pattern comparison.
   - Rank the basins by nitracline depth for each source, and compare each model's ordering with WOA23's.

5. Compare the models with WOA23 on a shared grid.
   - Use `regrid_horizontal` to put each model's nitracline map on the WOA23 1° grid.
   - Use `build_comparison_mask` over WOA23 and the three regridded model fields, with the default `ocean_only=True` and the WOA23 land-sea mask. Apply it with `.where(mask)` and report the cell counts it records.
   - Note [common-comparison-mask]: Use `compute_difference_metrics` on the shared mask for each model against WOA23. Report RMSE (1 decimal place, metres), `mean_data`, `mean_reference` and the bias. Use `analyze_correlation` for the spatial Pearson correlation (3 decimal places). If the native-grid means from step 4 are differenced, state that this difference is not the bias on the shared mask.
   - Decide which model reproduces the observed nitracline best from the RMSE and correlation together with the bias.

6. Figure 1.
   - Draw a 2×2 panel of global nitracline-depth maps for WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, all on one shared sequential colour scale where deeper is visually heavier.
   - Keep NaN cells blank.

7. Figure 2.
   - Take the nitrate section through the equatorial Pacific, averaged over 2S-2N, from 120E to 80W and the top 1000 m, with depth increasing downward. Use `compute_basic_stat` with the area-weighted latitude mean for the 2S-2N average.
   - Panels are in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, on one shared colour scale in one common concentration unit.
   - Overlay each panel's own nitracline depth along the section as a line.
   - Note [cmip6-bgc-units-per-volume]: Converting the models to µmol kg-1 needs the in-situ density from `compute_potential_density_anomaly` with `ref_pressure=None`. Use `transform_unit` for that conversion and state `attrs['density_assumption']`. Alternatively, show everything in mmol m-3 and state that WOA23 is per mass and the models per volume.

Other notes:
- State the order of operations: time-mean first, then the gradient, then the nitracline depth.
- State that the models are per volume (mol m-3) and WOA23 is per mass (µmol kg-1), and that the nitracline depth is unaffected.
- State that the basin means are on native grids and the RMSE and correlation are on the WOA23 1° shared mask, so the two sets of numbers are not the same quantity.
- State that GFDL-ESM4 is already on a 1° archive grid, and that MPI's coarse grid limits its equatorial Pacific.
- State that Southern and Arctic are excluded from the basin results.
