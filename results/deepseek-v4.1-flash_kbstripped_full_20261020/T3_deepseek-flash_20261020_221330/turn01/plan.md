# Plan

1. Load all source fields on their original grids. Use `load_regular_data` to read `cmip6_historical` variable `tos`, monthly for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 1971-01 to 2000-12, and to read each model’s fixed `areacello`. Use `load_climatology_data` for `woa23` variable `temperature`, climate normal 1971-2000, annual time slice, 0 m depth, on the 1-degree grid. Use `load_landsea_mask` for WOA23 `landsea_01` and treat cells with values greater than 1 as ocean. Use `load_areas` for the WOA23 1-degree cell areas. Keep all model data on its native grid at this stage.
   - Note [missing-values-and-land-cells]: Build each field’s valid mask from the field itself, mask the area weights to exactly the same valid cells, and treat values near 1e20 as unmasked fill values.
   - Note [cmip6-tos-bulk-top-cell]: `tos` is the bulk temperature of the top cell including under sea ice, so polar values near −1.8 °C are data and enter the native-grid means; make the model and WOA23 ocean masks consistent where WOA23 masks sea ice before comparing.
   - Note [cmip6-static-grid-fields]: Use each model’s published `areacello` as the only correct area weight on its grid, and report the total ocean area covered by each native mean; do not rescale total areas between models.

2. For each model, use `transform_unit` to convert `tos` to degrees Celsius if it is loaded in kelvin, then use `temporal_aggregate` to form the 1971-2000 annual-mean SST climatology from the monthly `tos` as a day-weighted long-term mean. On each native grid, compute the model’s global mean with `compute_basic_stat` using `areacello` as weights restricted to valid `tos`, and report it in degrees Celsius to 2 decimal places together with the total ocean area covered.
   - Note [omip-weighted-means-and-integrals]: Compute each global mean as Σ(tos·areacello)/Σ(areacello) over valid ocean cells with NaN-aware weights, never as an unweighted mean over the grid; use area-weighted RMSE and Pearson correlation, with `compute_difference_metrics` weighting by area by default and `analyze_correlation` using `weighting='area'`; state whether the WOA23 whole-grid mean is area-weighted.
   - Note [omip-native-grid-vs-regridded]: Compute the four global means on each source’s own grid with its own `areacello`; regrid only for maps, bias, RMSE and pattern correlation, and state the regridding method as conservative when going to the coarser 1-degree grid.

3. Compute the WOA23 0 m annual global mean over the whole 1-degree grid with `compute_basic_stat` and `load_areas` weights, restricted to valid WOA23 temperature cells with `landsea_01 > 1`; report it in degrees Celsius to 2 decimal places and state that it is an area-weighted mean over valid ocean cells, not an unweighted average over all 1-degree cells.

4. Use `regrid_horizontal` with the conservative method to regrid each model’s 1971-2000 annual SST climatology onto the WOA23 1-degree grid. Then use `build_comparison_mask` over the three regridded model fields and the WOA23 reference with its default `ocean_only=True`; apply `.where(mask)` to every regridded model field, to the WOA23 field, and to the area weights, and report the mask’s cell count.
   - Note [common-comparison-mask]: Regrid all three models to the WOA23 1-degree grid, build one mask over all fields with `build_comparison_mask`, keep its default `ocean_only=True` because WOA23 is the reference, apply the mask to every field before bias, RMSE and correlation, and report the mask’s cell count.

5. On the common comparison mask from [common-comparison-mask], use `compute_difference_metrics` for each model comparing `model_regridded.where(mask)` with `woa23.where(mask)`, with its default area weighting. Report each model’s `mean_data` and the WOA23 `mean_reference` as the shared-cell mean SSTs in degrees Celsius to 2 decimal places, alongside each model’s mean bias `delta_SST` in degrees Celsius to 2 decimal places and RMSE in degrees Celsius to 3 decimal places.

6. Use `analyze_correlation` with `weighting='area'`, as required by [omip-weighted-means-and-integrals], to compute the spatial Pearson correlation between each masked regridded model SST field and the masked WOA23 SST field, reported to 3 decimal places.

7. Rank CESM2, GFDL-ESM4 and MPI-ESM1-2-LR from smallest to largest absolute `delta_SST` using the values from Step 5.

8. Plot Figure 1 as a 2-by-2 panel using the masked regridded fields from Steps 4–5: the CESM2, GFDL-ESM4 and MPI-ESM1-2-LR bias fields on the common comparison cells with one shared diverging colour scale centred on zero, and the WOA23 0 m reference SST on the same cells with a sequential colour scale. Label every panel and give each colour scale a unit label in degrees Celsius, then use the maps to describe the spatial structure of the bias.

9. Plot Figure 2 as a bar chart of the four native-grid global mean SST values from Steps 2–3: CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23. Annotate each bar with its value and label the y axis in degrees Celsius.

Other notes:
- [cmip6-grid-labels] GFDL-ESM4 carries `gr`: its “native-grid” numbers are on its 1-degree regridded archive grid with a regridded land–sea mask, and its resolution differs from CESM2 (~100 km) and MPI (~250 km), so narrow features are not comparable in amplitude.
- Native-grid global means and shared-cell comparison means are different quantities; report and label them separately.
