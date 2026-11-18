# Plan

1. Load the data.
   - Use load_regular_data for cmip6_historical, variable `tos`, for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. The period is 1971-01 to 2000-12, monthly, on each model's own grid.
   - Use load_climatology_data for woa23 temperature, the 1971-2000 normal, annual time slice, at 0 m, on the 1° grid.
   - Use load_basin_mask only if needed. Use load_landsea_mask to get landsea_01 for the WOA23 1° ocean test (value > 1).
   - Use load_areas for each model's own areacello on its native grid. Use load_areas for the WOA23 1° grid as well.
   - Note [missing-values-and-land-cells]: use only CF-masked reads, so the 1e20 land and below-sea-floor fill values never enter any sum or mean. Check that the loaded fields contain no values near 1e20.
   - Note [cmip6-grid-labels]: CESM2 and MPI-ESM1-2-LR are `gn`. GFDL-ESM4 is `gr`, a 1° archive grid, so its "native" means are on that 1° grid with a regridded land-sea mask. Say this when reporting the native means.

2. Build each model's 1971-2000 SST climatology.
   - Use temporal_aggregate on each model's tos to get the DAY-WEIGHTED long-term mean over 1971-2000. This gives one 2D climatological field per model.
   - Note [cmip6-tos-bulk-top-cell]: tos is the top-cell bulk temperature. The top cell is 0–10 m for CESM2, 0–2 m for GFDL-ESM4 and 0–12 m for MPI-ESM1-2-LR. This is treated as like-for-like with WOA23 at 0 m, but the differing thicknesses must be mentioned. Keep sea-ice cells, which carry about −1.8 °C under-ice water, in all means.

3. Compute the native-grid global means (Figure 2).
   - Use compute_basic_stat with op='area_weighted_mean' on each model's climatology, using that model's own areacello as the weight. Restrict the weights to the same valid ocean cells as the data.
   - Compute the WOA23 0 m annual field global mean over the whole 1° grid with compute_basic_stat (`area_weighted_mean`), using 1° cell areas on valid (non-missing) cells.
   - Note [omip-weighted-means-and-integrals]: the mean is Σ T·area / Σ area over valid ocean cells only. Mask the weights to the same cells as the data. Do not use an unweighted mean. This step also applies to the 1° comparison statistics below.
   - Note [cmip6-static-grid-fields]: use each model's own areacello (load_areas) and record attrs['area_source']. Total ocean area differs between models by about 1%. State this and do not rescale the areas.
   - Note [omip-native-grid-vs-regridded]: these native means come from each model's own grid. State which numbers come from which grid.

4. Regrid to the WOA23 1° grid.
   - Use regrid_horizontal with the 'conservative' method to move each model's climatology onto the WOA23 1° grid. State the method used.

5. Build the common mask.
   - Use build_comparison_mask over all four fields: the three regridded models and WOA23. Use ocean_only=True with the landsea_01 ocean test (value > 1).
   - Apply the mask with `.where(mask)` before any metric. Report the cell counts the mask records.
   - Note [common-comparison-mask]: build ONE mask and score every model on that same sample. Do not use pairwise masks. Take the mean model, mean WOA23 and bias from the `compute_difference_metrics` outputs (`mean_data`, `mean_reference`) on this shared sample. State that these differ from the native-grid global means.

6. Compute the scores on the shared mask.
   - Use compute_difference_metrics (model vs WOA23, area weighted by default) to get, for each model, the mean model SST, mean WOA23 SST, mean bias, and RMSE. Also get the largest departure and where it falls.
   - Use analyze_correlation with weighting='area' for the spatial Pearson correlation of each model with WOA23.
   - Rank the models by smallest to largest absolute mean bias from the unrounded values.

7. Describe the spatial structure of the bias.
   - Use the masked bias maps for each model to identify regions of positive and negative bias. Use compute_difference_metrics' location of the largest departure, and locate_extreme if needed, to find the largest positive and negative biases and their locations.
   - Describe how the bias patterns differ between the models.

8. Produce the figures.
   - Figure 1 is a 2x2 panel of maps. Three panels show the bias maps (model minus WOA23) on the common mask, with one shared diverging colour scale centred on zero. The fourth panel shows WOA23 0 m SST on the same cells with a sequential colour scale. Label every panel and give each colorbar a °C label.
   - Figure 2 is a bar chart of the four native-grid global means (CESM2, GFDL-ESM4, MPI-ESM1-2-LR, WOA23). Annotate each bar with its value and label the y axis in °C.
   - Round only when reporting: 2 decimals for means and bias, 3 decimals for RMSE and correlation.

Other notes:
- The native-grid global means (Figure 2) are over each model's own ocean grid, and the WOA23 value is over the whole 1° grid. The 1° comparison means are over the shared mask, so they differ from the native means and are not interchangeable.
- GFDL-ESM4 "native" is a 1° regridded (`gr`) grid, not a tripolar grid.
- The tos top-cell thicknesses differ between models (CESM2 10 m, GFDL-ESM4 2 m, MPI-ESM1-2-LR 12 m), compared against WOA23 at 0 m.
- Regridding method is conservative; RMSE and correlation depend on this choice.
- Total ocean area differs between models by about 1% because each uses its own areacello.
- [common-comparison-mask] A single shared mask is used for all models and WOA23. Report its cell counts.
- [missing-values-and-land-cells] Ensure land and below-sea-floor 1e20 fill values are masked before any statistic.
