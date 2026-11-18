# Plan

1. Load CMIP6 historical monthly `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 1971-01 to 2000-12, together with each model’s fixed `areacello`, using `load_regular_data`; keep each model on its native grid.
2. For each model, use `temporal_aggregate` to compute the 1971-2000 day-weighted long-term mean SST climatology from the monthly `tos` series.
3. Load WOA23 `temperature` for the 1971-2000 climate normal, annual time slice, at 0 m depth on the 1-degree grid using `load_climatology_data`; load the WOA23 1-degree land-sea mask `landsea_01` using `load_landsea_mask`; load 1-degree cell areas using `load_areas`.
4. For each model, compute the global mean of its 1971-2000 SST climatology over its own native ocean grid using `compute_basic_stat` with spatial weighting by that model’s `areacello` and only native ocean cells; report to 2 decimal places.
5. Compute the global mean of the WOA23 0 m annual field over the whole 1-degree grid using `compute_basic_stat` with 1-degree cell-area weights, without applying the land-sea ocean mask; report to 2 decimal places.
6. Regrid each model’s 1971-2000 SST climatology horizontally to the WOA23 1-degree grid using `regrid_horizontal`.
7. Build the common comparison mask using `build_comparison_mask` from the three regridded model climatologies, the WOA23 0 m annual field, and the WOA23 ocean condition `landsea_01 > 1`.
8. On those common comparison cells, compute the mean SST of each model and of WOA23 using `compute_basic_stat` with 1-degree area weights; compute each model’s bias and RMSE against WOA23 using `compute_difference_metrics`; compute the spatial Pearson correlation with WOA23 using `analyze_correlation`; report means and bias to 2 decimal places, RMSE and correlation to 3 decimal places; rank models by smallest to largest absolute mean bias.
9. Produce Figure 1 as a 2-by-2 panel: bias maps (model minus WOA23) for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR on the common comparison cells with one shared diverging colour scale centred on zero, plus WOA23 0 m reference SST on the same cells with a sequential colour scale; label panels and colourbars in degrees Celsius. Produce Figure 2 as a bar chart of the four global mean SST values from steps 4-5 (CESM2, GFDL-ESM4, MPI-ESM1-2-LR, WOA23), annotate each bar with its value, and label the y axis in degrees Celsius.
Other notes:
- CMIP6 values are day-weighted 1971-2000 climatological means from monthly data, while WOA23 is the annual 1971-2000 normal at 0 m.
- Native-grid global means for the models use each model’s own `areacello` and native ocean cells, not the common WOA23 comparison mask.
- The WOA23 global mean uses all valid WOA23 0 m cells on the whole 1-degree grid and does not apply the `landsea_01` ocean mask.
- Comparison means, biases, RMSE and correlations are evaluated only on the common cells where all three regridded models and WOA23 are valid and `landsea_01 > 1`.
- Horizontal regridding to the WOA23 1-degree grid introduces interpolation, so comparison metrics are not native-grid metrics.
- The ranking is based on absolute mean bias on the common comparison cells and may differ from ranking by native-grid global means.
