# Plan

1. Load CMIP6 historical monthly `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 1971-01 to 2000-12 using load_regular_data; load each model's native published `areacello` cell areas using load_areas.
2. Use temporal_aggregate to reduce each model's monthly `tos` over 1971-2000 to a day-weighted long-term SST climatology.
3. For each model, use compute_basic_stat to compute the global mean of its 1971-2000 SST climatology over its own native ocean grid, weighted by that model's native `areacello`.
4. Load WOA23 temperature for the 1971-2000 climate normal, annual time slice, at 0 m depth on the 1-degree grid using load_climatology_data; load the WOA23 1-degree `landsea_01` mask using load_landsea_mask; load WOA23 1-degree cell areas using load_areas.
5. Use compute_basic_stat to compute the global mean of the WOA23 0 m annual field over the whole 1-degree grid, weighted by WOA23 cell areas and using valid WOA23 cells.
6. Use regrid_horizontal to regrid each model's 1971-2000 SST climatology onto the WOA23 1-degree grid.
7. Use build_comparison_mask on the WOA23 grid to select cells where WOA23 `landsea_01` is greater than 1, WOA23 SST is valid, and all three regridded model SST fields are valid.
8. On those comparison cells, use compute_basic_stat with WOA23 cell areas as spatial weights to compute the mean SST of each model and the mean SST of WOA23.
9. On the same comparison cells, form the per-cell model-minus-WOA23 difference for each model. Use compute_difference_metrics for each model against WOA23 to obtain RMSE, and use compute_basic_stat on the difference field with WOA23 cell areas to obtain the mean bias delta_SST. Use analyze_correlation for each model against WOA23 to obtain the spatial Pearson correlation.
10. Rank CESM2, GFDL-ESM4 and MPI-ESM1-2-LR from smallest to largest absolute mean bias delta_SST.
11. Produce Figure 1 as a 2 by 2 panel on the common comparison cells: bias maps for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR using one shared diverging colour scale centred on zero, and the WOA23 0 m reference SST using a sequential colour scale. Label all panels and colour scales in degrees Celsius.
12. Produce Figure 2 as a bar chart of the four global mean SST values from steps 3 and 5 for CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23, each over its own native grid, annotate each bar with its value, and label the y axis in degrees Celsius.
13. Use Figure 1 and the comparison statistics to describe the spatial structure of the bias, including sign, regional coherence and any large-scale zonal or regional patterns.

Other notes:
- Positive delta_SST means the model is warmer than WOA23; negative means the model is colder.
- Native-grid global means use each model's own areacello and are not directly cell-comparable across models until regridded.
- The WOA23 global mean is over the whole 1-degree grid using valid WOA23 cells; comparison means, biases, RMSE and correlations are only over the common comparison cells.
- CMIP6 monthly fields are reduced over 1971-2000 to a long-term climatology, while WOA23 is the 1971-2000 annual 0 m climate normal.
- Regridding to the WOA23 grid introduces interpolation, and the common comparison cells may omit parts of each model's native ocean domain.
- The model ranking uses absolute mean bias delta_SST, not RMSE or spatial correlation.
