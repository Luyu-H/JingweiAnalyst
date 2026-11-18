# Plan

1. Load CMIP6 historical monthly 1971-01 to 2000-12 data for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: `thetao`, `so` and `mlotst` via `load_regular_data`; load each model's published `areacello` via `load_areas`.
2. For each model, compute potential density anomaly from `thetao` and `so` with `compute_potential_density_anomaly` using TEOS-10, then compute the operator MLD monthly with `compute_mixed_layer_depth` using the density criterion defining `mlotst` (0.03 kg m^-3 potential-density increase relative to the reference/surface level), interpolating the crossing depth between bracketing levels.
3. For each model and separately for March and September, build the common valid-cell mask with `build_comparison_mask`, requiring both the operator MLD and `mlotst` to be defined in all 30 years of that calendar month (1971-2000); apply this per-calendar-month mask to both monthly fields before averaging.
4. Aggregate the masked operator MLD and masked `mlotst` to climatological-monthly values with `temporal_aggregate`, then select the March and September climatologies for each model.
5. For each model-month, compute operator MLD minus `mlotst`; use `compute_difference_metrics` for mean difference (bias), RMSE, and the largest single departure with its location, and use `analyze_correlation` for the spatial Pearson correlation. Use `compute_exceedance_fraction` on |difference| > 20 m, by area, over the common valid cells.
6. For each model-month, use `compute_basic_stat` on the binary common-valid mask over the model ocean cells with `areacello` area weights and the mean operator, yielding the valid area as a fraction of the model's total ocean area; report as a percentage.
7. Rank CESM2, GFDL-ESM4 and MPI-ESM1-2-LR from most to least internally consistent using the average of the March and September RMSE as the primary criterion (lower RMSE = more consistent), reporting bias and correlation alongside and noting where the largest differences and >20 m exceedances concentrate.
8. For Figure 1, regrid each native-grid difference map to a common regular longitude-latitude grid for display only using `regrid_horizontal`, then plot the 3-by-2 panel with one shared diverging colour scale centred on zero. For Figure 2, scatter operator MLD (x) against `mlotst` (y) for all common-valid ocean cells in each model-month, using identical logarithmic axes and ticks, a dashed y = x line, and annotated bias, RMSE and correlation.

Other notes:
- The operator MLD is an offline diagnostic from monthly `thetao` and `so`, while `mlotst` is the model's own native mixed-layer diagnostic; exact agreement is not expected because of differences in reference level, threshold implementation, vertical interpolation/resolution, and sea-ice or missing-data handling.
- The March and September climatologies are restricted to cells where both fields are defined in all 30 corresponding months; this validity mask differs by calendar month and by model.
- Area percentages use each model's published `areacello` and its own ocean mask; the denominator for the valid-area percentage is the model's total ocean area.
- Regridding is used only for the map display; all statistics are computed on the native model grids and common valid masks.
- Intermediate values are kept unrounded until final reporting at the requested precision.
