# Plan

1. Load the data with `load_regular_data` from cmip6_historical, for each of CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. Variables are `thetao`, `so` and `mlotst`, monthly, 1971-01 to 2000-12, global. Each model's own `areacello` is needed for the area weights (the question names it explicitly).
   - Note [cmip6-grid-labels]: GFDL-ESM4 is on the regridded `gr` grid, and CESM2 and MPI-ESM1-2-LR are on `gn`. Compare each model on its own grid with its own `areacello`. Use `regrid_horizontal` only to put curvilinear fields on a regular lon-lat grid for the Figure 1 maps, not for any metric.
   - Note [missing-values-and-land-cells]: Build each field's valid mask from the field itself, because the `thetao`, `so` and `mlotst` masks can differ within a model. Use NaN-aware reductions. A value near 1e20 means an unmasked fill value.

2. Use `load_depth_levels` for each model to get its shallowest level, which is the `ref_depth`. Use `load_areas` to obtain `areacello` for the area weights and for the total ocean area. Use `load_volumes_for_cmip6` only if a cell-volume field is needed for masking.

3. Compute the operator MLD for every monthly March and September profile (30 of each per model), using `compute_mixed_layer_depth` on `thetao` and `so`. The criterion is a density increase of 0.03 kg m-3, with the crossing depth interpolated between the bracketing levels.
   - Note [cmip6-mlotst-definition]: Use 0.03 kg m-3 and set `ref_depth` to each model's shallowest level, not the 10 m default. Compute the MLD of each monthly profile first, then average the MLDs over the 30 Marches or Septembers. Never average T and S first. Exclude columns where no crossing is found and report how much area this removes. The models encode no-crossing columns differently (GFDL-ESM4 and MPI give roughly the full column depth, CESM2 gives 0 m in many Arctic-shelf cells), so do not use a single fill value.
   - Note [cmip6-thetao-potential-temperature]: `thetao` is potential temperature, so use the potential-temperature route for density.
   - Note [cmip6-so-practical-salinity]: `so` is Practical Salinity, so do not rescale it. Convert S_P to S_A before any gsw routine that expects S_A.
   - Note [monthly-means-are-time-averages]: The operator MLD is computed from monthly-mean states, while `mlotst` is the time mean of instantaneous diagnostics, so say which order was used. Do the climatologies with `temporal_aggregate` (climatological-monthly, a plain mean over years), not with `compute_basic_stat(op='mean')`.

4. Build the 1971-2000 March and September climatologies of the operator MLD and of `mlotst`, per model, with `temporal_aggregate`.
   - Note [common-comparison-mask]: Use `build_comparison_mask` with `ocean_only=False` to select the cells where both fields are valid in all 30 instances of that month. Apply it to both fields with `.where(mask)` before any metric. Report the mask cell count, the mask area, and its share of the model's total ocean area.

5. Compute the metrics for each of the six model-month combinations on the shared mask.
   - Note [omip-weighted-means-and-integrals]: Use `compute_difference_metrics` (operator MLD minus `mlotst`, area-weighted by the model's own `areacello`) for the mean difference and the RMSE. Report `mean_data`/`mean_reference` and the bias. Use `analyze_correlation` with `weighting='area'` for the spatial Pearson correlation. The coverage is a percentage of the model's total ocean area.
   - Note [missing-values-and-land-cells]: Use `compute_exceedance_fraction` with an area basis and a threshold of |difference| > 20 m for the percentage of shared area. The denominator is the area where both fields are defined.
   - Report the mean difference and RMSE in m to 1 dp, the correlation to 3 dp, and both percentages to 1 dp.

6. Rank the models from most to least internally consistent. Use RMSE (with correlation and the >20 m area share as supporting criteria), and state the criterion. Locate where the disagreement concentrates, using `locate_extreme` and the difference maps (for example deep-convection and sea-ice regions such as the Southern Ocean, the North Atlantic and the Arctic shelves). Use `select_basin` with `load_basin_mask` if a basin breakdown is needed.

7. Make Figure 1, a 3×2 panel of global maps of the difference (rows CESM2, GFDL-ESM4, MPI-ESM1-2-LR, columns March and September). Use a common lon-lat frame, with curvilinear fields regridded for display only, and one shared diverging colour scale centred on zero.

8. Make Figure 2, a 3×2 scatter panel, with operator MLD on x and `mlotst` on y, one point per shared-mask ocean cell. Use log axes with identical range and ticks in every panel, a dashed y = x line, and an annotation of bias, RMSE and correlation in each panel.

Other notes:
- State that the operator MLD was computed from monthly-mean thetao and so, whereas `mlotst` is the time mean of instantaneous diagnostics, and that the operator is nonlinear, so exact agreement is not expected. Other reasons for disagreement are the vertical interpolation and the reference depth, the equation-of-state differences, and the encoding of no-crossing columns.
- Report the area removed by excluding no-crossing columns for each model.
- Report the mask size and coverage for each model-month, since a small mask means the products disagree about where data are defined.
