# Plan

1. Load the WOA23 reference: `load_climatology_data` for woa23 salinity, 1971-2000 normal, annual, levels 0, 5 and 10 m. Load the WOA23 landsea_01 mask (`load_landsea_mask`, ocean is value > 1) and basinmask_01 (`load_basin_mask`). Load the WOA23 grid cell areas (`load_areas`) if needed for weights. Compute the WOA23 SSS per cell as the mean of the 0, 5 and 10 m levels, each level weighted by its full cell volume (`compute_basic_stat`, volume-weighted form).

2. Load the model data with `load_regular_data`:
   - cmip6_historical `so` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, 1950-01 to 2014-12.
   - cmip6_ssp245 `so` for the same models, 2071-01 to 2100-12.
   - Each model's own `areacello` (`load_areas`) and `volcello` (`load_volumes_for_cmip6`). Use `load_depth_levels` to pick the levels with depths within 0-10 m, which is one or two levels per model.
   - Compute model SSS per cell as the mean of those levels, weighted by each level's full volcello (`compute_basic_stat`, volume-weighted). Do the same for all datasets.
   - Note [cmip6-static-grid-fields]: use volcello/areacello, not lev_bnds differences. Report `attrs['area_source']` and `attrs['volume_source']`. Handle the time axis of GFDL volcello. Do not rescale areas between models.

3. Mean-state bias, 1971-2000:
   - Take the model historical climatology for 1971-2000 with `temporal_aggregate` (day-weighted long-term mean).
   - Regrid the three models to the WOA23 1° grid with `regrid_horizontal`.
   - Build one mask over the three models and WOA23 with `build_comparison_mask` (default ocean_only=True), and apply it with `.where(mask)`.
   - On the shared sample, run `compute_difference_metrics` (area weighted) for bias, RMSE, `mean_data` and `mean_reference`. Run `analyze_correlation` with weighting='area' for the pattern correlation.
   - Report the global mean SSS of each model and of WOA23, and each model's difference from WOA23, to 2 decimals. Report RMSE to 2 decimals and Pearson r to 3 decimals.
   - Note [common-comparison-mask]: rank the models only on the shared mask. Report the cell counts the mask records. State that the reported means are shared-sample means, not each source's own global ocean mean. Do not apply this mask to trends or ΔSSS.
   - Figure 1: three maps in one row, with bias = model minus WOA23 on the shared-mask cells. Use one diverging colour scale centred on zero, shared by all panels, with the colour bar labelled in psu. Use a fixed colour per model.

4. Global-mean trends:
   - Compute the annual global mean SSS of each model with `temporal_aggregate` (annual) and `compute_basic_stat` (area_weighted_mean over valid cells).
   - Fit `fit_1D_trend` separately on historical 1950-2014 and on ssp245 2071-2100. Report the slope in psu/decade (3 significant figures), the 95% CI, R², and whether the CI excludes zero.
   - Compute the ratio of ssp245 slope to historical slope (1 decimal) and state whether it is meaningful.
   - Note [historical-ssp245-member-continuity]: fit each segment separately and never across 2014/2015. The ratio is not a continuous-trajectory quantity. For CESM2 never fit or draw across the join.
   - Note [anomaly-baselines-and-climate-indices]: trends use annual means of each model's own series, not anomalies against WOA23.

5. Cell-by-cell trend map: for each model, use `fit_trend_map` on annual-mean SSS over 1950-2014, with significance at p < 0.05. Put curvilinear fields on a regular lon-lat grid for display only (`regrid_horizontal`). Do not apply the comparison mask here.
   - Figure 2: three geographic maps in one row, with one shared diverging scale centred on zero, psu/decade, and significant cells stippled.

6. Figure 3: a single panel with the global-mean annual SSS of the three models. Draw the historical 1950-2014 and ssp245 2071-2100 segments separately, with the 2015-2070 gap empty. Overlay each segment's own fitted line. Label the y axis in psu and keep each model's colour consistent with Figures 1 and 2.

7. ΔSSS:
   - Compute the mean of ssp245 2071-2100 minus the mean of historical 1995-2014 per model (`temporal_aggregate`, day-weighted means, then difference).
   - Report it for the global ocean and for the Atlantic, Indian and Pacific basins, to 2 decimals.
   - Apply the WOA23 basinmask_01 with `select_basin` (include_marginal_seas=True; Southern and Arctic excluded) to each model. Take weighted means with `compute_basic_stat` (area_weighted_mean over valid cells).
   - Note [basin-and-region-definitions]: state this mask and these choices. Apply the same mask to models and observations, using surface labels. Marginal seas missing from a model have no counterpart.
   - Note [omip-weighted-means-and-integrals]: use NaN-aware weights, not plain mean. The volume-weighted form is for 0-10 m levels.

Other notes:
- [common-comparison-mask] Bias, RMSE and r are on the shared 1° ocean-only mask. Means are not each source's own global means.
- [historical-ssp245-member-continuity] CESM2 historical is r1i1p1f1 but ssp245 is r10i1p1f1, whose parent r10 is absent locally. Its ΔSSS therefore includes forced change plus the internal-variability difference between two members. GFDL-ESM4 and MPI-ESM1-2-LR use r1 for both, so their join is continuous.
- [anomaly-baselines-and-climate-indices] Baselines are the model climatology 1971-2000 and the WOA23 1971-2000 normal. The bias keeps each model's mean-state error.
- [cmip6-grid-labels] GFDL-ESM4 is already 1° (`gr`) with a regridded mask. CESM2 and MPI are `gn` (100 km and 250 km), so coastal and boundary-current features differ in amplitude between models. Moving vertical coordinates mean time-dependent volume weights.
