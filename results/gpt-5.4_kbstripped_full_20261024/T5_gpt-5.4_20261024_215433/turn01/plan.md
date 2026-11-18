# Plan

1. Load WOA23 salinity climatology and metadata needed to define observed SSS on the WOA23 grid. Datasets: woa23. Variables: salinity (1971-2000 climate normal, annual), depth levels. Tools: load_climatology_data, load_depth_levels, load_landsea_mask, load_basin_mask, load_areas.
   - Use load_climatology_data for WOA23 annual salinity on the full globe, load_depth_levels to identify the WOA23 levels within 0–10 m (0, 5, 10 m), load_landsea_mask for landsea_01 and load_basin_mask for basinmask_01, and load_areas for WOA23 cell areas.
   - Compute observed SSS as the volume-weighted mean over WOA23 levels whose depths lie within 0–10 m, with full-cell-volume weights as requested; keep values in psu.
   - Restrict ocean cells using the WOA23 land-sea mask condition ocean > 1.

2. Load CMIP6 historical and ssp245 salinity and fixed fields for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, and define each model’s SSS consistently with the requested 0–10 m full-cell-volume weighting. Datasets: cmip6_historical, cmip6_ssp245. Variables: so, areacello, volcello, depth levels. Models: CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Tools: load_regular_data, load_volumes_for_cmip6, load_depth_levels.
   - Load monthly historical so for 1950-01 to 2014-12 and monthly ssp245 so for 2071-01 to 2100-12 for each model.
   - Load each model’s volcello and areacello, and use load_depth_levels to identify the model levels whose depths lie within 0–10 m; compute SSS as the full-cell-volume-weighted mean over those one or two upper levels, in psu.
   - Note [historical-ssp245-member-continuity]: treat GFDL-ESM4 and MPI-ESM1-2-LR historical r1 + ssp245 r1 as continuous across experiments for interpretation, but do not join CESM2 historical and ssp245 across the experiment boundary because the local files use different members.

3. Compute 1971-2000 climatological mean SSS fields and global-mean SSS on each source’s own grid. Datasets: woa23, cmip6_historical. Variables: SSS derived in steps 1–2, areacello/WOA23 areas. Models: CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Tools: temporal_aggregate, compute_basic_stat.
   - For each model, aggregate monthly historical SSS to the 1971-2000 long-term mean using day-weighted temporal_aggregate.
   - Compute the global mean 1971-2000 SSS for each model on its own grid as an area-weighted surface mean with compute_basic_stat, and compute the WOA23 global mean on the WOA23 grid the same way.
   - Report each model’s global-mean difference from WOA23 using the model mean minus the WOA23 mean.
   - Note [omip-weighted-means-and-integrals]: use area-weighted surface means over valid top-field ocean cells, with weights restricted to the same valid cells as the SSS data.

4. Put model and observed 1971-2000 SSS climatologies on the WOA23 1° grid, build one shared comparison mask per model, and compute map-based mean-state metrics and bias fields for Figure 1. Datasets: woa23, cmip6_historical. Variables: 1971-2000 SSS climatologies. Models: CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Tools: regrid_horizontal, build_comparison_mask, compute_difference_metrics, analyze_correlation.
   - Regrid each model’s 1971-2000 SSS climatology horizontally to the WOA23 1° longitude-latitude grid for cell-by-cell comparison with WOA23.
   - Build the shared comparison mask from the regridded model field and WOA23 field, apply it to both, and then compute RMSE and the cellwise bias field (model minus WOA23) on that mask.
   - Compute the spatial Pearson correlation between each masked regridded model field and masked WOA23 field with area weighting.
   - Save the masked bias maps for the three Figure 1 panels, all using one shared diverging colour scale centred on zero.
   - Note [common-comparison-mask]: on the WOA23 1° grid, build one ocean-only shared mask from the fields being compared for each model and apply it before metrics, and report the mask cell counts with the metrics.
   - Note [omip-native-grid-vs-regridded]: do means and changes on each source’s own grid, and use the WOA23-grid regridded products only for the cell-by-cell comparison metrics and Figure 1 maps.
   - Note [cmip6-grid-labels]: treat GFDL-ESM4 calculations on its archive `gr` grid as that model’s own-grid results, while CESM2 and MPI-ESM1-2-LR are on `gn`; only the regridded comparison/display products are on a common lon-lat grid.

5. Compute annual global-mean SSS time series for historical 1950-2014 and ssp245 2071-2100, fit trends, and derive slope ratios and Figure 3 fitted segments. Datasets: cmip6_historical, cmip6_ssp245. Variables: monthly SSS, areacello. Models: CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Tools: compute_basic_stat, temporal_aggregate, fit_1D_trend.
   - For each model and experiment separately, compute the monthly global mean SSS on the model’s own grid as an area-weighted surface mean, then aggregate to annual means.
   - Fit a linear trend to each annual series: historical 1950-2014 and ssp245 2071-2100, and report slope in psu per decade, 95% confidence interval, R², and whether the interval excludes zero.
   - Compute the ratio of the ssp245 slope to the historical slope for each model, and state whether that ratio is meaningful in light of [historical-ssp245-member-continuity].
   - Save the annual global-mean series and fitted lines for Figure 3, plotting historical and ssp245 as separate segments with the 2015-2070 gap left empty.
   - Note [serial-correlation-and-significance]: for annual global-mean SSS trends, report slope with 95% confidence interval and R² from explicit linear fits, and call a trend significant only if the interval excludes zero using the tool’s default effective-sample-size correction.

6. Compute end-of-century mean change in SSS relative to late historical, globally and by Atlantic/Indian/Pacific basin on each model’s own grid where possible, and on the WOA23 basin definition requested. Datasets: cmip6_historical, cmip6_ssp245, woa23. Variables: SSS, areacello, basinmask_01. Models: CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Tools: temporal_aggregate, regrid_horizontal, select_basin, compute_basic_stat.
   - Compute each model’s historical 1995-2014 mean SSS and ssp245 2071-2100 mean SSS on its own grid, then form delta SSS = future mean minus historical mean.
   - Compute each model’s global-ocean mean delta SSS on its own grid as an area-weighted surface mean.
   - Regrid each model’s historical-mean and future-mean SSS fields to the WOA23 1° grid so the requested WOA23 basin mask can define Atlantic, Indian and Pacific consistently; then apply select_basin with the WOA23 basin mask and compute area-weighted basin means of delta SSS for Atlantic, Indian and Pacific, excluding Southern and Arctic.
   - Report the basin changes together with the experiment-member caveat required by [historical-ssp245-member-continuity].
   - Note [basin-and-region-definitions]: use the WOA23 basin mask with marginal seas included in Atlantic/Indian/Pacific, and exclude Southern and Arctic from the requested three-basin results.

7. Fit cell-by-cell historical SSS trends for Figure 2 and prepare display fields on a common longitude-latitude frame. Datasets: cmip6_historical. Variables: annual-mean SSS. Models: CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Tools: temporal_aggregate, fit_trend_map, regrid_horizontal.
   - For each model, aggregate monthly SSS to annual means over 1950-2014, fit a linear trend at every grid cell, and retain both slope in psu per decade and significance at p < 0.05.
   - Regrid the trend and significance products horizontally to a regular longitude-latitude grid only for display in Figure 2, preserving that the calculation itself was done on each model’s own grid.
   - Prepare three map panels on one shared diverging colour scale centred on zero, with significant cells stippled over the trend field.
   - Note [serial-correlation-and-significance]: for Figure 2, draw slopes everywhere and stipple cells significant at p < 0.05 using corrected significance rather than naïve yearly sample size.

Other notes:
- [historical-ssp245-member-continuity] For CESM2, the end-of-century minus historical change and the ssp245-to-historical slope ratio mix forced change with differences in internal variability between different members, so that ratio is not meaningfully interpreted as one member’s acceleration across experiments.
