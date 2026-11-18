# Plan

1. Load the WOA23 observed nitrate climatology and supporting masks/levels needed to derive the observed nitracline depth map.  
   - Datasets: woa23  
   - Variables/models: nitrate annual climatology (all-period composite), depth levels and layer bounds, landsea_01, basinmask_01  
   - Tools: load_climatology_data, load_depth_levels, load_landsea_mask, load_basin_mask, compute_vertical_gradient, compute_nutricline_depth  
   - Method: load the annual WOA23 nitrate profile on the 1° grid to at least the first level below 1000 m; use the WOA23 land-sea mask with ocean defined by landsea_01 > 1; compute the vertical nitrate gradient and diagnose nitracline depth as the depth of maximum vertical gradient searched over levels within the top 1000 m, with the deepest searched gradient allowed to use the level just below 1000 m.  
   - Note [missing-values-and-land-cells]: define valid WOA23 ocean cells from landsea_01 > 1 and keep NaN-aware masking consistent with the nitracline field itself.

2. Load and time-average the three CMIP6 historical nitrate fields, then convert them to the same concentration basis as WOA23 before deriving each model’s nitracline depth on its own grid.  
   - Datasets: cmip6_historical  
   - Variables/models: no3 for CESM2, GFDL-ESM4, MPI-ESM1-2-LR; monthly 1971-01 to 2000-12; each model’s native/archive grid  
   - Tools: load_regular_data, temporal_aggregate, load_depth_levels, compute_potential_density_anomaly, transform_unit, compute_vertical_gradient, compute_nutricline_depth  
   - Method: load monthly no3 for each model over 1971-01 to 2000-12 and form the long-term mean climatological profile before nitracline diagnosis; load enough depth levels to include the first level below 1000 m; convert CMIP6 no3 from per-volume to the same per-mass unit as WOA23 using in-situ density, then compute the vertical nitrate gradient and diagnose nitracline depth over the top 1000 m on each model grid.  
   - Note [cmip6-bgc-units-per-volume]: convert model no3 to a shared concentration unit with WOA23 before sections or vertical gradients, using density and recording the density assumption in metadata.

3. Compute global mean nitracline depth for WOA23 and each model on each source’s own grid with proper horizontal area weighting, then form model-minus-WOA23 differences.  
   - Datasets: woa23, cmip6_historical  
   - Variables/models: nitracline depth maps from Step 1 and Step 2; CESM2, GFDL-ESM4, MPI-ESM1-2-LR  
   - Tools: load_areas, compute_basic_stat  
   - Method: use each source’s own horizontal cell areas/weights to compute area-weighted global mean nitracline depth from its nitracline map, masking weights to valid nitracline cells; report each mean and each model’s difference from the WOA23 global mean.  
   - Note [omip-weighted-means-and-integrals]: use area-weighted horizontal means for global and basin nitracline depths, not unweighted means.

4. Compute Pacific, Atlantic and Indian mean nitracline depths for WOA23 and each model using the WOA23 basin definition on each grid, then form model-minus-WOA23 differences by basin.  
   - Datasets: woa23, cmip6_historical  
   - Variables/models: nitracline depth maps from Step 1 and Step 2; basinmask_01; CESM2, GFDL-ESM4, MPI-ESM1-2-LR  
   - Tools: select_basin, compute_basic_stat  
   - Method: apply the WOA23 basin mask to each source grid with Pacific, Atlantic and Indian selected separately; include marginal seas in the parent basin and exclude the Southern and Arctic Oceans; compute area-weighted mean nitracline depth within each selected basin on each source’s own grid; report basin means and each model’s difference from the corresponding WOA23 basin mean.  
   - Note [basin-and-region-definitions]: use the WOA23 basin mask with include_marginal_seas=True and depth_dependent=False so Pacific/Atlantic/Indian means exclude Southern and Arctic while surface basin labels apply through the full water column.

5. Regrid all three model nitracline depth maps to the WOA23 1° grid, build one shared comparison mask across WOA23 and all regridded models, and compute map skill metrics against WOA23 for each model.  
   - Datasets: woa23, cmip6_historical  
   - Variables/models: nitracline depth maps from Step 1 and Step 2; CESM2, GFDL-ESM4, MPI-ESM1-2-LR  
   - Tools: regrid_horizontal, build_comparison_mask, analyze_correlation, compute_difference_metrics  
   - Method: horizontally regrid each model nitracline map to the WOA23 1° grid; build one shared valid-data mask across WOA23 and all three regridded model maps; apply the common mask to WOA23 and each regridded model before computing area-weighted spatial Pearson correlation and RMSE for each model against WOA23; report the common-mask cell count used for the scored maps.  
   - Note [common-comparison-mask]: use one shared mask over WOA23 and all three regridded models for the cell-by-cell map metrics only, and report the mask cell counts with the metrics.

6. Identify which model best reproduces observed nitracline depth overall and whether each model reproduces the observed basin ordering.  
   - Datasets: woa23, cmip6_historical  
   - Variables/models: global and basin mean nitracline depths from Steps 3–4; spatial correlation and RMSE from Step 5; CESM2, GFDL-ESM4, MPI-ESM1-2-LR  
   - Tools: compute_difference_metrics  
   - Method: rank models primarily from the WOA23 map-comparison metrics on the shared mask and summarize consistency with the global and basin mean differences; compare each model’s Pacific/Atlantic/Indian basin-mean ordering with the WOA23 ordering and state whether it matches.

7. Prepare Figure 1: a 2×2 panel of global nitracline depth maps on one shared sequential colour scale.  
   - Datasets: woa23, cmip6_historical  
   - Variables/models: WOA23 nitracline depth; regridded nitracline depth maps for CESM2, GFDL-ESM4, MPI-ESM1-2-LR  
   - Tools: regrid_horizontal  
   - Method: use the WOA23 nitracline map and the three model nitracline maps regridded to the WOA23 1° grid; plot all four on one shared sequential scale with deeper depths mapped to visually heavier colours.

8. Prepare Figure 2: a 4-panel equatorial Pacific longitude–depth section of nitrate concentration with each source’s nitracline overlaid.  
   - Datasets: woa23, cmip6_historical  
   - Variables/models: nitrate concentration and nitracline depth for WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR  
   - Tools: load_regular_data, load_climatology_data, temporal_aggregate, transform_unit, compute_potential_density_anomaly, compute_basic_stat  
   - Method: extract nitrate over 120E–80W and 2S–2N for WOA23 and for each model’s 1971–2000 mean; ensure all four sections are in the same concentration unit as in Step 2; compute the latitude-band mean section for each source; restrict plots to the top 1000 m with depth increasing downward; overlay each source’s own nitracline depth line along longitude; use one shared colour scale across all four panels in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR.

Other notes:
- [omip-native-grid-vs-regridded] Compute nitracline means on each source’s own grid and use regridding only for the WOA23 map comparison metrics and Figure 1; note that nitracline depth is quantised by each source’s vertical grid.
- [cmip6-grid-labels] Note in the methods that GFDL-ESM4 no3 is on the provider’s gr archive grid, so its “native-grid” means are on that archive grid.
- [cmip6-static-grid-fields] Use each source’s own horizontal cell areas, with model areacello for CMIP6 curvilinear/archive grids.
