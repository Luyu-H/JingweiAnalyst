# Plan

1. Load WOA23 nitrate as the all-time-periods annual objectively analysed climatology on the 1-degree grid with `load_climatology_data`; load its depth levels with `load_depth_levels`, its land-sea mask `landsea_01` with `load_landsea_mask` (ocean where value > 1), and its basin mask `basinmask_01` with `load_basin_mask`. Load CMIP6 historical `no3`, `thetao` and `so` monthly for 1971-01 to 2000-12 for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR with `load_regular_data`; load each model’s native `areacello`/grid-cell areas with `load_areas` and depth levels with `load_depth_levels`.

2. For each CMIP6 model, aggregate `no3`, `thetao` and `so` over 1971-2000 into a day-weighted long-term mean using `temporal_aggregate` with each model’s calendar, before computing any diagnostic.
   - Note [monthly-means-are-time-averages]: use `temporal_aggregate` with day weights rather than an unweighted time mean; average the monthly fields first, then diagnose the nitracline, and state that this is not the time mean of monthly nitraclines.

3. Compute nitracline depth with `compute_nutricline_depth` for WOA23 and for each model from the annual/time-mean nitrate field on each source’s native vertical grid; search only the levels within the top 1000 m, but include the level just below 1000 m when evaluating the deepest gradient.
   - Note [cmip6-bgc-units-per-volume]: take the nitracline depth from the native `no3` field without density conversion because a gradient-extremum depth is scale-invariant; for the concentration section only, convert CMIP6 `no3` from mol m-3 to µmol kg-1 using in-situ density from `compute_potential_density_anomaly` with `ref_pressure=None`, and report `attrs['density_assumption']`.

4. Compute global mean nitracline depths: mask WOA23 with `landsea_01 > 1`; for each source use `compute_basic_stat` with `op='area_weighted_mean'` over the horizontal dimensions and that source’s own area weights, using NaN-aware valid cells. Form each model-minus-WOA23 difference.
   - Note [omip-weighted-means-and-integrals]: global and basin means must be area-weighted with the source’s own area weights; map RMSE/correlation must also be area-weighted, with `compute_difference_metrics` using its default area weighting and `analyze_correlation` using `weighting='area'`.
   - Note [omip-native-grid-vs-regridded]: take global and basin means on each source’s own grid with its own `areacello`; regrid to the WOA 1-degree grid only for the maps, the pattern correlation and the RMSE on the shared mask, using `regrid_horizontal` with `'conservative'` when going coarser; state this method and the native-vertical-grid quantisation caveat.

5. Compute basin mean nitracline depths using the WOA23 basin mask through `load_basin_mask`/`select_basin` for the Pacific, Atlantic and Indian Oceans for each source, then use `compute_basic_stat(op='area_weighted_mean')` with each source’s own area weights; exclude the Southern and Arctic Oceans and form each model-minus-WOA23 basin difference.
   - Note [basin-and-region-definitions]: use the WOA23 `basinmask_01` with `include_marginal_seas=True` and `depth_dependent=False` so marginal seas fold into their parent ocean and surface labels apply through the water column; state that this is the WOA labelling, not the OMIP 10-region one, and that the Southern and Arctic are separate basins.

6. Regrid each model nitracline map to the WOA 1-degree grid with `regrid_horizontal` using `'conservative'`; build one common mask over all four regridded nitracline maps with `build_comparison_mask` and apply it before the map comparison and Figure 1, recording the mask cell count.
   - Note [common-comparison-mask]: apply one shared comparison mask to every field before correlation/RMSE; the default `ocean_only=True` is correct because WOA23 is the reference, and report the cell count recorded by the mask.
   - Note [cmip6-grid-labels]: GFDL-ESM4 is already 1-degree `gr` output, while CESM2 is about 100 km and MPI-ESM1-2-LR about 250 km, so features narrower than the coarse grids are not comparable in amplitude.

7. On the shared masked WOA 1-degree grid, compute each model’s spatial Pearson correlation with WOA23 using `analyze_correlation(weighting='area')`, and its RMSE using `compute_difference_metrics`; report correlation to 3 decimal places, RMSE to 1 decimal place, and the common-mask cell count, referencing [common-comparison-mask] and [omip-weighted-means-and-integrals].

8. Produce Figure 1 as a 2-by-2 panel of global nitracline-depth maps for WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR on the shared WOA 1-degree grid and common mask, using one sequential colour scale on which deeper is visually heavier.

9. Produce Figure 2 from nitrate in the equatorial Pacific over 2S-2N and 120E-80W, top 1000 m plus the level below 1000 m needed for the nitracline gradient. For WOA23 use nitrate directly in µmol kg-1; for each model convert `no3` to µmol kg-1 using `transform_unit` for mol m-3 to µmol m-3 and divide by in-situ density from `compute_potential_density_anomaly(ref_pressure=None)` based on the time-mean `thetao`/`so`. Average over 2S-2N with `compute_basic_stat(dim='lat', op='area_weighted_mean')`; compute each panel’s nitracline line on its native vertical grid with `compute_nutricline_depth`; use `regrid_vertical` to WOA23 depth levels for a common plotted depth axis if needed while keeping the overlaid nitracline from the native diagnostic. Plot panels in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with depth increasing downward and one shared colour scale in µmol kg-1.

10. Compare each model with WOA23 using the global and basin mean differences, the spatial correlation/RMSE, and the ordering of mean nitracline depth across the Pacific, Atlantic and Indian basins; identify the best model and state whether each model reproduces the observed basin ordering.

Other notes:
- Report global and basin nitracline depths to the nearest metre, differences as model minus WOA23, correlations to 3 decimal places, and RMSE to 1 decimal place.
- State that the nitracline is the depth of the maximum vertical nitrate gradient within the top 1000 m, with the gradient at the deepest searched level still using the level just below 1000 m.
- State the panel order and shared colour scale for both figures, and that Figure 2 uses depth increasing downward.
