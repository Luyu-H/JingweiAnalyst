# Plan

1. Load WOA23 nitrate for the all-time-periods annual climatology on the 1-degree grid, plus WOA23 `landsea_01`, `basinmask_01`, depth levels and grid-cell areas. Use `load_climatology_data` for nitrate, `load_landsea_mask`, `load_basin_mask`, `load_depth_levels` and `load_areas`.

2. Load CMIP6 historical monthly `no3` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 1971-01 to 2000-12, each on its native grid, together with each model’s native `areacello` and depth levels. Use `load_regular_data`, `load_areas` and `load_depth_levels`.

3. For each CMIP6 model, aggregate monthly `no3` to a 1971-2000 day-weighted long-term mean with `temporal_aggregate`. No time aggregation is needed for the WOA23 annual climatology.

4. Compute the nitracline depth for WOA23 and for each model climatology as the depth of the maximum vertical gradient of nitrate, using `compute_nutricline_depth` on each native profile. Restrict the search to levels within the top 1000 m, retaining the level immediately below 1000 m so the gradient at the deepest top-1000 m level is defined.

5. Regrid each model’s native nitracline-depth map and nitrate climatology horizontally to the WOA23 1-degree grid with `regrid_horizontal`, using each model’s native `areacello` as the spatial weight. Keep WOA23 on its 1-degree grid.

6. Apply the WOA23 land-sea mask to define ocean cells as `landsea_01 > 1`. Build the common comparison mask across WOA23 and all three regridded model fields with `build_comparison_mask`, so all global statistics, correlations, RMSE values and maps use the same valid ocean cells.

7. Compute global mean nitracline depth for WOA23 and each model on the common WOA23 1-degree grid, area-weighted with WOA23 grid-cell areas, using `compute_basic_stat`. Compute each model’s difference from WOA23. Report global means and differences to the nearest whole metre.

8. Compute basin mean nitracline depths for the Pacific, Atlantic and Indian Oceans using the WOA23 `basinmask_01` and `select_basin`, with marginal seas folded into their parent ocean and the Southern and Arctic Oceans excluded. Use `compute_basic_stat` with WOA23 area weighting on the common mask. Compute each model’s difference from WOA23 in each basin and report to the nearest whole metre.

9. On the common WOA23 1-degree comparison mask, compute each model’s spatial Pearson correlation with the WOA23 nitracline depth using `analyze_correlation`, and compute RMSE using `compute_difference_metrics`. Report correlations to 3 decimal places and RMSE in metres to 1 decimal place.

10. Produce Figure 1 as a 2-by-2 panel of global nitracline-depth maps for WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, using one shared sequential colour scale with deeper values visually heavier/darker.

11. For Figure 2, convert WOA23 and model nitrate to a common concentration unit with `transform_unit`. For each source, average the nitrate climatology over 2S-2N, extract the 120E-80W longitude range and the top 1000 m, and plot the four panels in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR with depth increasing downward and one shared colour scale. Compute each panel’s own nitracline depth from its 2S-2N section using `compute_nutricline_depth` and overlay it as a line.

12. Compare the models’ global RMSE and correlation against WOA23 to identify the best reproduction of the observed nitracline, and compare the Pacific/Atlantic/Indian basin mean ordering of each model with the WOA23 ordering.

Other notes:
- Nitracline depth is defined as the depth of the maximum vertical nitrate gradient within the top 1000 m, with the gradient at the deepest top-1000 m level using the level immediately below 1000 m.
- Model values are 1971-2000 day-weighted long-term means of monthly `no3`; WOA23 is the all-time-periods annual nitrate climatology.
- Basin definitions use the WOA23 `basinmask_01`; marginal seas are folded into their parent ocean, and the Southern and Arctic Oceans are excluded.
- Spatial correlations and RMSE use the common WOA23 1-degree ocean mask (`landsea_01 > 1`) and only cells valid in all compared fields.
- Means and differences are area-weighted and rounded only for final reporting: whole metres for mean depths and differences, 3 decimal places for correlations, and 1 decimal place for RMSE.
- Figure 1 uses one shared sequential colour scale with deeper nitracline depths shown heavier/darker; Figure 2 uses one shared nitrate concentration unit and colour scale.
