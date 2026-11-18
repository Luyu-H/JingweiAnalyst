# Plan

1. Load the CMIP6 historical monthly `tos` fields for models CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 1971-01 to 2000-12 with `load_regular_data`, and load each model’s published `areacello` with `load_regular_data` from `cmip6_historical`.
- Note [cmip6-static-grid-fields]: use each model’s published `areacello` as the area weights for its own-grid global mean.
- Note [cmip6-tos-bulk-top-cell]: treat `tos` as the model SST to compare directly with WOA23 0 m, and state it is the bulk temperature of the top model cell.

2. Form each model’s 1971-2000 SST climatology by aggregating monthly `tos` with `temporal_aggregate` to a day-weighted long-term mean, then compute each model’s global mean SST on its own grid with `compute_basic_stat` as an area-weighted mean over valid ocean cells only.
- Note [omip-native-grid-vs-regridded]: these first global means must be computed on each source’s own grid with its own weights; regridding is only for the cell-by-cell comparison and maps.
- Note [missing-values-and-land-cells]: build validity from each SST field itself and ensure weights are masked to the same valid cells so land/missing cells do not enter the denominator.
- Note [cmip6-grid-labels]: state that GFDL-ESM4 `tos` is archived on `grid_label=gr`, so its “own-grid” mean is on that archive grid.
- Note [omip-weighted-means-and-integrals]: use area-weighted means, not plain latitude-longitude averages.

3. Load the WOA23 1971-2000 annual `temperature` climatology at 0 m depth on the 1° grid with `load_climatology_data`, and load the WOA23 1° land-sea mask with `load_landsea_mask`; apply the rule ocean where `landsea_01 > 1`, then compute the area-weighted global mean of the WOA23 0 m annual field over valid ocean cells on the 1° grid using `load_areas` and `compute_basic_stat`.
- Note [common-comparison-mask]: because WOA23 is part of the comparison, keep ocean-only treatment and report mask cell counts with the comparison metrics; use means on the shared sample for the comparison-stage means.
  
4. Regrid each model’s 1971-2000 SST climatology from its own grid to the WOA23 1° grid with `regrid_horizontal` using method `conservative`, then build one shared comparison mask across the three regridded model climatologies and WOA23 with `build_comparison_mask`; apply that mask to every regridded model field and to WOA23 before all comparison metrics and Figure 1.
- In later comparison steps, use the shared mask from [common-comparison-mask] so all models are scored on the same cells.
- In later reporting, identify that the common-grid results depend on the conservative regridding from [omip-native-grid-vs-regridded].

5. On the shared WOA23 1° comparison grid, compute for each masked regridded model and masked WOA23 the area-weighted mean SST with `compute_basic_stat`, and compute for each model against WOA23 the mean bias and RMSE with `compute_difference_metrics`.
- Use the shared masked sample from step 4 [common-comparison-mask].
- Use area weighting for these map-based statistics per [omip-weighted-means-and-integrals].

6. On the same shared masked 1° grid, compute the spatial Pearson correlation between each regridded model climatology and WOA23 with `analyze_correlation`, using area weighting.
- Use `weighting='area'` per [omip-weighted-means-and-integrals].

7. Rank CESM2, GFDL-ESM4 and MPI-ESM1-2-LR from smallest to largest absolute mean bias `delta_SST` using the step-5 shared-grid mean biases.

8. Prepare Figure 1 from the masked shared-grid fields of step 4: three global maps of bias (`model - WOA23`) for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, plus one global map of masked WOA23 0 m SST; use one shared diverging colour scale centred on zero for the three bias panels and one sequential colour scale for WOA23, and label every panel and colour bar in °C.

9. Prepare Figure 2 as a bar chart of the four own-grid global mean SST values from steps 2 and 3 (CESM2, GFDL-ESM4, MPI-ESM1-2-LR, WOA23), with each bar annotated by its value and the y-axis labelled in °C.

Other notes:
- Report, for each model comparison on the shared 1° grid, the masked cell count and the masked mean SST of both the model and WOA23 used by `compute_difference_metrics`, not means from different cell sets.
