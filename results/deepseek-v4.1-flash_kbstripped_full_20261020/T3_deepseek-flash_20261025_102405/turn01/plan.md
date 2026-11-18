# Plan

1. Load CMIP6 historical `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR with `load_regular_data`, monthly from 1971-01 to 2000-12, each on its own grid. Load each model’s `areacello` with `load_areas` and record `attrs['area_source']`. Convert each `tos` from K to degrees Celsius with `transform_unit`.
- Note [missing-values-and-land-cells]: treat fill values as NaN; derive each field’s valid-cell mask from the field itself and mask its weight array to those same cells so land remains missing rather than contributing zero; any output value near 1e20 or 9.97e36 is an unmasked fill, not a result.
- Note [cmip6-tos-bulk-top-cell]: `tos` is the top-cell bulk temperature including under sea ice, not a skin value; near-freezing polar values enter the mean, and if WOA23 masks sea ice the masks must be made consistent before comparison. Top-cell thickness differs by model (CESM2 0–10 m, GFDL-ESM4 0–2 m, MPI-ESM1-2-LR 0–12 m), so part of any inter-model mean difference is definitional.
- Note [cmip6-grid-labels]: GFDL-ESM4 `tos` sits on archive grid `gr` remapped to 1° with `conserve_order1`, so its “native” numbers are on the regridded grid with a regridded land–sea mask; CESM2 (~100 km) and MPI (~250 km) differ in resolution by about 2.5, so features narrower than the coarse grid are not amplitude-comparable when describing the bias structure.

2. Aggregate each model’s monthly `tos` to a single 1971–2000 SST climatology on its native grid with `temporal_aggregate`, using the day-weighted long-term mean.

3. Compute each model’s native-grid global mean SST for 1971–2000 with `compute_basic_stat` using `op='area_weighted_mean'` over latitude/longitude and that model’s own `areacello`, with the area weights masked to the valid `tos` cells. Record the ocean area covered by each mean for Figure 2.
- Note [omip-weighted-means-and-integrals]: native-grid global means must be Σ(tos·areacello)/Σ(areacello) over valid ocean cells via `compute_basic_stat` with `op='area_weighted_mean'` and masked weights, never an unweighted mean over y,x; `compute_difference_metrics` is area-weighted by default and `analyze_correlation` needs `weighting='area'`. Report the shared-sample `mean_data` and `mean_reference` it returns, and keep each model’s own-grid global mean separate from its shared-mask model-minus-WOA23 bias.
- Note [cmip6-static-grid-fields]: use each model’s own `areacello` for its native mean (`load_areas` records `attrs['area_source']`); total ocean area differs among models by about a per cent because coastlines are resolved differently, so state the area each mean was taken over and do not rescale one model to another.

4. Load WOA23 `temperature` for the 1971–2000 climate normal, annual time slice, 0 m depth, with `load_climatology_data`. Load WOA23 `landsea_01` with `load_landsea_mask` and define ocean cells as `landsea_01 > 1`. Load 1-degree cell areas with `load_areas`. Compute the WOA23 global mean over the whole 1-degree ocean grid with `compute_basic_stat` using `op='area_weighted_mean'`, with 1-degree areas masked to `landsea_01 > 1` and valid temperature cells; the missing-value handling from [missing-values-and-land-cells] applies here too.

5. Regrid each model’s 1971–2000 climatology horizontally to the WOA23 1-degree grid with `regrid_horizontal`; carry the GFDL-ESM4 grid and resolution caveat from [cmip6-grid-labels].

6. Build one common comparison mask with `build_comparison_mask` over the WOA23 0 m field and the three regridded model climatologies, using `ocean_only=True`, then apply `.where(mask)` to all four fields before any metric, ranking or map. Record the number of comparison cells.
- Note [common-comparison-mask]: put all three `tos` climatologies and WOA23 on the WOA23 1° grid, build one mask over all of them with `build_comparison_mask`, and apply it with `.where(mask)` to every field before any bias, RMSE or correlation; never score or rank models on per-model pairwise masks. `ocean_only=True` is right because WOA23 is the reference; report the mask cell counts next to the metrics, and note that a mask much smaller than the individual fields means the products disagree about where the ocean is.

7. On the shared comparison cells, compute the mean SST of each model and of WOA23 with `compute_basic_stat` using `op='area_weighted_mean'` and WOA23 1-degree cell areas. For each model compute `delta_SST` and RMSE against WOA23 with `compute_difference_metrics` (area-weighted by default), and compute the spatial Pearson correlation with `analyze_correlation` using `weighting='area'`. Retain unrounded values and report the shared-sample means with [omip-weighted-means-and-integrals].

8. Rank the three models from smallest to largest absolute `delta_SST` using the unrounded mean biases from step 7.

9. Produce Figure 1 as a 2-by-2 panel from the shared-mask fields: bias (model minus WOA23) for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR on one shared diverging colour scale centred on zero, and WOA23 reference SST on the same cells with a sequential colour scale. Label every panel and every colourbar in degrees Celsius.

10. Produce Figure 2 as a bar chart of the four native-grid global mean SST values from steps 3 and 4 (CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23), annotate each bar with its value, and label the y axis in degrees Celsius.

Other notes:
- [cmip6-experiments-historical-ssp245] These are coupled historical runs: statistics over the 1971–2000 window are comparable, but individual years and event timing are not, and with one member per model forced change cannot be separated from internal variability, so the bias and ranking results must carry that caveat.
