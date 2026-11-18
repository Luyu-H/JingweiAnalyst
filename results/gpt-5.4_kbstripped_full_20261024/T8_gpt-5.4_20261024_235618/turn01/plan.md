# Plan

1. Load CMIP6 historical monthly fields for models CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: variables `dissic`, `talk`, `ph`, `thetao`, `so` for 2002-01 to 2002-12, and `dissic`, `thetao`, `so` for 1971-01 to 2000-12; also load each model’s fixed-grid `areacello`/`volcello` support fields via the dedicated volume tool. Load GLODAP climatology variables `TCO2`, `TAlk`, `pHtsinsitutp`.
   - Tools: `load_regular_data`, `load_climatology_data`, `load_volumes_for_cmip6`.
   - Note [cmip6-carbonate-system-definitions]: compare CMIP6 `ph` only to GLODAP `pHtsinsitutp`, and note the CESM2 carbonate-constant caveat if interpretation depends on it.

2. Form the model time means needed for analysis: for each model, compute the 2002 annual mean of `dissic`, `talk`, `ph`, `thetao`, `so`, and the 1971-2000 mean of `dissic`, `thetao`, `so`.
   - Tools: `temporal_aggregate`.
   - Use annual and long-term means from the monthly mean fields without extra screening.

3. Convert CMIP6 `dissic` and `talk` from per-volume to per-mass units for both 2002 and 1971-2000 means by first deriving in-situ density from each model’s `thetao` and `so`, then applying the stated conversion to µmol kg-1.
   - Tools: `compute_potential_density_anomaly`, `transform_unit`.
   - Note [cmip6-bgc-units-per-volume]: use in-situ density from `thetao` and `so` with `ref_pressure=None` for the conversion `µmol kg-1 = mmol m-3 × 1000 / ρ`, and carry the density assumption in the reported caveat.

4. Compute the native-grid global 0-50 m means for each source: for each CMIP6 model, calculate the weighted mean over valid 0-50 m ocean cells of 2002 `dissic`, `talk`, and `ph`, and of 1971-2000 `dissic`; for GLODAP, calculate the weighted mean over valid 0-50 m cells of `TCO2`, `TAlk`, and `pHtsinsitutp`.
   - Tools: `load_areas`, `compute_basic_stat`.
   - Note [omip-weighted-means-and-integrals]: for 0-50 m means of `dissic`, `talk`, and `ph`, use each source’s own grid and valid-cell weights, with `volcello` weighting for the 3-D means.
   - Note [missing-values-and-land-cells]: build valid masks from each field itself and ensure weights are masked to the same valid cells so missing and land cells stay excluded from denominators.

5. Compute the requested source-vs-GLODAP scalar summaries from the native-grid 0-50 m means: for each model, difference from GLODAP for DIC, total alkalinity, and pH; for DIC also compute the period change (2002 minus 1971-2000) and the ratio of that change to the model’s 0-50 m DIC difference from GLODAP.
   - Tools: no extra tool beyond arithmetic on the means from step 4.
   - Interpret the ratio in the final answer as evidence about how much of the apparent 2002 bias may reflect model drift/forced change versus mean-state error.

6. Build each product’s 0-50 m mean DIC map on its own grid, then horizontally regrid the three model maps to the GLODAP 1° grid and construct one common comparison mask shared by GLODAP and all three regridded model maps; apply that same mask to every map before scoring and plotting.
   - Tools: `compute_basic_stat`, `regrid_horizontal`, `build_comparison_mask`.
   - Note [common-comparison-mask]: for the 1° map comparison, regrid to the GLODAP grid, build one shared mask over GLODAP plus all three models with `ocean_only=False`, apply it to every map, and report the common cell counts with the metrics.
   - Note [omip-native-grid-vs-regridded]: keep native-grid calculations for source-specific global 0-50 m means and period changes; use the regridded GLODAP-grid fields only for cell-by-cell map metrics and depthwise bias comparisons.

7. On the shared GLODAP-grid masked 0-50 m DIC maps, compute for each model the mean bias, RMSE, and area-weighted spatial correlation against GLODAP.
   - Tools: `compute_difference_metrics`, `analyze_correlation`.
   - Use area weighting for the correlation and the common masked sample from step 6 for all three model rankings and map statistics.

8. Prepare the depthwise DIC comparison on the GLODAP 1° grid: vertically interpolate each model’s 2002 DIC field and GLODAP `TCO2` to 0, 100, 500, 1000 and 2000 m, horizontally regrid the model depth slices to the GLODAP grid, and for each model-depth pair compute the global mean bias over the cells where that model and GLODAP are both defined at that depth.
   - Tools: `regrid_vertical`, `regrid_horizontal`, `compute_difference_metrics`, `compute_basic_stat`.
   - Use the exact requested target depths rather than nearest native levels [omip-native-grid-vs-regridded].

9. Produce Figure 1 from the shared-mask 0-50 m DIC bias maps (model minus GLODAP) on the GLODAP grid: one global panel per model, one shared diverging colour scale centred on zero, labelled in µmol kg-1.
   - Tools: no additional data tool; plotting from step-6 outputs.

10. Produce Figure 2 from the depthwise DIC bias results: one panel with bias on the horizontal axis, depth increasing downward on the vertical axis, one line per model, a vertical zero line, and markers at 0, 100, 500, 1000 and 2000 m.
   - Tools: no additional data tool; plotting from step-8 outputs.

11. Synthesize model performance across DIC, alkalinity, and pH using the requested metrics: identify which model best matches GLODAP overall using the 0-50 m mean-state differences plus DIC map bias/RMSE/correlation and the depthwise DIC bias structure, and assess whether DIC, alkalinity, and pH errors are consistent in sign within each model.
   - Tools: no additional data tool; comparative interpretation from prior results.

Other notes:
- [cf-cell-methods] Treat the CMIP6 monthly fields as monthly means and form the 2002 and 1971-2000 means from those time means before spatial analysis.
