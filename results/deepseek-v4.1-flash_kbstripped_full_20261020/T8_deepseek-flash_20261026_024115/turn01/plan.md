# Plan

1. Load CMIP6 historical monthly fields with `load_regular_data` from `cmip6_historical` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: `dissic`, `talk`, `ph`, `thetao`, `so` for calendar year 2002, and `dissic`, `thetao`, `so` for 1971-2000, each on the model's native grid.
   - Note [cmip6-carbonate-system-definitions]: treat `dissic` as total dissolved inorganic carbon, `talk` as total alkalinity (carbonate+borate+water+phosphate+silicate terms) and `ph` as in-situ pH on the TOTAL scale; all three models use OMIP best-practice constants, so they are directly comparable to GLODAP TCO2/TAlk/pH.
   - Note [cmip6-grid-labels]: GFDL-ESM4 is `gr` (already regridded to 1° with its own land–sea mask) and its `volcello` is monthly because 75 hybrid layers were remapped to 35 z levels; CESM2 (~100 km) and MPI (~250 km) differ in resolution, so features narrower than the coarse grid are not comparable in amplitude.

2. Load GLODAPv2.2016b with `load_climatology_data`: `TCO2`, `TAlk` and `pHtsinsitutp` (in-situ pH, total scale — the field matching the models' `ph`), plus its 1° grid and depth levels via `load_depth_levels` for GLODAP (and for CMIP6 for the vertical interpolation).

3. Convert the models' `dissic` and `talk` to micromoles per kilogram before any comparison: compute in-situ density from `thetao` and `so` with `compute_potential_density_anomaly(..., ref_pressure=None)` and divide, then report `attrs['density_assumption']`.
   - Note [cmip6-bgc-units-per-volume]: µmol kg-1 = mmol m-3 × 1000 / ρ using the in-situ density from `compute_potential_density_anomaly` with `ref_pressure=None`; do NOT use `transform_unit`'s constant default, and do not treat mmol m-3 as comparable to GLODAP's µmol kg-1.

4. Load each model's native-grid `areacello` with `load_areas` and `volcello` with `load_volumes_for_cmip6` (all three models, for both periods), and reduce the monthly 2002 fields to a calendar-year annual mean and the 1971-2000 fields to the long-term mean with `temporal_aggregate`.

5. Compute the 0–50 m means of DIC, total alkalinity and pH for each model and for GLODAP with `compute_basic_stat`: one mean per source over all its own 0–50 m cells, on each product's own grid; for the models use `op='volume_weighted_mean'` with `volcello` and NaN-aware weights. Report DIC/TAlk in µmol kg-1 to the nearest whole unit, pH to 3 decimals, and each model's difference from GLODAP.
   - Note [omip-weighted-means-and-integrals]: the "one mean over all its 0-50 m cells" is taken as a weighted mean — models weighted by `volcello` (NaN-aware weights); state that a plain cell mean would weight every 0-50 m cell equally instead.
   - Note [omip-native-grid-vs-regridded]: compute these means (and the period change in step 8) on each model's own grid with its own `areacello`/`volcello`; regrid to the GLODAP 1° grid only for maps and bias/RMSE/correlation, using `regrid_horizontal` with 'conservative', and state the method.

6. Build each model's 0–50 m mean DIC map (each cell's mean over its own 0–50 m levels) on its native grid and put it on the GLODAP 1° grid with `regrid_horizontal` ('conservative'); form the comparison set with `build_comparison_mask([...], ocean_only=False)` over model and GLODAP, apply `.where(mask)`, then get mean bias and RMSE with `compute_difference_metrics` and the spatial pattern correlation with `analyze_correlation(weighting='area')`; report bias/RMSE to the nearest whole µmol kg-1 and r to 3 decimals. Keep the gridded bias field of each model for Figure 1.
   - Note [common-comparison-mask]: regrid first, then build one mask over model and GLODAP with `build_comparison_mask(..., ocean_only=False)` (the reference is not WOA23, whose coastline would drop the hardest coastal cells) and apply `.where(mask)` before bias, RMSE and correlation.

7. For the vertical structure of the DIC bias: interpolate each model's annual-mean 2002 DIC and the GLODAP TCO2 onto 0, 100, 500, 1000 and 2000 m with `regrid_vertical` (not nearest level), regrid horizontally to the GLODAP 1° grid with `regrid_horizontal` ('conservative'), mask per depth with `build_comparison_mask(...ocean_only=False)` so that only cells where both the model and GLODAP are defined at that depth are used ([common-comparison-mask]), and compute the area-weighted global mean bias (model minus GLODAP) at each depth with `compute_difference_metrics`; report to the nearest whole µmol kg-1, together with the number of valid cells at each depth for each model. Keep the five biases per model for Figure 2.

8. For each model compute the 2002 minus 1971-2000 change of the global 0–50 m mean DIC on the model's own grid with `compute_basic_stat(op='volume_weighted_mean')` and `volcello` ([omip-native-grid-vs-regridded]), and divide it by that model's 0–50 m DIC difference from GLODAP obtained in step 5, reporting the ratio to 2 decimals.
   - Note [anomaly-baselines-and-climate-indices]: state the two periods/baselines used (2002 annual mean vs the model's own 1971-2000 long-term mean; GLODAP is the 1972-2013 mapped climatology); the period change uses the model's own fields so its mean-state bias cancels — do not read it as a model error.

9. Rank the three models from the bias/RMSE/correlation and depth-profile metrics and check whether the signs of their DIC, alkalinity and pH differences from GLODAP agree, stating which model reproduces the observed carbon system best.

10. Figures: Figure 1 — one row of three global maps of the 0–50 m DIC bias (model minus GLODAP) on the masked GLODAP 1° comparison cells, one shared diverging colour scale centred on zero, labelled µmol kg-1, one panel per model. Figure 2 — a single panel with depth (m) on the vertical axis increasing downward, global mean DIC bias (µmol kg-1) on the horizontal axis, one line per model, a vertical zero reference line and markers at 0, 100, 500, 1000 and 2000 m.

Other notes:
- [cmip6-experiments-historical-ssp245] `historical` is a coupled run with its own internal variability: a model year need not line up with the observed year, and a single member cannot separate forced change from internal variability — attach this caveat to the 2002-minus-1971-2000 DIC change.
- [omip-weighted-means-and-integrals] The question's wording "one mean over all its 0-50 m cells" could be read as an unweighted cell mean; the reported means are volume-weighted, so state this difference explicitly.
- [common-comparison-mask] Report the cell counts behind each depth in the bias profile, since coverage shrinks with depth.
- The GLODAP differences (DIC, TAlk, pH and the depth profile) carry the model's mean-state bias; the period change does not.
