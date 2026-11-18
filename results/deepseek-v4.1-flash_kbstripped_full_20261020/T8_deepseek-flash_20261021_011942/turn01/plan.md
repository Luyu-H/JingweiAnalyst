# Plan

1. Load the fields with `load_regular_data` (dataset `cmip6_historical`), for each of CESM2, GFDL-ESM4 and MPI-ESM1-2-LR: monthly `dissic`, `talk`, `ph`, `thetao`, `so` for calendar year 2002, and monthly `dissic`, `thetao`, `so` for 1971–2000; also load each model's fixed `areacello` and `volcello` (via `load_areas` and `load_volumes_for_cmip6`, noting the `attrs['area_source']`/`attrs['volume_source']`). Take the 2002 annual mean of each monthly field so the year is one annual-mean state.
   - Note [cmip6-carbonate-system-definitions]: compare model `ph` (total scale, in situ) only with GLODAP `pHtsinsitutp`, never `pHts25p0`; if the pH result leans on CESM2's carbonate constants, flag in the answer that the CESM2 paper does not document them (unverified assumption).

2. Load the reference with `load_climatology_data` (dataset `glodap`): `TCO2`, `TAlk` and `pHtsinsitutp` (33 standard levels, no time dimension) as the observational climatology for all comparisons.

3. Put the model carbon fields on a per-kilogram basis. Compute in-situ density for every model cell with `compute_potential_density_anomaly` using `thetao` + salinity and `ref_pressure=None`; use `transform_unit` only for the mol m-3 → mmol m-3 factor, then divide by that density in code to reach µmol kg-1.
   - Note [cmip6-so-practical-salinity]: `so` holds practical-salinity numbers (order 35) despite the units string "0.001" — never rescale it; convert S_P to S_A with the location-dependent gsw function before it enters the density calculation.
   - Note [cmip6-bgc-units-per-volume]: `dissic`/`talk` are mol m-3, so µmol kg-1 = mmol m-3 × 1000 / ρ with ρ from `compute_potential_density_anomaly` (`ref_pressure=None`); do not use `transform_unit`'s constant default (mean surface density), which errs by about a per cent in the deep ocean, and record the density used in `attrs['density_assumption']`.

4. Build the 0–50 m weights: for each model cell take the cell thickness as `volcello`/`areacello` (per the `cell_thickness_like` recipe) and restrict the vertical weighting to the 0–50 m levels each column actually has.
   - Note [cmip6-static-grid-fields]: weight with each model's own `areacello`/`volcello`; the 0–50 m thickness is `volcello`/`areacello` (`cell_thickness_like`), since `thkcello`/`masscello`/`rhozero` are absent and `lev_bnds` differences fail at partial bottom cells; report `attrs['area_source']` and `attrs['volume_source']`.

5. Compute, per source (three models and GLODAP), the single global 0–50 m mean of DIC, total alkalinity and pH over all of that source's own 0–50 m cells, and each model's difference from GLODAP; report DIC and alkalinity to the nearest µmol kg-1, pH to 3 decimals.
   - Note [omip-weighted-means-and-integrals]: use `compute_basic_stat` with `op='volume_weighted_mean'` for the means that span depth (all 0–50 m means) and `op='area_weighted_mean'` for the horizontal ones, with weights restricted to the same valid (NaN-free) cells; never take an unweighted `mean()` over (y, x) or over latitude.
   - Note [omip-native-grid-vs-regridded]: compute each source's global means, and the 2002-minus-1971–2000 change of step 8, on its own native grid with its own weights; regrid (`regrid_horizontal`, 'conservative') only for the map comparison, and state which numbers come from which grid.

6. Build the 0–50 m DIC map for each source: per column, the mean of DIC over that column's own 0–50 m levels (volume-weighted), then `regrid_horizontal` with 'conservative' onto the GLODAP 1° grid. Score the maps with `build_comparison_mask` (`ocean_only=False`, the reference is not WOA23), applied by `.where(mask)`, then `compute_difference_metrics` for mean bias and RMSE (µm kg-1, nearest whole unit) and `analyze_correlation` with `weighting='area'` for the spatial correlation (3 decimals).
   - Note [common-comparison-mask]: put all sources on the GLODAP 1° grid and take the DIC bias/RMSE/correlation on one mask from `build_comparison_mask` with `ocean_only=False`, applied by `.where(mask)`; report the cell counts and each model's `mean_data`/`mean_reference` beside the bias, and state that the per-source 0–50 m global means of step 5 each cover that source's own cells, so their model-minus-GLODAP difference is not the shared-mask bias.

7. Vertical-profile biases: interpolate every product onto 0, 100, 500, 1000 and 2000 m with `regrid_vertical` (never select the nearest level), regrid horizontally onto the GLODAP 1° grid as in step 6, and at each depth average model-minus-GLODAP DIC over every cell where both are defined at that depth, using `compute_difference_metrics` (area-weighted), reporting to the nearest µmol kg-1.

8. Period change: for each model, 2002 annual-mean 0–50 m global DIC minus the 1971–2000 mean of the same quantity, computed on that model's native grid with its own weights (`compute_basic_stat`, `op='volume_weighted_mean'`), reported to the nearest µmol kg-1; divide this change by the model's step-5 0–50 m DIC difference from GLODAP (2 decimals) and state what the ratio implies about reading the bias as a model error.

9. Synthesise: rank the models on DIC, alkalinity and pH difference from GLODAP and on the map bias/RMSE/correlation to say which reproduces the observed carbon system best, and check whether the three models' DIC, alkalinity and pH errors agree in sign.

10. Figures (own code, no tool covers plotting): Figure 1 — the three 0–50 m DIC-bias maps on the step-6 shared mask, one row, one diverging colour scale centred on zero, labelled µmol kg-1; Figure 2 — single panel, vertical DIC-bias profile from step 7, depth axis increasing downward, µmol kg-1 on the horizontal axis, one line per model, vertical zero line, the five sampled depths marked.

Other notes:
- All intermediate values are carried unrounded; rounding happens only where the question requests a precision (whole µmol kg-1, 3 decimals, 2 decimals).
- GLODAPv2 has no time dimension, so its single TCO2/TAlk/pH climatology is the reference state against which each model's 2002 annual mean is scored.
- The 0–50 m global means (step 5) and the map/depth metrics (steps 6–7) rest on different cell sets — native valid cells versus the shared GLODAP-grid mask — and the answer must keep them distinct.
- Report the per-model density assumption (`attrs['density_assumption']`) and grid-weight sources (`attrs['area_source']`, `attrs['volume_source']`) alongside the numbers.
