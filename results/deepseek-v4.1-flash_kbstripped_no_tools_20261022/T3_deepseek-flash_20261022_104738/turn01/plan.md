# Plan

1. Read the CMIP6 historical `tos` monthly fields for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR for 1971-01 to 2000-12, together with each model's own `areacello` fixed field, directly from the raw netCDF files with xarray (open the model/realization files, concatenate the time axis). Apply CF masking so the missing-value flag 1e20 (land and below-sea-floor cells) becomes NaN before any arithmetic.
   - Note [missing-values-and-land-cells]: mask every `tos`/`areacello` array with the CF `_FillValue`/`missing_value` (1e20) immediately after decode, before any scale/offset reuse, so land cells never enter sums or means.

2. Average each model's masked `tos` over all 360 months (1971-01 to 2000-12) to give its 1971-2000 SST climatology on its own native grid. Compute each model's global mean on that native grid as the area-weighted sum Σ tos·areacello / Σ areacello, with `areacello` masked to exactly the cells where `tos` is valid. Report these three native-grid means in °C to 2 d.p.
   - Note [omip-weighted-means-and-integrals]: every mean is area-weighted (Σ tos·areacello / Σ areacello) with the weights restricted to the valid cells; the difference statistics of step 5 use area weights as well.
   - Note [omip-native-grid-vs-regridded]: the native-grid global means above belong on each model's own grid and own `areacello`; regridding to the WOA 1° grid happens only for the maps, RMSE and pattern correlation, so native-grid and shared-mask means are not expected to agree.
   - Note [cmip6-grid-labels]: GFDL-ESM4 `tos`/`areacello` are already on a remapped 180×360 grid, so its "native" mean is a mean on the 1° archive grid with a regridded land-sea mask; CESM2 (~100 km) and MPI-ESM1-2-LR (~250 km) differ by ~2.5× in resolution, so features narrower than the coarse grid are not comparable in amplitude.

3. Read the WOA23 `temperature` field for the 1971-2000 climate normal, annual time slice, at 0 m depth, on the 1° grid, plus the WOA23 `landsea_01` field on the same grid. Define ocean cells as `landsea_01 > 1`. Compute the global mean of the WOA23 0 m annual field over the whole 1° grid as a latitude-area-weighted mean (weights ∝ cos(latitude)) over the valid temperature cells, and report it in °C to 2 d.p.
   - Note [cmip6-tos-bulk-top-cell]: `tos` is the top-cell bulk temperature including under-ice values near −1.8 °C rather than missing; keep this in mind when comparing against WOA23, whose sea-ice-covered cells are ice-free-water values.

4. Build the single shared comparison mask on the WOA23 1° grid across all three models and WOA23: require `landsea_01 > 1` (ocean, `ocean_only=True`) together with a valid `tos` value for CESM2, GFDL-ESM4, MPI-ESM1-2-LR and a valid WOA23 temperature cell, and apply it with `.where(mask)` to every field before computing any bias, RMSE or correlation.
   - Note [common-comparison-mask]: score all three models and WOA23 on this one shared mask (`build_comparison_mask` with `ocean_only=True`), never on each model's own pairwise cells; report the mask cell count and the shared-sample means together with the bias.

5. Regrid each model's 1971-2000 `tos` climatology from its native grid to the WOA23 1° grid with a conservative (coarsening) regridding method, state the method, then apply the step-4 mask. On those shared cells compute, for each model: the shared-sample mean SST of the model and of WOA23, the mean bias delta_SST = model − WOA23, the area-weighted RMSE against WOA23, and the area-weighted spatial Pearson correlation with WOA23. Report means and biases in °C to 2 d.p., RMSE and correlation to 3 d.p.

6. Rank the three models from smallest to largest absolute mean bias |delta_SST| using the unrounded biases from step 5, and report the ranking.

7. Figure 1: 2×2 panel of global maps on the step-4 comparison cells — the CESM2, GFDL-ESM4 and MPI-ESM1-2-LR bias fields (model minus WOA23) on one shared diverging colour scale centred on zero, and the WOA23 0 m reference SST on a sequential colour scale; label every panel and give each colour bar a unit label in °C.

8. Figure 2: bar chart of the four global mean SST values from steps 2 and 3 (CESM2, GFDL-ESM4, MPI-ESM1-2-LR on their own native ocean grids with their own `areacello`, and WOA23 on the whole 1° grid), each bar annotated with its value and the y axis labelled in °C. Use the unrounded values, rounding only the printed annotations.

Other notes:
- [cmip6-experiments-historical-ssp245] Coupled historical runs carry their own internal variability, so no single model year is expected to match an observed year; the comparison is made on the 1971-2000 climatological statistics only.
- The WOA23 global mean is quoted over the whole 1° grid as the question prescribes, which differs from the shared-mask WOA23 mean of item 3/step 5; both are reported and are different quantities.
- [cmip6-grid-labels] Bias-pattern amplitude differences between CESM2/GFDL-ESM4 and MPI-ESM1-2-LR partly reflect their ~2.5× resolution difference and must be stated when describing the spatial structure.
