# Plan

1. Load WOA23 annual climatology for `nitrate` and `phosphate` over the global ocean and all depth levels needed to cover 200–1000 m, and load the WOA23 basin mask `basinmask_01`; tools: `load_climatology_data`, `load_basin_mask`.
- Note [basin-and-region-definitions]: use the WOA23 basin mask with marginal seas included in parent oceans, and exclude Southern and Arctic Oceans from Pacific/Atlantic/Indian basin summaries.

2. Load CMIP6 historical monthly `no3` and `po4` for models `CESM2`, `GFDL-ESM4`, and `MPI-ESM1-2-LR` over 1971-01 to 2000-12 for the global ocean and all depth levels needed to cover 200–1000 m; load each model’s `volcello`; tools: `load_regular_data`, `load_volumes_for_cmip6`, `temporal_aggregate`.
- Note [cmip6-bgc-units-per-volume]: keep native per-volume units for volume-weighted 200–1000 m means and volume fractions, but convert model `no3` and `po4` to `µmol kg-1` for the requested regression/scatter using the stated constant density 1025 kg m-3.
- Note [cmip6-static-grid-fields]: use each source’s own `volcello` for all 200–1000 m volume weighting and exceedance fractions, relying on the published cell volumes rather than deriving thickness from level bounds.
- Note [cmip6-grid-labels]: treat GFDL-ESM4 as already on its archive grid label `gr`, while CESM2 and MPI-ESM1-2-LR are `gn` and only need horizontal regridding for Figure 2 display.

3. Form the 1971–2000 monthly climatological mean fields for CMIP6 `no3` and `po4`; tools: `temporal_aggregate`.

4. Restrict WOA23 and each model to the 200–1000 m layer on their own native vertical grids, then compute cellwise N:P ratios (`NO3/PO4`) on each native grid; tools: `compute_ratio`.
- Note [omip-native-grid-vs-regridded]: compute all reported statistics on each source’s own horizontal and vertical grid, and regrid only for Figure 2 display.
- Note [omip-weighted-means-and-integrals]: for each source, compute both the volume-weighted mean of the cellwise N:P ratio (headline mean) and the ratio of the separately volume-weighted mean nitrate and mean phosphate, using valid-cell `volcello` weights over 200–1000 m.

5. Compute, for each of the four sources globally over 200–1000 m: the volume-weighted mean nitrate, volume-weighted mean phosphate, the headline mean N:P, and the ratio of mean nitrate to mean phosphate; tools: `compute_basic_stat`.

6. Convert CMIP6 climatological `no3` and `po4` from `mol m-3` to `µmol kg-1` with the constant-density assumption, and place WOA23 nitrate and phosphate in matching `µmol kg-1` units if needed; then fit the equal-cell-weight ordinary least-squares regression `NO3 = a * PO4 + b` over all valid 200–1000 m cells for each source, and report slope, intercept and R-squared; tools: `transform_unit`, `fit_multivariate_relation`.

7. Compute the percentage of 200–1000 m water volume with N:P strictly between 14 and 18 for each source, globally; tools: `compute_exceedance_fraction`.

8. Apply the WOA23 basin mask to each source and compute the headline mean N:P separately for the Pacific, Atlantic and Indian Oceans over 200–1000 m; tools: `select_basin`, `compute_basic_stat`.

9. Compute the percentage of 200–1000 m water volume with N:P below 12 for each source: globally, within the Pacific, and within the Indian Ocean; tools: `compute_exceedance_fraction`, `select_basin`.

10. Prepare Figure 1 from the 200–1000 m valid-cell nitrate and phosphate values in `µmol kg-1`: a four-panel scatter in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with common axis limits/ticks, each panel showing its own fitted regression line plus the Redfield 16:1 line through the origin and annotation of `a`, `b`, and `R-squared`; use the regression outputs from step 6.

11. Compute each source’s 200–1000 m mean of the cellwise N:P ratio, subtract 16, and horizontally regrid only these map fields to a common regular longitude-latitude grid for display; then prepare Figure 2 as a 2×2 set of true geographic global maps on one shared diverging colour scale centred on zero; tools: `compute_basic_stat`, `regrid_horizontal`.

12. Synthesize the numerical and visual results to state, for WOA23 and each model, whether the intermediate-water stoichiometry lies at, above, or below Redfield, and whether each model reproduces the observed low-N:P water volume and its basin distribution.

Other notes:
- [bgc-model-structure-comparison] In interpretation, note that model spread in N:P can reflect both circulation and biogeochemical model structure, so departures from WOA23 should not be attributed to transport alone.
