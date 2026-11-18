# Plan

1. Load the WOA23 annual objectively analysed nitrate and phosphate climatology with `load_climatology_data` (`woa23`, variables `nitrate`, `phosphate`); load WOA23 depth levels and layer boundaries with `load_depth_levels`, WOA23 1-degree cell areas with `load_areas`, and the WOA23 `basinmask_01` with `load_basin_mask`. Restrict to the depth levels/layers within 200–1000 m. Derive WOA23 cell volumes from the 1-degree areas and the layer thicknesses from `load_depth_levels`, since no WOA23 volcello product is listed. Keep WOA23 nitrate and phosphate in micromoles per kilogram.

2. Load CMIP6 historical `no3` and `po4` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly 1971-01 to 2000-12, with `load_regular_data`; load each model's own `volcello` with `load_volumes_for_cmip6`; load each model's depth levels and layer boundaries with `load_depth_levels`, and restrict to the levels/layers within 200–1000 m.
- Note [cmip6-grid-labels]: GFDL-ESM4 is `gr` remapped to 1 degree, so its “native-grid” numbers are numbers on the archive grid; CESM2 is about 100 km and MPI about 250 km, so narrow-feature amplitudes are not directly comparable across models.

3. Average each model's monthly `no3` and `po4` over 1971-2000 with `temporal_aggregate` using day-weighted long-term means.
- Note [monthly-means-are-time-averages]: use `temporal_aggregate` with day weights, not `compute_basic_stat(op='mean', dim='time')`; state that N:P is formed from the 1971-2000 mean NO3 and PO4 fields (ratio of mean fields), so the headline “mean of the cell-by-cell N:P ratio” is the volume-weighted mean of that ratio field.

4. Convert model `no3` and `po4` from mol m-3 to micromoles per kilogram using the question's constant density of 1025 kg m-3: micromoles kg-1 = mmol m-3 times 1000 / 1025 = mol m-3 times 1e6 / 1025. WOA23 is already in micromoles per kilogram.
- Note [cmip6-bgc-units-per-volume]: convert before any comparison, quote `attrs['density_assumption']`, and state that the constant-density choice differs from in-situ density by about 1 % in deep and fresh water.

5. Compute the 3D N:P field for each source with `compute_ratio` (NO3 divided by PO4) on the depth-restricted, unit-converted fields [cmip6-bgc-units-per-volume].

6. Build the shared comparison sample with `build_comparison_mask` using WOA23 as the reference and the default `ocean_only=True`; apply `.where(mask)` to the WOA23 and model NO3, PO4 and N:P fields for the cell-by-cell scatter and regression sample; report the valid cell counts next to the metrics.
- Note [common-comparison-mask]: build one shared mask with `build_comparison_mask` and apply `.where(mask)` to all four sources before any cell-by-cell statistic or scatter/regression sample; with WOA23 as the reference the default `ocean_only=True` is right; report the cell counts next to the metrics.

7. Compute native-grid, volume-weighted layer statistics for each source over 200–1000 m using `compute_basic_stat(op='volume_weighted_mean')` with each source's own cell volumes (models: `load_volumes_for_cmip6`; WOA23: areas times layer thickness from step 1): (a) the headline mean N:P is the volume-weighted mean of the `compute_ratio` N:P field; (b) compute the volume-weighted mean NO3 and mean PO4 separately, then the ratio of means is mean NO3 divided by mean PO4. Label both ways and report to 2 decimal places.
- Note [omip-weighted-means-and-integrals]: layer means and exceedance percentages must be volume-weighted with NaN-aware weights on the same cells; use `compute_basic_stat(op='volume_weighted_mean')` and `compute_exceedance_fraction` (which returns the m3 behind the percentage, so report it); `compute_ratio` gives the cell-by-cell ratio (headline), while the ratio of means is the two volume-weighted means divided; the question's OLS with “each cell weighted equally” is a plain unweighted fit, unlike the area weighting `compute_difference_metrics` applies by default.
- Note [omip-native-grid-vs-regridded]: take the native-grid layer means, ratios and percentages on each source's own grid with its own cell volume; regrid only for the Figure 2 maps for display, and state the `regrid_horizontal` method used.

8. Fit the per-source OLS regression NO3 = a times PO4 + b with `fit_multivariate_relation`, using NO3 as the dependent variable and PO4 as the predictor, on the common-mask sample from step 6 [common-comparison-mask], with each cell weighted equally as a plain unweighted fit [omip-weighted-means-and-integrals]. Report slope a to 2 decimal places, intercept b in micromoles per kilogram to 2 decimal places, and R-squared to 3 decimal places.

9. Compute the volume percentage with N:P strictly between 14 and 18 for each source over 200–1000 m using `compute_exceedance_fraction` with volume weights on the native grid: fraction above 14 minus fraction at or above 18, or the tool's equivalent bounded call; report to 1 decimal place and retain the m3 behind the percentage [omip-weighted-means-and-integrals].

10. Compute basin means using the WOA23 `basinmask_01` with `load_basin_mask` and `select_basin` with `include_marginal_seas=True` to retain the Pacific, Atlantic and Indian Oceans, excluding the Southern and Arctic Oceans; for each source and each of these basins, compute the volume-weighted mean N:P over 200–1000 m with `compute_basic_stat(op='volume_weighted_mean')` and the source's own volumes, to 2 decimal places.
- Note [basin-and-region-definitions]: basin values use the WOA23 mask via `select_basin`/`load_basin_mask` with `include_marginal_seas=True`; the Ofx `basin` field is not archived; state which mask was used, whether the Southern Ocean is separate, and that `depth_dependent=False` surface labels are what a named ocean usually means.

11. Compute the volume percentage with N:P below 12 globally and within the Pacific and Indian Oceans for each source with `compute_exceedance_fraction` using volume weights on the native grid, applying the basin subset first for the basin values [basin-and-region-definitions]; report to 1 decimal place and retain the m3 behind each percentage [omip-weighted-means-and-integrals].

12. Figure 1: produce a 4-panel scatter of nitrate (y) against phosphate (x) in micromoles per kilogram over 200–1000 m from the common-mask sample, with panels ordered WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR; use identical axis ranges and ticks; show each panel's own fitted regression line from step 8 and the Redfield 16:1 line through the origin, made distinguishable by line style or colour; annotate a, b and R-squared.

13. Figure 2: compute the 200–1000 m vertical volume-weighted mean N:P for each source on its native grid (horizontal field), subtract 16, and, for display only, put the model fields on a common regular 1-degree longitude-latitude grid with `regrid_horizontal` (state the method, for example bilinear); WOA23 is already on the 1-degree grid. Plot the 2-by-2 maps on a common longitude-latitude frame with one shared diverging colour scale centred on zero [omip-native-grid-vs-regridded].

14. Interpret and report for each source whether the headline mean N:P is at, above or below Redfield 16; compare the global, Pacific and Indian volume percentages below 12 and the Pacific, Atlantic and Indian mean N:P values to judge whether each model reproduces the observed low-N:P volume and its basin distribution.
- Note [bgc-model-structure-comparison]: interpret the spread through model structure, with CESM2 having variable C:P, GFDL-ESM4 fixed N:P per plankton type, and MPI-ESM1-2-LR prescribed Redfield; “not detailed in paper” does not mean a feature is absent.

Other notes:
- WOA23 nitrate and phosphate are the single all-period annual objective analysis, not a 1971-2000 climate normal, according to the dataset description; the model fields are 1971-2000 day-weighted climatologies.
- Check per-model variable availability when loading CMIP6 `no3` and `po4`; if a model does not publish one of these variables, report it as missing rather than substituting another product.
- Keep intermediate values unrounded; round only where the question requests a reported precision.
- Report the common-mask cell counts alongside the scatter/regression results, and keep the native-grid volume-weighted layer and basin metrics clearly separate from the common-mask regression sample.
