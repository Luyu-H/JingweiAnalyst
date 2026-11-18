# Plan

1. **Load the model fields.** `load_regular_data` for `cmip6_historical`, models CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, variables `no3`, `po4` and `volcello`, monthly, 1971-01 to 2000-12. Use `load_volumes_for_cmip6` for each model's `volcello` (and `areacello`) on its archive grid, and `load_depth_levels` for each model's `lev` levels and layer boundaries.

2. **Load the observations.** `load_climatology_data` for `woa23`, variables `nitrate` and `phosphate`, annual objectively analysed climatology on the 1° grid; `load_depth_levels` for its 33 standard levels and layer boundaries; `load_areas` for the WOA23 grid-cell areas (volume weights for WOA23 = area × layer thickness from `load_depth_levels`); `load_basin_mask` for `basinmask_01`.
   - Note [cmip6-bgc-units-per-volume]: convert CMIP6 `no3`/`po4` from mol m-3 to µmol kg-1 as mol m-3 × 1e6/1025 with `transform_unit`, recording `attrs['density_assumption'] = 1025 kg/m3`; WOA23 nitrate and phosphate are already µmol kg-1 and need no conversion.

3. **Form the 1971–2000 annual means of the models.** `temporal_aggregate` with `freq='mean'` on each model's `no3` and `po4` (and on GFDL-ESM4's monthly `volcello`).
   - Note [monthly-means-are-time-averages]: use `temporal_aggregate(freq='mean')`'s day weights rather than `mean(dim='time')`, and average `no3` and `po4` first, forming N:P afterwards, so the diagnostic is of the mean state and like-for-like with the WOA23 climatology; state that ratio-first would instead give the mean monthly N:P.

4. **Restrict to the 200–1000 m layer.** Keep only levels whose layer lies in 200–1000 m for each of the four sources (all levels fully or partly within the layer; the interface levels 200 m and 1000 m are excluded as levels, the layer between them is kept). Carry each cell's volume from `volcello` (models) or area × layer thickness (WOA23).
   - Note [cmip6-static-grid-fields]: use the published `volcello`/`areacello`, with layer thickness `volcello`/`areacello` rather than `lev_bnds` differences because partial bottom cells matter; GFDL-ESM4 `volcello` is monthly and is annual-averaged in step 3; report `attrs['volume_source']`/`attrs['area_source']` wherever they feed an integral.

5. **Compute the cell-by-cell N:P ratio for all four sources.** `compute_ratio` of nitrate over phosphate at every cell of the 200–1000 m layer, giving one ratio field per source.

6. **Headline stoichiometry, two ways, per source.** With `compute_basic_stat` over the 200–1000 m layer: (i) the volume-weighted mean of the cell-by-cell N:P ratio (headline, 2 dp) and (ii) the ratio of the volume-weighted mean nitrate to the volume-weighted mean phosphate (2 dp).
   - Note [omip-weighted-means-and-integrals]: use Σ field·`volcello` / Σ `volcello` over valid ocean cells with NaN-aware weights (plain `mean()` is wrong); the headline is the volume-weighted mean of the cell-by-cell ratio, and the ratio-of-means is the ratio of the volume-weighted mean `no3` to the mean `po4` — label which is which; the regression in step 7 is equal-cell weighted; volume shares must be shares of volume with missing cells out of the denominator.

7. **Ordinary least squares NO3 = a·PO4 + b, one fit per source.** `fit_multivariate_relation` predicting `no3` from `po4` over all valid 200–1000 m cells in µmol kg-1, each cell weighted equally, reporting slope a (2 dp), intercept b (2 dp) and R² (3 dp), plus the number of cells used.

8. **Volume shares of N:P.** `compute_exceedance_fraction` by volume on each source's N:P field: the percentage of water volume with N:P strictly between 14 and 18 (1 dp), recorded together with the m3 behind the percentage and the valid-volume denominator.

9. **Basin breakdown of the headline mean N:P.** `select_basin` with the WOA23 `basinmask_01` on each source's N:P field, once per basin, then `compute_basic_stat` volume-weighted mean of the cell-by-cell ratio for Pacific, Atlantic and Indian Oceans (2 dp); repeat `compute_exceedance_fraction` by volume for N:P below 12 globally, in the Pacific and in the Indian Ocean (1 dp, with m3).
   - Note [basin-and-region-definitions]: use `load_basin_mask`/`select_basin` rather than any CMOR Ofx `basin` field, with `include_marginal_seas=True` folding marginal seas into their parent ocean and `depth_dependent=False` (surface labels applied through the water column); treat Southern and Arctic as basins of their own and exclude them from the Pacific/Atlantic/Indian breakdown (they remain in the global numbers).

10. **Quantify the model–observation agreement.** Put each source's 200–1000 m N:P column field (volume-weighted over levels) on the WOA23 1° grid with `regrid_horizontal`, then build one shared sample and score each model against WOA23 with `compute_difference_metrics` (bias, RMSE, largest departure and its location) in N:P units.
    - Note [common-comparison-mask]: build one mask with `build_comparison_mask` (`ocean_only=True`, WOA23 as reference) for the 200–1000 m layer rather than the surface, apply `.where(mask)` to all four fields before any model–observation metric, report the mask cell count, and take the compared means from `compute_difference_metrics` on that shared sample.

11. **Figures.** Figure 1: 4-panel scatter of nitrate (y) vs phosphate (x) in µmol kg-1 over 200–1000 m for WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, identical axis ranges/ticks, each with its own fitted line from step 7, the Redfield 16:1 line through the origin, and a/b/R² annotated. Figure 2: 2×2 maps of the 200–1000 m mean N:P minus 16 for the four sources on a common longitude–latitude frame with one shared diverging scale centred on zero.
    - Note [cmip6-grid-labels]: GFDL-ESM4 is `gr` (regridded 1°) while CESM2 and MPI-ESM1-2-LR are `gn`, with ~2.5× different resolution, so narrow features are not comparable; the regular longitude–latitude fields used in Figure 2 are for display only — all means, masks and volume shares in steps 6–10 are computed on each model's archive/native grid with its own `areacello`/`volcello`.

Other notes:
- The question prescribes a constant ρ = 1025 kg/m3 for the regression concentrations; state that in-situ density below the surface would change the µmol kg-1 concentrations (and hence a and b), while the ratio-based quantities are unaffected because the same ρ cancels.
- State that the headline mean N:P is the mean of cell-by-cell ratios and that ratio-of-means is reported separately, since the two differ because N:P is nonlinear.
- State whether each source's mean N:P sits at, above or below Redfield 16:1, and whether each model reproduces the observed low-N:P volume fraction (<12) and its Pacific/Indian distribution, including the comparison-mask cell counts behind those statements.
- State that the Pacific/Atlantic/Indian means exclude the Southern and Arctic basins, which are retained in the global numbers.
