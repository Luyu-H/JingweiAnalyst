# Plan

1. Load the data.
   - Use `load_regular_data` for cmip6_historical `no3` and `po4`, for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. Period is 1971-01 to 2000-12, monthly, all depth levels from 200 to 1000 m.
   - Load each model's own `volcello` with `load_volumes_for_cmip6`.
   - Use `load_depth_levels` to get the levels and layer boundaries for the 200-1000 m layer.
   - Use `load_climatology_data` for woa23 `nitrate` and `phosphate` (all-period annual composite).
   - Use `load_basin_mask` for the WOA23 basin mask (basinmask_01), and `load_landsea_mask` if needed.
   - Note [cmip6-static-grid-fields]: take layer thickness and volume weights from each model's own volcello, not from lev_bnds differences. Record and report `attrs['volume_source']`.
   - Note [cmip6-grid-labels]: GFDL-ESM4 is `gr` (1° archive grid, 35 levels, regridded land-sea mask). CESM2 and MPI are `gn`. State this when comparing the sources.

2. Build the model climatologies.
   - Use `temporal_aggregate` for the day-weighted long-term mean over 1971-2000 of model no3 and po4. This gives a 1971-2000 climatology per model.
   - Restrict all four sources to the 200-1000 m layer on their original levels. Integrate or average on the original levels, with no interpolation and no nearest-level selection.
   - Note [omip-native-grid-vs-regridded]: compute the means, volume shares and basin means on each source's native grid with its own weights. Regrid only for the maps and the common-mask regression and scatter. State which grid each number comes from.

3. Compute N:P and the per-source statistics on native grids.
   - Use `compute_ratio` (NO3/PO4) cell by cell for each source. Native fields are fine because N:P is scale-invariant. The mean of this cell-by-cell N:P is the headline value.
   - Use `compute_basic_stat` with `op='volume_weighted_mean'` for mean NO3 and mean PO4 over 200-1000 m. Use each model's own volcello and, for WOA23, volumes from its depth levels and areas. Divide the means to get the ratio of means, and label it as the secondary value.
   - Use `compute_basic_stat` volume-weighted mean on the cell-by-cell N:P field for the headline mean.
   - Note [omip-weighted-means-and-integrals]: use NaN-aware weights restricted to valid cells. Label the two means clearly, with the cell-by-cell mean as the headline. Use `compute_exceedance_fraction` by volume, with missing cells out of the denominator. Report the volume in m3 next to each percentage. The regression weights each cell equally, so it differs from the volume-weighted statistics; state this.

4. Compute the volume shares.
   - Use `compute_exceedance_fraction` (by volume) for N:P strictly between 14 and 18, globally for each source.
   - Use the same tool for N:P below 12, globally and within the Pacific and Indian basins (via step 5 masking).

5. Compute the basin means.
   - Use `select_basin` with the WOA23 basin mask to pick the Pacific, Atlantic and Indian basins. Marginal seas are folded into their parent ocean. Southern and Arctic are excluded.
   - Compute the volume-weighted mean of the cell-by-cell N:P per basin for each source.
   - Note [basin-and-region-definitions]: use `include_marginal_seas=True` and `depth_dependent=False` (surface labels over the whole layer). State both settings.

6. Regrid and build the common mask.
   - Use `regrid_horizontal` to put all four sources on one common 1° grid (WOA23 grid). Use the 'conservative' method before any area-weighted step.
   - Use `build_comparison_mask` with the default `ocean_only=True`. Build it from the 200-1000 m layer fields of all four sources, and apply it with `.where(mask)` before any metric. Report the cell counts.
   - Note [common-comparison-mask]: use this shared sample for the regression and the scatter and map comparisons. State that global means over each source's own ocean differ from the shared-sample means.

7. Convert units and run the regression.
   - Note [cmip6-bgc-units-per-volume]: convert model no3 and po4 from mol m-3 to µmol kg-1 with mmol m-3 × 1000 / 1025 (ρ = 1025 kg m-3), and state that density. WOA23 is already in µmol kg-1. N:P, slope and R² are unaffected by the scaling; the intercept b needs the conversion.
   - Use `fit_multivariate_relation` for NO3 = a·PO4 + b on the shared-mask cells, each cell weighted equally. Do this for each of the four sources. Report a (2 dp), b in µmol kg-1 (2 dp) and R² (3 dp).

8. Compare against Redfield and observations.
   - Judge each source's N:P as at, above or below 16 using the headline mean.
   - Judge whether each model reproduces the WOA23 low-N:P (<12) volume and its Pacific and Indian distribution. Use `compute_difference_metrics` against WOA23 on the shared mask for supporting bias and RMSE of N:P.

9. Figures.
   - Figure 1: a 4-panel scatter of NO3 (y) against PO4 (x) in µmol kg-1 over 200-1000 m. Panels in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with identical axis ranges and ticks. Each panel has its own fitted line and the Redfield 16:1 line through the origin, in distinct styles, annotated with a, b and R².
   - Figure 2: a 2×2 geographic map of the 200-1000 m mean N:P minus 16 on the common lon-lat frame. Use the regridded fields (display only) and one shared diverging colour scale centred on zero.

Other notes:
- Headline N:P is the mean of the cell-by-cell ratio. The ratio of means is reported separately and labelled.
- Native-grid statistics (means, volume shares, basin means) use each source's own volume weights. The regression and scatter use equal-weight cells on the shared mask, so they are not directly comparable to the volume-weighted values.
- Model concentrations are converted at constant ρ = 1025 kg m-3 for the regression, scatter and axes only.
- GFDL-ESM4 is on the 1° `gr` archive grid, so its "native" numbers are on that grid.
- The models cover 1971-2000 monthly means, while WOA23 is an all-period composite, so the periods are not identical.
- [cmip6-bgc-units-per-volume] Model no3 and po4 are per volume while WOA23 is per mass; only the intercept depends on the conversion.
- [common-comparison-mask] Shared-sample means differ from each source's own-ocean means.
