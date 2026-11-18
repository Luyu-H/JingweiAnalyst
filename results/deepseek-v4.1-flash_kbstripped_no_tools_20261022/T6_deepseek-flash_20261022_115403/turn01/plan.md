# Plan

1. Open the WOA23 `nitrate` field for period `all`, annual frequency, on its 1° grid, together with `landsea_01` and `basinmask_01`. Build the WOA ocean mask from `landsea_01 > 1` and keep nitrate fill values as NaN. This annual all-period field is the observational reference.

2. Open CMIP6 historical `no3` monthly for 1971-01 to 2000-12 from CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, each on its own native horizontal grid and native vertical levels. Also open each model’s `areacello`, and `thetao`/`so` where published for the density conversion required in Figure 2. Retain levels below 1000 m so the deepest allowed gradient can use the level just below 1000 m.
   - Note [missing-values-and-land-cells]: decode 1e20 to NaN, build each variable’s valid mask from that variable itself, mask weights to the same valid cells, and treat any near-1e20 output value as an unmasked fill value rather than a result.

3. Average each model’s monthly `no3` over 1971-2000 into one annual-mean climatological field using `temporal_aggregate` with day weights in that model’s calendar. Do not compute monthly nitraclines and then average them; average the nitrate field first and locate the gradient afterwards. WOA23 nitrate is already an annual climatology.
   - Note [monthly-means-are-time-averages]: state that the order used is time-average first, then locate the maximum gradient; an equal-weight average of monthly means is not a time average.

4. For each of the four nitrate fields, compute vertical gradients between adjacent native levels, considering only gradients whose shallower level is within the top 1000 m while retaining the level just below 1000 m for the deepest allowed gradient. Define nitracline depth as the depth of the maximum downward nitrate gradient, taken at the midpoint of the two levels forming that gradient. Keep each nitracline field on its source’s own horizontal and vertical grid; do not select a nearest common level.
   - Note [omip-native-grid-vs-regridded]: compute global and basin means on each source’s own grid with its own area weights; regrid only for the maps, pattern correlation and RMSE, using conservative horizontal regridding when going coarser.

5. Compute global and Pacific/Atlantic/Indian mean nitracline depths for each source on its native grid, using model `areacello` and WOA23 1° cell areas from the grid cell bounds. Obtain basin labels from WOA23 `basinmask_01` with `include_marginal_seas=True` and `depth_dependent=False`, carry those labels onto each model grid, and exclude the Southern and Arctic basins. Compute model-minus-WOA23 differences from the unrounded means.
   - Note [basin-and-region-definitions]: use the WOA23 basin labelling (not the OMIP 10-region mask), with marginal seas folded into their parent ocean and depth-independent surface labels; the Ofx `basin` field is unavailable, so basin labels must come from the WOA23 mask.
   - Note [omip-weighted-means-and-integrals]: every global and basin mean must be area-weighted, with weights restricted to the same valid cells; `compute_difference_metrics` is area-weighted by default and the correlation requires `weighting='area'`.

6. Regrid each model nitracline-depth field horizontally to the WOA23 1° grid with the conservative method, without vertical regridding. Build one shared valid-data mask from the four regridded nitracline-depth fields with `ocean_only=True` because WOA23 is the reference, apply it with `.where(mask)`, and report the mask cell counts.
   - Note [common-comparison-mask]: build the shared mask from the nitracline-depth fields themselves, not from a surface mask, and use the same masked cells for correlation and RMSE.
   - Note [cmip6-grid-labels]: GFDL-ESM4 is `gr` and already on a 1° archive grid; CESM2 (~100 km) and MPI-ESM1-2-LR (~250 km) have different resolutions, so narrow features are not comparably resolved.

7. On the shared mask and WOA 1° grid, compute each model’s spatial Pearson correlation with the WOA23 nitracline field and the RMSE in metres, using the area weighting specified in [omip-weighted-means-and-integrals].

8. Produce Figure 1: a 2×2 global map panel of nitracline depth for CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23 on the WOA 1° grid, using the shared mask and one common sequential colour scale chosen so that greater depth is visually darker/heavier (for example `YlGnBu`), with common vmin/vmax across all four panels.

9. Produce Figure 2: for each source, average nitrate over 2S–2N in the equatorial Pacific, subset 120E–80W and the top 1000 m, and plot longitude–depth panels in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with depth increasing downward and one shared colour scale in µmol kg⁻¹. Convert model `no3` from mol m⁻³ to µmol kg⁻¹ using in-situ density from the model `thetao`/`so` and `compute_potential_density_anomaly` with `ref_pressure=None`; WOA23 is already µmol kg⁻¹. Compute each panel’s nitracline depth per longitude from that same 2S–2N band-averaged profile and overlay it as a depth-versus-longitude line.
   - Note [cmip6-bgc-units-per-volume]: convert model nitrate with in-situ density for Figure 2 only, record `attrs['density_assumption']`, and state that nitracline depth itself is scale-invariant so the depth line is located on the native-unit profile without density conversion.

10. Report global mean nitracline depths and model-minus-WOA23 differences to the nearest metre; report Pacific, Atlantic and Indian means and differences to the nearest metre; report spatial correlation to 3 decimal places and RMSE to 1 decimal place. Identify the best model from the lowest area-weighted RMSE with the correlation as supporting evidence, and check whether each model reproduces the WOA23 ordering of nitracline depth across the three basins. Do not round intermediate values.

Other notes:
- WOA23 nitrate used here is the all-period annual composite rather than a 1971–2000 nutrient normal; the model fields are averaged over 1971–2000 as requested, so the comparison carries an unavoidable period mismatch that must be stated.
- All reported means, differences, correlations and RMSEs are computed from unrounded intermediate fields and rounded only at the stated reporting precision.
- “Best” model is judged by lowest area-weighted RMSE against WOA23, with spatial correlation as supporting evidence; basin-order agreement is assessed against the WOA23 Pacific/Atlantic/Indian nitracline-depth order.
