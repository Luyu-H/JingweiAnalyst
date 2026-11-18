# Plan

1. Load WOA23 climatology on the 1971-2000 climate normal: oxygen, percent_oxygen_saturation and apparent_oxygen_utilization with load_climatology_data; load landsea_01 with load_landsea_mask and basinmask_01 with load_basin_mask; load WOA23 depth levels/layer boundaries with load_depth_levels and cell areas with load_areas.

2. Load CMIP6 historical o2, thetao and so for GFDL-ESM4 and MPI-ESM1-2-LR, monthly 1971-01 to 2000-12, with load_regular_data; load each model's volcello with load_volumes_for_cmip6 and areacello/deptho with load_regular_data. Do not load CESM2 historical o2 because it is missing.

3. Load CMIP6 ssp245 o2, thetao and so for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly 2071-01 to 2100-12, with load_regular_data; load each model's volcello with load_volumes_for_cmip6 and areacello/deptho with load_regular_data.
- Note [historical-ssp245-member-continuity]: CESM2 historical r1i1p1f1 and ssp245 r10i1p1f1 are different realisations with no historical baseline, so show CESM2 ssp245 state only and no CESM2 change, and label its Figure 2 panel accordingly; GFDL-ESM4 and MPI-ESM1-2-LR are r1→r1, so their historical and ssp245 differences and percentage changes are valid.

4. For each CMIP6 model/experiment, compute in-situ density from thetao and so with compute_potential_density_anomaly using ref_pressure=None; convert o2 from mol m-3 to mmol m-3 with transform_unit, divide by the density field with compute_ratio, and convert to µmol kg-1 with transform_unit to obtain oxygen in µmol kg-1 for threshold comparisons.
- Note [cmip6-bgc-units-per-volume]: convert CMIP6 o2 (mol m-3) to µmol kg-1 as mmol m-3 × 1000/ρ using in-situ density from compute_potential_density_anomaly with ref_pressure=None, not transform_unit's constant default; take the dissolved oxygen inventory on the volume basis Σ o2·volcello (moles, no density), not from the converted per-mass field.

5. Average monthly converted o2, thetao and so over each period with temporal_aggregate (month-length weighted): historical 1971-2000 for GFDL-ESM4 and MPI-ESM1-2-LR, and ssp245 2071-2100 for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR; also average volcello over the same period. Separately, compute monthly dissolved oxygen content from raw o2 and volcello and temporal_aggregate it to the period mean for the inventory.
- Note [monthly-means-are-time-averages]: average monthly fields over each period with temporal_aggregate, which weights by month length; form threshold/OMZ diagnostics from the period-mean fields (average-first) to match WOA23, and state this rather than reporting monthly threshold frequencies.

6. Using the period-mean converted oxygen fields, call identify_oxygen_minimum_zone at thresholds 20 and 60 µmol kg-1 to obtain per-column OMZ thickness and total vertical thickness for each source/experiment. For CMIP6 supply cell thickness from volcello/areacello; for WOA23 use load_depth_levels layer boundaries. Apply WOA23 landsea_01 to WOA23 fields.
- Note [cmip6-static-grid-fields]: no thkcello/masscello locally; obtain column thickness from cell_thickness_like (= volcello/areacello), not lev_bnds differences, which are wrong at partial bottom cells; loaders handle GFDL-ESM4 time-varying volcello.

7. Use compute_exceedance_fraction by_volume at 20 and 60 µmol kg-1 on each period-mean oxygen field to get the total volume below each strict threshold and the fraction, using volcello for CMIP6 and WOA23 layer geometry from load_areas/load_depth_levels for WOA23; convert m3 to 1e15 m3 and report 3 significant figures for WOA23, GFDL-ESM4 and MPI-ESM1-2-LR historical, and CESM2, GFDL-ESM4 and MPI-ESM1-2-LR ssp245.
- Note [missing-values-and-land-cells]: 1e20 decodes to NaN via xarray; build each field's valid mask from that field and apply the same mask to its weights; use compute_exceedance_fraction rather than hand threshold comparisons so land and missing cells are excluded from numerator and denominator.
- Note [omip-weighted-means-and-integrals]: OMZ volume is a volume-weighted sum of volcello over sub-threshold cells; compute_exceedance_fraction omits missing cells and returns the m3 behind the percentage; mask weights to the same valid cells and state the Boussinesq assumption for models other than MPIOM and the basin/region mask used for each integral.

8. Call locate_extreme on each 60-threshold OMZ thickness field to find the maximum thickness and its latitude/longitude for the six combinations; report metres to 1 decimal place and coordinates to 2 decimal places.

9. Compute the global dissolved oxygen inventory for GFDL-ESM4 and MPI-ESM1-2-LR historical and CESM2, GFDL-ESM4 and MPI-ESM1-2-LR ssp245 from raw o2 (mol m-3) and volcello using convert_concentration_to_content, then temporal_aggregate to the period mean; sum over the ocean and report moles to 3 significant figures. Use the [cmip6-bgc-units-per-volume] volume-basis inventory Σ o2·volcello, with no density conversion.

10. Compute the percentage change in OMZ volume at 20 and 60 µmol kg-1 for GFDL-ESM4 and MPI-ESM1-2-LR as (ssp245 2071-2100 − historical 1971-2000)/historical ×100 using the unrounded step 7 volumes and compute_ratio; positive means an expanding OMZ. Do not form a CESM2 change ([historical-ssp245-member-continuity]).

11. Load WOA23 basinmask_01; select_basin for Pacific, Indian and Atlantic with include_marginal_seas=True and depth_dependent=False, excluding Southern and Arctic Oceans. For WOA23 and GFDL-ESM4/MPI-ESM1-2-LR historical, apply each selected basin to the 60-threshold oxygen field or its sub-threshold mask; where model and mask grids differ, align them with regrid_horizontal; compute basin OMZ volume with compute_exceedance_fraction by_volume and the [omip-weighted-means-and-integrals] volume weighting, using the model's volcello or WOA23 cell volumes. Report 1e15 m3 to 3 significant figures.
- Note [basin-and-region-definitions]: the CMOR Ofx basin field is absent; use the WOA23 basin mask with include_marginal_seas=True and depth_dependent=False, folding marginal seas into their parent oceans and treating Southern and Arctic as separate/excluded basins as requested; state the mask and these choices since they change basin means.

12. Using WOA23 oxygen, apparent_oxygen_utilization and percent_oxygen_saturation, restrict to cells with oxygen strictly below 60 µmol kg-1 and valid volume; compute volume-weighted means of AOU and percent saturation with compute_basic_stat; report each to 1 decimal place. Interpret high AOU/low saturation as respiration-dominated and low AOU/high saturation as ventilation-limited, based on the observed values.
- Note [cmip6-o2sat-and-derived-oxygen]: WOA23 AOU and percent saturation rest on a solubility/pressure convention that need not match the OMIP Garcia and Gordon recipe, so model-observation differences carry a definitional component; CESM2 has no historical o2/o2sat.

13. Compare historical GFDL-ESM4 and MPI-ESM1-2-LR OMZ volumes at 20 and 60 µmol kg-1 with WOA23 using compute_difference_metrics and relative differences via compute_ratio; identify which model is closer to the observed WOA23 volume. Compare the two models' step 10 percentage changes at both thresholds for sign and size, and state whether they agree.

14. Produce Figure 1 by mapping the 60-threshold OMZ thickness for GFDL-ESM4 historical, MPI-ESM1-2-LR historical and WOA23 in one row with one shared sequential colour scale in metres. Produce Figure 2 by mapping ssp245 minus historical 60-threshold thickness for GFDL-ESM4 and MPI-ESM1-2-LR using compute_difference_metrics, with one shared diverging scale centred on zero; place CESM2 ssp245 60-threshold thickness in the third panel on its own sequential scale, labelled as an ssp245 state with no historical baseline and no change. Use regrid_horizontal where a common display grid is needed.

Other notes:
- Keep all intermediate volumes, inventories, thicknesses and coordinates unrounded; round only when reporting the requested precision.
- The WOA23 comparison is a mean-state comparison between a climatology and CMIP6 period means, not a comparison of monthly threshold frequencies.
- Report the strict inequalities, the 1e15 m3 volume units, and the model/mask choices explicitly with the numbers.
