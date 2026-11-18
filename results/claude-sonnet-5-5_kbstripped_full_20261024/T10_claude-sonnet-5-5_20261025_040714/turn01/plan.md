# Plan

1. Load the data with load_regular_data (CMIP6 time series) and load_climatology_data (WOA23).
   - cmip6_historical: GFDL-ESM4 and MPI-ESM1-2-LR, `o2`, `thetao`, `so`, monthly 1971-01 to 2000-12. Also load each model's volcello, areacello and deptho.
   - cmip6_ssp245: CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, `o2`, `thetao`, `so`, monthly 2071-01 to 2100-12, with the same static fields.
   - woa23: oxygen, percent_oxygen_saturation and apparent_oxygen_utilization on the 1971-2000 normal (oxygen and the derived fields are the single all-period composite). Also load the landsea_01 and basinmask_01 masks with load_landsea_mask and load_basin_mask.
   - Note [historical-ssp245-member-continuity]: CESM2 ssp245 (r10) has no matching historical parent, so no CESM2 change is formed. Use the 2071-2100 minus 1971-2000 window differences for GFDL and MPI only, and do not fit a trend across the 2014/2015 join.

2. Take the 30-year long-term mean of each model field with temporal_aggregate (day-weighted long-term mean) over the stated windows. For WOA23, use the climatology as given.

3. Convert model `o2` to µmol kg-1.
   - Use load_areas and load_volumes_for_cmip6 for cell areas and volumes.
   - Compute in-situ density with compute_potential_density_anomaly from `thetao` and `so`.
   - Note [cmip6-bgc-units-per-volume]: ref_pressure=None, µmol kg-1 = mmol m-3 × 1000 / ρ. Do not use the constant default of transform_unit. State `attrs['density_assumption']`.
   - WOA23 is already in µmol kg-1.

4. Identify the OMZ at the 20 and 60 µmol kg-1 thresholds (strictly below) with identify_oxygen_minimum_zone, on each source's native grid and levels.
   - Note [cmip6-o2sat-and-derived-oxygen]: apply the thresholds only after the unit conversion. Use WOA23's own AOU and percent saturation for the OMZ mean, and do not mix them with model `o2`.
   - Note [omip-native-grid-vs-regridded]: compute volumes, inventories, basin totals and OMZ thickness on native grids and levels. Do not interpolate to the 33 WOA levels first.
   - Note [cmip6-static-grid-fields]: layer thickness is volcello / areacello (`cell_thickness_like`), not the lev_bnds difference. Report `attrs['area_source']` and `attrs['volume_source']`. Do not rescale models to a common ocean area.
   - Note [missing-values-and-land-cells]: build the valid mask from the `o2` field itself, use NaN-aware reductions, and treat values near 1e20 as fill.

5. Compute OMZ volume (1e15 m3) at both thresholds for all six source-experiment combinations.
   - Use compute_exceedance_fraction (strictly below, by volume) or compute_basic_stat with a volume-weighted sum on the thresholded cells.
   - Report the volume, the maximum OMZ thickness at 60, and its location. Locate the maximum thickness at 60 (summed thickness per column) with locate_extreme, giving latitude and longitude to 2 decimals.
   - Percentage change in OMZ volume, historical to ssp245, for GFDL-ESM4 and MPI-ESM1-2-LR at both thresholds: (ssp245 − hist) / hist × 100, positive for an expanding OMZ. Use compute_difference_metrics or direct arithmetic on the totals.

6. Compute the global dissolved oxygen inventory in moles with convert_concentration_to_content, using each model's own volcello and the µmol kg-1 → mol m-3 consistent conversion via density. This covers GFDL and MPI historical and all three models for ssp245. Sum with compute_basic_stat.

7. Compute basin OMZ volume at 60 for WOA23 and the two historical models.
   - Note [basin-and-region-definitions]: use select_basin on basinmask_01 with include_marginal_seas=True, depth_dependent=False, and keep the Pacific, Indian and Atlantic only. Exclude the Southern and Arctic Oceans. Report which option was used.
   - Sum the volumes on each source's own grid, per the volume note.

8. Compute the WOA23 OMZ mean AOU and percent saturation.
   - Mask the WOA23 cells where oxygen < 60.
   - Take the volume-weighted mean of apparent_oxygen_utilization and percent_oxygen_saturation with compute_basic_stat.
   - Interpret the result: high AOU and low saturation indicate respiration, while low AOU and high saturation indicate ventilation. Mention the solubility convention caveat.

9. Make the maps (Figures 1 and 2).
   - Note [common-comparison-mask]: regrid_horizontal (conservative) to the 1° grid. Build one shared mask with build_comparison_mask (ocean_only=True) and apply it with `.where(mask)`. Report the cell counts.
   - Figure 1: three maps in one row, historical OMZ thickness at 60 for GFDL-ESM4, MPI-ESM1-2-LR and WOA23, on one shared sequential colour scale in metres.
   - Figure 2: ssp245 minus historical thickness change for GFDL and MPI on one shared diverging scale centred on zero. The third panel is CESM2 ssp245 thickness on its own sequential scale, labelled "no change can be formed for CESM2".

10. Compare the models with the observations.
    - Compare each model's OMZ volumes at 20 and 60, and its basin volumes, against WOA23 with compute_difference_metrics. State which model is closer.
    - State whether GFDL and MPI agree on the sign and size of the scenario change.

Other notes:
- Maps use the shared 1° mask. All volumes, inventories, basin totals and thickness values are on each source's native grid. State this.
- CESM2 has no historical `o2` and no matching historical member, so it is shown only for ssp245. No change or percentage change is given for it.
- The OMZ thresholds are strict (< 20 and < 60 µmol kg-1).
- Report the density used for the unit conversion and the basin option used.
- [common-comparison-mask] The shared mask is used only for the maps and cell-by-cell comparisons, with cell counts reported.
- [historical-ssp245-member-continuity] The scenario change uses 30-year windows from the same member (r1) for GFDL and MPI.
