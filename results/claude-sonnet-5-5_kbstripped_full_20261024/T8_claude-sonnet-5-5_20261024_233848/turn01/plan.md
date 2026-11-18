# Plan

1. Load the GLODAPv2 climatology with load_climatology_data. Variables: TCO2, TAlk, pHtsinsitutp, plus temperature and salinity if needed for the density conversion. Also load its depth levels and layer boundaries with load_depth_levels, and the 1-degree area with load_areas.
   - Note [cmip6-carbonate-system-definitions]: compare model `ph` (total scale) only with GLODAP `pHtsinsitutp`, not `pHts25p0`. Convert model `dissic` and `talk` from mol m-3 to µmol kg-1 before comparing. State that the CESM2 carbonate constants are undocumented (unverified assumption). State that the models differ in their alkalinity terms (MPI-ESM1.2 added phosphate and silicate alkalinity).

2. Load the CMIP6 historical fields with load_regular_data for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR.
   - For calendar year 2002, load monthly `dissic`, `talk`, `ph`, `thetao` and `so`.
   - For 1971-2000, load monthly `dissic`, `thetao` and `so`.
   - Load each model's own areacello and volcello (load_areas, load_volumes_for_cmip6) and its depth levels (load_depth_levels).
   - Report a missing variable as missing for a model; do not substitute.
   - Note [cmip6-static-grid-fields]: layer thickness = volcello / areacello, not the lev_bnds difference. Use this thickness so partial cells are respected. State attrs['area_source'] and attrs['volume_source']. Do not rescale the models' differing ocean areas.

3. Aggregate the models in time with temporal_aggregate.
   - Take the 2002 annual mean (day-weighted) of dissic, talk, ph, thetao and so.
   - Take the 1971-2000 day-weighted long-term mean of dissic, thetao and so.

4. Convert units for DIC and TAlk.
   - Note [cmip6-bgc-units-per-volume]: µmol kg-1 = mmol m-3 × 1000 / ρ. Use the in-situ density from compute_potential_density_anomaly with ref_pressure=None, from thetao and so. Do not use the transform_unit default surface density. State attrs['density_assumption']. Treat the 2002 and 1971-2000 periods with the same treatment. `ph` needs no conversion. GLODAP is already in µmol kg-1.

5. Build the comparison masks with build_comparison_mask, using ocean_only=False.
   - Note [common-comparison-mask]: build ONE mask across the three models and GLODAP for the 0-50 m DIC map and apply it with .where(mask). Report the cell counts. At each of 0, 100, 500, 1000 and 2000 m build a per-model mask from the fields at that depth where that model and GLODAP are both defined, and report the cell counts.
   - Note [missing-values-and-land-cells]: build each field's mask from the field itself (dissic, talk and ph can differ). Use NaN-aware reductions with weights masked to the same valid cells. A value near 1e20 in the output is an unmasked fill value.

6. Compute the global 0-50 m means with compute_basic_stat.
   - Note [omip-weighted-means-and-integrals]: use op='volume_weighted_mean' for DIC, TAlk and pH, giving one mean over all of a source's 0-50 m cells. Use volcello for the models and the layer thicknesses for GLODAP, with weights masked to the valid cells. Never use a plain mean.
   - Report DIC and TAlk in µmol kg-1 to the nearest whole unit, and pH to 3 decimals. Give each model's difference from GLODAP (model minus GLODAP).

7. Evaluate the 0-50 m DIC map.
   - Compute each cell's mean over its own 0-50 m levels, and use the thickness weights from step 2.
   - Regrid each model to the GLODAP 1-degree grid with regrid_horizontal.
   - Apply the common mask from step 5.
   - Use compute_difference_metrics (area-weighted) for the mean bias and RMSE, to the nearest whole unit.
   - Use analyze_correlation with weighting='area' for the spatial correlation, to 3 decimals.
   - Report mean_data, mean_reference and the bias over the same sample.

8. Compute the DIC bias profile.
   - Use regrid_vertical if a model's levels do not match 0, 100, 500, 1000 and 2000 m, and regrid_horizontal to the GLODAP grid.
   - At each depth use the per-model mask from step 5.
   - Compute the area-weighted global mean bias with compute_difference_metrics, to the nearest whole unit.

9. Compute the period change.
   - Compute the global 0-50 m volume-weighted mean DIC (compute_basic_stat) for 2002 minus the 1971-2000 mean, using compute_anomaly_vs_own_climatology or a difference of the two aggregated means. Report it to the nearest whole unit.
   - Compute the ratio of this change to the model's 0-50 m DIC difference from GLODAP from step 6, to 2 decimals.
   - Note [anomaly-baselines-and-climate-indices]: state that the baseline is the model's own 1971-2000 climatology. This does not remove the model bias. Use the ratio to judge whether the offset is a model error or a period or drift effect.

10. Synthesise the results.
   - Rank the models by overall DIC, TAlk and pH errors, e.g. by the absolute 0-50 m differences and the map RMSE and correlation, and name the best model.
   - State whether the signs of the DIC, alkalinity and pH errors agree across the three models.

11. Make the figures.
   - Figure 1: three global maps in one row, one per model, of the 0-50 m DIC bias on the common-mask cells. Use one diverging colour scale centred on zero, shared by all panels, with a colourbar labelled in µmol kg-1.
   - Figure 2: a single panel of the global mean DIC bias profile. Depth is on the vertical axis increasing downward, bias in µmol kg-1 is on the horizontal axis, with one line per model and a vertical zero line. Mark the five sampled depths.

Other notes:
- Models are the 2002 annual mean (a single year), compared with a GLODAP climatology built from 1972-2013 data.
- DIC and TAlk are converted to µmol kg-1 with in-situ density, and the density used is stated.
- The CESM2 carbonate constants are undocumented (unverified assumption).
- The models differ in their alkalinity terms (MPI-ESM1.2 added phosphate and silicate alkalinity).
- pH is compared on the total scale against pHtsinsitutp.
- Each integral covers its own source's ocean area. The areas differ by about 1 % and are not rescaled.
- Coverage shrinks with depth, so the cell counts differ between depths.
- [cmip6-carbonate-system-definitions] Derived fields cover the whole water column.
- [common-comparison-mask] The mask uses ocean_only=False because GLODAP is the reference.
