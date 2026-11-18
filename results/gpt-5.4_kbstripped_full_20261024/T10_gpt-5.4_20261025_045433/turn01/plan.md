# Plan

1. Load the WOA23 1971-2000 climatology and masks needed for the observational OMZ baseline: datasets `woa23`; variables `oxygen`, `apparent_oxygen_utilization`, `percent_oxygen_saturation`; tools `load_climatology_data`, `load_landsea_mask`, `load_basin_mask`.
   - Note [cmip6-o2sat-and-derived-oxygen]: use WOA23 AOU and percent oxygen saturation directly for the WOA23 OMZ summary, and keep OMZ thresholding in the field’s own units.

2. Load CMIP6 historical monthly `o2`, `thetao`, `so` for models `GFDL-ESM4` and `MPI-ESM1-2-LR` over 1971-01 to 2000-12, and CMIP6 ssp245 monthly `o2`, `thetao`, `so` for models `CESM2`, `GFDL-ESM4`, `MPI-ESM1-2-LR` over 2071-01 to 2100-12; also load each model’s `volcello`, `areacello` and `deptho`; tools `load_regular_data`, `load_volumes_for_cmip6`.
   - Note [historical-ssp245-member-continuity]: form historical-to-ssp245 percentage changes only for GFDL-ESM4 and MPI-ESM1-2-LR; CESM2 has no historical `o2`, so no scenario change is available.

3. Form the period-mean model oxygen, temperature and salinity fields for 1971-2000 historical and 2071-2100 ssp245 with day-weighted averaging in each model calendar; tool `temporal_aggregate` with `freq='mean'`.
   - Note [monthly-means-are-time-averages]: use day-weighted long-term means from monthly data, and treat OMZ diagnostics from the period-mean oxygen field as mean-state diagnostics.

4. Convert each model’s period-mean `o2` from CMIP6 per-volume units to `µmol kg-1` using density from period-mean `thetao` and `so`; tools `compute_potential_density_anomaly`, `transform_unit`.
   - Note [cmip6-bgc-units-per-volume]: compute in-situ density with `ref_pressure=None`, convert model oxygen to `µmol kg-1` before any thresholding or comparison, use `Σ c·volcello` for oxygen inventories, and use layer thickness from `volcello/areacello` for OMZ thickness.

5. Diagnose global OMZs for each source-and-experiment combination at thresholds strictly below `20 µmol kg-1` and strictly below `60 µmol kg-1`: total OMZ volume, and at `60 µmol kg-1` the column-summed OMZ thickness field plus its global maximum and location; tools `identify_oxygen_minimum_zone`, `locate_extreme`.
   - Note [cmip6-static-grid-fields]: use each source’s own `volcello`, `areacello` and `deptho`, with layer thickness from `volcello/areacello`, and report totals on each source’s own ocean area/volume basis.

6. Compute the total dissolved oxygen inventory for GFDL-ESM4 and MPI-ESM1-2-LR historical and for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR ssp245 from the native-grid model oxygen concentrations and `volcello`; tool `convert_concentration_to_content`.

7. For GFDL-ESM4 and MPI-ESM1-2-LR, compute the percentage change in global OMZ volume from historical to ssp245 at both thresholds using the native-grid totals from step 5.

8. Compute WOA23 basin OMZ volumes at the `60 µmol kg-1` threshold for the Pacific, Indian and Atlantic, and compute the same basin totals for the two historical models after horizontal regridding of each model’s `60 µmol kg-1` OMZ mask/thickness support fields to the WOA23 1° grid and masking by the WOA23 basin definitions; tools `regrid_horizontal`, `select_basin`, `compute_basic_stat`.
   - Note [omip-native-grid-vs-regridded]: keep reported global integrals and inventories on native grids, and use regridded fields only for maps and WOA-grid basin-mask comparisons.

9. Within the WOA23 OMZ at the `60 µmol kg-1` threshold, compute the mean apparent oxygen utilization and mean percent oxygen saturation over the OMZ cells; tool `compute_basic_stat` with weighted mean.

10. Regrid the historical `60 µmol kg-1` OMZ thickness fields from GFDL-ESM4 and MPI-ESM1-2-LR to the WOA23 1° grid, and prepare Figure 1 together with WOA23 OMZ thickness on one shared sequential colour scale; tool `regrid_horizontal`.

11. Compute the `60 µmol kg-1` OMZ thickness change fields for GFDL-ESM4 and MPI-ESM1-2-LR as ssp245 minus historical on each model’s native grid, regrid those change fields to the WOA23 1° grid, and prepare Figure 2 with the two change maps on a shared diverging colour scale centred on zero plus the CESM2 ssp245 OMZ thickness map on its own sequential scale; tool `regrid_horizontal`.

12. Compare how well GFDL-ESM4 and MPI-ESM1-2-LR reproduce the observed global OMZ volume against WOA23 at both thresholds, and summarize whether they agree on the sign and magnitude of the scenario-driven OMZ volume change.

Other notes:
- [common-comparison-mask] Any cell-by-cell comparison or score after regridding must use one common valid-data mask over the fields being compared; do not impose the WOA23 mask on comparisons that do not include WOA23.
