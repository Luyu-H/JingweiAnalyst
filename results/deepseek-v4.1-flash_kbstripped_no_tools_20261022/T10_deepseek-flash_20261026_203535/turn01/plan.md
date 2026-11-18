# Plan

1. **WOA23 reference fields.** Open the `woa23` files for the 1971-2000 climate normal: `oxygen`, `percent_oxygen_saturation`, `apparent_oxygen_utilization`, plus the static `landsea_01` land-sea mask and `basinmask_01` basin mask. Build the ocean mask from `landsea_01`; build per-level layer thickness from the WOA23 depth coordinate bounds and cell areas from the 1° regular grid geometry, giving WOA23 cell volumes. Keep every field's own NaN pattern as its valid mask.

2. **CMIP6 fields and static geometry.** Open `cmip6_historical` `o2`, `thetao`, `so` for GFDL-ESM4 and MPI-ESM1-2-LR, monthly 1971-01 to 2000-12, and `cmip6_ssp245` `o2`, `thetao`, `so` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly 2071-01 to 2100-12. Open each model's `volcello`, `areacello` and `deptho`, keeping the model's own grid (GFDL-ESM4 `volcello` is time-varying).
   - Note [cmip6-bgc-units-per-volume]: convert `o2` from mol m-3 to µmol kg-1 using in-situ density from `compute_potential_density_anomaly` with `ref_pressure=None` (not the constant-density default); the inventory step later uses the native mol m-3 `o2` with `convert_concentration_to_content` (Σ c·volcello) and no density conversion.

3. **Time means and oxygen unit conversion.** Average the monthly `thetao` and `so` over each experiment window (1971-2000 historical; 2071-2100 ssp245), compute in-situ density from the averaged temperature and salinity, and convert the averaged `o2` to µmol kg-1 using that density ([cmip6-bgc-units-per-volume]). Do the same conversion for the single WOA23 normal where needed (WOA23 is already µmol kg-1).

4. **Threshold volumes at 20 and 60 µmol kg-1.** For each of the six source-and-experiment combinations (WOA23; GFDL-ESM4 and MPI-ESM1-2-LR historical; CESM2, GFDL-ESM4 and MPI-ESM1-2-LR ssp245) compute the volume of water with `o2` strictly below 20 and strictly below 60, using the exceedance-fraction operation that keeps missing cells out of both numerator and denominator and returns the amount in m3; WOA23 volumes use its layer thickness × cell areas, model volumes use each model's own `volcello`. Report in 1e15 m3 to 3 significant figures and record `area_source`/`volume_source`.
   - Note [missing-values-and-land-cells]: build each field's valid mask from that field itself (masks differ between variables) and mask the weights to the same cells; a hand threshold comparison treats NaN as False and counts land in the denominator, and outputs near 1e20 or 9.97e36 are unmasked fills.
   - Note [omip-weighted-means-and-integrals]: a threshold share is a share of ocean volume (m3), not of cells, and the amount behind the percentage is the returned volume; the AOU/percent-saturation means of step 7 must be volume-weighted with weights restricted to valid OMZ cells.

5. **Dissolved oxygen inventory.** For GFDL-ESM4 and MPI-ESM1-2-LR historical and for all three models in ssp245, sum native mol m-3 `o2` × `volcello` over the valid ocean cells ([cmip6-bgc-units-per-volume]) and report the global ocean inventory in moles to 3 significant figures.

6. **OMZ column thickness at the 60 threshold and its maximum.** For each model compute cell thickness as `volcello`/`areacello`; for WOA23 use the depth-level layer thickness from step 1. Sum the thickness over the cells where `o2` < 60 µmol kg-1 in each column, then locate the column maximum; report the maximum thickness in metres to 1 decimal place with its latitude and longitude to 2 decimal places.
   - Note [cmip6-static-grid-fields]: use `cell_thickness_like` = `volcello`/`areacello` (respects partial bottom cells; not `lev_bnds` differences), sum it over the sub-threshold cells, then locate the maximum; threshold volumes and inventories use each model's own `volcello`, and no model's area is rescaled to another's.

7. **WOA23 OMZ water-mass character.** Within the WOA23 cells where `o2` < 60 µmol kg-1, compute the volume-weighted mean `apparent_oxygen_utilization` and `percent_oxygen_saturation` with weights restricted to those valid OMZ cells ([omip-weighted-means-and-integrals]); report each to 1 decimal place and interpret whether the zone reflects respiration (high AOU, low saturation) or poor ventilation.
   - Note [cmip6-o2sat-and-derived-oxygen]: WOA23 AOU = `o2sat` − `o2` and percent saturation = 100·`o2`/`o2sat`; the volume basis cancels in percent saturation but not in AOU, and AOU is in µmol kg-1.

8. **Basin OMZ volumes at the 60 threshold.** Carry the WOA23 `basinmask_01` labels onto each model grid and combine with the ocean mask, folding marginal seas into their parent ocean (`include_marginal_seas=True`) and using the surface label through the column (`depth_dependent=False`); exclude the Southern and Arctic basins. For WOA23 and the two historical models compute Pacific, Indian and Atlantic OMZ volumes at the 60 threshold on each source's own grid, in 1e15 m3 to 3 significant figures.
   - Note [basin-and-region-definitions]: use the WOA23 basin mask carried onto each grid with `include_marginal_seas=True`, excluding Southern/Arctic as separate basins, and `depth_dependent=False` for named-ocean volumes; state the folding choice.

9. **Scenario change.** For GFDL-ESM4 and MPI-ESM1-2-LR only, compute the percentage change in OMZ volume from historical 1971-2000 to ssp245 2071-2100 at both the 20 and 60 thresholds, positive for expansion, to 2 decimal places; state whether the two models agree in sign and size.
   - Note [historical-ssp245-member-continuity]: GFDL-ESM4 and MPI-ESM1-2-LR historical r1 followed by ssp245 r1 are one continuous trajectory, so their historical-to-ssp245 changes are valid; CESM2 ssp245 is r10 while the local historical is r1, so no CESM2 change or join may be formed and its Fig. 2 panel must be labelled ssp245-only.

10. **Which model is closer to WOA23.** Compare the total OMZ volumes at 20 and 60 µmol kg-1 from step 4 against the WOA23 values and identify the closer historical model, reporting the differences in 1e15 m3.
    - Note [common-comparison-mask]: for any cell-by-cell model-vs-WOA metric, put the fields on one grid and apply one `build_comparison_mask` (`ocean_only=True`, since WOA23 is the reference) to every field before the metric, reporting the mask cell counts; the total OMZ volumes are each source's own-grid integrals, so a difference from WOA is not a shared-mask bias and must be described that way.

11. **Figure 1.** Plot three global maps in one row of the historical OMZ thickness at the 60 threshold (step 6) for GFDL-ESM4, MPI-ESM1-2-LR and WOA23, on one shared sequential colour scale in metres.

12. **Figure 2.** Plot three global maps in one row: the change in OMZ thickness at the 60 threshold (ssp245 2071-2100 minus historical 1971-2000, step 6 fields) for GFDL-ESM4 and MPI-ESM1-2-LR on one shared diverging colour scale centred on zero, and CESM2's ssp245 OMZ thickness on its own sequential scale, labelled to show that no change can be formed for CESM2 ([historical-ssp245-member-continuity]).

Other notes:
- CMIP6 `o2` is native mol m-3, converted to µmol kg-1 with in-situ density before applying the 20/60 µmol kg-1 thresholds; WOA23 oxygen is already µmol kg-1.
- WOA23 `apparent_oxygen_utilization` and `percent_oxygen_saturation` are used as published and equal `o2sat` − `o2` and 100·`o2`/`o2sat`; if any model saturation were compared with WOA23, the difference would include a solubility/pressure definitional component.
- OMZ thickness uses `volcello`/`areacello` for the models and depth-level thickness for WOA23, so the column-thickness definition differs between model and observed maps.
- Total OMZ volumes and inventories are each source's own-grid integrals; model–WOA volume differences are not shared-mask bias, and no area is rescaled between sources.
- CESM2 has no historical `o2` in this record and its ssp245 member differs from the local historical member, so CESM2 appears only as an ssp245 state and never in a change estimate.
