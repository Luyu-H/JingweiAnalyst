# Plan

1. Open the WOA23 files and read the 1971-2000 climate normal fields `oxygen`, `percent_oxygen_saturation` and `apparent_oxygen_utilization` (1° grid, 33-ish standard levels), together with the static masks `landsea_01` (land-sea) and `basinmask_01` (basin). Verify that `oxygen` is already in µmol kg⁻¹ (no density conversion) and that the vertical axis is depth in metres.

2. Open `cmip6_historical` and read `o2`, `thetao`, `so` for GFDL-ESM4 and MPI-ESM1-2-LR, monthly 1971-01 to 2000-12, together with each model's Ofx fixed fields `volcello`, `areacello` and `deptho`. Confirm that CESM2 publishes no `o2` in the historical experiment (record it as missing, do not substitute).

3. Open `cmip6_ssp245` and read `o2`, `thetao`, `so` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly 2071-01 to 2100-12, with the same `volcello`, `areacello`, `deptho`. Record the variant labels (GFDL-ESM4 and MPI-ESM1-2-LR ssp245 r1i1p1f1; CESM2 ssp245 r10i1p1f1).

4. For every model/experiment combination, average the monthly `o2`, `thetao` and `so` over its period with `temporal_aggregate(freq='mean')` so that month lengths are weighted, giving one period-mean 3-D field per model and experiment.
   - Note [monthly-means-are-time-averages]: the OMZ diagnostic is applied to the period-mean field (average-then-diagnostic), which is the like-for-like choice against the WOA23 climatology; do not threshold each month and then average, and do not use `compute_basic_stat(op='mean', dim='time')` (equal weight per calendar month).

5. Convert each model's `o2` from mol m⁻³ to µmol kg⁻¹: µmol kg⁻¹ = mmol m⁻³ × 1000/ρ, using the in-situ density computed from the period-mean `thetao` and `so` with `compute_potential_density_anomaly` and `ref_pressure=None`; keep and report `attrs['density_assumption']` with the density actually used.
   - Note [cmip6-bgc-units-per-volume]: the total O₂ inventory (step 9) is the exception — it stays in native mol m⁻³, integrated with `convert_concentration_to_content` and `volcello`, with no density applied.

6. Build the OMZ masks per depth level (coverage shrinks with depth) as `o2 < 20` and `o2 < 60` µmol kg⁻¹ and compute each OMZ volume as the sum of `volcello` over the cells below the threshold, separately for WOA23, GFDL-ESM4 historical, MPI-ESM1-2-LR historical, CESM2 ssp245, GFDL-ESM4 ssp245 and MPI-ESM1-2-LR ssp245; report in 1e15 m³ to 3 significant figures together with the cell counts.
   - Note [common-comparison-mask]: each volume is a mean-type sum taken over its own domain and grid — this is not a bias on a shared mask; use `ocean_only=True` because WOA23 is the reference and say that each source is evaluated on its own grid.
   - Note [cmip6-o2sat-and-derived-oxygen]: thresholding must be done on the µmol kg⁻¹ field (the 20 and 60 levels are compared in the field's own units); WOA23 AOU and percent saturation rest on WOA's own solubility and pressure conventions, so differences from the models' OMIP-derived oxygen are partly definitional.

7. At the 60 µmol kg⁻¹ threshold, compute the OMZ thickness of every water column as the sum over levels below the threshold of the cell thickness `volcello/areacello` (partial bottom cells respected), then take the global maximum and report it in metres to 1 decimal place with its latitude and longitude to 2 decimals, for each of the six combinations.
   - Note [cmip6-static-grid-fields]: column thickness comes from `volcello/areacello`, not from `lev_bnds` differences (partial steps differ substantially in GFDL-ESM4 and MPI-ESM1-2-LR); `deptho` is the model bathymetry and total ocean area differs by about a per cent between models, so report the domain each integral actually covers rather than rescaling.

8. Form the scenario change as ssp245 2071-2100 minus historical 1971-2000 for GFDL-ESM4 and MPI-ESM1-2-LR at both thresholds and express it as a percentage of the historical volume (positive = expansion), to 2 decimal places.
   - Note [historical-ssp245-member-continuity]: GFDL-ESM4 and MPI-ESM1-2-LR historical r1i1p1f1 followed by ssp245 r1i1p1f1 form one continuous trajectory, so their difference is valid; CESM2 ssp245 is r10i1p1f1 and historical r10i1p1f1 is absent, so no CESM2 change can be formed and it must be labelled as such.

9. Compute the total dissolved oxygen inventory of the global ocean as the volume integral of the native mol m⁻³ `o2` field over each model's own ocean domain using `convert_concentration_to_content` with `volcello` (no density conversion; see [cmip6-bgc-units-per-volume]), for GFDL-ESM4 and MPI-ESM1-2-LR historical and all three models in ssp245, in moles to 3 significant figures.

10. At the 60 µmol kg⁻¹ threshold, split each OMZ volume into Pacific, Indian and Atlantic using the WOA23 `basinmask_01` via `select_basin`/`load_basin_mask` with `include_marginal_seas=True` (marginal seas folded into the parent ocean) and apply surface labels through the column (`depth_dependent=False`); exclude the Southern and Arctic basins; report for WOA23 and the two historical models in 1e15 m³ to 3 significant figures.
    - Note [basin-and-region-definitions]: the CMIP6 Ofx `basin` field is not in the local archive, so basins come from the WOA23 mask; state the mask used and the `depth_dependent=False` (surface label applied through the column) choice.

11. For WOA23, average `apparent_oxygen_utilization` and `percent_oxygen_saturation` over the cells lying inside WOA23's own 60 µmol kg⁻¹ OMZ (each to 1 decimal place) and interpret the pair: high AOU with low percent saturation indicates the zone is respiration-driven rather than set by poor ventilation alone.

12. Compare the observed (WOA23) OMZ volumes at 20 and 60 µmol kg⁻¹ with the two historical models to say which reproduces the observed volume better, and compare the two models' percentage changes to say whether they agree on the sign and magnitude of the scenario change.

13. Figures, drawn with matplotlib from the fields already computed (no re-gridding to a common grid; each panel on its own model grid):
    (a) Figure 1 — one row of three global maps of historical OMZ thickness at 60 µmol kg⁻¹ (GFDL-ESM4, MPI-ESM1-2-LR, WOA23) sharing one sequential colour scale in metres;
    (b) Figure 2 — one row of three global maps: the thickness change (ssp245 2071-2100 minus historical 1971-2000) at 60 µmol kg⁻¹ for GFDL-ESM4 and MPI-ESM1-2-LR on one shared diverging scale centred on zero, plus CESM2's ssp245 OMZ thickness on its own sequential scale, annotated to state that no change can be formed for CESM2.

Other notes:
- Report each volume in 1e15 m³ to 3 significant figures, thickness to 0.1 m with latitude/longitude to 2 decimals, percentage changes to 2 decimals, inventories in moles to 3 significant figures, and the WOA23 AOU/percent-saturation means to 1 decimal place; carry unrounded values between steps.
- The six combinations compared are WOA23, GFDL-ESM4 and MPI-ESM1-2-LR historical 1971-2000, and CESM2, GFDL-ESM4 and MPI-ESM1-2-LR ssp245 2071-2100; the historical and ssp245 periods are different 30-year windows, so the change is a period difference, not a trend.
- CESM2 appears only in the ssp245 panels and in the ssp245 inventory; it has no historical counterpart and therefore no scenario change.
- OMZ volumes are integrals over each model's own wet domain and its own `volcello`; no interpolation onto a common grid is applied, so the numbers are domain differences as much as model differences.
