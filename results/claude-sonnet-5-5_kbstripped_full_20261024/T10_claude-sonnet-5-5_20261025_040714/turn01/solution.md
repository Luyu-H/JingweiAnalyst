## Direct answer

Both CMIP6 models overestimate the observed OMZ volume. GFDL-ESM4 is much closer to WOA23 than MPI-ESM1-2-LR. Under ssp245 (2071-2100 minus 1971-2000), the two models agree on the sign of the change in the <20 µmol/kg volume, a small decrease. They disagree on the sign of the change in the <60 µmol/kg volume, and all changes are small (about 0.6-4.4 %). CESM2 has no historical `o2`, so only its ssp245 state is reported and no change is formed for it.

## Methods as executed

- **Model data:** `load_regular_data` loaded monthly `o2`, `thetao` and `so` on each model's native grid. Values ≥1e19 were set to NaN as fill. `temporal_aggregate(freq="mean")` took the 30-year long-term mean over 1971-01 to 2000-12 (historical: GFDL, MPI) or 2071-01 to 2100-12 (ssp245: CESM2, GFDL, MPI). The OMZ criteria were therefore applied to 30-year mean fields, not to monthly values.
- **WOA23:** oxygen, apparent oxygen utilization (AOU) and percent saturation were used as given from `load_climatology_data`, on the native 1° grid with 102 depth levels. The landsea and basin masks came from `load_landsea_mask` and `load_basin_mask`.
- **Unit conversion:** the models' `o2` (mol m⁻³) was converted to µmol kg⁻¹ with `transform_unit`, using a cell-by-cell in-situ density from `compute_potential_density_anomaly(thetao, so, ref_pressure=None)`. The tool's `density_assumption` attribute reads: "density field, read as an anomaly (sigma, +1000) … each cell's own in-situ pressure; converted at each cell's own density". Density ranged from about 996 to 1057 kg m⁻³ for GFDL and from about 1003 to 1054 kg m⁻³ for MPI. The conversion is mmol m⁻³ × 1000 / ρ, not a constant density. Thresholds were applied only after conversion.
- **OMZ volumes:** `compute_exceedance_fraction` with strict "less than" at 20 and 60 µmol/kg and volume weighting. Weights were `volcello` for the models and the WOA23 `lat_bnds`/`lon_bnds` × `depth_bnds` volume for WOA23.
- **Grid sources:** the model volume source was "volcello published by cmip6_historical" (or `cmip6_ssp245`), and the area source was `areacello`. For GFDL ssp245 the area source is reported as "areacello published by cmip6_ssp245 (borrowed: this experiment publishes none)", exactly as the tool printed it.
- **Native grids:** volumes, thickness, basin totals and inventories were computed on each source's native grid and levels. No interpolation to the WOA levels and no rescaling to a common ocean area.
- **Thickness and location:** `identify_oxygen_minimum_zone(threshold=60)` gave the per-column OMZ thickness. `locate_extreme(op="max")` gave its location.
- **Inventory:** `convert_concentration_to_content` on the native `o2` (mol m⁻³) using each model's own `volcello`, summed with `compute_basic_stat(op="sum")`.
- **Basins:** `select_basin` on `basinmask_01` with `include_marginal_seas=True` and `depth_dependent=False`. Only Pacific, Indian and Atlantic were kept; Southern and Arctic were excluded.
- **Scenario change:** (ssp245 − hist)/hist × 100 on the volume totals, for GFDL and MPI only. No trend was fitted across the 2014/2015 join.

## Key results

**OMZ volume (10¹⁵ m³, strictly below the threshold)**

| Source | <20 | <60 |
|---|---|---|
| WOA23 | 14.2 | 80.3 |
| GFDL-ESM4 hist | 43.1 | 105 |
| MPI-ESM1-2-LR hist | 92.0 | 212 |
| CESM2 ssp245 | 60.9 | 165 |
| GFDL-ESM4 ssp245 | 41.3 | 110 |
| MPI-ESM1-2-LR ssp245 | 89.2 | 211 |

**Maximum OMZ thickness at 60 µmol/kg**

| Source | Max thickness (m) | Latitude | Longitude |
|---|---|---|---|
| WOA23 | 3125.0 | 7.50 | 121.50 |
| GFDL hist | 4343.2 | 13.50 | −93.50 |
| MPI hist | 4960.0 | −9.39 | −81.22 |
| CESM2 ssp245 | 3631.9 | 7.45 | 121.44 |
| GFDL ssp245 | 4318.2 | 13.50 | −93.50 |
| MPI ssp245 | 4950.0 | −9.39 | −81.22 |

The MPI latitude and longitude are cell-centre coordinates on its curvilinear grid (j=129, i=78). Each maximum was a single cell, with no ties.

**Percentage change in OMZ volume, historical to ssp245**

| Model | <20 | <60 |
|---|---|---|
| GFDL-ESM4 | −4.06 % | +4.41 % |
| MPI-ESM1-2-LR | −3.03 % | −0.63 % |

**Global dissolved O₂ inventory (mol)**

- GFDL hist: 2.20×10¹⁷
- MPI hist: 2.23×10¹⁷
- CESM2 ssp245: 2.02×10¹⁷
- GFDL ssp245: 2.15×10¹⁷
- MPI ssp245: 2.17×10¹⁷

By my arithmetic from these values, the inventory falls by about 2.4 % (GFDL) and 2.3 % (MPI) between the two periods.

**Basin OMZ volume at <60 (10¹⁵ m³; include_marginal_seas=True, depth_dependent=False)**

| Source | Pacific | Indian | Atlantic |
|---|---|---|---|
| WOA23 | 67.4 | 11.6 | 1.29 |
| GFDL hist | 92.8 | 5.43 | 6.39 |
| MPI hist | 180 | 19.2 | 11.3 |

**WOA23 OMZ interior (cells with O₂ < 60 µmol/kg, volume-weighted mean)**

- Mean AOU is 271.3 µmol/kg.
- Mean percent saturation is 11.3 %.
- This is a very large oxygen debt, with only about 11 % saturation remaining. The water has accumulated a strong respiration (remineralisation) signal in poorly ventilated water. The zone is therefore set by respiration acting on weakly ventilated water, not by ventilation alone.
- Caveat: AOU and %sat both depend on the O₂ solubility formulation (WOA23 uses Garcia & Gordon), so they are not independent. AOU also assumes surface saturation at water-mass formation, so preformed disequilibrium biases it.

## Model versus observations

Errors relative to WOA23, from the code's bias and relative-error calculation on native-grid totals:

| Metric | GFDL-ESM4 | MPI-ESM1-2-LR |
|---|---|---|
| <20 volume | +204 % | +549 % |
| <60 volume | +30.9 % | +164 % |
| Pacific <60 | +37.7 % | +167 % |
| Indian <60 | −53.2 % | +65.4 % |
| Atlantic <60 | +397 % | +778 % |
| Mean absolute relative error, 5 volumes | ≈145 % | ≈345 % |

- **Closer model:** GFDL-ESM4, on both the mean error over all five volumes and the <60 volume alone.
- **Basin pattern:** GFDL underestimates the Indian OMZ but overestimates the Pacific and Atlantic. MPI overestimates all three basins, and the Atlantic is by far the worst case in relative terms for both models.
- **Map comparison on the shared 1° mask** (40,893 cells, thickness at 60):

  | Metric | GFDL | MPI |
  |---|---|---|
  | Bias | +68.6 m | +373 m |
  | RMSE | 493 m | 853 m |
  | Centred RMSE | 488 m | 768 m |
  | Std ratio | 1.56 | 2.29 |

- **Scenario change:** both models give a small decrease in the <20 volume, so the sign agrees and the sizes are comparable (−4.1 % and −3.0 %). At <60, GFDL expands by 4.4 % while MPI contracts by 0.6 %, so the sign disagrees. Both changes are small compared with the inter-model and model-observation spread.
- **CESM2:** its ssp245 volumes (60.9 and 165) lie between GFDL and MPI. Its inventory is the lowest of the three. Because it has no historical `o2` or matching historical member, none of this can be read as a scenario change.

## Visualizations

- **Figure 1:** three maps in one row of historical OMZ thickness at <60 µmol/kg (GFDL-ESM4 and MPI-ESM1-2-LR for 1971-2000, and WOA23). The model thickness fields were conservatively regridded to the WOA23 1° grid (`regrid_horizontal`). One shared ocean mask from `build_comparison_mask(ocean_only=True)` was applied with `.where`, leaving 40,893 cells. All three panels share one sequential colour scale in metres, from 0 to the pooled 98th percentile, with the top end clipped.
- **Figure 2:** three maps in one row.
  - The first two panels show the ssp245 (2071-2100) minus historical (1971-2000) thickness change for GFDL and MPI. The difference was formed on the native grids and then regridded to 1°.
  - They share one diverging colour scale centred on zero (±98th percentile of the pooled absolute change).
  - The third panel shows CESM2 ssp245 thickness on its own sequential scale and is labelled "no change can be formed for CESM2".
  - The shared mask covers 40,829 cells.

## Caveats

- Maps and cell-by-cell metrics use the shared 1° mask. All volumes, inventories, basin totals and thickness maxima are on each source's native grid.
- OMZ thresholds were applied to 30-year mean fields. Overestimation by the models may partly reflect this, but the code did not test it.
- The ssp245 change comes from different 30-year windows of the same member (r1), so it includes internal variability.
