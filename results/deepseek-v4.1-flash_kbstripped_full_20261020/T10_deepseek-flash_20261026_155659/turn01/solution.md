**Assessment**

The CMIP6 models overestimate the WOA23 oxygen-minimum-zone volume at both thresholds. By integrated OMZ volume, **GFDL-ESM4 is closer to WOA23 than MPI-ESM1-2-LR** at both <20 and <60 µmol kg⁻¹. For the scenario change, the two models agree that the <20 µmol kg⁻¹ volume contracts, but they disagree at the 60 threshold: GFDL-ESM4 expands, while MPI-ESM1-2-LR slightly contracts. CESM2 has no historical `o2`, so no CESM2 historical value or CESM2 scenario change can be formed.

**Methodology used**

- **WOA23:** loaded annual `oxygen`, `apparent_oxygen_utilization`, and `percent_oxygen_saturation` with `load_climatology_data(dataset="woa23", ..., time_slice="annual")`, plus `load_landsea_mask()`, `load_basin_mask()`, `load_depth_levels(dataset="woa23", variable="oxygen", time_slice="annual")`, and `load_areas(dataset="woa23")`.
  - The WOA23 oxygen/AOU/saturation fields in this archive are an **all-period composite**, although the question labels them 1971–2000. This product difference is recorded.
- **CMIP6:** loaded monthly `o2`, `thetao`, and `so` with `load_regular_data`:
  - historical GFDL-ESM4 and MPI-ESM1-2-LR, 1971-01 to 2000-12;
  - ssp245 CESM2, GFDL-ESM4, and MPI-ESM1-2-LR, 2071-01 to 2100-12.
- **Average-first:** each monthly field was reduced to a period mean with `temporal_aggregate(freq="mean", op="mean")`. Threshold, thickness, and AOU diagnostics were computed on these period-mean fields, not averaged after thresholding or integration.
- **Oxygen units for thresholds:** CMIP6 `o2` was converted from mol m⁻³ to µmol kg⁻¹ using cell-by-cell in-situ density from `compute_potential_density_anomaly(ref_pressure=None, temperature_type="in_situ")`. The returned “In-situ density anomaly” was read as sigma and +1000 kg m⁻³ for the conversion. WOA23 oxygen was already in µmol kg⁻¹.
- **OMZ volume:** `compute_exceedance_fraction(..., comparison="less", weighting="volume")` was used at thresholds 20 and 60 µmol kg⁻¹; absolute volumes are the returned `exceeding_amount`. This is strictly below the threshold.
- **OMZ thickness:** `identify_oxygen_minimum_zone(o2_data=..., threshold=60.0)` produced column thickness and the 60-threshold mask; `locate_extreme(..., op="max")` gave the maximum thickness and location.
- **Common-grid comparisons:** CMIP6 oxygen was regridded conservatively to the WOA 1° grid with `regrid_horizontal(target="woa23", method="conservative")`. Native-grid OMZ thickness maps were regridded bilinearly to WOA 1° for the common-grid thickness comparison and figures. A single `build_comparison_mask(..., ocean_only=True, reduce_other_dims="all")` was built from WOA23, GFDL-ESM4 historical, MPI-ESM1-2-LR historical, GFDL-ESM4 ssp245, and MPI-ESM1-2-LR ssp245, then applied to all comparison fields, including CESM2.
  - The common comparison mask retained only **80 of 64,800 horizontal cells** (`N_MASK_CELLS=80`, `N_MASK_TOTAL=64800`), so the spatial metrics and shared-scale maps are based on this restricted intersection.
- **Basin volumes:** `select_basin(..., include_marginal_seas=True, depth_dependent=False)` was applied to each native field for Pacific, Indian, and Atlantic, followed by `compute_exceedance_fraction(threshold=60, comparison="less", weighting="volume")`.
- **WOA within-OMZ means:** AOU and percent oxygen saturation were masked to the WOA 60-threshold OMZ and averaged with `compute_basic_stat(op="volume_weighted_mean", dim=["depth","lat","lon"])`.
- **Scenario change:** \(100 \times (V_{\rm ssp245} - V_{\rm historical}) / V_{\rm historical}\) for GFDL-ESM4 and MPI-ESM1-2-LR at 20 and 60 µmol kg⁻¹; positive means expansion.
- **Inventory:** computed from native mol m⁻³ `o2` with `convert_concentration_to_content(concentration=..., density=1025.0)` and summed to moles.

**Key numerical results**

OMZ volumes and maximum 60-threshold thickness:

| Source | V(<20) [1e15 m³] | V(<60) [1e15 m³] | Max thickness at <60 [m] | Location |
|---|---:|---:|---:|---|
| WOA23 | 14.2 | 80.3 | 3125.0 | 7.50°N, 121.50°E |
| GFDL-ESM4 historical | 43.1 | 105 | 4343.2 | 13.50°N, 93.50°W |
| MPI-ESM1-2-LR historical | 92.0 | 212 | 4960.0 | 9.39°S, 81.22°W |
| CESM2 ssp245 | 60.9 | 165 | 3631.9 | 7.45°N, 121.44°E |
| GFDL-ESM4 ssp245 | 41.3 | 110 | 4318.2 | 13.50°N, 93.50°W |
| MPI-ESM1-2-LR ssp245 | 89.2 | 211 | 4950.0 | 9.39°S, 81.22°W |

Scenario percentage changes in OMZ volume, positive = expansion:

| Model | <20 µmol kg⁻¹ | <60 µmol kg⁻¹ |
|---|---:|---:|
| GFDL-ESM4 | -4.06% | +4.41% |
| MPI-ESM1-2-LR | -3.03% | -0.62% |

Global dissolved oxygen inventory [moles, 3 s.f.]:

| Source | Inventory |
|---|---:|
| GFDL-ESM4 historical | 2.20×10¹⁷ |
| MPI-ESM1-2-LR historical | 2.23×10¹⁷ |
| CESM2 ssp245 | 2.02×10¹⁷ |
| GFDL-ESM4 ssp245 | 2.15×10¹⁷ |
| MPI-ESM1-2-LR ssp245 | 2.17×10¹⁷ |

60-threshold OMZ volume by WOA basin [1e15 m³]:

| Source | Pacific | Indian | Atlantic |
|---|---:|---:|---:|
| WOA23 | 67.4 | 11.6 | 1.29 |
| GFDL-ESM4 historical | 92.8 | 5.43 | 6.39 |
| MPI-ESM1-2-LR historical | 180 | 19.2 | 11.3 |

WOA23 volume-weighted means within the 60-threshold OMZ:
- Mean apparent oxygen utilization: **271.3 µmol kg⁻¹**
- Mean percent oxygen saturation: **11.3%**

These values indicate strongly oxygen-depleted, high-AOU water inside the observed OMZ. That is consistent with substantial respiration-driven oxygen consumption, rather than the zone being set only by poor ventilation.

**Model–WOA comparison and scenario interpretation**

- At <20 µmol kg⁻¹: GFDL-ESM4 historical is high by 28.9×10¹⁵ m³ (+204.1% relative to WOA23); MPI-ESM1-2-LR historical is high by 77.8×10¹⁵ m³ (+549.0%).
- At <60 µmol kg⁻¹: GFDL-ESM4 historical is high by 24.8×10¹⁵ m³ (+30.9%); MPI-ESM1-2-LR historical is high by 131.9×10¹⁵ m³ (+164.2%).
- Therefore, **by integrated OMZ volume, GFDL-ESM4 reproduces WOA23 better** at both thresholds.
- On the restricted common-grid 60-thickness comparison (only 80 cells):
  - GFDL-ESM4: bias = -220.52 m, MAE = 220.52 m, RMSE = 373.70 m, maximum absolute error = 887.54 m at 52.50°N, 162.50°E.
  - MPI-ESM1-2-LR: bias = -64.03 m, MAE = 115.70 m, RMSE = 240.90 m, maximum absolute error = 750.00 m at 29.50°N, 147.50°E.
  - Thus, on this limited common mask the thickness-pattern metrics favour MPI-ESM1-2-LR, while the total OMZ volumes favour GFDL-ESM4. The mask is very small, so the spatial metrics are not robust.
- The two models agree on the sign at <20 µmol kg⁻¹: both show contraction under ssp245. They disagree at <60 µmol kg⁻¹: GFDL-ESM4 shows a modest expansion (+4.41%), while MPI-ESM1-2-LR shows a small contraction (-0.62%). The scenario changes are small (within about 5%) and threshold-dependent.

**Visualizations**

- **Figure 1** shows three global maps in one row of historical 60-threshold OMZ thickness for GFDL-ESM4, MPI-ESM1-2-LR, and WOA23 on the common WOA 1° grid, using one shared sequential scale in metres (pooled 98th-percentile upper limit). It is intended to compare the spatial pattern and column thickness of the historical OMZ on the same colour scale.
- **Figure 2** shows three global maps in one row. The first two are ssp245-minus-historical 60-threshold OMZ thickness changes for GFDL-ESM4 and MPI-ESM1-2-LR on a shared diverging scale centred on zero. The third panel shows CESM2’s ssp245 OMZ thickness alone on its own sequential scale and is labelled to indicate that no historical CESM2 `o2` exists, so no CESM2 change can be formed. All panels use the common comparison grid/mask.

**Important caveats**

- CESM2 contributes only ssp245 OMZ volume, thickness, and inventory; it has no historical `o2` in the available archive here, so no CESM2 historical value or scenario-change percentage is reported.
- The WOA23 fields are the archive’s all-period composite, not a separately constructed 1971–2000 normal.
- The common comparison mask retained only 80 horizontal cells, so Figure 1/2 and the spatial bias/RMSE metrics represent that restricted intersection, not the full global ocean.
