## Methodological summary

The executed code defined OMZs using strict threshold comparisons (`oxygen < 20` and `< 60 µmol kg⁻¹`) via `identify_oxygen_minimum_zone`. For WOA23, the annual 1971–2000 climatologies of oxygen, apparent oxygen utilization (AOU), percent oxygen saturation, and the WOA23 basin mask were loaded. WOA23 oxygen was converted to µmol kg⁻¹ and volume-weighted exceedance fractions were computed with `compute_exceedance_fraction(..., comparison="less", weighting="volume")`; the WOA volume weighting reported by the tool used WOA23 latitude/longitude and depth bounds. No separate `landsea_01` load appears in the executed code.

For CMIP6, monthly `o2`, `thetao`, and `so` were loaded with `load_regular_data` for:
- historical 1971-01 to 2000-12 for GFDL-ESM4 and MPI-ESM1-2-LR;
- ssp245 2071-01 to 2100-12 for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR.

The monthly fields were reduced to time means with `temporal_aggregate(freq="mean", op="mean")`. Potential density anomaly was computed from time-mean `thetao` and `so` with `compute_potential_density_anomaly(..., ref_pressure=None, temperature_type="potential")` and used to convert oxygen to µmol kg⁻¹. OMZ volumes were obtained from the volume-weighted `exceeding_amount` returned by `compute_exceedance_fraction` and are reported below in 10¹⁵ m³. Maximum thickness and its location were obtained with `locate_extreme` on the 60-threshold column-thickness maps.

Oxygen inventories were computed by applying `convert_concentration_to_content` to the time-mean `o2` field and summing globally in moles. Scenario changes are `(ssp245 − historical) / historical × 100`, positive for expansion. Basin totals used the WOA23 basin mask, with basin-code grouping from name matching plus centroid heuristics for marginal seas, excluding Southern and Arctic codes; for models the basin mask was mapped by nearest selection onto each model grid and depth, then used with each model’s own volume weighting.

## OMZ volumes

Values are in units of 10¹⁵ m³, displayed to about three significant figures.

| Source | Experiment | <20 µmol kg⁻¹ | <60 µmol kg⁻¹ |
|---|---:|---:|---:|
| WOA23 | 1971–2000 climatology | 14.2 | 80.3 |
| GFDL-ESM4 | historical | 43.1 | 105 |
| MPI-ESM1-2-LR | historical | 92.0 | 212 |
| CESM2 | ssp245 | 60.9 | 165 |
| GFDL-ESM4 | ssp245 | 41.3 | 110 |
| MPI-ESM1-2-LR | ssp245 | 89.2 | 211 |

## Maximum 60-threshold OMZ thickness

| Source | Experiment | Max thickness (m) | Latitude | Longitude |
|---|---:|---:|---:|---:|
| WOA23 | climatology | 3125.0 | 7.50 | 121.50 |
| GFDL-ESM4 | historical | 4343.2 | 13.50 | −93.50 |
| MPI-ESM1-2-LR | historical | 4960.0 | −9.39 | −81.22 |
| CESM2 | ssp245 | 3631.9 | 7.45 | 121.44 |
| GFDL-ESM4 | ssp245 | 4318.2 | 13.50 | −93.50 |
| MPI-ESM1-2-LR | ssp245 | 4950.0 | −9.39 | −81.22 |

## Scenario changes for the two models with historical and ssp245 data

| Model | Threshold (µmol kg⁻¹) | Percent change |
|---|---:|---:|
| GFDL-ESM4 | 20 | −4.06% |
| GFDL-ESM4 | 60 | +4.41% |
| MPI-ESM1-2-LR | 20 | −3.03% |
| MPI-ESM1-2-LR | 60 | −0.63% |

At the 20 µmol kg⁻¹ threshold both models show a small decrease in OMZ volume. At 60 µmol kg⁻¹ the models disagree in sign: GFDL-ESM4 expands by +4.41%, while MPI-ESM1-2-LR contracts slightly by −0.63%. Thus the two models do not agree on both the sign and size of the 60-threshold scenario change.

## Global dissolved oxygen inventories

| Model | Experiment | Inventory (mol) |
|---|---:|---:|
| GFDL-ESM4 | historical | 2.20 × 10¹⁷ |
| MPI-ESM1-2-LR | historical | 2.23 × 10¹⁷ |
| CESM2 | ssp245 | 2.02 × 10¹⁷ |
| GFDL-ESM4 | ssp245 | 2.15 × 10¹⁷ |
| MPI-ESM1-2-LR | ssp245 | 2.17 × 10¹⁷ |

CESM2 has no historical oxygen inventory in the executed workflow, so no CESM2 historical-to-ssp245 inventory change is available.

## Basin OMZ volumes at the 60 µmol kg⁻¹ threshold

Units: 10¹⁵ m³. Southern and Arctic basins are excluded; marginal seas are folded into parent ocean groups where the code’s basin-code mapping assigned them.

| Source | Pacific | Indian | Atlantic |
|---|---:|---:|---:|
| WOA23 | 67.4 | 11.6 | 1.29 |
| GFDL-ESM4 historical | 91.7 | 5.34 | 6.80 |
| MPI-ESM1-2-LR historical | 179 | 18.9 | 11.3 |

For WOA23, the Pacific dominates the observed 60-threshold OMZ volume, followed by the Indian and Atlantic. Both models place the largest 60-threshold volume in the Pacific, with MPI-ESM1-2-LR producing much larger Pacific, Indian and Atlantic totals than GFDL-ESM4.

## WOA23 AOU and oxygen saturation inside the 60-threshold OMZ

Volume-weighted mean AOU within the WOA23 60-threshold OMZ is **271.3 µmol kg⁻¹**, and mean percent oxygen saturation is **11.3%**. Under the code’s diagnostic, high AOU together with low oxygen saturation indicates that **respiration-driven oxygen consumption dominates within the OMZ** rather than a low-AOU, poorly ventilated signature.

## Comparison with WOA23 observed OMZ volume

The code found **GFDL-ESM4 closer to WOA23 at both thresholds**. At 60 µmol kg⁻¹, GFDL-ESM4 historical is 105 × 10¹⁵ m³ versus WOA23 80.3 × 10¹⁵ m³, an absolute difference of 24.8 × 10¹⁵ m³ (~30.9%). MPI-ESM1-2-LR historical is 212 × 10¹⁵ m³, an absolute difference of 132 × 10¹⁵ m³ (~164%). At 20 µmol kg⁻¹, both models overestimate strongly, but GFDL-ESM4 remains closer: ~204% high versus ~549% high for MPI-ESM1-2-LR.

## Figures

**Figure 1** shows three global maps in one row of the historical 60-threshold OMZ thickness for GFDL-ESM4, MPI-ESM1-2-LR and WOA23. The maps were regridded bilinearly to the WOA23 plotting grid and plotted on one shared sequential `viridis` colour scale in metres, with the colour maximum clipped at the 98th percentile. The figure is designed to compare the spatial pattern and column-integrated thickness of water below 60 µmol kg⁻¹ among the two historical models and the WOA23 climatology.

**Figure 2** shows three global maps in one row. The first two are the ssp245-minus-historical change in 60-threshold OMZ thickness for GFDL-ESM4 and MPI-ESM1-2-LR, regridded bilinearly to the WOA23 plotting grid and plotted on a shared diverging `RdBu_r` scale centred on zero and clipped at the 98th percentile of absolute change. The third panel shows CESM2’s ssp245 60-threshold OMZ thickness on its own sequential `viridis` scale and is labelled to indicate that no CESM2 historical change can be formed. The figure is designed to show where OMZ thickness expands or contracts in the two change-capable models and to display CESM2’s ssp245 spatial pattern only.

## Important caveats and patterns

- CESM2 contributes only ssp245 results: no CESM2 historical OMZ, no CESM2 historical-to-ssp245 change, and no CESM2 observed-volume comparison can be formed.
- WOA23 is a 1971–2000 climatological normal, whereas the CMIP6 results are 30-year time means (historical 1971–2000; ssp245 2071–2100).
- The GFDL-ESM4 converted oxygen fields include negative minima (−27.17 µmol kg⁻¹ historical and −21.07 µmol kg⁻¹ ssp245). Because the threshold comparison is strictly “less than”, those negative values are counted in the below-threshold volumes and are a quality caveat for the GFDL-ESM4 OMZ volumes.
- The basin grouping uses the WOA23 basin mask and code-mapping heuristics for marginal seas; this is not a purely geometric regrouping and could affect small marginal-sea contributions.
