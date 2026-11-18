**Methodology as executed.** WOA23 oxygen, apparent oxygen utilization (AOU) and percent oxygen saturation were loaded as annual 1971–2000 climate normals with `load_climatology_data`; WOA23 land–sea mask, basin mask, depth levels/bounds and cell areas were also loaded. WOA23 fields were masked with `landsea > 1`. CMIP6 monthly `o2`, `thetao` and `so` were loaded for GFDL-ESM4 and MPI-ESM1-2-LR historical (1971-01 to 2000-12) and for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR ssp245 (2071-01 to 2100-12); CESM2 historical `o2` was not used. Model `volcello` was loaded with `load_volumes_for_cmip6` and time-averaged where it had a time dimension. In-situ density was computed with `compute_potential_density_anomaly(..., ref_pressure=None, temperature_type="auto")`, and oxygen was converted from mol m⁻³ to µmol kg⁻¹ with `transform_unit(..., density=rho)`. Monthly converted oxygen, temperature, salinity and volume were averaged over each period with `temporal_aggregate(freq="mean")`, and OMZ diagnostics were formed from those period-mean fields. OMZ volumes used `compute_exceedance_fraction(..., comparison="less", weighting="volume")` at strict thresholds of 20 and 60 µmol kg⁻¹ and are reported in 10¹⁵ m³. Thickness diagnostics used `identify_oxygen_minimum_zone(..., threshold=...)`; maxima were located with `locate_extreme(..., op="max")`. The executed thickness calls did not explicitly pass cell thickness/volume arrays. Basin volumes used `select_basin(... include_marginal_seas=True, depth_dependent=False)` for Pacific, Indian and Atlantic only, then the same volume exceedance calculation at 60 µmol kg⁻¹. Scenario and model–WOA percentage differences were computed directly from unrounded volumes; `compute_difference_metrics` was used for the Figure 2 display differences after horizontal regridding to WOA23.

**Total water volume with O₂ strictly below the thresholds** (10¹⁵ m³, 3 s.f.):

| Source / experiment | < 20 µmol kg⁻¹ | < 60 µmol kg⁻¹ |
|---|---:|---:|
| WOA23 1971–2000 | 14.2 | 80.3 |
| GFDL-ESM4 historical 1971–2000 | 43.1 | 105 |
| MPI-ESM1-2-LR historical 1971–2000 | 92.0 | 212 |
| CESM2 ssp245 2071–2100 | 60.9 | 165 |
| GFDL-ESM4 ssp245 2071–2100 | 41.3 | 110 |
| MPI-ESM1-2-LR ssp245 2071–2100 | 89.2 | 211 |

**Maximum OMZ thickness at < 60 µmol kg⁻¹**:

| Source / experiment | Max thickness (m) | Latitude | Longitude |
|---|---:|---:|---:|
| WOA23 1971–2000 | 3125.0 | 7.50 | 121.50 |
| GFDL-ESM4 historical | 4343.2 | 13.50 | −93.50 |
| MPI-ESM1-2-LR historical | 4960.0 | −9.39 | −81.22 |
| CESM2 ssp245 | 3631.9 | 7.45 | 121.44 |
| GFDL-ESM4 ssp245 | 4318.2 | 13.50 | −93.50 |
| MPI-ESM1-2-LR ssp245 | 4950.0 | −9.39 | −81.22 |

**Historical-to-ssp245 OMZ volume change** (positive = expanding):

| Model | < 20 µmol kg⁻¹ | < 60 µmol kg⁻¹ |
|---|---:|---:|
| GFDL-ESM4 | −4.06% | +4.41% |
| MPI-ESM1-2-LR | −3.03% | −0.63% |

Thus both models shrink the strict <20 µmol kg⁻¹ volume, but they disagree at <60 µmol kg⁻¹: GFDL-ESM4 expands by 4.41%, while MPI-ESM1-2-LR contracts slightly by 0.63%.

**Global dissolved oxygen inventory** (moles, 3 s.f., as computed by the executed content conversion and time averaging):

| Source / experiment | O₂ inventory (mol) |
|---|---:|
| GFDL-ESM4 historical | 2.20 × 10¹⁷ |
| MPI-ESM1-2-LR historical | 2.23 × 10¹⁷ |
| CESM2 ssp245 | 2.02 × 10¹⁷ |
| GFDL-ESM4 ssp245 | 2.15 × 10¹⁷ |
| MPI-ESM1-2-LR ssp245 | 2.17 × 10¹⁷ |

The executed inventory calculation used `convert_concentration_to_content(concentration=o2, density=1025.0)` on monthly raw oxygen, then `temporal_aggregate(freq="mean")` and a global sum; it did not use an explicit `o2 × volcello` expression in the executed code.

**OMZ volume at < 60 µmol kg⁻¹ by basin** (10¹⁵ m³, 3 s.f.; WOA23 basin definition, marginal seas included in parent oceans, Southern and Arctic excluded):

| Source / experiment | Pacific | Indian | Atlantic |
|---|---:|---:|---:|
| WOA23 1971–2000 | 67.4 | 11.6 | 1.29 |
| GFDL-ESM4 historical | 92.8 | 5.43 | 6.39 |
| MPI-ESM1-2-LR historical | 180 | 19.2 | 11.3 |

**WOA23 oxygen-status within the observed <60 µmol kg⁻¹ OMZ:** volume-weighted mean AOU = **271.3 µmol kg⁻¹**; volume-weighted mean percent oxygen saturation = **11.3%**. The very high AOU and very low saturation imply substantial oxygen consumption by respiration, rather than a zone merely defined by low saturation with little respiratory drawdown. Poor ventilation can still allow the respired signature to accumulate, but the observed AOU strongly supports respiration as the dominant oxygen sink within the OMZ.

**Model–observation comparison and scenario agreement.** Both historical models overestimate the WOA23 OMZ volume at both thresholds. Relative to WOA23, GFDL-ESM4 is high by **204.1%** at <20 and **30.9%** at <60; MPI-ESM1-2-LR is high by **549.0%** and **164.2%**, respectively. GFDL-ESM4 therefore reproduces the observed OMZ volume better than MPI-ESM1-2-LR (mean absolute relative difference across thresholds 117.5% versus 356.6%). The two models agree on the sign of the change only at <20 µmol kg⁻¹ (both shrink); at <60 µmol kg⁻¹ they disagree, and the magnitudes also differ. CESM2 is shown only as an ssp245 state because its historical `o2` is unavailable/no matching historical baseline exists; no CESM2 change is formed.

**Visualizations.** Figure 1 shows three global maps of historical <60 µmol kg⁻¹ OMZ thickness, one for GFDL-ESM4, one for MPI-ESM1-2-LR and one for WOA23, using one shared sequential viridis colour scale in metres, with the CMIP6 thickness fields bilinearly regridded to WOA23 for display and the scale clipped at the pooled 98th percentile. Figure 2 contains three panels: GFDL-ESM4 and MPI-ESM1-2-LR maps of ssp245 minus historical <60 µmol kg⁻¹ thickness, both regridded to WOA23 and plotted on a shared diverging `RdBu_r` scale centred on zero and clipped at the pooled 98th percentile of absolute change; the third panel shows CESM2’s ssp245 <60 µmol kg⁻¹ thickness on its own sequential viridis scale and is labelled as an ssp245 state with no historical baseline and no formable change.
