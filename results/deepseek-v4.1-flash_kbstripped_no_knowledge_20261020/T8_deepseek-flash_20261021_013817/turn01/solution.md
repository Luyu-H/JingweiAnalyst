## Direct answer

The three CMIP6 models were evaluated against the GLODAPv2.2016b annual mapped climatology for 0–50 m DIC, total alkalinity and pH. On the executed comparison, **CESM2 reproduces the observed carbon system best overall**, with the lowest composite error score, the highest 0–50 m DIC map spatial correlation and the lowest DIC map RMSE. **GFDL-ESM4** has smaller absolute 0–50 m DIC and alkalinity mean biases than CESM2, but its DIC/alkalinity/pH biases are all positive, its pH bias is larger, and its spatial correlation is lower. **MPI-ESM1-2-LR** is clearly the poorest of the three.

Within each model, the DIC, total-alkalinity and pH biases are **consistent in sign**:
- CESM2: all negative.
- GFDL-ESM4: all positive.
- MPI-ESM1-2-LR: all negative.

Across models, however, there is no single common sign: GFDL-ESM4 is positive for all three variables, while CESM2 and MPI-ESM1-2-LR are negative for all three.

## Method as executed

- CMIP6 historical monthly `dissic`, `talk`, `ph`, `thetao` and `so` were loaded for 2002, and monthly `dissic`, `thetao` and `so` for 1971–2000, for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. All requested variables were present in the execution; none was reported missing.
- 2002 annual means and 1971–2000 long-term means were formed with `temporal_aggregate(freq='mean', op='mean')`.
- Model DIC and total alkalinity were converted from mol m⁻³ to µmol kg⁻¹ using density calculated from `thetao` and `so` with `compute_potential_density_anomaly(ref_pressure=None)`, then `transform_unit(..., to_unit='umol/kg', density=...)`. This was applied to 2002 `dissic`, 2002 `talk`, and 1971–2000 `dissic`.
- GLODAP annual climatology fields were loaded as `TCO2`, `TAlk` and `pHtsinsitutp`. GLODAP DIC and alkalinity were used in µmol kg⁻¹. pH was compared as supplied against GLODAP `pHtsinsitutp`, with no pH-scale conversion.
- Global 0–50 m means were computed with `compute_basic_stat(..., op='volume_weighted_mean', dim=None)` over the depth slice 0–50 m. The 0–50 m maps were computed with the same operation but along depth: `dim='depth'`.
- Horizontal model fields were regridded conservatively to the GLODAP 1° grid with `regrid_horizontal(..., method='conservative')`. Map and depth comparisons used `build_comparison_mask(..., ocean_only=True)` and only common defined cells.
- Map and depth bias/RMSE statistics used `compute_difference_metrics(..., weighting='area')`; spatial correlation used `analyze_correlation(..., method='pearson', weighting='area')`.
- Vertical comparison used `regrid_vertical(..., method='linear', fill_surface=True)` to 0, 100, 500, 1000 and 2000 m, followed by horizontal regridding to the GLODAP grid.
- The period change is the 2002 global 0–50 m DIC mean minus the 1971–2000 mean; its ratio to the 2002 DIC difference from GLODAP was computed from unrounded values.

## Global 0–50 m means and differences

DIC and total alkalinity in µmol kg⁻¹; pH on the total scale. Differences are model minus GLODAP, computed from unrounded means and rounded only for reporting.

| Model/source | DIC mean | DIC diff | TA mean | TA diff | pH mean | pH diff |
|---|---:|---:|---:|---:|---:|---:|
| CESM2 | 1997 | -39 | 2264 | -40 | 8.080 | -0.001 |
| GFDL-ESM4 | 2059 | +23 | 2336 | +31 | 8.088 | +0.008 |
| MPI-ESM1-2-LR | 1974 | -61 | 2218 | -87 | 8.059 | -0.022 |
| GLODAP | 2036 | — | 2304 | — | 8.080 | — |

CESM2 and MPI-ESM1-2-LR are low in DIC and alkalinity; GFDL-ESM4 is high in both. CESM2 has a near-zero pH bias, GFDL-ESM4 is slightly high in pH, and MPI-ESM1-2-LR is low in pH.

## 0–50 m DIC map comparison

Comparison is on the GLODAP 1° grid, over common defined cells.

| Model | Mean bias (µmol kg⁻¹) | RMSE (µmol kg⁻¹) | Spatial correlation |
|---|---:|---:|---:|
| CESM2 | -36 | 58 | 0.838 |
| GFDL-ESM4 | +26 | 61 | 0.738 |
| MPI-ESM1-2-LR | -58 | 79 | 0.768 |

CESM2 has the highest spatial correlation and lowest RMSE. GFDL-ESM4 has a smaller absolute map-mean bias than CESM2 but a lower correlation. MPI-ESM1-2-LR has the largest absolute map bias and RMSE.

## Vertical profile of global mean DIC bias

DIC bias is model minus GLODAP, in µmol kg⁻¹, at the requested depths, over common cells on the GLODAP grid.

| Depth (m) | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---:|---:|---:|---:|
| 0 | -28 | +31 | -57 |
| 100 | -46 | +20 | -59 |
| 500 | -13 | -11 | -78 |
| 1000 | +2 | -1 | -67 |
| 2000 | -2 | -11 | -61 |

CESM2 has near-zero DIC bias at 1000 and 2000 m and a negative bias in the upper ocean. GFDL-ESM4 is positive near the surface and slightly negative below 500 m. MPI-ESM1-2-LR is strongly negative at all sampled depths, with the largest bias at 500 m.

## Period change in 0–50 m DIC and interpretation

| Model | 2002 minus 1971–2000 DIC (µmol kg⁻¹) | Ratio to 2002 DIC–GLODAP bias |
|---|---:|---:|
| CESM2 | +9 | -0.24 |
| GFDL-ESM4 | +12 | +0.53 |
| MPI-ESM1-2-LR | +11 | -0.18 |

All absolute ratios are below one. Therefore the 2002-minus-1971–2000 DIC change is smaller than the model–GLODAP DIC bias in every case. The biases should not be interpreted as dominated by the model’s transient DIC change. For GFDL-ESM4 the ratio is +0.53, so roughly half of the positive 2002 DIC bias is of the same sign and scale as its period change, but it still does not explain the full bias. For CESM2 and MPI-ESM1-2-LR the ratios are negative, meaning their DIC increase from 1971–2000 to 2002 works opposite to their negative 2002 DIC bias and cannot explain it. Overall, the biases are more consistent with model mean-state/climatology mismatch than with a transient period change alone.

## Which model is best and sign consistency

Using the executed composite score—based on normalized absolute DIC, alkalinity and pH differences, DIC map RMSE, and 1 minus spatial correlation—the ranking is:

1. **CESM2**, score 2.81 — best overall; best DIC map correlation, lowest RMSE, near-zero pH bias.
2. **GFDL-ESM4**, score 3.17 — smaller absolute DIC/alkalinity biases than CESM2 but all positive and lower spatial correlation.
3. **MPI-ESM1-2-LR**, score 6.68 — largest DIC, alkalinity and pH errors and largest vertical DIC biases.

Sign consistency:
- CESM2: DIC negative, TA negative, pH negative — consistent.
- GFDL-ESM4: DIC positive, TA positive, pH positive — consistent.
- MPI-ESM1-2-LR: DIC negative, TA negative, pH negative — consistent.

Thus, **each model has internally consistent sign errors across DIC, alkalinity and pH**, but the sign is not common across all three models.

## Figures generated

- **Figure 1**: three global maps in one row, one per model, showing the 0–50 m DIC bias against GLODAP on the common GLODAP 1° comparison cells. The panels share one diverging colour scale centred on zero, labelled in µmol kg⁻¹, with the scale clipped at the 98th percentile of pooled absolute bias.
- **Figure 2**: a single vertical-profile panel with depth increasing downward, DIC bias in µmol kg⁻¹ on the horizontal axis, one line per model, a vertical zero reference line, and markers/ticks at 0, 100, 500, 1000 and 2000 m.

The map and vertical biases differ slightly from the global-mean DIC differences because the map/depth comparisons are restricted to common GLODAP-grid cells and involve horizontal/vertical regridding.
