## Direct answer

Using the executed workflow, the three CMIP6 models do not all reproduce the GLODAPv2.2016b carbon system equally well. **CESM2 is the best overall choice on the pattern/RMSE and pH metrics, while GFDL-ESM4 has the smallest absolute 0–50 m mean DIC and total-alkalinity biases.** MPI-ESM1-2-LR is the weakest of the three on nearly all metrics. Within each model, the DIC, total-alkalinity and pH offsets have the same sign; across models, GFDL-ESM4 has the opposite sign to CESM2 and MPI-ESM1-2-LR. No requested Omon variable was reported missing for any model; `volcello` was not opened or used.

## Methodology actually used by the code

- **GLODAPv2.2016b mapped climatology:** opened TCO2, TAlk, `pHtsinsitutp`, salinity and temperature. The pH comparison uses the in-situ, total-scale field `pHtsinsitutp`, not `pHts25p0`. GLODAP has no time dimension and 33 standard depth levels. GLODAP 1° cell areas were computed from the lat/lon axes with \(R_E=6{,}371{,}000\) m and constant 1° `dlat`/`dlon` in radians.
- **CMIP6:** monthly `dissic`, `talk`, `ph`, `thetao`, `so` for 2002 and monthly `dissic`, `thetao`, `so` for 1971–2000, plus `areacello`, were used for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. The executed code **did not use `volcello`**; all model weights are `areacello`. CESM2 `areacello` was additionally masked by `sftof > 0`.
- **Density and units:** for each month, in-situ density was computed from that month’s `thetao`/`so` and depth using the UNESCO 1983 EOS-80 formula, with pressure approximated as depth/10 dbar then expressed in bar. Monthly DIC and TAlk were converted to µmol kg⁻¹ as `value * 1e6 / rho`. Only after that conversion were months averaged with `time.dt.days_in_month` weights. For 1971–2000, all selected months were weighted by calendar length. Sub-monthly variance cannot be recovered and is not reported.
- **Global 0–50 m means:** computed on each source’s own grid as one area-weighted mean over all level-cells whose depth ≤50 m. Models use `areacello`; GLODAP uses its computed 1° areas. These are not unweighted lat/lon means.
- **Map comparison:** each horizontal cell’s 0–50 m DIC is the simple mean over that cell’s own levels with depth ≤50 m. Those 2-D fields were put on the GLODAP 1° grid: GFDL-ESM4 was treated as a regular lat/lon grid and linearly interpolated; CESM2 and MPI-ESM1-2-LR were treated as curvilinear and interpolated with `scipy.interpolate.griddata(..., method="linear")`. A single finite-value common mask across GLODAP plus all three models was applied, then area-weighted bias and RMSE and area-weighted spatial correlation were computed.
- **Fixed-depth biases:** model annual-mean DIC and GLODAP TCO2 were vertically interpolated/extrapolated to 0, 100, 500, 1000 and 2000 m with xarray `interp(..., fill_value="extrapolate")`, then horizontally placed on the GLODAP 1° grid as above. At each depth the model-minus-GLODAP bias is the area-weighted mean over exactly the cells where both fields are finite at that depth.
- **Period change:** each model’s 2002 native-grid 0–50 m DIC mean minus its 1971–2000 days-weighted native-grid 0–50 m DIC mean. The ratio is this change divided by the same model’s 0–50 m DIC difference from GLODAP.
- **Caveats:** GLODAP is a 1972–2013 climatology with no time dimension, while the model comparison uses the 2002 annual mean; they are not synchronous. CESM2’s carbonate-system constants are not documented in the model paper, so the CESM2 pH comparison rests on an unverified assumption. No lagged correlation or lag range was calculated.

## Key numerical results

All DIC and TAlk values below are in µmol kg⁻¹ after the per-month in-situ-density conversion. Differences are model minus GLODAP.

### 0–50 m global means and differences

| Source | DIC | DIC − GLODAP | TAlk | TAlk − GLODAP | pH | pH − GLODAP |
|---|---:|---:|---:|---:|---:|---:|
| GLODAP | 2031 | — | 2302 | — | 8.083 | — |
| CESM2 | 1997 | −34 | 2264 | −38 | 8.080 | −0.003 |
| GFDL-ESM4 | 2056 | +25 | 2333 | +31 | 8.088 | +0.006 |
| MPI-ESM1-2-LR | 1975 | −56 | 2218 | −85 | 8.058 | −0.024 |

### 0–50 m DIC map comparison on the GLODAP 1° common mask

| Model | Mean bias (µmol kg⁻¹) | RMSE (µmol kg⁻¹) | Spatial correlation |
|---|---:|---:|---:|
| CESM2 | −33 | 57 | 0.833 |
| GFDL-ESM4 | +28 | 63 | 0.736 |
| MPI-ESM1-2-LR | −54 | 77 | 0.760 |

### Fixed-depth global mean DIC bias (model − GLODAP)

Values are µmol kg⁻¹, with the number of common finite cells in parentheses.

| Model | 0 m | 100 m | 500 m | 1000 m | 2000 m |
|---|---:|---:|---:|---:|---:|
| CESM2 | −28 (40641) | −46 (38407) | −12 (36459) | +2 (35863) | −2 (33672) |
| GFDL-ESM4 | +31 (40034) | +20 (37848) | −11 (35903) | −1 (35173) | −11 (33041) |
| MPI-ESM1-2-LR | −60 (40658) | −59 (38427) | −78 (36477) | −67 (35872) | −61 (33680) |

### Period change of 0–50 m DIC and ratio to the GLODAP offset

| Model | DIC 2002 | DIC 1971–2000 | Change | DIC − GLODAP | Ratio |
|---|---:|---:|---:|---:|---:|
| CESM2 | 1997 | 1987 | +9 | −34 | −0.27 |
| GFDL-ESM4 | 2056 | 2044 | +13 | +25 | +0.50 |
| MPI-ESM1-2-LR | 1975 | 1964 | +11 | −56 | −0.19 |

## Interpretation

- **GFDL-ESM4** has the smallest absolute 0–50 m DIC and TAlk biases (+25 and +31 µmol kg⁻¹) and the smallest absolute map-mean DIC bias (+28 µmol kg⁻¹), but it has the lowest spatial correlation (0.736), a larger RMSE (63 µmol kg⁻¹) than CESM2, and a positive pH offset (+0.006).
- **CESM2** has the best spatial pattern of 0–50 m DIC (correlation 0.833), the lowest RMSE (57 µmol kg⁻¹), and the smallest pH offset (−0.003), but its 0–50 m DIC and TAlk means are low by −34 and −38 µmol kg⁻¹. Its DIC bias is largest at 100 m (−46 µmol kg⁻¹) and becomes near zero at 1000–2000 m.
- **MPI-ESM1-2-LR** is strongly and persistently low in DIC and TAlk: −56 µmol kg⁻¹ DIC and −85 µmol kg⁻¹ TAlk at 0–50 m, with depth biases between about −60 and −78 µmol kg⁻¹. It has the largest RMSE (77 µmol kg⁻¹) and a large negative pH offset (−0.024).
- **Ranking:** CESM2 is the best overall by the combined pattern, RMSE and pH metrics, narrowly ahead of GFDL-ESM4; GFDL-ESM4 is best on absolute mean bias; MPI-ESM1-2-LR is worst overall. No single model wins every metric.
- **Period-change ratio:** the 2002 minus 1971–2000 change is +9, +13 and +11 µmol kg⁻¹ for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. The ratio shows how large that period change is relative to the apparent GLODAP offset. For GFDL-ESM4 the ratio is +0.50: the period change has the same sign and is about half the apparent offset, so about half of its bias could be consistent with period/drift mismatch rather than a purely structural model error. For CESM2 (−0.27) and MPI-ESM1-2-LR (−0.19), the period change has the opposite sign to the offset, so drift/period mismatch does not reduce the apparent bias; structural or other errors dominate the apparent offset. These ratios quantify consistency, not causation, and GLODAP is not synchronous with the 2002 model year.
- **Sign consistency:** within each model, DIC, TAlk and pH differences all share the same sign: CESM2 all negative, GFDL-ESM4 all positive, MPI-ESM1-2-LR all negative. Across the three models the signs are not consistent, because GFDL-ESM4 is positive while CESM2 and MPI-ESM1-2-LR are negative.
- **Missing data:** no requested `dissic`, `talk`, `ph`, `thetao`, `so` or `areacello` field was reported missing for any model. `volcello` was not used.

## Visualization description

- **Figure 1** shows three side-by-side global maps, one for each model, of 0–50 m DIC bias (model minus GLODAP) on the common GLODAP 1° masked grid. All panels share one diverging colour scale centred on zero and labelled in µmol kg⁻¹; the scale is clipped at the 98th percentile of pooled absolute bias.
- **Figure 2** is a single vertical-profile panel: depth on the vertical axis increasing downward, area-weighted global mean DIC bias on the horizontal axis, one line per model, a vertical zero reference line, markers and y-ticks at 0, 100, 500, 1000 and 2000 m.
