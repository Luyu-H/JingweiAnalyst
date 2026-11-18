**Overall answer.** All three CMIP6 models reproduce the broad structure of the observed surface ocean carbon system, but each has systematic biases. Using the combined criteria specified in the code (absolute 0–50 m DIC/TAlk/pH biases, DIC-map RMSE and correlation, and vertical-profile DIC bias), **CESM2 is the best overall**, followed by GFDL-ESM4, with MPI-ESM1-2-LR clearly the worst. However, GFDL-ESM4 has slightly smaller absolute 0–50 m mean DIC and TAlk biases than CESM2; CESM2 wins mainly through better pH, lower DIC-map RMSE, higher spatial correlation, and a better combined score.

**Methodology used (from the executed code).** Monthly CMIP6 `dissic`, `talk`, `ph`, `thetao` and `so` were used for 2002, and `dissic`, `thetao` and `so` for 1971–2000. DIC and TAlk were converted to µmol kg−1 using in-situ density from EOS-80 with `ref_pressure=None` and pressure estimated from depth; no constant-density approximation was used. The 2002 annual means and 1971–2000 mean were formed before the requested comparisons. Global 0–50 m means were volume-weighted using each model’s own `volcello` on its native grid; for GFDL-ESM4 the monthly `volcello` was averaged onto the 2002 annual-mean grid and aligned to the selected months. GLODAP has no `volcello`, so its 0–50 m volume weights were built from its 1° cell area times standard-level thickness, truncated at 50 m. Map comparisons were made on the GLODAP 1° grid after conservative horizontal regridding, using one common map mask with `ocean_only=False`; the mask has **22,210 cells**. Profile comparisons used vertical regridding to 0, 100, 500, 1000 and 2000 m and per-depth shared finite masks, with area weighting by GLODAP 1° cell area.

**0–50 m means and model-minus-GLODAP differences.**

| Source | DIC (µmol kg−1) | ΔDIC | TAlk (µmol kg−1) | ΔTAlk | pH (total scale) | ΔpH |
|---|---:|---:|---:|---:|---:|---:|
| GLODAPv2.2016b | 2032 | — | 2303 | — | 8.082 | — |
| CESM2 | 1997 | −36 | 2264 | −39 | 8.080 | −0.003 |
| GFDL-ESM4 | 2059 | +27 | 2336 | +32 | 8.088 | +0.006 |
| MPI-ESM1-2-LR | 1974 | −58 | 2218 | −86 | 8.059 | −0.024 |

Thus CESM2 and MPI-ESM1-2-LR are low in DIC, TAlk and pH, while GFDL-ESM4 is high in all three. MPI-ESM1-2-LR has by far the largest TAlk and pH biases.

**0–50 m DIC map comparison against GLODAP on the GLODAP 1° grid.**

| Model | Mean bias (µmol kg−1) | RMSE (µmol kg−1) | Spatial correlation |
|---|---:|---:|---:|
| CESM2 | −32 | 60 | 0.799 |
| GFDL-ESM4 | +26 | 79 | 0.611 |
| MPI-ESM1-2-LR | −56 | 88 | 0.660 |

The map bias values differ slightly from the native-grid 0–50 m mean differences because they are computed only on the shared 22,210-cell GLODAP comparison grid and are area-weighted there. CESM2 has the lowest DIC-map RMSE and highest spatial correlation; MPI-ESM1-2-LR has the largest negative map bias and largest RMSE.

**Vertical profile of global mean DIC bias (model minus GLODAP), µmol kg−1.**

| Model | 0 m | 100 m | 500 m | 1000 m | 2000 m |
|---|---:|---:|---:|---:|---:|
| CESM2 | −27 | −46 | −12 | +2 | −2 |
| GFDL-ESM4 | +31 | +19 | −11 | −1 | −11 |
| MPI-ESM1-2-LR | −57 | −60 | −75 | −74 | −77 |

CESM2 is negative near the surface and at 100 m but nearly unbiased at 1000–2000 m. GFDL-ESM4 is positive at the surface and 100 m, slightly negative at 500–1000 m, and negative at 2000 m. MPI-ESM1-2-LR is persistently and strongly negative at all five sampled depths, with biases around −57 to −77 µmol kg−1.

**2002 minus 1971–2000 0–50 m DIC change and ratio to the GLODAP DIC bias.**

| Model | ΔDIC 2002−1971/2000 (µmol kg−1) | Ratio to 0–50 m DIC bias |
|---|---:|---:|
| CESM2 | +9 | −0.26 |
| GFDL-ESM4 | +12 | +0.46 |
| MPI-ESM1-2-LR | +11 | −0.19 |

All models show a positive DIC increase over the period, about 9–12 µmol kg−1. The ratios are small in magnitude. This means the 2002-minus-1971–2000 change is small compared with the standing model-minus-GLODAP DIC bias; because the temporal change is internal to each model, its ratio to the GLODAP bias does not convert the standing GLODAP bias into a model error.

**Sign consistency.** Within each model the DIC, TAlk and pH biases have the same sign: CESM2 is negative for all three, GFDL-ESM4 positive for all three, and MPI-ESM1-2-LR negative for all three. Across models, however, the signs are not consistent for any of DIC, TAlk or pH, because GFDL-ESM4 is opposite to CESM2 and MPI-ESM1-2-LR.

**Best model.** The computed combined score is lowest for CESM2 (3.40), then GFDL-ESM4 (3.73), then MPI-ESM1-2-LR (8.80), so CESM2 is judged best overall. This ranking should still be read with the grid caveat: CESM2 is about 100 km, MPI-ESM1-2-LR about 250 km, and GFDL-ESM4 is on a 1° archive grid with a remapped land–sea mask and only approximately conservative budgets, so narrow features are not equally resolved.

**Figures.** Figure 1 shows three side-by-side global maps of the 0–50 m DIC bias (model minus GLODAP) on the common GLODAP 1° comparison cells, one panel per model. The panels share a diverging colour scale centred on zero and labelled in µmol kg−1; the scale is symmetric and clipped at the pooled 98th percentile of absolute bias for display. Figure 2 shows the vertical profile of global mean DIC bias: depth increases downward, DIC bias is on the horizontal axis, one line is drawn per model, a vertical zero reference line is included, and the five sampled depths (0, 100, 500, 1000 and 2000 m) are marked.
