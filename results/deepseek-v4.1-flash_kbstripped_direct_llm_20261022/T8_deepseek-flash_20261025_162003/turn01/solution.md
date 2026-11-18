The request was executed and the figures were produced, but the GLODAP 0–50 m calculation contains a serious land-mask artifact: all-NaN land columns become zero after the depth integration, so the GLODAP surface means are biased low. This contaminates the global means, map biases, and sign comparison. The table below therefore reports the numbers returned by the run, followed by a caution about their interpretation.

### 1. Global 0–50 m means and model–GLODAP differences

All DIC and total alkalinity values are in µmol kg⁻¹; pH is on the total scale.

| Source | DIC mean | DIC diff. | TAlk mean | TAlk diff. | pH mean | pH diff. |
|---|---:|---:|---:|---:|---:|---:|
| GLODAPv2.2016b | 1384 | — | 1569 | — | 5.478 | — |
| CESM2 | 1983 | +599 | 2249 | +680 | 8.025 | +2.547 |
| GFDL-ESM4 | 2059 | +675 | 2336 | +767 | 8.088 | +2.611 |
| MPI-ESM1-2-LR | 1974 | +590 | 2217 | +649 | 8.059 | +2.581 |

These large positive “biases” are not physically credible; they mainly reflect the zero-valued land cells in the GLODAP surface map.

### 2. 0–50 m DIC map comparison on the GLODAP 1° grid

| Model | Mean bias (µmol kg⁻¹) | RMSE (µmol kg⁻¹) | Spatial correlation |
|---|---:|---:|---:|
| CESM2 | +628 | 1098 | 0.391 |
| GFDL-ESM4 | +710 | 1166 | 0.323 |
| MPI-ESM1-2-LR | +675 | 1172 | 0.196 |

Again, these metrics are strongly affected by the land-zero artifact in the GLODAP map.

### 3. Global mean DIC bias versus depth

Bias = model − GLODAP, in µmol kg⁻¹, averaged over cells where both model and GLODAP are defined.

| Depth | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---:|---:|---:|---:|
| 0 m | −28 | +25 | −68 |
| 100 m | −53 | +15 | −68 |
| 500 m | −21 | −6 | −88 |
| 1000 m | −9 | 0 | −84 |
| 2000 m | −13 | −12 | −77 |

This depth-profile comparison is less contaminated because it excludes GLODAP land/missing cells. It shows GFDL-ESM4 generally closest to GLODAP, while MPI-ESM1-2-LR has a large persistent negative DIC bias.

### 4. DIC change from 1971–2000 to 2002 and ratio to surface DIC bias

| Model | 2002 − 1971–2000 mean 0–50 m DIC | Ratio to 0–50 m DIC difference from GLODAP |
|---|---:|---:|
| CESM2 | +9 µmol kg⁻¹ | 0.02 |
| GFDL-ESM4 | +12 µmol kg⁻¹ | 0.02 |
| MPI-ESM1-2-LR | +11 µmol kg⁻¹ | 0.02 |

The exact ratios are about 0.015, 0.018, and 0.018, respectively, all rounding to 0.02. Thus the 2002–(1971–2000) change is only about 1.5–1.8% of the reported surface DIC “bias.” If that bias were interpreted as a real model error, it would be roughly 50–65 times larger than the recent anthropogenic DIC increase, so the mismatch could not be explained by a small temporal offset between model year and climatology. In this run, however, the bias is clearly dominated by the GLODAP land-zero processing error.

### 5. Best model and sign consistency

- The code’s combined score selected **CESM2** as the best model, mainly because it has the lowest DIC map RMSE (1098 µmol kg⁻¹) and highest spatial correlation (0.391). The scores were CESM2 = 4924.5, GFDL-ESM4 = 5219.0, and MPI-ESM1-2-LR = 4991.8.
- This ranking is not scientifically trustworthy because the map/global metrics are corrupted by the GLODAP land-zero problem.
- On the cleaner vertical-profile comparison, **GFDL-ESM4** is closest to GLODAP at most depths, with small DIC biases near the surface and around zero at 500–1000 m; MPI-ESM1-2-LR is the poorest.
- In the contaminated 0–50 m global means, all three models have positive DIC, positive TAlk, and positive pH differences from GLODAP. However, the depth profile shows that near the surface CESM2 and MPI-ESM1-2-LR actually have negative DIC biases while GFDL-ESM4 has a positive bias. Therefore the apparent sign consistency in DIC is an artifact of the surface-mean land contamination.

### Figures

- **Figure 1** shows the three maps of 0–50 m DIC bias with a shared diverging colour scale centred on zero. The panels contain large positive and negative regional features, but the overall positive bias pattern is dominated by comparison against the zero-valued GLODAP land cells.
- **Figure 2** shows the vertical profile of global mean DIC bias. CESM2 is negative at all sampled depths, GFDL-ESM4 is near zero with a small positive surface bias and negative deep bias, and MPI-ESM1-2-LR is strongly negative throughout the water column.

**Bottom line:** the executed calculation reports very large positive surface DIC/TAlk/pH biases and selects CESM2 as “best,” but these results are not a valid evaluation of the model carbon systems because of the GLODAP land-masking error. After discounting that artifact, the depth-profile comparison suggests GFDL-ESM4 reproduces the observed DIC profile best among the three, while MPI-ESM1-2-LR is the worst; the apparent same-sign error structure is not robust.
