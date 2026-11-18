Using 2002 annual means for the CMIP6 models and the GLODAPv2.2016b mapped climatology, with model DIC and total alkalinity expressed in μmol kg⁻¹ and pH on the total scale:

## 1. Global 0–50 m mean carbon-system state

| Source | DIC (μmol kg⁻¹) | TAlk (μmol kg⁻¹) | pH |
|---|---:|---:|---:|
| CESM2 | 1997 | 2264 | 8.080 |
| GFDL-ESM4 | 2059 | 2336 | 8.088 |
| MPI-ESM1-2-LR | 1974 | 2217 | 8.059 |
| GLODAP | 2032 | 2303 | 8.082 |

Model minus GLODAP differences:

| Model | ΔDIC (μmol kg⁻¹) | ΔTAlk (μmol kg⁻¹) | ΔpH |
|---|---:|---:|---:|
| CESM2 | −36 | −39 | −0.003 |
| GFDL-ESM4 | +27 | +32 | +0.006 |
| MPI-ESM1-2-LR | −58 | −86 | −0.024 |

Thus CESM2 and MPI-ESM1-2-LR are lower than GLODAP in both DIC and alkalinity, while GFDL-ESM4 is higher. The pH bias is very small in CESM2 and GFDL-ESM4, but MPI-ESM1-2-LR is substantially too acidic in this 0–50 m mean.

## 2. 0–50 m DIC map comparison on the GLODAP 1° grid

| Model | Mean bias (μmol kg⁻¹) | RMSE (μmol kg⁻¹) | Spatial correlation |
|---|---:|---:|---:|
| CESM2 | −35 | 66 | 0.770 |
| GFDL-ESM4 | +21 | 83 | 0.637 |
| MPI-ESM1-2-LR | −68 | 101 | 0.660 |

CESM2 has the best spatial agreement, with the highest correlation and the smallest RMSE. GFDL-ESM4 has the smallest absolute mean bias on the GLODAP grid, but a poorer spatial pattern than CESM2. MPI-ESM1-2-LR is clearly the worst, with a large negative bias, largest RMSE, and only moderate correlation.

The vertical-profile comparison of global-mean DIC bias gives:

| Depth | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---:|---:|---:|---:|
| 0 m | −28 | +25 | −68 |
| 100 m | −49 | +15 | −68 |
| 500 m | −18 | −6 | −87 |
| 1000 m | −8 | 0 | −81 |
| 2000 m | −12 | −12 | −77 |

(all values μmol kg⁻¹, model minus GLODAP). CESM2 is negative at all sampled depths, GFDL-ESM4 changes from positive near the surface to near-zero/negative at depth, and MPI-ESM1-2-LR is strongly negative throughout.

## 3. Temporal DIC change and interpretation of the bias

Global 0–50 m DIC change, 2002 minus the 1971–2000 average:

| Model | DIC 2002 | DIC 1971–2000 mean | Change 2002−1971–2000 | 2002 DIC bias vs GLODAP | Ratio change / bias |
|---|---:|---:|---:|---:|---:|
| CESM2 | 1997 | 1987 | +9 | −36 | −0.26 |
| GFDL-ESM4 | 2059 | 2047 | +12 | +27 | +0.47 |
| MPI-ESM1-2-LR | 1974 | 1963 | +11 | −58 | −0.19 |

The 1971–2000 to 2002 DIC increase is about 9–12 μmol kg⁻¹, i.e. roughly 19–47% of the magnitude of the model-minus-GLODAP DIC bias. Because GLODAPv2.2016b is a 1972–2013 climatology while the model field is for 2002, the comparison contains some time-window mismatch. Therefore the model-minus-GLODAP DIC offset should not be read as a pure model error. Nevertheless, the offset is generally larger than the recent 1971–2000 to 2002 change, so temporal mismatch alone is unlikely to explain all of it; structural model biases remain important, especially for MPI-ESM1-2-LR.

## 4. Best model and sign consistency

Overall, **CESM2 reproduces the observed carbon system best**: it has the highest 0–50 m DIC spatial correlation, the smallest RMSE, and a small pH bias. GFDL-ESM4 is close and actually has the smallest absolute global-mean DIC and TAlk biases, but its spatial DIC pattern is poorer. MPI-ESM1-2-LR performs worst by a large margin, with strongly negative DIC and alkalinity biases and the largest pH deficit.

The errors are **not consistent in sign across the three models**. GFDL-ESM4 is positive for DIC, TAlk and pH, whereas CESM2 and MPI-ESM1-2-LR are negative for all three. Within each individual model, the DIC, TAlk and pH biases happen to share the same sign in this calculation, but the magnitudes differ strongly, especially for MPI-ESM1-2-LR, whose alkalinity and pH deficits are much larger than those of the other models.

Figure 1 shows the three 0–50 m DIC bias maps on a common diverging colour scale centred on zero, and Figure 2 shows the corresponding global-mean DIC bias vertical profiles with the five sampled depths marked.
