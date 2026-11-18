The executed analysis compares the CMIP6 historical carbon system from CESM2, GFDL-ESM4 and MPI-ESM1-2-LR with the GLODAPv2.2016b mapped climatology. CMIP6 monthly fields were averaged to a 2002 annual mean and to a 1971–2000 period mean. GLODAP was used as its time-independent 1° mapped climatology with 33 depth levels.

**Methodology actually used**
- `dissic` and `talk` were converted from mol m−3 to µmol kg−1 using `×1.0e6/ρ`, where ρ is GSW in-situ density computed from `thetao` and `so` via `p_from_z`, `SA_from_SP`, `CT_from_pt(SA, thetao)` and `rho(SA, CT, p)`. Thus `thetao` was treated as potential temperature, and `so` as practical salinity converted to absolute salinity before the density calculation.
- `ph` was used directly as the total-scale in-situ pH field (`pHtsinsitutp` for GLODAP); no H+ transformation was applied.
- The 0–50 m column mean was computed over all model/GLODAP levels with depth ≤50 m, with no partial-layer interpolation. Global 0–50 m means were volume-weighted: by `volcello` for models, and by level thickness × cos(latitude) area for GLODAP.
- Map comparisons and depth-level biases used conservative regridding to the GLODAP 1° grid. The map comparison used the single GLODAP 0–50 m DIC mask built with `ocean_only=False` and applied before statistics. Area weighting used cos(latitude). Global means and period changes stayed on each model’s native grid.
- For the depth profile, model DIC was linearly interpolated vertically to 0, 100, 500, 1000 and 2000 m and conservatively regridded. GLODAP values were taken at its exact standard levels; all five requested depths exist in the GLODAP 33-level set. Depth-level biases were averaged over pairwise valid model–GLODAP cells at each depth, so the sample differs by model and is not the same as the map mask.

**Global 0–50 m means and model-minus-GLODAP offsets**

| Source | DIC (µmol kg−1) | ΔDIC | TAlk (µmol kg−1) | ΔTAlk | pH | ΔpH |
|---|---:|---:|---:|---:|---:|---:|
| GLODAPv2.2016b | 2036 | — | 2305 | — | 8.080 | — |
| CESM2 | 1997 | −39 | 2264 | −41 | 8.080 | −0.0005 (≈0.000 at 3 d.p.) |
| GFDL-ESM4 | 2059 | +24 | 2336 | +31 | 8.088 | +0.008 |
| MPI-ESM1-2-LR | 1974 | −61 | 2217 | −87 | 8.059 | −0.021 |

**2002 0–50 m DIC map comparison on the GLODAP 1° grid**

| Model | Mean bias (µmol kg−1) | RMSE (µmol kg−1) | Spatial correlation | Valid cells |
|---|---:|---:|---:|---:|
| CESM2 | −37 | 64 | 0.805 | 41,450 |
| GFDL-ESM4 | +24 | 74 | 0.654 | 41,958 |
| MPI-ESM1-2-LR | −60 | 87 | 0.710 | 41,555 |

**Global mean DIC bias versus depth (model − GLODAP, µmol kg−1)**

| Depth | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---:|---:|---:|---:|
| 0 m | −29 | +28 | −61 |
| 100 m | −46 | +19 | −59 |
| 500 m | −12 | −11 | −78 |
| 1000 m | +2 | −1 | −67 |
| 2000 m | −2 | −11 | −61 |

**2002 minus 1971–2000 change in global 0–50 m DIC and ratio to the 2002 DIC offset**

| Model | ΔDIC 2002−(1971–2000) (µmol kg−1) | Ratio ΔDIC / 2002 DIC offset |
|---|---:|---:|
| CESM2 | +9 | −0.24 |
| GFDL-ESM4 | +12 | +0.53 |
| MPI-ESM1-2-LR | +11 | −0.18 |

The ratio is the model’s 30-year DIC change divided by its 2002 0–50 m DIC offset from GLODAP. A ratio of magnitude ~0.5, as for GFDL-ESM4, means the 1971–2000-to-2002 model trend is about half the size of the 2002 offset, so period mismatch between the model year and the GLODAP climatology can account for a non-negligible part of the apparent bias. For CESM2 and MPI-ESM1-2-LR the 30-year trend is much smaller than, and opposite in sign to, their negative offsets, so the offset is not plausibly explained as a simple period drift. In all cases, the offset should not be read as pure model error without accounting for the climatology period mismatch.

**Best model and sign consistency**
By the combined error score used in the code—sum of absolute DIC offset, absolute TAlk offset, 1000× absolute pH offset, map RMSE, and sum of absolute depth-profile biases—**GFDL-ESM4 is best** (score ≈208.4, versus CESM2 ≈235.3 and MPI-ESM1-2-LR ≈582.7). However, judged only by the 0–50 m DIC map comparison, **CESM2 has the lowest RMSE (64 µmol kg−1) and highest spatial correlation (0.805)**; GFDL-ESM4 has a smaller mean bias but a higher RMSE and lower spatial correlation.

Within each model, the DIC, TAlk and pH offsets have consistent signs: CESM2 is negative for all three, GFDL-ESM4 is positive for all three, and MPI-ESM1-2-LR is negative for all three. Across the three models, the signs are not all the same: CESM2 and MPI-ESM1-2-LR are negative for DIC/TAlk/pH, while GFDL-ESM4 is positive for all three.

**Visualizations**
- **Figure 1** shows three panels, one per model, of the 2002 0–50 m DIC bias (model − GLODAP) on the GLODAP 1° comparison cells. All panels share one diverging colour scale centred on zero, clipped at the pooled 98th percentile of |bias|, with the colour bar labelled in µmol kg−1. It is intended to show the spatial structure of the DIC map bias consistently across models.
- **Figure 2** shows a single vertical profile of the global mean DIC bias at 0, 100, 500, 1000 and 2000 m, with depth increasing downward, bias in µmol kg−1 on the horizontal axis, one line per model, a vertical zero reference line, and the five sampled depths marked. It is intended to show how each model’s DIC bias varies with depth.

Key caveats from the executed implementation are that the 0–50 m pH mean is a volume-weighted arithmetic mean of pH values, not a mean of H+ concentration, and that the 1971–2000 period change is a model trend while GLODAP is a multi-decadal climatology. Only the map comparison and depth-level biases were regridded; global means and period changes used each model’s native grid and own weights.
