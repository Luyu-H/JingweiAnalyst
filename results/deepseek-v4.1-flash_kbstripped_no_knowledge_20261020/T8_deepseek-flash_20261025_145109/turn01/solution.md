The computation and visualizations completed successfully. The CMIP6 ocean carbon fields were evaluated against the GLODAPv2.2016b mapped climatology for DIC (TCO2), total alkalinity (TAlk), and total-scale in-situ-temperature pH (`pHtsinsitutp`). Model DIC and TAlk were converted to µmol kg−1 using TEOS-10 potential density anomaly computed from each model’s `thetao` and `so`; GLODAP TCO2/TAlk are already in µmol kg−1. Model 2002 fields are annual means from monthly data, and the 1971–2000 fields are long-term means from monthly data. For global 0–50 m means, the code used the products’ own depth levels: GLODAP was thickness-weighted over 0–50 m and then area-weighted; models were volume-weighted over their own 0–50 m levels. Map and vertical-profile comparisons were performed on the GLODAP 1° grid using common defined cells.

### 1. Global 0–50 m means and differences from GLODAP

DIC and TAlk are in µmol kg−1; pH differences are on the total scale.

| Source | DIC | DIC bias | TAlk | TAlk bias | pH | pH bias |
|---|---:|---:|---:|---:|---:|---:|
| GLODAP | 2034 | — | 2302 | — | 8.080 | — |
| CESM2 | 1997 | −37 | 2264 | −38 | 8.080 | −0.001 |
| GFDL-ESM4 | 2059 | +26 | 2336 | +34 | 8.088 | +0.008 |
| MPI-ESM1-2-LR | 1975 | −59 | 2218 | −84 | 8.059 | −0.022 |

CESM2 and MPI-ESM1-2-LR are low in DIC, TAlk, and pH. GFDL-ESM4 is high in all three. Within each model, the signs of the DIC, TAlk, and pH errors are therefore consistent.

### 2. 0–50 m DIC map comparison on the GLODAP 1° grid

| Model | Mean bias (µmol kg−1) | RMSE (µmol kg−1) | Spatial correlation |
|---|---:|---:|---:|
| CESM2 | −32 | 55 | 0.841 |
| GFDL-ESM4 | +29 | 62 | 0.744 |
| MPI-ESM1-2-LR | −54 | 75 | 0.772 |

CESM2 has the highest spatial correlation and the lowest RMSE for the 0–50 m DIC map. GFDL-ESM4 has a smaller absolute global-mean DIC/TAlk bias than CESM2, but poorer spatial pattern skill. MPI-ESM1-2-LR has the largest map bias and RMSE and the largest global-mean carbon-system biases.

### 3. Vertical profile of global-mean DIC bias

Model minus GLODAP, µmol kg−1, on common cells at each depth:

| Depth | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---:|---:|---:|---:|
| 0 m | −28 | +31 | −57 |
| 100 m | −45 | +20 | −58 |
| 500 m | −8 | −6 | −73 |
| 1000 m | +12 | +9 | −57 |
| 2000 m | +18 | +9 | −41 |

MPI-ESM1-2-LR is strongly and persistently low in DIC throughout the sampled water column. CESM2 is low in the upper water column and high at 1000–2000 m. GFDL-ESM4 has relatively modest biases below the surface, with positive near-surface DIC and small positive deep biases.

### 4. 2002 minus 1971–2000 change in global 0–50 m DIC

| Model | DIC change (µmol kg−1) | Ratio to GLODAP 0–50 m DIC bias |
|---|---:|---:|
| CESM2 | +9 | −0.25 |
| GFDL-ESM4 | +12 | +0.49 |
| MPI-ESM1-2-LR | +11 | −0.18 |

A ratio near one would mean that the apparent model–GLODAP DIC bias is comparable to the model’s own recent historical DIC change, so period choice or model drift could be a substantial part of the apparent bias. None of the ratios is near one. GFDL-ESM4’s ratio of +0.49 is same-signed as its positive bias, so roughly half of its positive DIC bias is of the same size as its 2002 minus 1971–2000 change. CESM2 and MPI-ESM1-2-LR have negative ratios because their recent DIC increase is opposite in sign to their negative GLODAP bias; for these two models the 2002-versus-1971–2000 period choice cannot explain the negative DIC bias, which therefore reflects model error or other structural differences.

### 5. Which model performs best?

No model is best in every metric, but CESM2 is the best overall compromise: it has the lowest 0–50 m DIC map RMSE (55 µmol kg−1), the highest spatial correlation (0.841), and the pH mean closest to GLODAP (bias −0.001). GFDL-ESM4 has the smallest absolute global-mean DIC and TAlk biases (+26 and +34 µmol kg−1) and a relatively modest deep vertical DIC bias, but its map correlation is lower (0.744) and its RMSE is larger (62 µmol kg−1). MPI-ESM1-2-LR performs worst overall, with the largest DIC and TAlk biases, the largest pH bias, and the largest map RMSE.

The DIC, TAlk, and pH errors are consistent in sign within each model: CESM2 is uniformly low, GFDL-ESM4 uniformly high, and MPI-ESM1-2-LR uniformly low. The comparison is against a GLODAP climatology rather than a year-2002 observational field, so some temporal representativeness difference remains, but the large MPI biases and the depth-dependent CESM2/GFDL biases are not removed by the model’s own 2002–1971/2000 DIC change.

### 6. Figures

- **Figure 1** shows three side-by-side global maps of 0–50 m DIC bias (model minus GLODAP) on the common GLODAP 1° grid. Each panel is one model, and all panels share a diverging colour scale centred on zero and labelled in µmol kg−1. The scale limit is set from the 98th percentile of the absolute bias across all three model maps. The maps are based on the 0–50 m per-cell DIC mean, conservative horizontal regridding to GLODAP, and a common ocean mask.
- **Figure 2** shows the vertical profile of global-mean DIC bias versus depth. Depth increases downward, DIC bias is on the horizontal axis in µmol kg−1, one line is drawn per model, a vertical zero reference line is included, and the five sampled depths (0, 100, 500, 1000, and 2000 m) are marked. The profile is based on vertical regridding to GLODAP depths, horizontal conservative regridding to the GLODAP grid, and area-weighted model-minus-GLODAP differences over common cells at each depth.
