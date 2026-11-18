## Scientific summary

All salinities are reported in psu. The WOA23 reference is the 1971–2000 annual climatology. Model mean-state fields are 1971–2000 means from the historical experiment. SSS was formed by selecting the model/observation levels within 0–10 m and applying `compute_basic_stat` with `op="volume_weighted_mean"` along depth; no partial-layer overlap weighting was applied. The levels used were: WOA23 0, 5, 10 m; CESM2 5 m; GFDL-ESM4 2.5 and 10 m; MPI-ESM1-2-LR 6 m. Global and basin means used `compute_basic_stat` with `op="area_weighted_mean"` over latitude/longitude.

### 1. Mean-state SSS, bias against WOA23, RMSE and correlation

WOA23 global mean 1971–2000 SSS: **34.63 psu**.

| Model | 1971–2000 global mean SSS | Difference from WOA23 | RMSE vs WOA23 | Spatial Pearson r |
|---|---:|---:|---:|---:|
| CESM2 | 34.32 | **−0.32** | **0.78** | **0.896** |
| GFDL-ESM4 | 34.56 | **−0.07** | **0.63** | **0.920** |
| MPI-ESM1-2-LR | 34.36 | **−0.27** | **0.77** | **0.902** |

On the common WOA23 comparison cells, the mean biases were also negative: CESM2 −0.31, GFDL-ESM4 −0.04, MPI-ESM1-2-LR −0.27 psu. Thus all three models are globally slightly fresher than WOA23. GFDL-ESM4 has the smallest global-mean bias, the lowest RMSE, and the highest spatial correlation; CESM2 and MPI-ESM1-2-LR are fresher by about 0.27–0.32 psu and have larger RMSE.

### 2. Historical and ssp245 global-mean SSS trends

Trends were fitted with `fit_1D_trend(method="linear")` to annual global-mean SSS series over 1950–2014 and 2071–2100 separately. Values are psu decade⁻¹.

| Model | Historical 1950–2014 slope (95% CI), R², significant? | ssp245 2071–2100 slope (95% CI), R², significant? | ssp/historical ratio | Ratio meaningful? |
|---|---|---|---:|---|
| CESM2 | **−0.00268** (−0.00693, 0.00157), R²=0.245, **no** | **−0.0209** (−0.0255, −0.0164), R²=0.910, **yes** | **7.8** | **No** |
| GFDL-ESM4 | **−0.00110** (−0.00361, 0.00141), R²=0.0639, **no** | **−0.00534** (−0.0106, −0.0000386), R²=0.324, **yes** | **4.9** | **No** |
| MPI-ESM1-2-LR | **−0.000551** (−0.00504, 0.00393), R²=0.00718, **no** | **−0.0231** (−0.0513, 0.00504), R²=0.701, **no** | **42.0** | **No** |

All historical 1950–2014 global-mean trends are small and have 95% confidence intervals that include zero. The ssp245 trends are negative (freshening); they are significant for CESM2 and GFDL-ESM4, but not for MPI-ESM1-2-LR. Because the historical slopes are not significantly nonzero, the ssp245/historical slope ratios are not meaningful; CESM2 additionally compares different historical and ssp245 ensemble members.

### 3. End-of-century SSS change: ssp245 2071–2100 minus historical 1995–2014

| Model | Global | Atlantic | Indian | Pacific |
|---|---:|---:|---:|---:|
| CESM2 | **−0.12** | **0.11** | **0.06** | **−0.24** |
| GFDL-ESM4 | **−0.01** | **0.19** | **0.00** (raw −0.004) | **−0.04** |
| MPI-ESM1-2-LR | **−0.09** | **−0.04** | **−0.05** | **−0.11** |

All three models show a globally fresher end-of-century SSS under ssp245 relative to their late-historical state. Basin responses are mixed: CESM2 freshens most strongly in the Pacific but salinifies the Atlantic and Indian; GFDL-ESM4 salinifies the Atlantic with near-zero Indian and slight Pacific freshening; MPI-ESM1-2-LR freshens all three basins. These basin values were obtained with `select_basin(..., include_marginal_seas=True)` for the Atlantic, Indian and Pacific, excluding the Southern and Arctic Oceans.

**Ensemble caveat:** CESM2 uses historical `r1i1p1f1` but ssp245 `r10i1p1f1`, so its delta and slope ratio include differences between realizations as well as the forced response. GFDL-ESM4 and MPI-ESM1-2-LR use `r1i1p1f1` for both experiments. All results are single-member estimates, not ensemble means.

### 4. Figures

- **Figure 1:** Three global maps in one row of model-minus-WOA23 1971–2000 SSS bias on the common WOA23 1-degree comparison cells. The panels share one zero-centred diverging `RdBu_r` colour scale, clipped at the 98th percentile of absolute bias, with the colour bar labelled in psu. Model names are coloured consistently.
- **Figure 2:** Three longitude-latitude map panels of the cell-by-cell linear trend of annual-mean SSS over historical 1950–2014 in psu decade⁻¹. Trends were computed on each model’s native grid with `fit_trend_map(..., per_decade=True)` and regridded to the WOA23 1-degree grid for display only: slope by conservative regridding, p-value by nearest-source-to-destination regridding. The panels share one zero-centred diverging colour scale, clipped at the 98th percentile, and cells with p < 0.05 are stippled.
- **Figure 3:** A single time-series panel of annual global-mean SSS for all three models. The historical 1950–2014 and ssp245 2071–2100 segments are drawn separately with the 2015–2070 gap left empty. Each segment is overlaid with its own fitted linear line from `fit_1D_trend`, the y-axis is labelled in psu, and model colours are used consistently with the other figures.

### 5. Key interpretation

The models are all slightly fresh in their global mean relative to WOA23, with GFDL-ESM4 clearly closest and best correlated. Historical global-mean SSS trends over 1950–2014 are weak and not statistically distinguishable from zero in these single members. Under ssp245, the late-century global-mean SSS trends are negative; they are significant for CESM2 and GFDL-ESM4 but not for MPI-ESM1-2-LR. The slope ratios should not be interpreted quantitatively because the historical trends are non-significant, and CESM2 mixes different ensemble members between experiments. The end-of-century delta is negative globally for all models, but the basin pattern is model-dependent.
