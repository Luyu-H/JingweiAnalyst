The offline density-threshold MLD does not reproduce the published CMIP6 `mlotst` exactly, but it closely tracks its large-scale spatial pattern. Across all six model-month comparisons the Pearson correlation is 0.950–0.998, while mean differences are all negative, meaning the offline operator MLD is generally shallower than the published `mlotst`. The best internal consistency is for CESM2, followed by MPI-ESM1-2-LR, with GFDL-ESM4 least consistent by the stated RMSE criterion.

Methodology used by the executed code:
- Loaded CMIP6 historical monthly `thetao` and `so` over 1971-01 to 2000-12 with depth range 0–1000 m, and monthly `mlotst` plus fixed native-grid `areacello`, using `load_regular_data` for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR.
- Computed potential density anomaly with TEOS-10 via `compute_potential_density_anomaly`, then computed the operator MLD with `compute_mixed_layer_depth` using `method='density'`, `ref_depth=10.0`, and `density_threshold=0.03`.
- For March and September separately, built a common validity mask requiring both operator MLD and published `mlotst` to be non-null in all 30 instances of that calendar month. Applied that mask and formed climatologies with `temporal_aggregate(..., freq='mean')` over the 30 selected months.
- Computed mean difference and RMSE with `compute_difference_metrics`, spatial Pearson correlation with `analyze_correlation(dim=['lat','lon'])` on native-grid common cells, exceedance of |difference| > 20 m by area with `compute_exceedance_fraction`, and valid-area percentages using `areacello` sums. The ranking used the mean of the March and September RMSEs, with mean >20 m area fraction as tie-breaker.
- Figure 1 used `regrid_horizontal(..., target='1deg', method='bilinear')` for display only, then plotted the six difference maps on a common longitude-latitude frame with one shared diverging colour scale centred on zero and clipped symmetrically at the pooled 98th percentile of |difference|.
- Figure 2 used native-grid climatology pairs, plotted on identical logarithmic axes with a dashed y = x line and annotations of bias, RMSE, and correlation. Because the axes are logarithmic, the scatter panels show positive finite values only; the reported metrics were computed on the common valid cells.

Key results:

| Model | Month | Bias, operator − mlotst (m) | RMSE (m) | Pearson r | Area with \|Δ\| > 20 m (%) | Common valid area (% of total ocean) |
|---|---:|---:|---:|---:|---:|---:|
| CESM2 | March | -4.4 | 6.8 | 0.996 | 2.5 | 95.9 |
| CESM2 | September | -2.6 | 5.1 | 0.998 | 0.2 | 97.9 |
| GFDL-ESM4 | March | -7.0 | 14.7 | 0.950 | 5.5 | 96.7 |
| GFDL-ESM4 | September | -5.8 | 12.9 | 0.986 | 6.8 | 98.3 |
| MPI-ESM1-2-LR | March | -2.6 | 5.8 | 0.995 | 1.3 | 96.4 |
| MPI-ESM1-2-LR | September | -2.2 | 12.4 | 0.977 | 3.3 | 97.4 |

Largest absolute departures reported by `locate_extreme`:
- CESM2 March: -205.4 m near 58.1°N, 58.1°W.
- CESM2 September: -51.2 m near 75.5°S, 56.3°W.
- GFDL-ESM4 March: +528.7 m near 60.5°N, 169.5°E.
- GFDL-ESM4 September: +429.7 m near 37.5°S, 139.5°E.
- MPI-ESM1-2-LR March: -433.7 m near 54.0°N, 51.8°W.
- MPI-ESM1-2-LR September: -454.0 m near 64.9°S, 39.5°W.

Ranking from most to least internally consistent, using mean March/September RMSE as the primary criterion:
1. CESM2 — mean RMSE 5.9 m; mean >20 m area fraction 1.3%.
2. MPI-ESM1-2-LR — mean RMSE 9.1 m; mean >20 m area fraction 2.3%.
3. GFDL-ESM4 — mean RMSE 13.8 m; mean >20 m area fraction 6.2%.

Interpretation: the operator and published `mlotst` agree well in spatial pattern, with very high correlations and relatively small mean biases compared with typical MLD magnitudes. However, localized discrepancies can be large, especially at high latitudes and in regions associated with deep convection, marginal seas, or sea-ice influence. These differences are expected because the operator is an offline calculation from monthly mean temperature and salinity, whereas `mlotst` is an online model diagnostic; they can differ because of threshold/reference-depth implementation, vertical discretization and interpolation, sea-ice and brine effects, short-lived convection, temporal averaging before MLD calculation, and masking/missing-value details.

The generated Figure 1 presents the six operator-minus-`mlotst` difference maps in a 3×2 layout, with rows CESM2, GFDL-ESM4, MPI-ESM1-2-LR and columns March and September, using the common regular grid and shared diverging scale. Figure 2 presents the corresponding 3×2 native-grid scatter comparisons on identical logarithmic axes, with the dashed y = x line and panel annotations of bias, RMSE, and correlation. These visualizations are intended to show where and how the offline and online MLD diagnostics diverge, while the quantitative metrics are computed on the native model grids.
