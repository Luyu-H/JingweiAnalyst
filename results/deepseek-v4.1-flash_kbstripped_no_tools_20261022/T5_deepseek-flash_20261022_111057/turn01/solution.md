The analysis was completed successfully. SSS was defined as the full-cell-volume-weighted mean salinity over model levels whose nominal depths lie in 0–10 m, with no partial-layer weighting and no conversion to Absolute Salinity; model `so` and WOA23 salinity were treated directly as practical salinity in psu. For WOA23, the objectively analysed 1971–2000 annual climatology (`s_an`) was used with the 0, 5 and 10 m levels, weighted by 1° cell area times the level thickness implied by the depth-coordinate gradient. For models, annual-mean SSS was formed from monthly fields, and global/basin means used native-grid `areacello` weights. Mean-state bias, RMSE and correlation were computed after conservative area-weighted regridding of each model’s 1971–2000 climatology to the WOA23 1° grid and applying a single common mask requiring finite WOA23 and all three model fields over WOA23 ocean cells (`landsea_01 > 1`). Trends were fitted separately over historical 1950–2014 and ssp245 2071–2100 with an effective-sample-size correction based on residual autocorrelation lags 1 to `min(n/4, 10)` (up to 10 years for the 65-year historical segment and up to 7 years for the 30-year ssp245 segment). No trend was fitted across the 2014/2015 join.

### 1. Mean-state SSS vs WOA23, 1971–2000

The common WOA23 1° comparison mask contains **21,961 cells** (of 64,800), notably fewer than the individual regridded model fields (CESM2 37,450; GFDL-ESM4 44,564; MPI-ESM1-2-LR 26,313), indicating coastal/marginal-cell disagreement among products. On these shared cells, WOA23 has a mean SSS of **34.70 psu**.

| Model | Model mean SSS (psu) | WOA23 mean SSS (psu) | Bias (psu) | RMSE (psu) | Spatial Pearson *r* |
|---|---:|---:|---:|---:|---:|
| CESM2 | 34.48 | 34.70 | **−0.23** | **1.04** | **0.866** |
| GFDL-ESM4 | 34.70 | 34.70 | **−0.01** | **0.91** | **0.894** |
| MPI-ESM1-2-LR | 34.45 | 34.70 | **−0.25** | **1.03** | **0.872** |

RMSE is area-weighted on the WOA23 grid; the Pearson correlation is the unweighted spatial correlation over the shared cells. GFDL-ESM4 has the smallest mean bias, lowest RMSE and highest pattern correlation. CESM2 and MPI-ESM1-2-LR are both fresher than WOA23 by about 0.23–0.25 psu on the common sample. GFDL-ESM4 is natively supplied on a 1° grid in this archive; CESM2 and MPI were conservatively regridded to the WOA23 grid for this comparison.

### 2. Global-mean SSS trends and ssp245/historical ratios

All values are psu per decade. Historical trends (1950–2014) are not statistically distinguishable from zero in any model: the 95% confidence interval includes zero. The ssp245 2071–2100 trends are all negative and the 95% intervals exclude zero.

| Model | Historical 1950–2014 slope [95% CI], R² | ssp245 2071–2100 slope [95% CI], R² | Ratio ssp/hist |
|---|---|---:|---:|
| CESM2 | **−0.00268** [−0.00626, 0.000906], R²=0.245, not significant | **−0.0209** [−0.0235, −0.0184], R²=0.910, significant | **7.8** |
| GFDL-ESM4 | **−0.00110** [−0.00369, 0.00150], R²=0.064, not significant | **−0.00534** [−0.00877, −0.00191], R²=0.325, significant | **4.9** |
| MPI-ESM1-2-LR | **−0.000687** [−0.00316, 0.00179], R²=0.0107, not significant | **−0.0229** [−0.0335, −0.0122], R²=0.684, significant | **33.3** |

The negative global-mean slopes indicate a freshening of the area-weighted global mean in these models, but only the ssp245 segments have trends distinguishable from zero. The ssp/historical slope ratios are descriptive only and **not meaningful as like-for-like forced-amplification measures**, because they compare a 30-year period with a 65-year period, different forcing levels, and for CESM2 different realisations. Because the historical slopes are non-significant, the ratios should not be interpreted physically.

### 3. End-of-century delta SSS

Delta is ssp245 2071–2100 mean minus historical 1995–2014 mean, area-weighted on each model’s native grid. Basin values use the WOA23 basin mask applied as surface labels through the whole water column, with marginal seas folded into the nearest major ocean label and the Southern and Arctic Oceans excluded. Counts are given as Atlantic / Pacific / Indian.

| Model | Global ΔSSS (psu) | Atlantic ΔSSS (psu) | Pacific ΔSSS (psu) | Indian ΔSSS (psu) | Basin cell counts (A/P/I) |
|---|---:|---:|---:|---:|---:|
| CESM2 | **−0.12** | **−0.01** | **−0.24** | **+0.04** | 25,036 / 40,657 / 16,420 |
| GFDL-ESM4 | **−0.01** | **+0.10** | **−0.07** | **−0.01** | 12,492 / 20,625 / 7,677 |
| MPI-ESM1-2-LR | **−0.09** | **−0.05** | **−0.12** | **−0.05** | 17,326 / 10,590 / 6,239 |

All three models show a negative global-mean end-of-century SSS change under ssp245, strongest in CESM2 and MPI-ESM1-2-LR and weak in GFDL-ESM4. CESM2 and MPI show the strongest Pacific freshening; GFDL-ESM4 instead shows weak Atlantic salinification and weak Pacific/Indian freshening. Basin values are not expected to average exactly to the global value because the Southern and Arctic Oceans are excluded from the basin statistics.

CESM2 uses historical **r1i1p1f1** but ssp245 **r10i1p1f1** (historical r10 is not archived), so its delta and slope ratio include an inter-member internal-variability difference and are not a single continuous realisation. GFDL-ESM4 and MPI-ESM1-2-LR use **r1i1p1f1** for both experiments.

### 4. Figures generated

- **Figure 1** shows three global maps of 1971–2000 SSS bias (model minus WOA23) on the shared comparison cells, one panel per model. The panels share a diverging `RdBu_r` colour scale centred on zero, with the shared limit set from the pooled 98th percentile of absolute bias and the colour bar labelled in psu. The comparison is made on the common WOA23 1° mask, so coastal/marginal cells where the products disagree are excluded.
- **Figure 2** shows three global maps of the cell-by-cell 1950–2014 annual-mean SSS trend in psu per decade, on a regular longitude–latitude frame after conservative regridding for display only. The three panels share a diverging scale centred on zero, with the limit set from the pooled 98th percentile of absolute trend. Cells with p < 0.05 from the effective-sample-size-corrected trend fits are stippled over the trend field.
- **Figure 3** shows the area-weighted global-mean SSS time series for all three models. Historical 1950–2014 and ssp245 2071–2100 are drawn as separate segments with the 2015–2070 gap left empty; each segment has its own fitted line. The y axis is labelled in psu, and each model uses a consistent colour across all three figures.

### 5. Key interpretation

The mean-state comparison shows small global-mean biases but spatially structured errors. GFDL-ESM4 performs best on the common 1° mask by bias, RMSE and pattern correlation; CESM2 and MPI-ESM1-2-LR are both about 0.23–0.25 psu fresher than WOA23 globally. Historical global-mean SSS trends over 1950–2014 are weak and not significant in any model. Under ssp245 over 2071–2100, all models show significant global-mean freshening, but the amplitude differs: MPI-ESM1-2-LR and CESM2 have steeper ssp245 trends than GFDL-ESM4. The end-of-century delta is negative globally in all models, with the strongest global freshening in CESM2 (−0.12 psu) and MPI (−0.09 psu). The basin response is model-dependent, and the CESM2 delta/ratio carries the additional caveat of mixed historical and ssp245 ensemble members.
