## 1. What was done (as executed)

**SSS definition.** SSS is the volume-weighted mean over the 0–10 m levels, using `compute_basic_stat(op="volume_weighted_mean", dim="depth")` on each dataset's own full-cell volumes. Units are psu.
- WOA23 uses the 0, 5 and 10 m levels of the 1971–2000 annual climatology.
- CESM2 uses one level, centred at 5 m (bounds 0–10 m).
- GFDL-ESM4 uses two levels, centred at 2.5 m (0–5 m) and 10 m (5–15 m).
- MPI-ESM1-2-LR uses one level, centred at 6 m (bounds 0–12 m).
- All three models report `area_source = areacello published by cmip6_historical` and `volume_source = volcello published by cmip6_historical`.
- GFDL's volcello carries a time dimension, while CESM2's and MPI's do not.
- No areas were rescaled between models.

**Mean-state bias (1971–2000).**
- Each model's historical SSS was averaged over 1971–2000 with `temporal_aggregate(freq="mean")`.
- It was then conservatively regridded to the WOA23 1° grid with `regrid_horizontal`.
- `build_comparison_mask(ocean_only=True)` was built over the three models and WOA23 and applied with `.where`.
- On that shared sample, `compute_difference_metrics` gave bias, RMSE and the means, and `analyze_correlation(weighting="area")` gave Pearson r.
- The mask keeps **40,829 cells**. The smallest input was WOA23 with 41,088 valid cells. The model fields had 42,576 (CESM2), 44,564 (GFDL) and 42,713 (MPI) valid cells. The land-sea mask removed 0 further cells.
- The reported means are shared-sample (area-weighted) means, not each source's own global ocean mean.

**Global-mean trends.**
- Annual means came from `temporal_aggregate(freq="year")`, then an area-weighted mean over lat/lon.
- `fit_1D_trend(method="linear")` was fitted separately to historical 1950–2014 (65 years) and ssp245 2071–2100 (30 years). Nothing was fitted across 2014/2015.
- Slopes are per year ×10 to give psu per decade.
- The 95% CI and the "significant" flag are the tool's autocorrelation-adjusted values. The effective sample size is far below n: about 9–15 for the historical fits and 5–12 for ssp245.
- The series are the models' own annual means, not anomalies against WOA23.

**Trend maps.**
- `fit_trend_map` was applied on each model's native grid to annual-mean SSS for 1950–2014.
- The slope and the significance field (p < 0.05) were conservatively regridded to the 1° grid for display.
- The significance field was regridded as a 0/1 fraction and thresholded at 0.5.
- The comparison mask was not applied.

**ΔSSS.**
- ΔSSS is the ssp245 2071–2100 mean minus the historical 1995–2014 mean, each from `temporal_aggregate(freq="mean")`.
- The global value is an area-weighted mean over valid cells.
- Basin values use `select_basin(include_marginal_seas=True)`, which applies the WOA23 basin mask, and then `compute_basic_stat(area_weighted_mean)`. The Southern and Arctic Oceans are excluded from the three basins.

## 2. Results

### Mean state vs WOA23 (shared 1° ocean mask, 1971–2000)

| Model | Model mean (psu) | WOA23 mean (psu) | Difference (psu) | RMSE (psu) | Pearson r |
|---|---|---|---|---|---|
| CESM2 | 34.37 | 34.68 | −0.31 | 0.77 | 0.903 |
| GFDL-ESM4 | 34.64 | 34.68 | −0.03 | 0.60 | 0.926 |
| MPI-ESM1-2-LR | 34.41 | 34.68 | −0.27 | 0.76 | 0.893 |

- GFDL's difference is −0.03 rather than −0.04 because the difference was rounded from unrounded values.
- CESM2 and MPI are fresher than WOA23 by about 0.3 psu on average. GFDL has almost no global-mean bias and the lowest RMSE and highest correlation, so on this shared sample it ranks best.
- GFDL's small mean bias is not a sign of small errors. Its centred RMSE (0.60) equals its RMSE, and its bias maps range from about −12 to +11 psu in individual cells, so positive and negative errors cancel.
- Individual-cell biases reach roughly ±10–12 psu in all three models, in a small number of cells.

### Global-mean SSS trends

| Model | Historical slope (psu/dec) [95% CI], R² | CI excludes 0? | ssp245 slope (psu/dec) [95% CI], R² | CI excludes 0? | Ratio ssp/hist |
|---|---|---|---|---|---|
| CESM2 | −0.00268 [−0.00693, 0.00157], 0.245 | No | −0.0209 [−0.0255, −0.0164], 0.910 | Yes | 7.8 |
| GFDL-ESM4 | −0.00110 [−0.00361, 0.00141], 0.064 | No | −0.00534 [−0.0106, −3.86e-05], 0.324 | Yes (marginal) | 4.9 |
| MPI-ESM1-2-LR | −0.000551 [−0.00504, 0.00393], 0.007 | No | −0.0231 [−0.0513, 0.00504], 0.701 | No | 42.0 |

- No model has a historical global-mean SSS trend distinguishable from zero.
- In ssp245 2071–2100, CESM2 and GFDL have significantly negative slopes. GFDL's CI only just excludes zero. MPI's slope is larger in magnitude but its CI includes zero because n_eff is only about 5.
- **The ratios are not meaningful in any model.** The historical slope is statistically indistinguishable from zero in all three, so the denominator is not a robust quantity. The ratios (7.8, 4.9, 42.0) compare two separately fitted segments, not a continuous trajectory.
- For CESM2 the segments also come from different members, so the ratio is further compromised.

### ΔSSS (ssp245 2071–2100 minus historical 1995–2014, psu)

| Model | Global | Atlantic | Indian | Pacific |
|---|---|---|---|---|
| CESM2 | −0.12 | +0.11 | +0.06 | −0.24 |
| GFDL-ESM4 | −0.01 | +0.19 | −0.00 | −0.04 |
| MPI-ESM1-2-LR | −0.09 | −0.04 | −0.05 | −0.11 |

**Ensemble-member caveat.**
- CESM2's historical run is r1i1p1f1 but its ssp245 run is r10i1p1f1, and the parent r10 run is absent locally. Its ΔSSS is therefore the forced change plus the internal-variability difference between two members.
- GFDL-ESM4 and MPI-ESM1-2-LR use r1 in both experiments, so their join is continuous.
- Marginal seas absent from a model have no counterpart.

## 3. Figures

- **Figure 1:** three maps in one row (one per model) of SSS bias, model 1971–2000 minus WOA23, on the shared 40,829-cell 1° mask.
  - All panels use one RdBu_r colour scale centred on zero, with limits set at the 98th percentile of the pooled |bias|.
  - The colour bar is labelled in psu.
  - Panel titles are coloured by model and give each model's mean difference and RMSE.
- **Figure 2:** three maps in one row of the cell-by-cell linear trend of annual-mean SSS for 1950–2014, in psu per decade.
  - Native-grid trend fields were regridded to a regular lon-lat grid for display only.
  - One shared diverging scale is centred on zero, clipped at the 98th percentile of the pooled |trend|.
  - Cells with p < 0.05 are hatched with dots on top of the trend field.
  - The maps are drawn with `pcolormesh` on lon/lat axes with equal aspect, not with a cartographic projection library.
  - No comparison mask was applied.
- **Figure 3:** one panel of global-mean annual SSS (psu) for the three models.
  - Historical 1950–2014 and ssp245 2071–2100 are drawn as separate line segments, so the 2015–2070 gap has no line.
  - Each segment has its own dashed linear fit.
  - Colours are consistent across all three figures: CESM2 blue, GFDL-ESM4 orange, MPI-ESM1-2-LR green.

## 4. Caveats and insights

- The mean-state numbers describe the shared-sample ocean only. They differ from each model's own global ocean mean. For example, the 1950–2014 mean of each model's own annual global-mean SSS is 34.31 (CESM2), 34.56 (GFDL) and 34.36 (MPI).
- Bias, RMSE and r were not used for trends or ΔSSS, which use no comparison mask.
- The historical global-mean changes are small compared with interannual noise and strong autocorrelation. The 30-year ssp245 fits have few effective degrees of freedom, so their CIs are wide.
- Basin ΔSSS signs differ between models (for example, Atlantic +0.11 and +0.19 in CESM2 and GFDL but −0.04 in MPI), so the basin-scale pattern is model-dependent. The CESM2 Pacific and global values include a member-difference contribution.
- Differing resolutions (CESM2 and MPI `gn`, GFDL `gr` at 1°) affect coastal and boundary-current features.
