## 1. Direct answer

Using the WOA23 1971–2000 annual objectively-analysed salinity climatology as the reference, sea-surface salinity (defined as the full-cell-volume-weighted mean of model/observed levels within 0–10 m) was assessed in CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. All three models are slightly **fresher than WOA23** in the global mean (by −0.32, −0.07 and −0.27 psu respectively), with GFDL-ESM4 clearly the closest to observations (RMSE 0.63 psu, r = 0.924). The **historical (1950–2014) global-mean trends are all weak and statistically indistinguishable from zero** once autocorrelation is accounted for; the **ssp245 (2071–2100) trends are larger and negative (freshening)** for all three models, but only CESM2 and GFDL-ESM4 have 95 % intervals excluding zero. Because none of the historical slopes is significantly different from zero, the ssp245/historical slope ratios (7.8, 4.9, 42.0) are **not meaningful**. End-of-century ΔSSS (ssp245 2071–2100 minus historical 1995–2014) is small and negative globally (−0.11, −0.01, −0.09 psu), with basin-dependent structure (Atlantic salinification in CESM2 and GFDL-ESM4; Pacific freshening in all three).

## 2. Methodology as actually executed

**Data.** CMIP6 historical `so` was loaded with `load_regular_data(dataset='cmip6_historical', variable='so', model=…, time_range=['1950-01','2014-12'], depth_range=[0,10])`; ssp245 `so` with the analogous call over `['2071-01','2100-12']`. WOA23 salinity came from `load_climatology_data(dataset='woa23', variable='salinity', time_slice='annual', period='1971-2000', depth_range=[0,10])`, giving three levels (0, 5, 10 m). The WOA23 1° `landsea_01` mask, `basinmask_01`, WOA cell areas and WOA depth levels/boundaries were loaded (`load_landsea_mask`, `load_basin_mask`, `load_areas`, `load_depth_levels`). No explicit `load_volumes_for_cmip6` or CMIP6 `load_depth_levels` call was made; model areacello/volcello/depth information was handled through the data provenance carried by `load_regular_data` and `compute_basic_stat`.

**SSS definition.** For both models and WOA23, SSS was formed with `compute_basic_stat(op='volume_weighted_mean', dim='depth')` on the 0–10 m level subset — i.e. all selected levels weighted by their full cell volume, not by their overlap with 0–10 m. For models this was applied to every monthly time step; for WOA23 to the annual climatology. Units are psu (CMIP6 salinity stored in 1e-3, numerically psu).

**Mean state.** 1971–2000 climatologies were formed with `temporal_aggregate(freq='mean', op='mean')` on monthly fields from 1971-01 to 2000-12. Native-grid global means used `compute_basic_stat(op='area_weighted_mean', dim=['lat','lon'])`; WOA23's global mean used the same operator with WOA area weights. Differences were model minus WOA23.

**Bias/RMSE/correlation (Fig. 1).** Model climatologies were regridded to the WOA23 1° grid with `regrid_horizontal(target='woa23', method='conservative')`; a common valid-cell set was built with `build_comparison_mask(ocean_only=True)`; bias and RMSE came from `compute_difference_metrics(weighting='area')` and the spatial Pearson correlation from `analyze_correlation(weighting='area')`. The cell-by-cell difference field was retained for plotting.

**Trends.** Annual means were formed with `temporal_aggregate(freq='year', op='mean')`; global means with area weighting; trends from `fit_1D_trend(method='linear')`, whose confidence intervals account for lag-1 autocorrelation (effective sample size n_eff). Slopes and CIs were multiplied by 10 to report psu per decade; R² from `fit_1D_trend['metrics']['r2']`. Cell-by-cell trends used `fit_trend_map(dim='year', per_decade=True)` on native grids, regridded conservatively (trend) and nearest-neighbour (significance) to the WOA23 grid for display only.

**ΔSSS.** ssp245 2071–2100 mean minus historical 1995–2014 mean was computed with `temporal_aggregate(freq='mean', op='mean')`, regridded to WOA23, masked, and area-averaged globally; basin averages used `select_basin(basin='Atlantic'/'Indian'/'Pacific', include_marginal_seas=True)` followed by area-weighted means.

## 3. Key numerical results

**Global mean SSS, 1971–2000 (psu)**

| | global mean | difference vs WOA23 | bias (WOA grid) | RMSE | spatial r |
|---|---|---|---|---|---|
| WOA23 | 34.63 | — | — | — | — |
| CESM2 | 34.32 | −0.32 | −0.31 | 0.78 | 0.901 |
| GFDL-ESM4 | 34.56 | −0.07 | −0.04 | 0.63 | 0.924 |
| MPI-ESM1-2-LR | 34.36 | −0.27 | −0.27 | 0.77 | 0.891 |

**Linear trends of annual global-mean SSS (psu per decade; autocorrelation-adjusted 95 % CI)**

| | historical 1950–2014 | R² | excludes 0 | ssp245 2071–2100 | R² | excludes 0 |
|---|---|---|---|---|---|---|
| CESM2 | −0.00268 [−0.00693, 0.00157] | 0.245 | no (p=0.179) | −0.0209 [−0.0255, −0.0164] | 0.910 | yes (p=1.0e−6) |
| GFDL-ESM4 | −0.00110 [−0.00361, 0.00141] | 0.064 | no (p=0.362) | −0.00534 [−0.0106, −3.9e−5] | 0.324 | yes (p=0.049) |
| MPI-ESM1-2-LR | −0.000551 [−0.00504, 0.00393] | 0.007 | no (p=0.791) | −0.0231 [−0.0513, 0.00504] | 0.701 | no (p=0.079) |

Effective sample sizes were small (n_eff ≈ 9, 15 and 12 of 65 historical years; ≈ 12, 12 and 5 of 30 ssp245 years), which is why several unadjusted p-values were formally significant while the autocorrelation-adjusted intervals were not.

**ssp245/historical slope ratios:** CESM2 7.8, GFDL-ESM4 4.9, MPI-ESM1-2-LR 42.0 — all reported as **not meaningful**, because each historical slope's confidence interval includes zero and the two experiment segments need not use the same ensemble member.

**ΔSSS, ssp245 2071–2100 minus historical 1995–2014 (psu)**

| | global | Atlantic | Indian | Pacific |
|---|---|---|---|---|
| CESM2 | −0.11 | +0.10 | +0.06 | −0.24 |
| GFDL-ESM4 | −0.01 | +0.19 | −0.00 | −0.04 |
| MPI-ESM1-2-LR | −0.09 | −0.04 | −0.05 | −0.11 |

**Caveat.** The ΔSSS values combine the ssp245 and historical experiments, which may use different ensemble members or variant labels for a given model; they should be interpreted with that in mind. Southern and Arctic Oceans were excluded from the basin breakdown, and marginal seas were folded into their parent ocean.

## 4. Figures (what they show)

- **Figure 1** — one row of three global maps (CESM2, GFDL-ESM4, MPI-ESM1-2-LR) of the 1971–2000 SSS bias (regridded model climatology minus WOA23) on the shared WOA23 1° comparison cells. All panels use a single diverging RdBu_r scale centred on zero, clipped at the 98th percentile of |bias|, with a colour bar labelled in psu. The figure conveys where each model is saltier or fresher than the WOA23 normal.
- **Figure 2** — one row of three global maps, on a common −180…180°/ −90…90° frame, of the cell-by-cell linear trend of annual-mean SSS over 1950–2014 in psu per decade (native-grid `fit_trend_map` results regridded for display only). A shared diverging scale centred on zero is used, and cells with p < 0.05 are stippled.
- **Figure 3** — a single time-series panel of annual global-mean SSS for all three models, with the historical 1950–2014 and ssp245 2071–2100 segments drawn separately and the 2015–2070 gap left empty. Each segment is overlaid with its own `fit_1D_trend` line, the y-axis is labelled in psu, and each model keeps the colour used in Figures 1–2.

## 5. Insights and patterns

1. **Systematic fresh bias.** All three models underestimate global SSS relative to WOA23, and the native-grid and WOA-grid differences agree closely, so the bias is not a regridding artefact. GFDL-ESM4 is the best mean-state match on every metric (smallest difference, lowest RMSE, highest correlation); CESM2 and MPI-ESM1-2-LR are similarly biased (≈ −0.3 psu) and have similarly degraded pattern statistics.
2. **Pattern skill is high despite the offset.** Spatial correlations of 0.89–0.92 show that the large-scale SSS pattern is well captured; the residual RMSE (0.63–0.78 psu) is dominated by regional features (marginal seas, high-latitude and river-influenced regimes) and by the uniform fresh offset.
3. **Historical trends are not detectable.** Integrated global-mean SSS over 1950–2014 changes by only a few thousandths of a psu per decade in each model, with R² ≤ 0.25; after accounting for strong lag-1 autocorrelation (0.62–0.76) the intervals easily include zero. No model shows a robust historical global SSS trend.
4. **End-of-century freshening is stronger.** Under ssp245 the global mean freshens by −0.021, −0.005 and −0.023 psu/decade in CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, i.e. an order of magnitude larger than the historical rates (though over a shorter 30-year window). CESM2 and GFDL-ESM4 have intervals excluding zero; MPI-ESM1-2-LR's interval is wider because its effective sample size is only ~5 years, so its (largest-magnitude) freshening is not statistically resolved here. Because the historical slopes are not distinguishable from zero, the ssp245/historical ratios should not be interpreted as an acceleration factor.
5. **Basin contrast.** The global ΔSSS is small because regional changes partly cancel: the Atlantic salinifies in CESM2 (+0.10) and GFDL-ESM4 (+0.19), while the Pacific freshens in all three models (−0.24, −0.04, −0.11). MPI-ESM1-2-LR is the only model that freshens in all three basins. This Atlantic/Pacific contrast is consistent with the familiar "wet-gets-wetter / salty-gets-saltier" amplification pattern and is a more robust signature than the global-mean change.
