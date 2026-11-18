# Plan

1. Load monthly SST for 1950-01 to 2014-12 from `ersst` (`sst`) and from `cmip6_historical` (`tos`) for models CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, and load the published monthly PDO series from `cpc_indices` (product `pdo`) with `load_regular_data`.
- Note [anomaly-baselines-and-climate-indices]: build each source’s monthly anomalies against its own 1950–2014 monthly climatology, and compare the ERSSTv5-derived index with the published PDO on the same stated baseline.

2. Regrid ERSSTv5 SST and each model SST by bilinear interpolation onto the common 1-degree grid with `regrid_horizontal`, then subset the PDO domain 20N–70N, 110E–100W for the EOF analysis while retaining the full field needed to remove the global-mean SST anomaly.
- Note [omip-native-grid-vs-regridded]: do the regridding before computing PDO metrics and patterns, use bilinear interpolation for all sources, and state that explained variance is comparable only on this common grid and period.

3. For ERSSTv5 and each model, compute the 1950–2014 monthly climatology and subtract it with `compute_anomaly_vs_own_climatology`; compute the monthly global-mean SST anomaly from the regridded full field with `compute_basic_stat` using area weighting; subtract that global-mean anomaly from each source’s monthly SST anomaly field; then compute the PDO with `compute_pdo_index` over the North Pacific domain, obtaining the leading-mode standardized monthly index, explained variance, and regression pattern amplitude field.
- Note [eof-modes-and-pdo-index]: form PDO from monthly SST anomalies after removing both the monthly climatology and the global-mean SST anomaly, use area weighting with sqrt(cos latitude) in the EOF, fix the sign so positive PDO matches cool central/western North Pacific, and take the mode amplitude from the regression pattern in °C per standard deviation.

4. Build one common comparison mask over the ERSSTv5 regression pattern and the three model regression patterns with `build_comparison_mask`, then apply it to all four fields and compute area-weighted spatial Pearson correlations between each model regression pattern and the ERSSTv5 regression pattern with `analyze_correlation`.
- Note [common-comparison-mask]: compare cell by cell only on the cells shared by all fields on the common grid, using one common mask for every model-vs-observation pattern score.

5. From each standardized monthly PDO index, compute the 12-month lag autocorrelation with `analyze_correlation` by correlating each series with itself at a 12-month lag; for ERSSTv5 and the published PDO series, compute their Pearson correlation over 1950-01 to 2014-12 with `analyze_correlation`.
- Note [serial-correlation-and-significance]: for ERSSTv5-derived vs published PDO, correlate the undetrended monthly series and, if significance is mentioned, base it on effective sample size rather than the raw monthly count; for the 12-month autocorrelation, avoid over-interpreting nominal significance.

6. Use `temporal_aggregate` to compute the centred 108-month running mean of the ERSSTv5 standardized PDO index with month-length weighting, then identify the sign-change years in that smoothed series; merge any crossing occurring less than 3 years after the previous counted transition by keeping the earlier year. Report the resulting observed phase-transition years and note which correspond to documented Pacific climate shifts.
- Note [monthly-means-are-time-averages]: use `temporal_aggregate` for the 108-month running mean because the source is monthly means with varying month lengths, while keeping the monthly climatology baseline itself as a plain month-by-month mean.

7. Summarize and report, for ERSSTv5 and each model, the leading-mode explained variance (percent, 1 decimal), mode amplitude from the regression pattern (°C per standard deviation, 3 decimals), and 12-month autocorrelation (3 decimals); report the ERSSTv5-vs-published PDO correlation (3 decimals) and explain what it does and does not demonstrate; identify which model has the highest pattern correlation with ERSSTv5 and assess whether each model’s amplitude and persistence are consistent with observations.

8. Prepare Figure 1 from the ERSSTv5 and model regression patterns on one diverging colour scale centered on zero, ordered ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, annotating each panel with explained variance and each model panel also with pattern correlation versus ERSSTv5. Prepare Figure 2 as a single panel of the standardized monthly PDO indices for ERSSTv5, the three models, and the published PDO series, adding the 9-year smoothed ERSSTv5 curve and markers for the transition years.

Other notes:
- [anomaly-baselines-and-climate-indices] Any offset between the ERSSTv5-derived and published PDO series can arise from baseline choice and should not be treated as reconstruction error.
