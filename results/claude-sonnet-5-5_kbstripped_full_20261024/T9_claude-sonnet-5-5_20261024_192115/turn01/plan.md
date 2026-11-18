# Plan

1. Load the data with `load_regular_data` for the PDO domain (20N-70N, 110E-100W), 1950-01 to 2014-12.
   - ERSSTv5: dataset `ersst`, variable `sst`.
   - CMIP6 historical `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, from `cmip6_historical`.
   - Published PDO: dataset `cpc_indices`, product `pdo`, loaded through `load_regular_data`.

2. Regrid the four SST fields to a common 1-degree grid with `regrid_horizontal`, method 'bilinear'. Use the same 1-degree grid and the same period for all four sources.
   - Note [omip-native-grid-vs-regridded]: the EOF, explained variance and regression pattern are grid-dependent, so take all numbers on the 1° grid. State that bilinear was used, as the question prescribes, rather than conservative.
   - Note [cmip6-grid-labels]: GFDL-ESM4 is already on a 1° `gr` grid, whereas CESM2 and MPI-ESM1-2-LR are on native `gn` grids, so bilinear interpolation affects them differently. The coarse MPI grid cannot reproduce narrow features in amplitude.
   - Note [missing-values-and-land-cells]: build each field's valid mask from the field itself. Use NaN-aware reductions with weights masked to the same valid cells. Any value near 1e20 is an unmasked fill value.

3. Build the anomalies with `compute_anomaly_vs_own_climatology`, using each field's own 1950-2014 monthly climatology, not WOA23. State the baseline.
   - Note [anomaly-baselines-and-climate-indices]: model bias cancels with each field's own climatology. The published `pdo` index may sit on another baseline or construction, so offsets against it are not reconstruction errors.
   - Remove the global-mean SST anomaly as well. Compute it with `compute_basic_stat` using an area-weighted mean on the 1° field, over the global SST field if it is available, otherwise state the domain used.

4. Compute the leading EOF (the PDO mode) for each source on the 1° grid, with the anomalies weighted by sqrt(cos(latitude)). Use `compute_pdo_index`, with code only for what the tool does not cover.
   - Note [eof-modes-and-pdo-index]: fix the sign so the positive phase is cool in the central and western North Pacific, for every pattern and index. Standardise each index to zero mean and unit variance. Report the percentage of variance explained, to 1 decimal place, only across fields on the same grid and period.
   - Amplitude comes from the regression pattern (anomalies regressed on the standardised index, °C per standard deviation), using `fit_multivariate_relation` or `fit_trend_map` as appropriate. Its area-weighted RMS, computed with `compute_basic_stat`, is the mode strength, reported to 3 decimal places.

5. Compute the lag-12-month autocorrelation of each monthly standardised index with `analyze_correlation` (lead-lag analysis), to 3 decimal places.

6. Correlate the ERSSTv5-derived index with the published PDO series using `analyze_correlation`, as it is and without detrending, to 3 decimal places.
   - Note [serial-correlation-and-significance]: quote the effective-sample-size-corrected significance, not the raw one. Monthly autocorrelation makes the number of independent samples far smaller than the number of months. Explain what the correlation shows (the two series are estimates of the same mode) and what it does not show (identical construction or baseline, or independence of the samples).

7. Compute the spatial Pearson correlation between each model's regression pattern and the ERSSTv5 regression pattern with `analyze_correlation` (`weighting='area'`), to 3 decimal places.
   - Note [common-comparison-mask]: do not apply a shared mask. Correlate over the cells the two patterns share, which `analyze_correlation` handles by leaving out cells either one lacks. Report the cell count for each model, and compare the models on these correlations. Do not give a p-value for the spatial correlations.
   - Name the model with the highest pattern correlation as the best at reproducing the observed pattern.

8. Find the PDO phase transitions of the ERSSTv5 index.
   - Smooth the ERSSTv5 index with a 108-month centred running mean using `temporal_aggregate`.
   - Detect the sign changes, treating them as persistence-based events, with `detect_threshold_events` at threshold 0 where useful.
   - Merge any crossing less than 3 years after the last counted transition into it, keeping the earlier year.
   - Report the transition years and say which of them match documented Pacific climate shifts (for example 1976/77, and the late-1990s shift if one is found).
   - Compare each model's amplitude (step 4) and 12-month autocorrelation (step 5) with the observed values and say whether they are consistent.

9. Figures.
   - Figure 1: a 2x2 panel of the regression patterns over the PDO domain, in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Use one diverging colour scale centred on zero, labelled °C per standard deviation. Annotate each panel with its explained variance and, for the models, its pattern correlation with ERSSTv5.
   - Figure 2: one time-series panel with the standardised indices of the four sources and the published PDO as a distinguishable fifth line. Add the 9-year smoothed ERSSTv5 curve, mark the transition years, and label the y axis as a standardized index.

Other notes:
- The baseline is each field's own 1950-2014 monthly climatology, and the bilinear regridding to 1° is the question's prescribed method.
- Transitions are defined on the 108-month running mean, with crossings less than 3 years apart merged into the earlier one.
- Explained variance is comparable only across fields on the same grid and period.
- The published PDO is a NOAA NCEI index and may differ in baseline and construction, so offsets against it are not reconstruction errors.
- The sign convention is positive phase cool in the central and western North Pacific, so a negative pattern correlation is a real disagreement.
