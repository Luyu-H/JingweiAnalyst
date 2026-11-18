# Plan

1. Load the data with `load_regular_data`: CMIP6 historical `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly 1950-01 to 2014-12; ERSSTv5 `sst`, monthly 1950-01 to 2014-12; and the published PDO series from `cpc_indices` using product `pdo` over the same window.
- Note [cmip6-tos-bulk-top-cell]: `tos` is the bulk temperature of the top model cell, including sea-ice-covered cells, which remain valid near-freezing water rather than missing data.
- Note [cmip6-experiments-historical-ssp245]: the historical runs are coupled and contain internal variability, and with one member per model the forced change cannot be separated from internal variability; compare statistics and patterns, not phase timing.

2. Restrict all SST fields to the PDO domain 20N-70N, 110E-100W, handling the dateline, and put the model `tos` fields and the ERSSTv5 `sst` field on the 1-degree grid with `regrid_horizontal` using bilinear interpolation.
- Note [omip-native-grid-vs-regridded]: bilinear interpolation to 1 degree is used as the question specifies, although conservative regridding would be preferred going to a coarser grid and before an area-weighted step, so state that results may differ from a conservative-regridded analysis.
- Note [cmip6-grid-labels]: GFDL-ESM4 `tos` is already on a 1-degree `gr` grid regridded from its tripolar grid with `conserve_order1`; bilinear interpolation to 1 degree does not restore native resolution, and its land-sea mask is the regridded one.

3. For each of the four sources (ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR), compute the PDO with `compute_pdo_index` over 20N-70N, 110E-100W for 1950-2014, using anomalies relative to that source's own 1950-2014 monthly climatology, and obtain the standardized index, the regression pattern in degrees Celsius per standard deviation, and the leading-mode explained variance.
- Note [eof-modes-and-pdo-index]: remove each source's monthly climatology and global-mean SST anomaly, area-weight anomalies by `sqrt(cos(latitude))`, take the leading EOF, fix the sign so positive index means negative SST anomalies in the central/western North Pacific, standardize the index, and report amplitude as the area-weighted RMS of the regression pattern; compare patterns with `analyze_correlation(..., weighting='area')`, and treat explained variance as comparable only on the same grid and period.
- Note [anomaly-baselines-and-climate-indices]: remove the stated 1950-2014 monthly climatology and use each model's own climatology so mean-state bias cancels; state the climatology and baseline, and note that a different published-PDO `base_period` produces a baseline offset, not an error.

4. Use `compute_basic_stat` with area weighting on each regression pattern to obtain its area-weighted root-mean-square amplitude in degrees Celsius per standard deviation; report the explained variance from step 3 to 1 decimal place and the amplitude to 3 decimal places.

5. Standardize each monthly index and compute the lag-12-month autocorrelation directly from the standardized monthly series, reporting it to 3 decimal places; no listed tool computes autocorrelation.

6. Use `analyze_correlation` to correlate the ERSSTv5-derived standardized index with the published PDO series as-is, without detrending, and report the Pearson coefficient to 3 decimal places with effective-sample-size-corrected significance.
- Note [serial-correlation-and-significance]: correlate the two estimates of the same quantity without detrending because their shared trend is part of the agreement, and quote effective-sample-size-corrected significance rather than naive monthly significance.

7. For each model, use `analyze_correlation` with area weighting to correlate its regression pattern with the ERSSTv5 regression pattern over the cells the two products share, reporting the spatial Pearson coefficient to 3 decimal places.
- Note [common-comparison-mask]: do not apply a shared comparison mask to the regression patterns; each mode uses its own source's cells, and `analyze_correlation` omits missing cells, so correlate each model pattern with ERSSTv5 only over their shared cells.
- Use [serial-correlation-and-significance] for spatial patterns by reading the coefficient rather than its p-value, since no spatial effective-sample correction is available.

8. Use `temporal_aggregate` to form the 9-year (108-month) centred running mean of the ERSSTv5 standardized index, then identify the years in which that smoothed series changes sign; apply the stated merge rule that a crossing less than 3 years after the last counted transition is merged into it, keeping the earlier year. Compare the resulting observed PDO transition years over 1950-2014 with documented Pacific climate shifts such as the mid-1970s, late-1980s and late-1990s shifts.

9. Build Figure 1 as a 2-by-2 panel of the regression patterns in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, using one diverging colour scale centred on zero labelled in degrees Celsius per standard deviation, and annotate each panel with its explained variance and, for models, its pattern correlation with observations.

10. Build Figure 2 as a single time-series panel of the standardized indices of all four sources plus the published PDO series as a distinguishable fifth line, overlay the 9-year-smoothed ERSSTv5 curve, mark the transition years, and label the y axis as a standardized index.

11. Identify the model with the highest spatial pattern correlation against ERSSTv5 as the best pattern match, and compare each model's area-weighted RMS amplitude and lag-12-month persistence with the observed values.

Other notes:
- Figure and table results should state that the published PDO correlation is not an independent validation; it measures agreement between two estimates of the same quantity and includes their shared low-frequency signal and trend.
- Model-versus-observed differences in amplitude and persistence should be interpreted with the caveats that one member per model cannot separate forced change from internal variability and that model resolutions differ substantially.
- The phase-transition comparison should use the 1950-2014 ERSSTv5 index crossings, not model phase timing, because the coupled historical runs are not expected to reproduce observed event timing.
