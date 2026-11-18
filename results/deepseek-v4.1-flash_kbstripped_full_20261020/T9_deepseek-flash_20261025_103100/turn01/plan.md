# Plan

1. Load the four SST/tos sources over 1950-01 to 2014-12 and the PDO domain 20°N-70°N, 110°E-100°W (use 110°E-260°E to avoid a dateline split): `cmip6_historical` `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR; `ersst` `sst`; and `cpc_indices` product `pdo`. Use `load_regular_data` for each.
   - Note [cmip6-tos-bulk-top-cell]: `tos` is the bulk temperature of the top model cell, whose thickness differs between models (0-10 m CESM2, 0-2 m GFDL-ESM4, 0-12 m MPI-ESM1-2-LR); sea-ice-covered cells carry under-ice water near freezing and are not missing, which limits exact amplitude comparability.

2. Put every source, including ERSSTv5, on the common 1-degree grid with `regrid_horizontal` using bilinear interpolation, restricted to the PDO domain, before computing anomalies or EOFs.
   - Note [omip-native-grid-vs-regridded]: the question prescribes bilinear regridding onto the 1-degree grid, whereas the fragment prefers conservative regridding before area-weighted steps; follow the question, state bilinear as the method used, and state that regridded explained variance and pattern correlations depend on this.
   - Note [cmip6-grid-labels]: GFDL-ESM4 `thetao`/`tos` is on the `gr` 180x360 1-degree archive grid, so its "native-grid" numbers are archive-grid numbers; CESM2 (~100 km) and MPI (~250 km) differ in resolution, so features narrower than the coarse grid are not comparable in amplitude.

3. For each of the four regridded fields, compute monthly anomalies relative to its own 1950-2014 monthly climatology with `compute_anomaly_vs_own_climatology`, using `base_period` 1950-2014 and removing each calendar month's mean. Use the published PDO series as provided, restricted to 1950-2014.
   - Note [anomaly-baselines-and-climate-indices]: build each model and ERSSTv5 index against its own 1950-2014 monthly climatology and state the climatology and baseline; the published PDO has its own baseline and a constant per-calendar-month offset is not an error, so only variability/phase is comparable.

4. For ERSSTv5 and each model, compute the leading EOF/PDO mode over 20°N-70°N, 110°E-100°W from the monthly anomalies with `compute_pdo_index`, returning the leading-mode explained variance, the standardized monthly index, and the regression pattern.
   - Note [eof-modes-and-pdo-index]: fix the arbitrary EOF sign by the convention that a positive index corresponds to negative SST anomalies in the central/western North Pacific; a negative pattern correlation is then a real disagreement, not a sign flip. The EOF is dimensionless: take amplitude from the regression pattern, and compare explained variance only on the same grid and period.

5. Compute each mode's amplitude as the area-weighted RMS of its regression pattern over the PDO domain with `compute_basic_stat` using the weighted RMS operator, giving degrees Celsius per standard deviation [eof-modes-and-pdo-index]. Report to 3 decimal places.

6. Compute the lag-1-year (12-month) autocorrelation of each standardized monthly index with `analyze_correlation`, correlating each index with itself at a 12-month lag over the common valid period. Report to 3 decimal places.

7. Compute the Pearson correlation between the ERSSTv5-derived standardized index and the published PDO series over 1950-2014 with `analyze_correlation`.
   - Note [serial-correlation-and-significance]: comparing two estimates of the same quantity needs no detrending; the shared trend is part of the agreement being checked, not a spurious link. The correlation demonstrates agreement in temporal variability/phase, not absolute amplitude agreement, independent validation, or causality.

8. For each model, compute the spatial Pearson correlation between its regression pattern and the ERSSTv5 regression pattern with `analyze_correlation`. Build the pairwise shared-cell mask with `build_comparison_mask` for each model-versus-ERSSTv5 pair, and report the shared cell count for each correlation.
   - Note [common-comparison-mask]: do not mask the variability-mode regression pattern to all sources' coverage; correlate each model pattern with the reference pattern over only the cells the two share and report the shared cell counts. Read the pattern-correlation coefficient, not its p-value, because no effective-sample-size correction exists across space [serial-correlation-and-significance].

9. Compute the 9-year (108-month) centred running mean of the ERSSTv5 standardized index with `temporal_aggregate` using a running window. Identify the years in which this smoothed series changes sign, applying the stated merge rule: a crossing less than 3 years after the last counted transition is merged into it, keeping the earlier year. Report the observed PDO phase-transition years over 1950-2014 and compare them with documented Pacific climate shifts such as 1976/77, 1988/89 and 1998/99.

10. Compare the four sources: rank the models by spatial pattern correlation with ERSSTv5, and assess whether each model's amplitude (area-weighted RMS of the regression pattern) and persistence (12-month autocorrelation) are consistent with the ERSSTv5 values. State which model reproduces the observed pattern best and whether amplitude and persistence are consistent.

11. Produce Figure 1: a 2-by-2 panel of the regression patterns over the PDO domain in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, on one diverging colour scale centred on zero and labelled in degrees Celsius per standard deviation; annotate each panel with its explained variance and, for models, its pattern correlation with observations.

12. Produce Figure 2: a single time-series panel showing the standardized index of all four sources plus the published PDO series as a distinguishable fifth line, with the 9-year smoothed ERSSTv5 curve and the detected transition years marked, and the y axis labelled as a standardized index.

Other notes:
- The 108-month centred running mean requires a full window, so sign-change detection is limited near the 1950 and 2014 ends; state this with the transition years.
- The published PDO is a dimensionless, independently baselined index, so the ERSSTv5-published correlation checks phase/variability agreement rather than absolute amplitude.
- Report the explained variance as a percentage to 1 decimal place and amplitudes/autocorrelations/correlations to 3 decimal places, using unrounded intermediates.
