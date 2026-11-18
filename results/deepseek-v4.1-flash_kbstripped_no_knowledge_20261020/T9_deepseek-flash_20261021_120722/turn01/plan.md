# Plan

1. Load monthly CMIP6 historical `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over 1950-01 to 2014-12, and monthly ERSSTv5 `sst` over the same window, using `load_regular_data`; select the North Pacific PDO domain 20N-70N, 110E-100W, with a slightly larger loading box if needed for later bilinear regridding.
2. Load the published PDO monthly series from `cpc_indices` using `load_regular_data` with product `pdo` and variable `pdo`, restricted to 1950-01 to 2014-12.
3. Bilinearly regrid each CMIP6 `tos` field onto the 1-degree grid using `regrid_horizontal`, so all four SST sources are compared cell by cell over the PDO domain.
4. For each source, compute monthly anomalies relative to its own 1950-2014 monthly climatology using `compute_anomaly_vs_own_climatology`.
5. For each source, derive the leading North Pacific SST mode and standardized PDO index over 20N-70N, 110E-100W using `compute_pdo_index`; where the tool returns only the index, use the same anomaly fields and standardized index to obtain the leading-mode explained variance, amplitude in degrees Celsius per standard deviation, and the regression pattern by custom EOF/regression calculation, without rounding intermediate results.
6. Compute the lag-1-year (12-month) autocorrelation of each monthly standardized index using `analyze_correlation` with a 12-month lag; compare CESM2, GFDL-ESM4 and MPI-ESM1-2-LR with ERSSTv5.
7. Compute the Pearson correlation between the ERSSTv5-derived standardized PDO index and the published `pdo` series using `analyze_correlation`, and report what the correlation demonstrates and does not demonstrate.
8. Compute the spatial Pearson correlation between each model’s PDO regression pattern and the ERSSTv5 regression pattern using `analyze_correlation`; use the highest pattern correlation to identify the model reproducing the observed pattern best, while comparing amplitudes and lag-12 autocorrelations separately for consistency.
9. Compute the 108-month centred running mean of the ERSSTv5 standardized PDO index using `temporal_aggregate` with a running window, then identify sign changes in the smoothed ERSSTv5 index over 1950-2014; merge any crossing less than 3 years after the last counted transition into that transition, keeping the earlier year. Match the resulting phase-transition years to documented Pacific climate shifts, including 1976/77, 1988/89 and 1998/99.
10. Create Figure 1: a 2-by-2 panel of regression patterns over the PDO domain in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, using one diverging colour scale centred on zero and labelled in degrees Celsius per standard deviation; annotate each panel with explained variance and, for the models, pattern correlation with ERSSTv5.
11. Create Figure 2: a single time-series panel of the standardized PDO index for all four sources plus the published PDO series as a distinguishable fifth line, with the 9-year smoothed ERSSTv5 curve and observed transition years marked, and the y axis labelled as a standardized index.

Other notes:
- The Pearson correlation with the published PDO series measures temporal agreement with an observation-derived index and does not establish independence, causality or a physical mechanism.
- Pattern correlations assess spatial similarity of regression patterns only, not amplitude, timing or persistence.
- Amplitude and lag-12 autocorrelation are separate consistency checks from spatial pattern correlation.
- The PDO domain is 20N-70N, 110E-100W, and anomalies are relative to each source’s 1950-2014 monthly climatology.
- Transition years use the 108-month centred running mean of the ERSSTv5 index and merge crossings less than 3 years after the previous counted transition, keeping the earlier year.
- Documented Pacific climate shifts relevant to 1950-2014 include 1976/77, 1988/89 and 1998/99.
