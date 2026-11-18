# Plan

1. Load the data with load_regular_data for the window 1950-01 to 2014-12.
   - cmip6_historical: variable `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. Each model keeps its own native grid and its own published areacello. Do not apply a common mask to the four fields.
   - ersst: variable `sst`, native 2-degree grid, same window.
   - Note [cmip6-tos-bulk-top-cell]: use `tos` directly and do not substitute `thetao`. The differences in top-cell thickness between models are model definitions, not errors.
   - Note [common-comparison-mask]: the index is computed from each source's own cells, so no shared mask is applied before the index.
   - Note [cmip6-grid-labels]: use each model's own grid and its own areacello (GFDL is on `gr`, CESM2 and MPI on `gn`).

2. Compute the monthly Nino 3.4 anomaly for each of the four sources (box 5S-5N, 170W-120W) over the 780 months, using compute_enso_nino_index.
   - For the models, use the model's own areacello and mask the weights to valid cells. For ERSSTv5, use its own grid.
   - Note [anomaly-baselines-and-climatology]: this is not a fragment id, so see the next line.
   - Note [anomaly-baselines-and-climate-indices]: take each source's anomaly from its own 1950-2014 monthly climatology (calendar-month matching removes the seasonal cycle). Do not use `climatology='woa23'` on models. If the index tool is not used, use compute_anomaly_vs_own_climatology.
   - Note [omip-weighted-means-and-integrals]: the box mean is area-weighted over valid cells only (compute_basic_stat with op='area_weighted_mean'), with weights masked to the same cells as the data.
   - Note [missing-values-and-land-cells]: use NaN-aware reductions and temporal_aggregate for any weighted time aggregation. A value near 1e20 in an index is an unmasked fill value, not a result.

3. Overall amplitude. For each source, compute the population standard deviation (ddof=0) of the Nino 3.4 anomaly over all 780 months, using compute_basic_stat. Report to 2 decimals.

4. Amplitude bias and ratio for each model.
   - Bias = model std minus ERSSTv5 std. Ratio = model std divided by ERSSTv5 std. Report both to 2 decimals.
   - Rank the three models from smallest to largest absolute bias, using unrounded values.
   - Use compute_difference_metrics only if it helps. Otherwise do the arithmetic on the four std values.

5. Seasonal cycle of amplitude.
   - For each source, compute the std (ddof=0) of the anomaly for each calendar month, Jan-Dec. This gives 48 values, reported to 2 decimals.
   - For each source, identify the peak month and its value, the month of smallest std, and the largest-to-smallest ratio (2 decimals).
   - State whether each model reproduces the observed peak month and the strength of phase locking, judged from the ratio compared with ERSSTv5.

6. Exceedance counts.
   - For each source, count the months with anomaly strictly > +1.0 or < -1.0 °C in absolute value (|anomaly| > 1.0), out of 780. Report the number and the percentage to 1 decimal.
   - Note [omip-weighted-means-and-integrals]: this is a plain count of months, not an area share. Do not use compute_exceedance_fraction.
   - Count individual months, not persistent events, so do not use detect_threshold_events with min_duration.

7. Figures.
   - Figure 1: line plot of calendar-month std (Jan-Dec on x, std in °C on y), four curves on one axis. ERSSTv5 is visually distinct (for example black and thick), and each curve's peak month is marked.
   - Figure 2: four stacked panels with a shared time axis over 1950-2014, ordered ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR.
     - Each panel shows that source's Nino 3.4 anomaly series.
     - All panels use an identical y-range that contains every series.
     - Dashed lines mark +1.0 and -1.0 °C.
     - Each panel is annotated with its std and exceedance percentage.

Other notes:
- All anomalies are relative to each source's own 1950-2014 monthly climatology, so model mean-state biases are excluded from the index. Std values are population std (ddof=0) over 780 months.
- The Nino 3.4 box is 5S-5N, 170W-120W.
- [common-comparison-mask] No shared mask is used, because the index is computed per source on its own cells.
- [anomaly-baselines-and-climate-indices] Exceedances are counted as individual months with strict >1.0 °C, as the question defines them.
- Resolution differs between sources: ERSST is 2°, MPI is about 250 km, and GFDL is regridded to 1°. The top-cell thickness of `tos` also differs between models. These may affect the equatorial detail and amplitude.
- Context only: the literature suggests GFDL ENSO is somewhat stronger than observed. Check this against the computed numbers before stating it, and do not use it as a target.
