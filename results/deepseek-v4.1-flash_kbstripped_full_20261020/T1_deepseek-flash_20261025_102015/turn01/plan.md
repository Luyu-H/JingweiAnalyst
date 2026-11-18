# Plan

1. Load ERSSTv5 monthly SST for 1950-01 to 2014-12 on its native 2° grid with `load_regular_data` (dataset `ersst`, variable `sst`).
2. Load CMIP6 historical `tos` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR for 1950-01 to 2014-12 with `load_regular_data` (dataset `cmip6_historical`), each on its own native grid and with its own published cell areas.
   - Note [missing-values-and-land-cells]: make sure the read path applies CF masking so the 1e20 land/below-sea-floor fill values never enter sums or means.
3. Build the ERSSTv5 Niño 3.4 index with `compute_enso_nino_index`: area-weighted mean SST anomaly over 5S–5N, 170W–120W, monthly series 1950-01 to 2014-12, anomaly against ERSSTv5's own 1950–2014 monthly climatology.
   - Note [anomaly-baselines-and-climate-indices]: each source's index uses its own 1950–2014 monthly climatology (each calendar month compared with its 1950–2014 mean), never WOA23 or another source's climatology; all four share the same 1950–2014 baseline, so no cross-baseline offset arises.
   - Note [omip-weighted-means-and-integrals]: the regional mean is Σ T·dA / Σ dA with the area weights restricted to the same valid cells as the data, not a plain mean over (y, x).
   - Note [monthly-means-are-time-averages]: the removed climatology is the plain mean of each calendar month over 1950–2014 (no day weights); the resulting standard deviations measure interannual/inter-month variability only, not sub-monthly variance.
4. Build the three model Niño 3.4 indices with `compute_enso_nino_index` on each model's own grid and published cell areas over the identical 5S–5N, 170W–120W box and window, each against its own 1950–2014 monthly climatology ([anomaly-baselines-and-climate-indices], [omip-weighted-means-and-integrals], [missing-values-and-land-cells] apply as in step 3).
   - Note [cmip6-grid-labels]: GFDL-ESM4 is `grid_label='gr'`, its `tos` already remapped to 1°, so compute its index on that archive grid with the regridded mask, not the native tripolar grid; resolution differs (CESM2 ~100 km, MPI ~250 km), so features narrower than the coarse grid are not comparable in amplitude.
   - Note [cmip6-experiments-historical-ssp245]: historical is the all-forcing 1850–2014 coupled run branched from piControl, with free internal variability, so its years are not synchronized with the observed years.
5. For each of the four index series, compute the population standard deviation (ddof = 0) over all 780 months with `compute_basic_stat` (operator `std`, dimension time), giving four amplitude values in °C to 2 dp.
6. From the step-5 values, compute for each model the amplitude bias (model std − ERSSTv5 std) and the amplitude ratio (model std / ERSSTv5 std), to 2 dp, and rank the three models by smallest to largest absolute bias (own arithmetic on the tool outputs; no auxiliary masking or adjustment).
7. For each source, compute the standard deviation (ddof = 0) of its Niño 3.4 anomaly separately for each of the 12 calendar months with `compute_basic_stat`, grouped by calendar month, yielding 48 values in °C to 2 dp.
8. From step 7, report for each source the calendar month of largest and of smallest standard deviation, their values, and the largest-to-smallest ratio to 2 dp; then compare each model's peak month and peak strength with ERSSTv5's to judge whether the observed seasonal phase locking is reproduced.
9. For each source, apply `detect_threshold_events` to |Niño 3.4 anomaly| with threshold 1.0 °C over 1950-01 to 2014-12; count the exceedances as the sum of the reported episode durations (strictly greater than 1.0 °C), and take the percentage as that count divided by 780 months × 100, to 1 dp.
10. Figure 1: line plot of the 12 calendar-month standard deviations (January–December x axis, °C y axis, four curves on one axis, ERSSTv5 visually distinct from the three models) with each curve's peak month marked.
11. Figure 2: four vertically stacked panels sharing one 1950–2014 time axis in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, each plotting that source's monthly Niño 3.4 anomaly on an identical y-axis range spanning all four series, with +1.0 and −1.0 °C levels marked and each panel annotated with its standard deviation and exceedance percentage.

Other notes:
- All standard deviations are population values (ddof = 0), as specified in the question.
- The model-to-observation comparison is statistical (amplitude, seasonal phase locking, exceedance frequency); because the coupled historical runs are not year-synchronized with the observed record (see [cmip6-experiments-historical-ssp245] in step 4), no event-by-event matching is implied.
- Exceedances use a strict |anomaly| > 1.0 °C criterion, evaluated on the monthly index of each source over the full 780-month window.
