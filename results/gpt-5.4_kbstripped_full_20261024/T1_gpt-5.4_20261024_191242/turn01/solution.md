Using the executed workflow, I assessed ENSO amplitude from monthly Niño 3.4 SST anomalies over 1950-01 to 2014-12 for ERSSTv5 and three CMIP6 historical runs: CESM2, GFDL-ESM4, and MPI-ESM1-2-LR.

Method used
- Data loaded:
  - CMIP6 historical `tos` for CESM2, GFDL-ESM4, MPI-ESM1-2-LR
  - ERSST `sst`
- Spatial domain loaded directly: 5°S–5°N, 170°W–120°W, implemented as `lat_range=[-5,5]`, `lon_range=[-170,-120]`
- Time window: `1950-01` to `2014-12` inclusive, giving 780 monthly values for each source
- Each source was analyzed on its own native grid; no regridding was performed
- Monthly Niño 3.4 SST was computed with `compute_basic_stat(..., op="area_weighted_mean", dim=["lat","lon"])`
- Anomalies were computed with `compute_anomaly_vs_own_climatology(..., freq="month", base_period=["1950-01","2014-12"], method="subtract")`, so each source used its own monthly climatology over 1950-2014
- Amplitudes are population standard deviations over time (`ddof=0` as requested in the task statement, and consistent with the reported outputs)
- Calendar-month amplitude was computed from within-month standard deviations across years using `temporal_aggregate(..., freq="climatology", op="std")`
- Exceedances were counted for months with `|anomaly| > 1.0 °C` using `compute_exceedance_fraction(..., comparison="abs_greater", threshold=1.0, dim="time", weighting="none", reference=0.0)`

1. Direct answer: ENSO amplitude against ERSSTv5
Full-period Niño 3.4 anomaly standard deviation, 1950-01 to 2014-12:
- ERSSTv5: 0.85 °C
- CESM2: 1.29 °C
- GFDL-ESM4: 1.00 °C
- MPI-ESM1-2-LR: 0.85 °C

Model amplitude bias relative to ERSSTv5 and amplitude ratio:
- CESM2: +0.44 °C, ratio 1.52
- GFDL-ESM4: +0.15 °C, ratio 1.18
- MPI-ESM1-2-LR: +0.01 °C, ratio 1.01

Ranking from smallest to largest absolute amplitude bias:
1. MPI-ESM1-2-LR
2. GFDL-ESM4
3. CESM2

Interpretation:
- CESM2 substantially overestimates Niño 3.4 amplitude.
- GFDL-ESM4 also overestimates amplitude, but more moderately.
- MPI-ESM1-2-LR matches observed overall amplitude very closely.

2. Seasonal cycle of Niño 3.4 amplitude
Calendar-month standard deviation of Niño 3.4 anomalies (°C):

- ERSSTv5:
  - Jan 1.08, Feb 0.88, Mar 0.67, Apr 0.60, May 0.61, Jun 0.61,
    Jul 0.68, Aug 0.76, Sep 0.85, Oct 0.99, Nov 1.08, Dec 1.12

- CESM2:
  - Jan 1.62, Feb 1.43, Mar 1.20, Apr 1.07, May 1.01, Jun 0.94,
    Jul 0.98, Aug 1.09, Sep 1.24, Oct 1.40, Nov 1.56, Dec 1.62

- GFDL-ESM4:
  - Jan 1.13, Feb 1.06, Mar 0.96, Apr 0.91, May 0.84, Jun 0.78,
    Jul 0.75, Aug 0.85, Sep 0.98, Oct 1.11, Nov 1.21, Dec 1.25

- MPI-ESM1-2-LR:
  - Jan 0.96, Feb 0.91, Mar 0.78, Apr 0.82, May 0.86, Jun 0.78,
    Jul 0.71, Aug 0.74, Sep 0.77, Oct 0.86, Nov 0.97, Dec 1.03

Peak month, minimum month, and max/min ratio:
- ERSSTv5: peak Dec 1.12 °C; minimum Apr 0.60 °C; ratio 1.88
- CESM2: peak Jan 1.62 °C; minimum Jun 0.94 °C; ratio 1.72
- GFDL-ESM4: peak Dec 1.25 °C; minimum Jul 0.75 °C; ratio 1.67
- MPI-ESM1-2-LR: peak Dec 1.03 °C; minimum Jul 0.71 °C; ratio 1.45

Does each model reproduce the observed peak month and seasonal phase locking?
- Observed peak month is December.
- CESM2: No; its largest monthly standard deviation occurs in January, not December.
- GFDL-ESM4: Yes; peak in December.
- MPI-ESM1-2-LR: Yes; peak in December.

Strength of seasonal phase locking relative to observations:
- Observed max/min ratio is 1.88.
- CESM2: 1.72, somewhat weaker than observed
- GFDL-ESM4: 1.67, weaker than observed
- MPI-ESM1-2-LR: 1.45, much weaker than observed

Interpretation:
- All three models show weaker seasonal phase locking than ERSSTv5.
- GFDL-ESM4 and MPI-ESM1-2-LR capture the observed peak timing better than CESM2.
- MPI-ESM1-2-LR gets the overall amplitude right but underplays the seasonal contrast most strongly.

3. Months exceeding |Niño 3.4 anomaly| > 1.0 °C
- ERSSTv5: 170 months, 21.8%
- CESM2: 360 months, 46.2%
- GFDL-ESM4: 237 months, 30.4%
- MPI-ESM1-2-LR: 199 months, 25.5%

Interpretation:
- CESM2 produces strong ENSO-amplitude months far too often.
- GFDL-ESM4 is also too active relative to observations.
- MPI-ESM1-2-LR is closest to observed exceedance frequency, though still slightly high.

4. Overall scientific interpretation
- ERSSTv5 gives a full-period Niño 3.4 anomaly standard deviation of 0.85 °C over 1950-2014.
- Among the three models, MPI-ESM1-2-LR best matches the observed overall ENSO amplitude.
- GFDL-ESM4 is moderately too energetic.
- CESM2 is clearly too strong, both in total variance and in the frequency of >1 °C events.
- In seasonal phase locking, ERSSTv5 peaks in December and is weakest in April.
- GFDL-ESM4 and MPI-ESM1-2-LR reproduce the observed December peak month, while CESM2 peaks in January.
- However, all three models have weaker peak-to-minimum seasonal contrast than observed, especially MPI-ESM1-2-LR.

5. Visualization summary
Two figures were generated.

Figure 1
- A single-axis line plot of the seasonal cycle of Niño 3.4 anomaly standard deviation
- X-axis: calendar month from January to December
- Y-axis: standard deviation in °C
- Four curves are shown together: ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR
- ERSSTv5 is visually distinguished from the models
- Each source’s peak month is marked on its curve

Figure 2
- Four stacked time-series panels sharing the same time axis over 1950-2014
- Panel order: ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR
- Each panel shows that source’s Niño 3.4 anomaly series
- All panels use the same y-axis range spanning all four series
- Horizontal reference lines at +1.0 °C and -1.0 °C are included
- Each panel is annotated with the source’s full-period standard deviation and exceedance percentage

These figures are designed to show, respectively, the seasonal modulation of Niño 3.4 variability and the temporal evolution and event-frequency structure of the anomaly series on a directly comparable scale.
