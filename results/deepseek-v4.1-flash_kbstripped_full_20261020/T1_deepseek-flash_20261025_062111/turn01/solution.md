The code successfully computed and visualized the ENSO amplitude assessment for ERSSTv5 and three CMIP6 historical runs over 1950-01 to 2014-12.

**Methodology used by the executed code**
- Time range: 1950-01 to 2014-12, giving 780 monthly values per source.
- Region: Nino 3.4, 5°S–5°N, 170°W–120°W. Fields were loaded on their native grids over a containing latitude–longitude box `[-10, 10]` × `[-180, -110]`; no horizontal regridding or common spatial mask was applied.
- Anomalies: each source’s Nino 3.4 index was computed with `compute_enso_nino_index(..., nino_region="nino3.4", depth_level=None, anomaly_freq="month", climatology="self", base_period=["1950-01","2014-12"])`, i.e. relative to that source’s own 1950–2014 monthly climatology. No WOA23 climatology, no `ssp245`, and no `cpc_indices` were used.
- Cell areas: `load_areas` was called for each source. For ERSSTv5 the area source is “computed from ersst cell centres (the product publishes no cell bounds)”; for the CMIP6 models it is the published `areacello` from `cmip6_historical`. The recorded valid area within 5°S–5°N, 170°W–120°W is approximately: ERSSTv5 6.39×10¹² m², CESM2 6.35×10¹² m², GFDL-ESM4 6.19×10¹² m², MPI-ESM1-2-LR 6.21×10¹² m². The index call itself did not include a separate explicit `weights` argument in the executed code.
- Statistics: population standard deviation with `ddof=0`. Values are rounded only for reporting; ratios and rankings used unrounded intermediate values.
- No lagged analysis was performed.

**Overall Nino 3.4 amplitude and model bias/ratio**

| Source | Nino 3.4 std (°C) | Bias vs ERSSTv5 (°C) | Ratio vs ERSSTv5 |
|---|---:|---:|---:|
| ERSSTv5 | 0.85 | — | 1.00 |
| CESM2 | 1.29 | +0.44 | 1.52 |
| GFDL-ESM4 | 1.00 | +0.15 | 1.18 |
| MPI-ESM1-2-LR | 0.85 | +0.01 | 1.01 |

The MPI-ESM1-2-LR bias is +0.00514°C before rounding, hence reported as +0.01°C. Ranking of models from smallest to largest absolute amplitude bias: **MPI-ESM1-2-LR < GFDL-ESM4 < CESM2**. Thus MPI-ESM1-2-LR has the closest total ENSO amplitude to ERSSTv5, GFDL-ESM4 is moderately over-amplified, and CESM2 is strongly over-amplified.

**Calendar-month standard deviations, Jan–Dec (°C)**

ERSSTv5:  
Jan 1.08, Feb 0.88, Mar 0.67, Apr 0.60, May 0.61, Jun 0.61, Jul 0.68, Aug 0.76, Sep 0.85, Oct 0.99, Nov 1.08, Dec 1.12

CESM2:  
Jan 1.62, Feb 1.43, Mar 1.20, Apr 1.07, May 1.01, Jun 0.94, Jul 0.98, Aug 1.09, Sep 1.24, Oct 1.40, Nov 1.56, Dec 1.62

GFDL-ESM4:  
Jan 1.13, Feb 1.06, Mar 0.96, Apr 0.91, May 0.84, Jun 0.78, Jul 0.75, Aug 0.85, Sep 0.98, Oct 1.11, Nov 1.21, Dec 1.25

MPI-ESM1-2-LR:  
Jan 0.96, Feb 0.91, Mar 0.78, Apr 0.82, May 0.86, Jun 0.78, Jul 0.71, Aug 0.74, Sep 0.77, Oct 0.86, Nov 0.97, Dec 1.03

**Peak/minimum seasonal amplitude and phase locking**

| Source | Peak month | Peak std (°C) | Minimum month | Minimum std (°C) | Peak/min ratio |
|---|---:|---:|---:|---:|---:|
| ERSSTv5 | Dec | 1.12 | Apr | 0.60 | 1.88 |
| CESM2 | Jan | 1.62 | Jun | 0.94 | 1.72 |
| GFDL-ESM4 | Dec | 1.25 | Jul | 0.75 | 1.67 |
| MPI-ESM1-2-LR | Dec | 1.03 | Jul | 0.71 | 1.45 |

The observed ERSSTv5 peak phase-locking month is December, with a maximum-to-minimum ratio of 1.88. GFDL-ESM4 and MPI-ESM1-2-LR also peak in December, but both have weaker phase-locking ratios. CESM2 peaks in January, not December, and although its peak amplitude is large (1.62°C), its peak/min ratio (1.72) is weaker than observed.

**Months with |Nino 3.4 anomaly| strictly greater than 1.0°C**

| Source | Count / 780 | Percentage |
|---|---:|---:|
| ERSSTv5 | 170 | 21.8% |
| CESM2 | 360 | 46.2% |
| GFDL-ESM4 | 237 | 30.4% |
| MPI-ESM1-2-LR | 199 | 25.5% |

All three models exceed the ±1.0°C threshold more often than ERSSTv5. The ordering follows the total amplitude: CESM2 has the most frequent large anomalies, GFDL-ESM4 is intermediate, and MPI-ESM1-2-LR is closest to the observed exceedance rate.

**Visualizations**
- **Figure 1** is a line plot of the 12 calendar-month Nino 3.4 anomaly standard deviations from January to December. It puts all four sources on one axis, with ERSSTv5 visually distinguished from the three models, and marks each curve’s peak month with a star. The y-axis is standard deviation in °C.
- **Figure 2** consists of four vertically stacked panels in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. The panels share one 1950–2014 time axis and an identical y-axis range covering all anomaly series. Each panel shows that source’s monthly Nino 3.4 anomaly series, marks 0 and ±1.0°C levels, shades anomalies beyond ±1.0°C, and annotates that source’s all-month standard deviation and ±1.0°C exceedance count/percentage.

**Overall interpretation**
Total ENSO amplitude is best matched by MPI-ESM1-2-LR, which has essentially the observed standard deviation (ratio 1.01). GFDL-ESM4 is moderately stronger (ratio 1.18), while CESM2 is substantially over-amplified (ratio 1.52). The exceedance statistics show the same ranking. For seasonal phase locking, GFDL-ESM4 and MPI-ESM1-2-LR reproduce the December peak month but with weaker maximum-to-minimum ratios than ERSSTv5; CESM2 peaks in January and does not reproduce the observed peak month.
