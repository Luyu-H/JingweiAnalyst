**Answer in brief**
Over 1950-2014, MPI-ESM1-2-LR reproduces the ERSSTv5 Niño 3.4 amplitude almost exactly. GFDL-ESM4 is about 18% too strong, and CESM2 is about 52% too strong. All three models place the seasonal peak in boreal winter. Only MPI and GFDL peak in December like the observations. CESM2's January peak is only 0.003 °C above its December value. MPI has the weakest phase locking of the three.

**Method, as executed**
- **Data:**
  - ERSSTv5 `sst` and CMIP6 historical `tos` (CESM2, GFDL-ESM4, MPI-ESM1-2-LR) were loaded with `load_regular_data` for 1950-01 to 2014-12.
  - Each source was subset to lat −10 to 10 and lon −180 to −110, which is slightly larger than the Niño 3.4 box.
  - Each model stayed on its own native grid, and no common mask was applied.
  - All four sources were already in °C, so the Kelvin conversion helper never changed anything.
- **Index:**
  - `compute_enso_nino_index` was called with `nino_region="nino3.4"`, `anomaly_freq="month"`, `climatology="self"` and `base_period=["1950-01","2014-12"]`.
  - Each source's anomaly is therefore relative to its own 1950-2014 monthly climatology.
  - The 5S-5N, 170W-120W region masking and the area weighting are done inside the tool. The code did not pass `areacello` or weights explicitly.
  - Every index has 780 months, no NaNs and no fill-value outliers (checked by assertions).
- **Overall amplitude:** `compute_basic_stat(op="std", dim="time")` on each 780-month index. The code did not set `ddof` for this call.
- **Bias, ratio and ranking:** plain arithmetic on the unrounded standard deviations. Bias is the model minus ERSSTv5, and ratio is the model divided by ERSSTv5.
- **Seasonal cycle:** NumPy `std(ddof=0)` of the anomalies in each calendar month. Peak month, minimum month and the max/min ratio were taken from these 12 values.
- **Phase-locking judgement:** each model's max/min ratio was divided by the ERSSTv5 ratio. A value above 1.15 was labelled stronger, below 0.85 weaker, and anything between similar. These thresholds were chosen in the code and are not standard.
- **Exceedance:** a plain count of months with |anomaly| > 1.0 °C out of 780, with no event or duration logic.

**Overall amplitude, bias and ratio**

| Source | Std (°C) | Bias (°C) | Ratio |
|---|---|---|---|
| ERSSTv5 | 0.85 | – | – |
| CESM2 | 1.29 | +0.44 | 1.52 |
| GFDL-ESM4 | 1.00 | +0.15 | 1.18 |
| MPI-ESM1-2-LR | 0.85 | +0.01 | 1.01 |

Ranking from smallest to largest absolute bias: MPI-ESM1-2-LR, GFDL-ESM4, CESM2. All three models are at or above the observed amplitude. This agrees with the literature remark that GFDL is somewhat too strong, and the computed numbers support it.

**Calendar-month standard deviation (°C), Jan to Dec**
- ERSSTv5: 1.08, 0.88, 0.67, 0.60, 0.61, 0.61, 0.68, 0.76, 0.85, 0.99, 1.08, 1.12
- CESM2: 1.62, 1.43, 1.20, 1.07, 1.01, 0.94, 0.98, 1.09, 1.24, 1.40, 1.56, 1.62
- GFDL-ESM4: 1.13, 1.06, 0.96, 0.91, 0.84, 0.78, 0.75, 0.85, 0.98, 1.11, 1.21, 1.25
- MPI-ESM1-2-LR: 0.96, 0.91, 0.78, 0.82, 0.86, 0.78, 0.71, 0.74, 0.77, 0.86, 0.97, 1.03

**Seasonal summary**

| Source | Peak month (value) | Smallest-std month (value) | Max/min ratio |
|---|---|---|---|
| ERSSTv5 | Dec (1.12) | Apr (0.60) | 1.88 |
| CESM2 | Jan (1.62) | Jun (0.94) | 1.72 |
| GFDL-ESM4 | Dec (1.25) | Jul (0.75) | 1.67 |
| MPI-ESM1-2-LR | Dec (1.03) | Jul (0.71) | 1.45 |

- **CESM2:** it does not formally reproduce the December peak, since its maximum falls in January. The December value (1.62) is almost identical to January's (1.6236 vs 1.6208). The ratio of ratios is 0.92, so its phase-locking strength is judged similar to observed.
- **GFDL-ESM4:** it reproduces the December peak, with a ratio of ratios of 0.89, which is judged similar to observed.
- **MPI-ESM1-2-LR:** it reproduces the December peak, but its ratio of ratios is 0.77, so its phase locking is judged weaker than observed. Its seasonal range is the flattest of the four.
- All models have their minimum in boreal summer (June or July), while ERSSTv5's falls in April.

**Months with |Niño 3.4 anomaly| > 1.0 °C (out of 780)**

| Source | Count | Percent |
|---|---|---|
| ERSSTv5 | 170 | 21.8% |
| CESM2 | 360 | 46.2% |
| GFDL-ESM4 | 237 | 30.4% |
| MPI-ESM1-2-LR | 199 | 25.5% |

The exceedance rates follow the amplitude ordering. Even MPI, whose overall std matches the observations, exceeds the threshold somewhat more often, at 25.5% against 21.8%.

**Figures**
- **Figure 1:** a line plot of calendar-month standard deviation (Jan to Dec on x, °C on y) with four curves on one axis.
  - ERSSTv5 is drawn thick and black, and the models are in distinct colours.
  - Each curve's peak month is marked with a star and a label giving the month and value.
- **Figure 2:** four stacked panels sharing a time axis (1950-2014) and a common y-axis, ordered ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR.
  - Each panel shows that source's monthly Niño 3.4 anomaly series.
  - The y-limits are symmetric and set from the largest absolute anomaly across all four series, plus a margin, so every series fits.
  - Dashed lines mark +1.0 and −1.0 °C, and a thin line marks zero.
  - Each panel has a text box with its standard deviation and the percentage of months with |anomaly| > 1.0 °C.

**Insights and caveats**
- Because each source is anomalised against its own climatology, mean-state SST biases are excluded. The differences reflect interannual variability only.
- MPI's close match in overall std coexists with a flatter seasonal cycle. Matching the overall amplitude therefore does not mean matching the seasonal phase locking.
- CESM2 is too strong in every calendar month, with a minimum (0.94) that already exceeds the observed December peak's neighbours.
- Differences in grid resolution (ERSST 2°, MPI coarse and curvilinear, GFDL at 1°) and in `tos` top-cell definition may contribute to the amplitude differences. The analysis did not test this.
- Each source is a single realisation of 65 years. Sampling variability of ENSO could account for part of the differences, especially in the single calendar-month values. No uncertainty estimate was computed.
