## 1. Answer in brief

Over 1950–2014 on the 1° bilinear grid, all three CMIP6 models produce a leading North Pacific SST mode that resembles the ERSSTv5 PDO. **CESM2 reproduces the observed spatial pattern best (r = 0.925)**, followed by GFDL-ESM4 (0.891) and MPI-ESM1-2-LR (0.750). Amplitudes are within about ±25% of observed for all three models. The 12-month autocorrelation is within ±0.1 of observed for GFDL-ESM4 and MPI-ESM1-2-LR, but CESM2 is 0.118 lower. The ±25% and ±0.1 tolerances were my own choices in the code, not values given in the question, and were not tested against sampling uncertainty.

## 2. What was actually done

- **Data.**
  - `load_regular_data` was called with only a time range of 1950-01 to 2014-12, so global fields were loaded. The PDO box was applied later, inside `compute_pdo_index` and in the masking of the regression patterns.
  - Sources were ERSSTv5 `sst` (native 89×180, about 2°), CMIP6 historical `tos` for the three models on their native grids (CESM2 nlat/nlon, GFDL-ESM4 lat/lon, MPI-ESM1-2-LR j/i), and the published `pdo` series from `cpc_indices`.
  - Values with |x| ≥ 1e10 (unmasked fill values) were set to NaN. Units were already °C, so no conversion was needed.
- **Regridding.** Every global field was regridded with `regrid_horizontal(target="1deg", method="bilinear")`. ERSSTv5 is therefore interpolated up from about 2° and carries no finer detail than its native grid.
- **Anomalies.** `compute_anomaly_vs_own_climatology(freq="month", base_period=["1950-01","2014-12"])` was applied to each field, so the baseline is each field's own 1950–2014 monthly climatology. No WOA23 climatology was used. The `source_dataset: woa23` label in the output metadata is an attribute inherited by the tool outputs and does not change the baseline the code specified.
- **Global-mean anomaly.**
  - The area-weighted global-mean SST anomaly on the 1° field was computed with `compute_basic_stat` and stored for reporting.
  - The removal itself is done inside `compute_pdo_index`, which received the global anomaly field. The code comment says the tool removes it. I did not verify this independently.
- **EOF and index.**
  - `compute_pdo_index(an, standardize=True)` returned the standardized index, EOF, explained variance and regression pattern.
  - The sqrt(cos lat) weighting, the 20–70°N / 110°E–100°W box and the sign convention are delegated to that tool. The script applied no sign flip of its own.
  - All patterns have negative domain means, and all pattern correlations with ERSSTv5 are positive. This is consistent with a common sign convention but is not a separate check of the "cool central/western North Pacific in the positive phase" requirement.
- **Amplitude.**
  - It was taken from the regression pattern returned by `compute_pdo_index` (SST anomaly regressed on the standardized index, °C per SD). `fit_multivariate_relation` and `fit_trend_map` were not used.
  - The pattern was masked to the PDO box. The amplitude is sqrt of the area-weighted mean of the squared pattern, using `compute_basic_stat(op="area_weighted_mean")`.
- **Autocorrelation.** `analyze_correlation(idx, idx, dim="time", lag=12)` gives the single 12-month-lag value. No other lags were examined.
- **Published PDO.** `analyze_correlation(..., detrend=False)` was applied to the ERSSTv5 index and the published series, after aligning both on month-start timestamps.
- **Pattern correlation.** `analyze_correlation(..., dim=["lat","lon"], weighting="area")` was applied between each model pattern and the ERSSTv5 pattern. No common mask was imposed across models, and no p-values were computed for these correlations.
- **Transitions.**
  - The ERSSTv5 index was smoothed with `temporal_aggregate(freq="rolling", window=108, center=True)`, which leaves 673 valid months (the first and last 54 are NaN).
  - `detect_threshold_events` at 0 (above and below, `min_duration=1`) gave the runs. The first run is dropped because it begins at the edge of the valid record.
  - Crossings less than 3 years after the last counted transition are merged into it, keeping the earlier one.
  - Each remaining year is compared with ±1 year of 1977 and ±2 years of 1998.

## 3. Results (1950–2014, 1° grid)

| Source | Explained variance (%) | Amplitude (°C/SD) | Lag-12 autocorr. | Pattern r vs ERSSTv5 | Cells |
|---|---|---|---|---|---|
| ERSSTv5 | 23.0 | 0.330 | 0.342 | – | – |
| CESM2 | 26.1 | 0.364 | 0.223 | 0.925 | 4445 |
| GFDL-ESM4 | 14.2 | 0.291 | 0.305 | 0.891 | 4624 |
| MPI-ESM1-2-LR | 22.9 | 0.409 | 0.425 | 0.750 | 4368 |

The explained variances are comparable because all four fields share the same grid, period and procedure. The cell counts are the cells the two patterns share, with cells missing in either excluded.

**ERSSTv5 index vs published NCEI PDO.** r = 0.962 (n = 780 months, effective n ≈ 110.6, no detrending). The corrected and uncorrected p-values both print as 0.0, which means below numerical precision rather than exactly zero. The effective sample size is about one-seventh of the nominal one.
- **What it shows:** the ERSSTv5-derived index and the NCEI series are two estimates of essentially the same North Pacific mode.
- **What it does not show:**
  - It does not show identical construction or baseline. The published series has mean 0.003 and SD 1.02 over this window, so it is not on exactly the same footing as my standardized index.
  - It does not show that the 780 months are independent samples.
  - It does not validate the models.

**Pattern fidelity.** CESM2 is highest (0.925), then GFDL-ESM4 (0.891), then MPI-ESM1-2-LR (0.750). The MPI pattern also has the strongest extreme regression values (minimum about −2.35 °C/SD against −0.76 for ERSSTv5), which is a localized feature of that model.

**Phase transitions of the ERSSTv5 index.** The smoothed index changes sign at 1961-02, 1962-07, 1976-08 and 2002-05. The 1962 crossing falls within 3 years of 1961 and is merged. The counted transitions are therefore **1961, 1976 and 2002**.
- **1976 (August)** matches the documented 1976/77 shift, from negative to positive.
- **1961** matches no documented shift. It is a brief (17-month) and very weak (peak about 0.08) positive excursion of the smoothed index. The merge keeps the earlier year, so 1961 is counted even though the smoothed series returned to negative in 1962. I would not treat it as a real regime shift.
- **2002 (May)** is the positive-to-negative crossing. Under the code's ±2-year rule around 1998 it is not matched to the late-1990s shift. It lags 1998/99 by about 3–4 years, and the 9-year smoothing may well account for that lag. This is an interpretation, not something the code tested.
- Because 54 months are lost at each end of the smoothed series, any transition in the first or last 4.5 years cannot be detected.

**Model consistency with observations.**

| Model | Amplitude ratio to obs | Lag-12 difference from obs | Amplitude within ±25% | Persistence within ±0.1 |
|---|---|---|---|---|
| CESM2 | 1.10 | −0.118 | yes | no |
| GFDL-ESM4 | 0.88 | −0.037 | yes | yes |
| MPI-ESM1-2-LR | 1.24 | +0.084 | yes | yes |

MPI-ESM1-2-LR's amplitude is close to the edge of the tolerance. Each model is represented by a single historical simulation. A 65-year record gives a noisy lag-12 estimate, and differences of about 0.1 were not tested for significance, so I cannot say they are statistically distinguishable.

GFDL-ESM4 has the smallest explained variance (14.2%), yet its pattern is the second best.

## 4. Figures

- **Figure 1** is a 2×2 map of regression patterns (SST anomaly per SD of the standardized index) over 20–70°N, 110–260°E, in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR.
  - All panels use one `RdBu_r` colour scale centred on zero, with limits ±(98th percentile of |values| pooled over the four patterns). The colourbar is labelled "°C per standard deviation" and has extended ends.
  - Values beyond the limits therefore saturate.
  - Each panel is annotated with its explained variance and, for the models, its pattern correlation with ERSSTv5.
- **Figure 2** is a single time-series panel. It shows the standardized monthly indices of ERSSTv5, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, and the published NCEI PDO re-standardized as a red dashed fifth line.
  - A thick purple line shows the 108-month (9-year) centred running mean of the ERSSTv5 index.
  - Dotted vertical lines with year labels mark the counted transitions (1961, 1976, 2002).
  - The y axis is labelled "Standardized index (unitless)".

## 5. Caveats

- The EOF sign convention, the sqrt(cos lat) weighting and the global-mean removal are performed inside `compute_pdo_index`; the script does not independently verify them.
- The results use bilinear regridding, as prescribed. The native-grid differences (CESM2 and MPI-ESM1-2-LR are native `gn`, GFDL-ESM4 is 1° `gr`, ERSSTv5 is about 2°) mean the regridding affects the sources differently. Coarse grids cannot reproduce narrow features in amplitude.
- Offsets against the published PDO reflect baseline and construction differences rather than reconstruction errors.
- Because the pattern correlations are over about 4400–4600 spatially autocorrelated cells, no significance is attached to them.
