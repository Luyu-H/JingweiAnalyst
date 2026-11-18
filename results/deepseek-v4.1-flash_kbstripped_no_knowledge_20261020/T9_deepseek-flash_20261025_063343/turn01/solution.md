**Overall assessment.** The executed code compared the leading North Pacific SST mode from ERSSTv5 and three CMIP6 historical simulations (CESM2, GFDL-ESM4, MPI-ESM1-2-LR) over 1950-01 to 2014-12. All four sources produce a leading mode, but the models differ in variance explained, amplitude, persistence, and spatial pattern agreement. By the signed area-weighted pattern correlation computed by the code, **MPI-ESM1-2-LR is the best reproduction of the ERSSTv5 regression pattern (r = 0.765)**. However, CESM2 and GFDL-ESM4 have larger-magnitude but negative pattern correlations (−0.929 and −0.891), which is important because the sign of an EOF/PC is arbitrary; a sign-insensitive comparison would favor CESM2 and GFDL-ESM4. The ERSSTv5-derived index is strongly anti-correlated with the published PDO series (r = −0.957), indicating strong temporal correspondence up to sign, not causality or independent model validation.

**Methodology actually applied in the executed code**
- Time range: 1950-01 to 2014-12 for both loading and anomaly baseline.
- Data loaded with `load_regular_data`: ERSSTv5 `sst`; CMIP6 historical `tos` for CESM2, GFDL-ESM4, MPI-ESM1-2-LR; CPC `pdo` product.
- Units: temperature fields were checked and converted from K to °C if their units indicated Kelvin; otherwise used as °C. Regression and amplitude are in °C per standard deviation.
- Regridding: every SST/tos field was bilinearly regridded with `regrid_horizontal(..., target="1deg", method="bilinear")`.
- Masking: a 1° land/sea mask was loaded with `load_landsea_mask() > 1`, reindexed to the reference 1° grid by nearest latitude/longitude, filled with `False`, and applied so only ocean cells were retained.
- Anomalies: each source was converted to monthly anomalies relative to **its own** 1950–2014 monthly climatology using `compute_anomaly_vs_own_climatology(..., freq="month", base_period=["1950-01","2014-12"])`. No additional detrending, screening, or filtering was applied beyond the anomaly calculation and the later 108-month running mean.
- PDO mode: for each source, `compute_pdo_index(anoms[k], standardize=True)` was used. Per the code comment, this tool internally selects the 20°N–70°N, 110°E–100°W PDO box and removes the global-mean anomaly. The returned leading-mode variance-explained value was taken directly from `compute_pdo_index`.
- Sign orientation: the code attempted to fix the arbitrary EOF sign by requiring the mean regression loading in the eastern North Pacific (20°N–45°N, 200–232°E) to be positive; if negative, the index, regression pattern, and EOF were multiplied by −1.
- Amplitude: the spatial standard deviation of each regression pattern was computed with `compute_basic_stat(..., op="std", dim=["lat","lon"])`, giving °C per standard deviation of the standardized index.
- Lag-1-year persistence: `analyze_correlation(index, index, dim="time", lag=12)` was used, i.e. a 12-month lag autocorrelation of each standardized monthly index.
- ERSSTv5-derived index vs published PDO: `analyze_correlation(..., dim="time")` gave the temporal correlation.
- Pattern correlation: area-weighted spatial Pearson correlation via `analyze_correlation(model_regression, ERSSTv5_regression, dim=["lat","lon"], method="pearson", weighting="area")`.
- Smoothing and transitions: the ERSSTv5 standardized index was smoothed with `temporal_aggregate(..., freq="rolling", window=108, center=True)`. Zero crossings of the smoothed series were found, and crossings less than 3 years after the last counted transition were merged by keeping the earlier transition year.

**Key numerical results (rounded as requested)**

| Source | Variance explained (%) | Amplitude (°C per SD) | Lag-12 autocorrelation | Pattern correlation vs ERSSTv5 |
|---|---:|---:|---:|---:|
| ERSSTv5 | 24.9 | 0.321 | 0.355 | not computed in code |
| CESM2 | 26.5 | 0.354 | 0.230 | −0.929 |
| GFDL-ESM4 | 15.1 | 0.287 | 0.306 | −0.891 |
| MPI-ESM1-2-LR | 23.1 | 0.396 | 0.427 | 0.765 |

- **ERSSTv5-derived index vs published PDO:** computed correlation = **−0.957**. This is a strong linear temporal relationship with opposite sign. It demonstrates temporal phase/amplitude agreement up to the arbitrary sign convention; it does **not** by itself demonstrate causality, independently validate the published PDO index, or establish CMIP6 model skill.
- **Observed amplitude:** 0.321 °C per SD. CESM2 (0.354) is slightly high, GFDL-ESM4 (0.287) is lower, and MPI-ESM1-2-LR (0.396) is the highest.
- **Observed lag-12 persistence:** 0.355. GFDL-ESM4 (0.306) is closest, CESM2 (0.230) has weaker persistence, and MPI-ESM1-2-LR (0.427) has stronger persistence.
- **Variance explained:** observed 24.9%; CESM2 (26.5%) and MPI-ESM1-2-LR (23.1%) are close, while GFDL-ESM4 (15.1%) is substantially lower.
- **Pattern comparison:** signed maximum is MPI-ESM1-2-LR (0.765), so the code identifies it as the best signed match. But CESM2 (−0.929) and GFDL-ESM4 (−0.891) have larger absolute correlations with opposite sign, showing that spatial similarity is present but the EOF sign convention was not aligned for those models. Because EOF/PC sign is arbitrary, the negative signs should not be interpreted as poorer spatial structure without sign correction.

**Observed PDO transition years**
- Zero crossings before merging: 1961, 1962, 1976, 2002.
- After applying the rule that crossings less than 3 years after the last counted transition are merged into the earlier one: **1961, 1976, 2002**.
- Comparison with documented Pacific climate shifts:
  - **1976** matches the well-known 1976/77 shift.
  - **2002** is nearest to the 1998/99 shift in the code’s nearest-match table, but it is about 4 years later and not an exact match.
  - **1961** has no listed documented shift in [1976, 1988, 1998, 2014].
  - The documented **1988/89** and **2014** shifts were not detected by this 108-month smoothed zero-crossing procedure.
- These transition results are sensitive to smoothing and to the EOF/index sign convention.

**Visualizations generated**
- **Figure 1** is a 2×2 panel of PDO regression patterns in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. It uses one shared diverging colour scale centred on zero, labelled in °C per standard deviation. The colour limits were set from the 98th percentile of absolute regression values pooled across all four sources. Each panel is annotated with its variance explained, and model panels also show their pattern correlation with ERSSTv5. The ERSSTv5 panel does not show a model-style pattern correlation.
- **Figure 2** is a single time-series panel of the standardized PDO indices for all four sources plus the published PDO series as a dashed fifth line. It also shows the thick 9-year centred running mean of the ERSSTv5 index and marks the detected transition years with vertical dotted lines. The y-axis is labelled “Standardized index.”

**Main insight.** The models capture a leading North Pacific SST mode, but fidelity is mixed. MPI-ESM1-2-LR has the best signed pattern correlation, yet its amplitude and persistence are somewhat too strong. CESM2 and GFDL-ESM4 have high absolute spatial similarity with ERSSTv5 but opposite sign under the stored orientation, so their signed correlations are negative. The observed PDO’s 9-year smoothed transitions over this window are detected as 1961, 1976, and 2002; only 1976 clearly corresponds to the canonical 1976/77 Pacific climate shift, while 1988/89 and 2014 are not detected by this smoothing/merging procedure.
