## Methodology actually used

The analysis covered **1950-01 through 2014-12** with baseline/climatology period **1950-01 to 2014-12**. Global monthly SST fields were loaded with `load_regular_data`: ERSSTv5 `sst` from `ersst`, and CMIP6 historical `tos` for **CESM2, GFDL-ESM4, and MPI-ESM1-2-LR**. Each source was bilinearly regridded to a common **1° grid** using `regrid_horizontal(target="1deg", method="bilinear")`. Monthly anomalies were computed relative to **each source’s own 1950–2014 monthly climatology** with `compute_anomaly_vs_own_climatology(freq="month", base_period=["1950-01","2014-12"])`.

The leading North Pacific SST mode and standardized PDO index were obtained with `compute_pdo_index` applied to the global anomaly fields; the code notes that this tool removes the global-mean anomaly itself and cuts the PDO domain **20°N–70°N, 110°E–100°W**. The returned explained variance, EOF, and regression pattern were used. Amplitude was computed from the regression pattern as the **maximum absolute regression value** in °C per standard deviation of the standardized index; the code also computed a domain-RMS amplitude as a supplementary diagnostic. Lag-1-year persistence was computed as the **12-month autocorrelation** of each standardized monthly index with `analyze_correlation(..., dim="time", lag=12)`. Pattern correlations between model regression patterns and the ERSSTv5 regression pattern were computed with `analyze_correlation(..., dim=["lat","lon"], method="pearson", weighting="area")`, i.e. **area-weighted spatial Pearson correlation**.

The published PDO series was loaded from `cpc_indices` with product `pdo`, restricted to the same time range, and correlated with the ERSSTv5-derived standardized index over time using `analyze_correlation(..., dim="time")` with no lag. The 9-year smoothed ERSSTv5 index was formed with `temporal_aggregate(freq="rolling", op="mean", window=108, center=True)`. Sign changes were identified from the 108-month centred running mean; crossings less than 3 years after the last counted transition were merged into the earlier transition. Documented shifts checked by the code were **1976/77, 1988/89, and 1998/99**, with a match allowed within 1 year.

## Key numerical results

| Source | Explained variance (%) | Peak amplitude (max |regression|) (°C SD⁻¹) | Domain RMS amplitude (°C SD⁻¹) | Lag-12 autocorrelation |
|---|---:|---:|---:|---:|
| ERSSTv5 | 23.0 | 0.762 | 0.326 | 0.342 |
| CESM2 | 26.1 | 1.042 | 0.366 | 0.223 |
| GFDL-ESM4 | 14.2 | 1.114 | 0.293 | 0.305 |
| MPI-ESM1-2-LR | 22.9 | 2.352 | 0.419 | 0.425 |

The ERSSTv5-derived standardized PDO index correlates with the published PDO series at **r = 0.962** (p reported as 0.0, n = 780). This demonstrates strong temporal agreement between the ERSSTv5-derived index and the published observation-derived PDO index over 1950–2014. It does **not** demonstrate independence, causality, or a physical mechanism, because both series are observation-derived and are not independent validation datasets.

Spatial pattern correlations of each model’s PDO regression pattern with the ERSSTv5 regression pattern:

| Model | Pattern correlation vs ERSSTv5 |
|---|---:|
| CESM2 | 0.925 |
| GFDL-ESM4 | 0.891 |
| MPI-ESM1-2-LR | 0.750 |

**CESM2 reproduces the observed spatial PDO regression pattern best** by this metric. Pattern correlation assesses spatial similarity only; it does not assess amplitude, timing, or persistence. Using the code’s peak-amplitude definition, all three models overestimate the observed peak amplitude, especially MPI-ESM1-2-LR (2.352 vs 0.762 °C SD⁻¹). Using the supplementary domain-RMS amplitude, CESM2 and MPI-ESM1-2-LR are somewhat larger than observed, while GFDL-ESM4 is smaller. For lag-12 persistence, GFDL-ESM4 is closest to ERSSTv5 (0.305 vs 0.342), CESM2 is lower (0.223), and MPI-ESM1-2-LR is higher (0.425). Thus the models’ amplitude and persistence are only partly consistent with observations and do not uniformly match the observed spatial pattern ranking.

## Observed PDO phase transitions

The sign changes of the 108-month centred running mean of the ERSSTv5 standardized PDO index occurred at **1961, 1962, 1976, and 2002**. After merging a crossing less than 3 years after the previous counted transition, keeping the earlier year, the retained transitions are:

**1961, 1976, 2002**

Among the documented Pacific climate shifts checked by the code—1976/77, 1988/89, and 1998/99—only **1976** is matched, corresponding to the **1976/77** shift. Under this smoothed-index sign-change criterion, no transitions were matched to 1988/89 or 1998/99. The 1961 and 2002 transitions were detected but were not matches to the listed documented shifts.

## Figures generated

**Figure 1** is a 2×2 panel of the PDO regression patterns over 20°N–70°N, 110°E–100°W in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. It uses a shared diverging colour scale centred on zero, labelled in °C per standard deviation, with a zero contour. Each panel is annotated with its explained variance; model panels are also annotated with their area-weighted spatial pattern correlation against the ERSSTv5 regression pattern. The colour scale is clipped at the 98th percentile of the pooled absolute regression values across the four panels.

**Figure 2** is a single time-series panel showing the standardized PDO index for ERSSTv5, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR, plus the published PDO series as a distinguishable fifth line. It overlays the 108-month centred running mean of the ERSSTv5 index and marks the detected transition years with vertical dotted lines. The y-axis is labelled as the standardized PDO index.

## Main insights

- ERSSTv5’s leading North Pacific mode explains **23.0%** of monthly SST anomaly variance over the PDO domain. CESM2 and MPI-ESM1-2-LR are close to this (26.1% and 22.9%), while GFDL-ESM4 explains substantially less (14.2%).
- The ERSSTv5-derived index is highly correlated with the published PDO series (**r = 0.962**), confirming close temporal agreement with an observation-derived PDO product, but not independent validation or causality.
- CESM2 has the highest spatial pattern correlation with the observed regression pattern (**0.925**), followed by GFDL-ESM4 (**0.891**) and MPI-ESM1-2-LR (**0.750**).
- Amplitude and persistence checks separate from pattern correlation show mixed consistency: all models overestimate peak regression amplitude, and lag-12 autocorrelation ranges from too low in CESM2 to too high in MPI-ESM1-2-LR, with GFDL-ESM4 closest for persistence.
- The smoothed ERSSTv5 index identifies retained phase transitions in **1961, 1976, and 2002**; only **1976** corresponds to the documented 1976/77 Pacific climate shift within the code’s matching criterion.
