**Direct answer**

Using the executed analysis, the three CMIP6 models reproduce the observed PDO spatial pattern with the following area-weighted pattern correlations against the ERSSTv5-derived PDO regression pattern: **CESM2 r = 0.925**, **GFDL-ESM4 r = 0.891**, and **MPI-ESM1-2-LR r = 0.750**. Therefore, **CESM2 best reproduces the observed PDO pattern**. Its amplitude is close to observed, but its lag-12 persistence is lower. GFDL-ESM4 has the closest lag-12 persistence and a near-observed amplitude, but lower explained variance and slightly lower pattern correlation. MPI-ESM1-2-LR has the largest amplitude and highest persistence but the weakest pattern match.

**Methodology actually executed**

- Time range: **1950-01 through 2014-12** (780 monthly samples).
- SST sources: ERSSTv5 (`ersst`, variable `sst`), and CMIP6 historical `tos` for **CESM2**, **GFDL-ESM4**, and **MPI-ESM1-2-LR**.
- CMIP6 `tos` was converted from K to °C if its units metadata indicated kelvin, bilinearly regridded to a common **1° grid** with `regrid_horizontal(..., target="1deg", method="bilinear")`, and if a depth dimension was present the first depth index was selected and dropped. `compute_pdo_index` was then called with `depth_level=None`.
- ERSSTv5 was also bilinearly regridded to the same 1° grid.
- Seasonal cycles were removed using `compute_anomaly_vs_own_climatology(..., freq="month", base_period=["1950-01","2014-12"], method="subtract")`. Thus each SST source used **its own 1950–2014 monthly climatology**, not an observed climatology applied to models.
- PDO modes were obtained with `compute_pdo_index(..., depth_level=None, standardize=True)`. The returned regression patterns have dimensions **lat = 50, lon = 150**, displayed over the PDO sector **20–70°N, 110–260°E** with longitudes continuous across the dateline.
- Mode amplitude was computed as the **area-weighted RMS** of each regression pattern using `compute_basic_stat(..., op="area_weighted_mean", dim=["lat","lon"])` followed by a square root. Units are **°C per standard deviation**.
- Lag-12 persistence was computed with `analyze_correlation(var1=index, var2=index, dim="time", method="pearson", lag=12)`, i.e. a single lag of **12 months** on the standardized monthly index.
- Spatial pattern correlations used `analyze_correlation` on each model regression pattern versus the ERSSTv5 regression pattern, with `dim=["lat","lon"]`, Pearson correlation, and `weighting="area"`, over shared cells.
- The published PDO was loaded from `cpc_indices`, product `pdo`, over 1950–2014, with no region argument. Its loaded units are dimensionless (`"1"`). The checked metadata did not provide a base period; the code recorded **“see product metadata.”** For plotting only, the published PDO was standardized over 1950–2014; for correlation it was used as loaded, without detrending.
- The ERSSTv5-derived standardized index was correlated with the published PDO using Pearson correlation after inner time alignment. Significance was taken from `analyze_correlation`’s effective sample size, not the raw 780 months.
- The ERSSTv5 standardized index was smoothed with a **108-month (9-year) centred running mean** using `temporal_aggregate(..., freq="rolling", op="mean", window=108, center=True)`. Sign changes were identified from the smoothed series; crossings less than 3 years after the last counted transition were merged into that earlier transition.

**Key numerical results**

| Source | Variance explained by leading mode (%) | Amplitude (°C per SD) | Lag-12 autocorrelation | Pattern r vs ERSSTv5 | Shared cells |
|---|---:|---:|---:|---:|---:|
| ERSSTv5 | 23.0 | 0.330 | 0.342 | — | — |
| CESM2 | 26.1 | 0.364 | 0.223 | 0.925 | 4445 |
| GFDL-ESM4 | 14.2 | 0.291 | 0.305 | 0.891 | 4624 |
| MPI-ESM1-2-LR | 22.9 | 0.409 | 0.425 | 0.750 | 4368 |

- ERSSTv5-derived standardized PDO index versus published PDO: **Pearson r = 0.962**, p ≈ 0.0, **n = 780**, **effective n = 110.6**. This demonstrates temporal agreement between the two PDO estimates; it does **not** demonstrate equal variance or equal amplitude.
- The lag-12 autocorrelation p-values from the executed analysis were: ERSSTv5 0.000471, CESM2 0.0690, GFDL-ESM4 0.000454, MPI-ESM1-2-LR 0.00591.
- All model spatial pattern correlations with ERSSTv5 had p ≈ 0.0 in the reported output.

**Observed PDO phase transitions**

The raw sign changes of the 108-month centred running mean of the ERSSTv5 standardized PDO index were **1961, 1962, 1976, and 2002**. Merging crossings less than 3 years after the last counted transition into the earlier transition gave:

**1961, 1976, 2002**

Comparison with the documented shifts listed in the code:

- **1976** corresponds to the **1976/77 PDO regime shift**.
- **1961** and **2002** were not among the commonly cited documented shifts explicitly listed in the code.
- The code’s documented-shift list also included 1989, 1998/99, and 2014, but those were not returned as detected transitions from the smoothed ERSSTv5 index. The centred 108-month smoother has missing endpoints (673 valid smoothed months), so near-end transitions such as 2014 may be unresolved.

**Interpretation**

- **Best pattern reproduction:** **CESM2**, with spatial pattern correlation **0.925**.
- **Amplitude comparison without an imposed threshold:** ERSSTv5 amplitude is **0.330 °C per SD**. CESM2 is slightly higher at **0.364**, GFDL-ESM4 slightly lower at **0.291**, and MPI-ESM1-2-LR higher at **0.409**. CESM2 has the closest amplitude in absolute difference, followed by GFDL-ESM4, then MPI-ESM1-2-LR.
- **Persistence comparison:** ERSSTv5 lag-12 autocorrelation is **0.342**. GFDL-ESM4 is closest at **0.305**, CESM2 is lower at **0.223**, and MPI-ESM1-2-LR is higher at **0.425**.
- **Variance explained:** ERSSTv5 gives **23.0%**; CESM2 **26.1%**, GFDL-ESM4 **14.2%**, and MPI-ESM1-2-LR **22.9%**. GFDL-ESM4’s leading mode explains substantially less variance over this common grid and period.
- Bilinear regridding to 1° was used for all SST fields because that was the prescribed comparison method in the executed code. The code’s printed note states that conservative regridding would be preferred when coarsening before an area-weighted step. The executed code did not quantify native-grid resolution differences among the models.

**Visualizations generated**

- **Figure 1** is a 2×2 panel of PDO regression patterns in the order **ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR**. It uses one shared diverging `RdBu_r` colour scale centred on zero, with limits set from the 98th percentile of the absolute pooled regression values. The colourbar is labelled **“SST regression (°C per SD).”** Each panel is annotated with its variance explained; model panels also show the spatial pattern correlation with ERSSTv5. The panels display the PDO sector over **110–260°E** and **20–70°N**.
- **Figure 2** is a single time-series panel showing the standardized indices of all four SST sources plus the **published PDO as a distinct fifth dashed black line**. It also overlays the **ERSSTv5 9-year centred running mean** and marks the detected transition years with vertical dotted lines and labels. The y-axis is labelled **“Standardized index.”**
