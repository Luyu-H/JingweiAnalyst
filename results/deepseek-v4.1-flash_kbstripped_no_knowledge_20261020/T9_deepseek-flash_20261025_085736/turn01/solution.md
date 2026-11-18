**Assessment.** The executed workflow computed the leading North Pacific SST mode for ERSSTv5 and for the three CMIP6 historical `tos` simulations over 1950-01 to 2014-12. The input monthly fields were loaded with `load_regular_data`, bilinearly regridded to the 1° grid with `regrid_horizontal`, and then passed directly to `compute_pdo_index` with `depth_level=None` and `standardize=True`. The code did **not** include a separate explicit call to `compute_anomaly_vs_own_climatology`; any anomaly/climatology handling is internal to the PDO tool. Regression patterns were compared over a common ocean mask built with `build_comparison_mask(ocean_only=True)`, using Pearson correlation with area weighting. Lag-12-month autocorrelations were computed on each standardized monthly index with `analyze_correlation(..., dim="time", lag=12, method="pearson")`. The published CPC/NCEI PDO series was loaded from `cpc_indices` product `pdo` and correlated with the ERSSTv5-derived index over the same time axis. The 9-year smoothed series was formed with `temporal_aggregate(freq="rolling", window=108, center=True, op="mean")`; zero crossings were detected from that smoothed ERSSTv5 index and crossings less than 3 years after the last counted transition were skipped.

### Key numerical results

| Source | Variance explained (%) | RMS regression amplitude (°C SD⁻¹) | Peak |regression| (°C SD⁻¹) | Lag-12 autocorrelation | Pattern correlation vs ERSSTv5 |
|---|---:|---:|---:|---:|---:|
| ERSSTv5 | 23.0 | 0.330 | 0.762 | 0.342 | 1.000 |
| CESM2 | 26.1 | 0.364 | 1.042 | 0.223 | 0.926 |
| GFDL-ESM4 | 14.2 | 0.291 | 1.114 | 0.305 | 0.895 |
| MPI-ESM1-2-LR | 22.9 | 0.409 | 2.352 | 0.425 | 0.752 |

- The ERSSTv5-derived standardized PDO index correlates with the published PDO series at **r = 0.962** over 1950–2014. This indicates strong temporal/phase agreement with an existing PDO index. It does **not** independently validate the EOF/regression calculation, prove a physical mechanism, or demonstrate model skill beyond agreement with that shared index.
- Pattern correlations measure only spatial similarity of the regression patterns. Amplitude and persistence must be judged separately from the regression amplitude and lag-12 autocorrelation.
- **CESM2** has the highest spatial pattern correlation with ERSSTv5 (**0.926**), followed by GFDL-ESM4 (**0.895**) and MPI-ESM1-2-LR (**0.752**).
- CESM2 captures the pattern best, but its lag-12 persistence is lower than observed (**0.223** vs **0.342**). Its RMS amplitude is slightly higher than observed (**0.364** vs **0.330**).
- GFDL-ESM4 has a high pattern correlation and relatively close persistence (**0.305**), but it explains less variance (**14.2%** vs **23.0%**) and has a slightly lower RMS amplitude (**0.291**).
- MPI-ESM1-2-LR has the weakest pattern correlation and is over-amplified and over-persistent: RMS amplitude **0.409** vs **0.330**, peak |regression| **2.352** vs **0.762**, and lag-12 autocorrelation **0.425** vs **0.342**.

### Observed transition years

The smoothed ERSSTv5 index produced raw zero crossings at **1961, 1962, 1976, and 2002**. Applying the less-than-3-year merging rule kept the earlier year, yielding transition years:

**1961, 1976, 2002**

Among the documented Pacific climate shifts checked in the code (1976/77, 1988/89, 1998/99), only **1976** matches a detected transition, corresponding to the **1976/77** shift. The **1988/89** and **1998/99** shifts were not detected by this smoothed zero-crossing method. The years **1961** and **2002** are detected transitions in this analysis but are not mapped to those documented shifts in the executed code.

### Visualization description

- **Figure 1** is a 2×2 panel of PDO regression patterns in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. It uses one shared diverging `RdBu_r` colour scale centred on zero, with limits set symmetrically from the pooled 98th percentile of absolute regression values. The colour bar is labelled in °C per standard deviation. Each panel is annotated with its variance explained, and model panels are also annotated with their spatial pattern correlation against ERSSTv5. The displayed domain is 20°N–70°N and 110°E–260°E.
- **Figure 2** is a single time-series panel showing the standardized monthly indices from all four sources plus the standardized published PDO series as a fifth, distinguishable line. It also overlays the 9-year centred running mean of the ERSSTv5 standardized index and marks the detected transition years with vertical dotted lines and year labels. The y-axis is labelled “Standardized index.”

### Main interpretation

CESM2 reproduces the observed PDO spatial pattern best, but none of the models is fully consistent with observations in all respects. CESM2 has the strongest pattern match but too little year-to-year persistence. GFDL-ESM4 has a reasonably high pattern correlation and near-observed persistence, but its explained variance and RMS amplitude are lower. MPI-ESM1-2-LR has the weakest spatial pattern and is too strong and too persistent. The correlation with the published PDO confirms close temporal agreement between the ERSSTv5-derived index and an existing PDO product, but it is not independent validation of the EOF calculation or of physical mechanisms.
