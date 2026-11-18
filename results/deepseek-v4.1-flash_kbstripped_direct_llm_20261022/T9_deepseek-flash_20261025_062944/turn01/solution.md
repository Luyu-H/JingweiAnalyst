Using the specified PDO definition (leading EOF of monthly SST anomalies over 20°N–70°N, 110°E–100°W, relative to a 1950–2014 monthly climatology, bilinearly regridded to 1°), the diagnostics are:

| Source | Variance explained (%) | Peak amplitude (°C sd⁻¹) | Area-weighted RMS amplitude (°C sd⁻¹) | Lag-12 autocorrelation | Pattern correlation with ERSSTv5 |
|---|---:|---:|---:|---:|---:|
| ERSSTv5 | 23.0 | 0.679 | 0.327 | 0.341 | — |
| CESM2 | 26.4 | 1.167 | 0.390 | 0.524 | −0.793 |
| GFDL-ESM4 | 14.3 | 0.879 | 0.295 | 0.222 | 0.903 |
| MPI-ESM1-2-LR | 22.7 | 2.451 | 0.428 | 0.591 | 0.699 |

The Pearson correlation between the ERSSTv5-derived standardized index and the published PDO series is **−0.962**. This very high absolute correlation shows that the leading EOF captures essentially the same temporal evolution as the published PDO index, but with the **opposite sign**. The sign of an EOF is arbitrary, and the sign convention used in the code (forcing warm western North Pacific and cool central/eastern North Pacific) is evidently reversed relative to the official PDO convention. Thus, after sign reversal, the temporal agreement is about +0.962. This correlation demonstrates temporal/methodological consistency with the published PDO; it does **not** demonstrate correctness of the regression amplitude, the spatial pattern, model skill, or physical causality.

For the models:
- **GFDL-ESM4** reproduces the observed ERSSTv5 regression pattern best, with a spatial pattern correlation of **0.903**.
- **MPI-ESM1-2-LR** has a moderate pattern correlation of **0.699**.
- **CESM2** has a pattern correlation of **−0.793**, i.e. its leading pattern is strongly reversed relative to the ERSSTv5 pattern under the chosen sign convention. If sign is treated as arbitrary, its absolute pattern correlation is 0.793, but it is still not the best match.

Amplitude and persistence compare as follows. ERSSTv5 has peak amplitude **0.679°C sd⁻¹** and lag-12 autocorrelation **0.341**.
- **GFDL-ESM4** is closest in amplitude: peak amplitude is 1.29× observed, RMS amplitude is 0.90× observed. Its persistence is slightly too weak (lag-12 autocorrelation 0.222, i.e. −0.120 vs observed).
- **CESM2** overestimates the peak amplitude (1.72× observed; RMS 1.19× observed) and is too persistent (lag-12 autocorrelation 0.524, i.e. +0.183 vs observed).
- **MPI-ESM1-2-LR** strongly overestimates the peak amplitude (3.61× observed; RMS 1.31× observed) and is also too persistent (lag-12 autocorrelation 0.591, i.e. +0.250 vs observed).

The observed PDO phase-transition years obtained from sign changes of the 108-month centred running mean of the ERSSTv5-derived index are **1960, 1963, 1976, and 2001**. The documented Pacific climate shifts commonly cited for this period are **1976/77, 1988/89, 1998/99, and 2014/15**. Under the stated matching criterion (±1 year), only **1976** corresponds to a documented shift. The 2014 transition cannot be detected here because the 108-month centred running mean is undefined near the end of the record; 1988 and 1998 are not sign changes in this particular ERSSTv5-derived series. The 1960 and 1963 crossings are close together (the merging rule is strict “less than 3 years”, so 3 years apart are retained separately) and are not among the canonical listed shifts; 2001 is also not one of the four documented shifts.

Overall, **GFDL-ESM4** gives the best spatial pattern and the most reasonable amplitude and persistence among the three models. CESM2 and MPI-ESM1-2-LR produce PDO-like variability but are too strong and/or too persistent, with MPI-ESM1-2-LR especially over-amplified. The figures generated show the four regression patterns on a common diverging colour scale and the standardized time series together with the published PDO and the smoothed ERSSTv5 index with transition years marked.
