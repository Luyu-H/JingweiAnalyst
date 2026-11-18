I assessed the PDO as the leading EOF of monthly North Pacific SST anomalies over 20°N–70°N, 110°E–100°W, relative to a 1950–2014 monthly climatology. The target grid was 1°, with area weighting by \(\sqrt{\cos\phi}\). The leading PC was standardized, and the regression pattern is in °C per one standard deviation.

| Source | Variance explained (%) | Amplitude, max |reg| (°C/std) [RMS] | Lag-12 autocorrelation | Pattern correlation vs ERSSTv5 |
|---|---:|---:|---:|---:|
| ERSSTv5 | 22.8 | 0.679 [0.328] | 0.342 | 1.000 |
| CESM2 | 25.8 | 1.073 [0.384] | 0.530 | 0.790 |
| GFDL-ESM4 | 15.1 | 0.822 [0.300] | 0.219 | 0.900 |
| MPI-ESM1-2-LR | 21.6 | 2.405 [0.415] | 0.572 | 0.700 |

The ERSSTv5-derived PDO index correlates with the published PDO series at **r = 0.957**. This demonstrates that the leading EOF over the chosen domain closely reproduces the canonical temporal PDO variability. It does **not** demonstrate that the two products are identical or interchangeable: the published PDO is based on ERSSTv6 and uses its own preprocessing, baseline, and global-mean treatment. It also does not by itself validate the simulated spatial pattern or mechanisms.

Using the 9-year centered running mean of the ERSSTv5-derived index, the detected sign-change years over 1950–2014 are:

- **1960** (negative→positive)
- **1963** (positive→negative)
- **1976** (negative→positive)
- **2001** (positive→negative)

Applying the stated rule that crossings less than 3 years after the last counted transition are merged, these remain the detected transitions. Of these, only **1976** corresponds directly to a documented Pacific climate shift in the supplied list, namely the **1976/77 Pacific climate shift**. The 1960, 1963, and 2001 crossings are not among the listed canonical shift years; 2001 is near the 1998/99 shift in time but is not the same year. The method as configured did not flag 1989, 1998, or 2014.

**Model assessment:** GFDL-ESM4 reproduces the observed spatial regression pattern best, with pattern correlation **0.900**, compared with **0.790** for CESM2 and **0.700** for MPI-ESM1-2-LR. However, no model is fully consistent with the observations in both amplitude and persistence. All models overestimate the maximum regression amplitude relative to ERSSTv5 (0.679 °C/std): GFDL-ESM4 is closest at 0.822, CESM2 is 1.073, and MPI-ESM1-2-LR is much too strong at 2.405. In RMS amplitude, GFDL-ESM4 (0.300) is closest to observed (0.328), while CESM2 (0.384) and MPI-ESM1-2-LR (0.415) are high. For persistence, observed lag-12 autocorrelation is 0.342: CESM2 (0.530) and MPI-ESM1-2-LR (0.572) are too persistent, whereas GFDL-ESM4 (0.219) is too weak. Thus, GFDL-ESM4 best captures the pattern, but not the persistence; CESM2 has reasonable variance but too much amplitude and memory; MPI-ESM1-2-LR has the poorest pattern, strongly excessive amplitude, and excessive persistence.

Figure 1 shows the four regression patterns on one common diverging color scale, annotated with variance explained and model–observation pattern correlations. Figure 2 shows the standardized ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, and published PDO series, with the 9-year smoothed ERSSTv5 index and detected transition years marked.
