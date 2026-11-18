## Scientific summary

The PDO was evaluated over 1950-01 to 2014-12 for ERSSTv5 and the CMIP6 historical `tos` simulations CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. The executed code bilinearly regridded each SST field to the 1° grid with `regrid_horizontal(..., target="1deg", method="bilinear")`, then called `compute_pdo_index(field, depth_level=None, standardize=True)` for each source. The code did not separately call a detrending, domain-selection, climatology or global-mean-removal routine; those steps, if applied, are internal to `compute_pdo_index`. The returned PDO pattern fields were used on their 50×150 grid; the plots select 20–70°N and 110–260°E, i.e. 110°E–100°W after the dateline shift. Bilinear regridding was used as specified; results could differ from a conservative-regridded analysis, and bilinear regridding does not restore GFDL-ESM4’s native resolution.

Amplitude was computed as the area-weighted RMS of the regression pattern: the pattern was squared, passed to `compute_basic_stat(op="area_weighted_mean", dim=["lat","lon"])`, and square-rooted, giving °C per standard deviation. Lag-12-month autocorrelation was the direct Pearson correlation between the standardized monthly index and the same index shifted by 12 months. Pattern correlations were computed with `analyze_correlation(..., dim=["lat","lon"], weighting="area")` against the ERSSTv5 regression pattern, omitting missing cells separately for each pair. The published PDO correlation used `analyze_correlation(var1=ERSSTv5_index, var2=published_PDO, dim="time", detrend=False)`.

### Key numerical results

| Source | Leading-mode variance explained (%) | Amplitude (°C per std) | Lag-12-month autocorrelation | Pattern correlation with ERSSTv5 |
|---|---:|---:|---:|---:|
| ERSSTv5 | 23.0 | 0.330 | 0.342 | 1.000 |
| CESM2 | 26.1 | 0.364 | 0.223 | 0.925 |
| GFDL-ESM4 | 14.2 | 0.291 | 0.305 | 0.891 |
| MPI-ESM1-2-LR | 22.9 | 0.409 | 0.425 | 0.750 |

The ERSSTv5-derived standardized PDO index correlates with the published PDO series at **r = 0.962** (effective sample size **n_eff = 110.6**, p ≈ 0.0). This high correlation demonstrates strong agreement between two estimates of the same PDO quantity, including shared low-frequency variability and trend. It is **not** an independent validation and should not be interpreted as independent evidence of model skill or of a unique physical mechanism.

### Observed PDO phase transitions

The 9-year (108-month) centred running mean of the ERSSTv5 standardized index was computed with `temporal_aggregate(..., freq="rolling", window=108, center=True, op="mean")`. Raw sign crossings occurred at **1961, 1962, 1976 and 2002**. Applying the merge rule that crossings less than 3 years after the last counted transition are merged into the earlier year removed 1962, leaving observed transition years:

**1961, 1976, 2002**

Among these, **1976** corresponds to the documented mid-1970s Pacific climate shift. The late-1980s and late-1990s documented shifts were not detected by this 9-year running-mean sign-crossing criterion over 1950–2014. The 1961 and 2002 crossings are not among the listed documented shifts. This comparison uses the ERSSTv5 index crossings, not model phase timing, because the coupled historical runs are not expected to reproduce observed event timing.

### Model comparison

The best spatial pattern match is **CESM2**, with a pattern correlation of **0.925** against the ERSSTv5 PDO regression pattern.

Amplitude and persistence are mixed relative to ERSSTv5 (0.330 °C per std; lag-12 = 0.342):
- **CESM2**: amplitude 0.364 (about 11% higher), lag-12 0.223 (less persistent than observed), but best spatial pattern correlation.
- **GFDL-ESM4**: amplitude 0.291 (about 12% lower), lag-12 0.305 (close to observed, slightly less persistent), and pattern correlation 0.891; its leading-mode variance explained is notably lower (14.2%).
- **MPI-ESM1-2-LR**: amplitude 0.409 (about 24% higher), lag-12 0.425 (more persistent than observed), and the lowest pattern correlation (0.750).

Thus, no model is uniformly closest in all metrics: CESM2 best reproduces the spatial pattern, GFDL-ESM4 is closest in amplitude and persistence but has low explained variance, and MPI-ESM1-2-LR has near-observed explained variance but too-large amplitude and too-high persistence. These differences should be interpreted with the caveat that one member per model cannot separate forced change from internal variability, and model resolutions differ substantially.

### Visualizations

**Figure 1** is a 2×2 panel of PDO regression patterns in the order ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Each panel is a `pcolormesh` of the regression pattern over 20–70°N and 110–260°E, with a shared diverging `RdBu_r` colour scale centred on zero and clipped at the 98th percentile of pooled absolute values. The colourbar is labelled in °C per standard deviation. Each panel title gives the leading-mode explained variance, and model panels also give the spatial pattern correlation with ERSSTv5. The figure is designed to compare spatial structure and relative amplitude on a common scale.

**Figure 2** is a single time-series panel of the standardized monthly PDO indices from ERSSTv5, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, plus the published PDO series as a distinguishable dashed black line. The ERSSTv5 9-year centred running mean is overlaid, and the observed transition years 1961, 1976 and 2002 are marked with vertical dotted lines. The y-axis is labelled “Standardized index.” The figure is intended to compare temporal variability and phase behaviour, while the transition marking highlights the observed ERSSTv5 phase changes rather than model event timing.
