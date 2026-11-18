## Key numerical results

The PDO was defined as the leading area-weighted EOF of monthly SST anomalies over 20°N–70°N, 110°E–100°W, relative to the 1950–2014 monthly climatology. All fields were bilinearly regridded to a 1° grid. The index was standardized; the regression patterns are in °C per standard deviation of the standardized index.

| Source | Variance explained by leading mode (%) | Mode amplitude (°C SD⁻¹) | Lag-12-month autocorrelation | Spatial pattern correlation with ERSSTv5 |
|---|---:|---:|---:|---:|
| ERSSTv5 | 22.8 | 0.326 | 0.342 | 1.000 |
| CESM2 | 25.9 | 0.383 | 0.536 | 0.776 |
| GFDL-ESM4 | 14.5 | 0.285 | 0.218 | 0.903 |
| MPI-ESM1-2-LR | 22.0 | 0.416 | 0.592 | 0.689 |

The correlation between the ERSSTv5-derived PDO index and the published PDO series over 1950–2014 is **0.957**.

## Interpretation of the published-PDO correlation

The value **r = 0.957** shows that the ERSSTv5-derived leading EOF index has very similar month-to-month temporal variability and phase to the published PDO index after sign orientation. It supports the interpretation that the leading North Pacific SST EOF in ERSSTv5 captures the canonical PDO time series.

However, this correlation **does not** demonstrate that:

- the spatial PDO pattern is identical to the published product;
- the EOF is a unique physical mode or a causal mechanism;
- the CMIP6 models are validated by this correlation;
- the agreement is independent evidence, because both the published PDO and the ERSSTv5-derived EOF are based on similar SST data and EOF-type methodology and may share data/method lineage.

It is a temporal consistency check, not a full spatial or process-level validation.

## Observed PDO phase transitions and documented climate shifts

Using the specified 9-year (108-month) centred running mean of the ERSSTv5 index, with crossings less than 3 years after the last counted transition merged into the earlier transition, the detected transition years are:

**1960, 1963, 1976, 2001**

Of these, only **1976** corresponds to the documented Pacific climate-shift list used here (the canonical 1976/77 shift). The algorithm does **not** identify the 1988/89 or late-1990s documented shifts as separate sign-change transitions under this smoothing/merging rule. The year 2001 is near the late-1990s period but is not within the ±1-year match criterion used for the documented-shift list here. Thus, under this specific transition definition, the clear documented match is **1976/77**; the result is sensitive to the smoothing window and merging rule.

## Which model performs best?

**GFDL-ESM4 reproduces the observed spatial PDO pattern best**, with a pattern correlation of **0.903** against ERSSTv5. The area-weighted correlation is similarly high at about **0.899**. CESM2 is intermediate (**0.776**), and MPI-ESM1-2-LR is poorest (**0.689**).

## Are model amplitude and persistence consistent with observations?

The observed mode amplitude is **0.326 °C SD⁻¹** and the observed lag-12 autocorrelation is **0.342**.

- **CESM2:** amplitude **0.383** (+0.057, slightly too strong) and lag-12 autocorrelation **0.536** (+0.194, much too persistent). Its amplitude is close, but its persistence is substantially too high.
- **GFDL-ESM4:** amplitude **0.285** (−0.041, slightly too weak) and lag-12 autocorrelation **0.218** (−0.123, insufficiently persistent). It has the best spatial pattern and a reasonably close amplitude, but it underestimates PDO persistence and explains much less variance (14.5%) than observed (22.8%).
- **MPI-ESM1-2-LR:** amplitude **0.416** (+0.090, too strong) and lag-12 autocorrelation **0.592** (+0.251, much too persistent). It also has the lowest spatial pattern correlation.

Overall, the model amplitudes are all within about **0.09 °C SD⁻¹** of the observed value, so they are broadly of the right order, with CESM2 and GFDL-ESM4 closest. Persistence is less consistent: CESM2 and MPI-ESM1-2-LR are too persistent, while GFDL-ESM4 is under-persistent. No single model reproduces the observed pattern, amplitude, and persistence simultaneously. **GFDL-ESM4 is the best for spatial pattern; CESM2 is closest for amplitude; none captures the observed persistence accurately.**

## Figures

The generated **Figure 1** shows the 2 × 2 PDO regression patterns in the requested order (ERSSTv5, CESM2, GFDL-ESM4, MPI-ESM1-2-LR) on a common diverging colour scale centred on zero, with explained variance and pattern correlations annotated. Visually, GFDL-ESM4 is the closest to the ERSSTv5 pattern, consistent with its highest pattern correlation.

The generated **Figure 2** shows the standardized indices from all four sources plus the published PDO series, together with the 9-year smoothed ERSSTv5 curve and the marked transition years. It illustrates the high temporal agreement between the ERSSTv5-derived and published PDO indices, and the differing model persistence and phase behaviour.
