**1. What was done (as executed)**

Each CMIP6 model (CESM2, GFDL‑ESM4, MPI‑ESM1‑2‑LR) was loaded from `cmip6_historical` with `load_regular_data`: `dissic`, `talk`, `ph`, `thetao`, `so` monthly for calendar year 2002, plus `dissic`, `thetao`, `so` monthly for 1971‑01 to 2000‑12. GLODAPv2.2016b mapped climatology fields `TCO2`, `TAlk` and `pHtsinsitutp` were loaded with `load_climatology_data`. No separate `areacello`, `volcello` or depth‑level files were loaded in the executed code; the weighting was supplied internally by `compute_basic_stat(op="volume_weighted_mean")` (and by `compute_difference_metrics`/`analyze_correlation` area weighting), so all means are volume weighted in the vertical and area weighted in the horizontal.

Time aggregation used `temporal_aggregate(freq="mean", op="mean")` for both the 2002 fields and the 1971–2000 series (the long‑term mean is therefore an unweighted monthly mean over the 360 months, not an explicitly day‑weighted one). Model `dissic` and `talk` (mol m⁻³) were converted to µmol kg⁻¹ with `transform_unit`, using an in‑situ density from `compute_potential_density_anomaly(..., ref_pressure=None)` built from the matching `thetao`/`so`; the density assumption was recorded in `attrs['density_assumption']`. GLODAP `TCO2`/`TAlk` are already µmol kg⁻¹ and GLODAP `pHtsinsitutp` is total‑scale in‑situ pH, so neither was unit‑converted. Fills near 1e20 are handled as NaN by the xarray reading path and masked per field.

**2. 0–50 m global means (each source over its own valid cells)**

| Quantity | CESM2 | GFDL‑ESM4 | MPI‑ESM1‑2‑LR | GLODAP |
|---|---|---|---|---|
| DIC (µmol kg⁻¹) | 1997 | 2059 | 1974 | 2036 |
| DIC bias | −39 | +23 | −61 | — |
| TAlk (µmol kg⁻¹) | 2264 | 2336 | 2218 | 2304 |
| TAlk bias | −40 | +31 | −87 | — |
| pH (total scale) | 8.080 | 8.088 | 8.059 | 8.080 |
| pH bias | −0.001 | +0.008 | −0.022 | — |

GFDL‑ESM4 has the smallest absolute DIC and TAlk biases at 0–50 m; CESM2 has the smallest pH bias (essentially zero) and MPI‑ESM1‑2‑LR is the poorest on all three.

**3. 0–50 m DIC map comparison (GLODAP 1° grid, one shared mask of 41,227 cells out of 64,800)**

Each per‑cell 0–50 m DIC was formed as a depth volume‑weighted mean, conservatively regridded with `regrid_horizontal(method="conservative")` onto the GLODAP 1° grid, and a single common mask was built with `build_comparison_mask(ocean_only=False)`.

| Model | Mean bias (µmol kg⁻¹) | RMSE (µmol kg⁻¹) | Spatial r |
|---|---|---|---|
| CESM2 | −36 | 60 | 0.828 |
| GFDL‑ESM4 | +25 | 65 | 0.709 |
| MPI‑ESM1‑2‑LR | −59 | 83 | 0.736 |

These shared‑mask biases differ slightly from the separately taken global‑mean differences (e.g. CESM2 −36 vs −39), as expected because the two averages cover different cells and weightings. Note that GFDL‑ESM4’s “native” values arrive on a 180×360 1° archive grid with a regridded land–sea mask, so its budgets are only approximately conservative.

**4. Vertical DIC bias profile (linear `regrid_vertical` to 0, 100, 500, 1000, 2000 m, then conservative horizontal regrid; pairwise valid cells, ~33,000–41,300 per depth)**

| Depth (m) | CESM2 | GFDL‑ESM4 | MPI‑ESM1‑2‑LR |
|---|---|---|---|
| 0 | −29 | +29 | −58 |
| 100 | −46 | +19 | −59 |
| 500 | −12 | −11 | −78 |
| 1000 | +2 | −1 | −67 |
| 2000 | −2 | −11 | −61 |

CESM2 and GFDL‑ESM4 are strongly biased only in the upper water column and are within ~10 µmol kg⁻¹ below 500 m; MPI‑ESM1‑2‑LR carries a near‑uniform −60 to −80 µmol kg⁻¹ deficit through the whole column.

**5. DIC period change and its ratio to the 2002 bias (native grids, own cell measures)**

| Model | 2002 | 1971–2000 | Change | Bias | Ratio |
|---|---|---|---|---|---|
| CESM2 | 1997 | 1987 | +9 | −39 | −0.24 |
| GFDL‑ESM4 | 2059 | 2047 | +12 | +23 | +0.53 |
| MPI‑ESM1‑2‑LR | 1974 | 1964 | +11 | −61 | −0.18 |

All three models have a plausible positive anthropogenic DIC increase (~9–12 µmol kg⁻¹ over ~30 years). The ratios (via `compute_ratio`) are all well below one in magnitude, so the 2002‑minus‑1971/2000 change is only ~18–53% of the model–GLODAP bias. The bias therefore cannot be attributed to period/background drift; it remains the more likely model‑error signal. (The ratio is a diagnostic, not a causal proof.)

**6. Which model is best, and are the error signs consistent?**

A rank score summed over five error measures (|DIC| bias, |TAlk| bias, |pH| bias, map RMSE, map correlation) gives CESM2 = 7, GFDL‑ESM4 = 9, MPI‑ESM1‑2‑LR = 14, so **CESM2 reproduces the observed carbon system best overall** — highest spatial correlation (0.828), lowest RMSE (60 µmol kg⁻¹), and a near‑zero pH bias. GFDL‑ESM4 is, however, the best on the plain 0–50 m DIC and TAlk mean biases (−39/+23 vs +23/+31 … GFDL smaller), so the ranking is not unanimous; MPI‑ESM1‑2‑LR is consistently the worst.

**Signs are not consistent across the three models.** For DIC, TAlk and pH alike, CESM2 and MPI‑ESM1‑2‑LR are negative while GFDL‑ESM4 is positive. The three models also show consistent internal sign agreement between DIC and TAlk biases within each model (both negative for CESM2 and MPI, both positive for GFDL), which is the expected carbonate‑system coupling, but the inter‑model sign disagreement means the bias is not a common systematic offset of the CMIP6 generation.

**7. Figures**

- **Figure 1** shows three global maps (one panel per model) of the 0–50 m DIC bias (model − GLODAP) on the shared 1° comparison mask, using a single diverging `RdBu_r` scale centred on zero and clipped at the 98th percentile of pooled |bias|, labelled in µmol kg⁻¹; each panel title reports the shared‑mask bias, RMSE and spatial correlation, and the common‑mask cell count (41,227) is given in the figure title. The intended logic is to reveal where each model’s near‑surface DIC error is geographically concentrated (e.g. high‑latitude and upwelling regimes), given that the global means alone differ by tens of µmol kg⁻¹.
- **Figure 2** shows the vertical profile of the global mean DIC bias: depth on the vertical axis increasing downward, bias on the horizontal axis, one line per model with distinct markers, a dashed vertical zero reference, and the five sampled depths marked on the axis. Its logic is to separate surface‑intensified errors (CESM2, GFDL‑ESM4) from a depth‑uniform offset (MPI‑ESM1‑2‑LR), which the single 0–50 m number cannot show.

**8. Key caveats observed in the results**

- The 1971–2000 baseline is a plain monthly mean over the 360 loaded months rather than a day‑weighted climatology, and the period change uses 2002 (a single year) versus that 30‑year mean, so a small interannual component is folded into the “change”.
- DIC and alkalinity conversion depends on the in‑situ density field; GFDL‑ESM4’s density anomaly reaches small negative values in the very fresh surface cells, and GFDL‑ESM4’s DIC minimum (347 µmol kg⁻¹) reflects those near‑zero‑salinity archive cells, which inflate its spread and RMSE relative to its mean bias.
- The 0–50 m global means (Section 2) and the map bias (Section 3) are deliberately different estimators (own cells vs a shared mask) and should not be interchanged.
