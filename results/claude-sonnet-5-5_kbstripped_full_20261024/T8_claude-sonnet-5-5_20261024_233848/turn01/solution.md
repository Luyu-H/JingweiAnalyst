## Summary

**Answer in brief.** All three models were evaluated against GLODAPv2.2016b in the upper 50 m and down to 2000 m. None is free of bias.

- **CESM2** ranks best overall. It has the lowest map RMSE (60 µmol kg⁻¹), the highest spatial correlation (0.828) and the smallest 0–50 m pH difference.
- **GFDL-ESM4** has the smallest absolute 0–50 m global DIC and alkalinity differences, but a weaker spatial pattern (r = 0.709).
- **MPI-ESM1-2-LR** is worst on every metric. It is too low in DIC, alkalinity and pH.

No variable was missing for any model, so nothing was substituted.

## Methods as executed

- **Data.** The models are CESM2, GFDL-ESM4 and MPI-ESM1-2-LR (`cmip6_historical`). The 2002 fields were monthly `dissic`, `talk`, `ph`, `thetao` and `so`, averaged with `temporal_aggregate(freq="mean")`. For 1971–2000 the code loaded `dissic`, `thetao` and `so` and used the same aggregation. GLODAP provided `TCO2`, `TAlk` and `pHtsinsitutp`, so pH is compared on the in-situ total scale.
- **Unit conversion.** The code used `compute_potential_density_anomaly(thetao, so, ref_pressure=None)` and `transform_unit` to convert `dissic` and `talk` from mol m⁻³ to µmol kg⁻¹, cell by cell. The same conversion was applied to the 2002 and 1971–2000 periods. `ph` was not converted.
  - The recorded `density_assumption` describes the density as "In-situ density anomaly … converted at each cell's own density".
  - The reported density ranges were 1005–1031 (CESM2), 996–1045 (GFDL-ESM4) and 1004–1036 kg m⁻³ (MPI-ESM1-2-LR).
  - The GFDL-ESM4 maximum is high for the upper 50 m, which suggests some odd cells.
- **Grid provenance.**
  - `area_source` is "areacello published by cmip6_historical" for each model, and "1-degree WGS84 grid" for GLODAP.
  - `volume_source` is "volcello published by cmip6_historical" for each model.
  - Layer thickness was computed as volcello / areacello. Areas were not rescaled.
- **Global 0–50 m means.** The code called `compute_basic_stat(op="volume_weighted_mean", dim=["depth","lat","lon"])` once per source and variable. The weights were not passed explicitly in the code, so they were selected inside the tool.
- **GLODAP layers.** GLODAP's 0–50 m selection has five levels (0, 10, 20, 30, 50 m) with bounds reaching 62.5 m. The GFDL-ESM4 0–50 m selection has four levels, and the other two models have five.
- **DIC map.**
  - Each cell's mean over its own 0–50 m levels was weighted by layer thickness (valid levels only). For GLODAP the thickness came from the layer bounds, and for the models from volcello / areacello.
  - Models were regridded to the GLODAP grid with `regrid_horizontal(method="conservative")`.
  - One common mask was built with `build_comparison_mask(..., ocean_only=False)` over the three models and GLODAP.
  - Metrics came from `compute_difference_metrics` and from `analyze_correlation(weighting="area")`.
- **Profile.** Models and GLODAP were interpolated to 0, 100, 500, 1000 and 2000 m with `regrid_vertical`. Models were then regridded horizontally to the GLODAP grid. A per-model mask (model and GLODAP both valid) was built at each depth, and `compute_difference_metrics` gave the area-weighted mean bias.
- **Period change.** The change is the model's 2002 0–50 m volume-weighted mean DIC minus its own 1971–2000 volume-weighted mean. The code differenced the two aggregated means rather than calling `compute_anomaly_vs_own_climatology`. The ratio is this change divided by the model's 0–50 m DIC difference from GLODAP.

## Results

### Global 0–50 m means (one mean over all of each source's cells)

| | DIC (µmol kg⁻¹) | TAlk (µmol kg⁻¹) | pH (total) |
|---|---|---|---|
| GLODAP | 2036 | 2304 | 8.080 |
| CESM2 | 1997 (−39) | 2264 (−40) | 8.080 (−0.001) |
| GFDL-ESM4 | 2059 (+23) | 2336 (+31) | 8.088 (+0.008) |
| MPI-ESM1-2-LR | 1974 (−61) | 2218 (−87) | 8.059 (−0.022) |

Values in parentheses are model minus GLODAP. The CESM2 pH difference is −0.0005 before rounding.

### 0–50 m DIC map on the GLODAP 1° grid

- **Common mask.** There are 41,227 common cells, out of 41,977 GLODAP cells. The input fields had 42,576 (CESM2), 44,564 (GFDL-ESM4), 42,713 (MPI-ESM1-2-LR) and 41,977 (GLODAP) valid cells.
- **GLODAP mean on this sample.** 2034.0 µmol kg⁻¹.

| Model | Mean bias (µmol kg⁻¹) | RMSE (µmol kg⁻¹) | Spatial r |
|---|---|---|---|
| CESM2 | −37 | 60 | 0.828 |
| GFDL-ESM4 | +25 | 65 | 0.709 |
| MPI-ESM1-2-LR | −59 | 83 | 0.736 |

The model means on this sample were 1997.5, 2059.3 and 1974.7 µmol kg⁻¹. The bias equals the model mean minus the GLODAP mean over the same cells. Some very large single-cell differences exist (about −1770 µmol kg⁻¹ in GFDL-ESM4 and −1729 µmol kg⁻¹ in MPI-ESM1-2-LR), probably coastal or regridding artefacts. They inflate the RMSE.

### Global mean DIC bias profile (µmol kg⁻¹)

| Depth (m) | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---|---|---|---|
| 0 | −29 | +29 | −58 |
| 100 | −46 | +19 | −59 |
| 500 | −12 | −11 | −78 |
| 1000 | +2 | −1 | −67 |
| 2000 | −2 | −11 | −61 |

The number of common cells shrinks with depth (CESM2 40,921 → 33,041; GFDL-ESM4 41,301 → 33,795; MPI-ESM1-2-LR 40,828 → 33,109, from 0 to 2000 m).

### Period change of global 0–50 m DIC (2002 minus the model's own 1971–2000 mean)

| Model | Change (µmol kg⁻¹) | Ratio to GLODAP difference |
|---|---|---|
| CESM2 | +9 | −0.24 |
| GFDL-ESM4 | +12 | 0.53 |
| MPI-ESM1-2-LR | +11 | −0.18 |

**Interpretation.**
- The baseline is each model's own climatology, so this change does not remove the model bias.
- For CESM2 and MPI-ESM1-2-LR the change is positive while the bias is negative, which gives a negative ratio. Its magnitude is small (|ratio| ≤ 0.24), so a period or drift effect cannot account for their offsets. These offsets are mean-state model errors.
- For GFDL-ESM4 the change is about half the bias and has the same sign. A time-offset effect cannot be fully excluded there, but most of its +23 µmol kg⁻¹ offset is still not explained by it.

### Synthesis

- **Ranking.** The mean rank over five criteria (|ΔDIC|, |ΔTAlk|, |ΔpH|, map RMSE, map correlation) was CESM2 1.4, GFDL-ESM4 1.8 and MPI-ESM1-2-LR 2.8, so CESM2 is best.
  - GFDL-ESM4 wins on the global DIC and alkalinity differences. CESM2 wins on pH, RMSE and correlation. With only three models and equal weights, the CESM2–GFDL-ESM4 order is not very robust.
  - MPI-ESM1-2-LR is last on every metric.
- **Sign consistency.** Within each model, the DIC, alkalinity and pH differences share a sign: CESM2 and MPI-ESM1-2-LR are all negative, and GFDL-ESM4 is all positive. Across models the signs are not consistent for any of the three variables.

## Figures

- **Figure 1.** Three side-by-side global maps (one per model) of the 0–50 m DIC bias (model 2002 minus GLODAP) on the common-mask cells. All panels share one diverging `RdBu_r` scale centred on zero. Its limits are ± the 98th percentile of the pooled |bias|, with extended end arrows, and the colourbar is labelled "DIC bias (µmol kg-1)". Each panel title gives that model's mean bias and RMSE.
- **Figure 2.** A single panel of the global mean DIC bias profile. Depth is on the vertical axis, increasing downward. Bias in µmol kg⁻¹ is on the horizontal axis. There is one line per model with markers at 0, 100, 500, 1000 and 2000 m, a vertical zero line and light horizontal guides at the sampled depths.

## Caveats

- **Averaging period.** Models are a single year (2002 annual mean), while GLODAP is a climatology built from 1972–2013 data.
- **Carbonate chemistry.** The CESM2 carbonate constants are undocumented, which is an unverified assumption. The models also differ in their alkalinity terms, since MPI-ESM1.2 added phosphate and silicate alkalinity. Alkalinity and pH differences therefore partly reflect these definitions and not only the physical or biological state.
- **Ocean areas.** Each model's integrals cover its own ocean area, which differs by about 1 % between models, and the areas were not rescaled.
- **Mask changes with depth.** Coverage shrinks with depth, so the cell counts, and hence the sample, differ between depths.
- **Odd values.** Some GFDL-ESM4 and MPI-ESM1-2-LR minima in the shallow DIC fields (about 310–350 µmol kg⁻¹) and the high GFDL-ESM4 density value suggest a few odd cells that were not screened out. They have little effect on the volume-weighted means but a larger effect on RMSE.
