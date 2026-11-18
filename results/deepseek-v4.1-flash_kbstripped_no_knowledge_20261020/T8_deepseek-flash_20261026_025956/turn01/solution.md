## 1. What was done (as executed)

**Reference.** GLODAPv2.2016b mapped climatology (`TCO2`, `TAlk`, `pHtsinsitutp`) was loaded with `load_climatology_data` on its native 33-level / 1°×1° grid (depth levels 0, 10, 20, 30, 50, 75, …, 5500 m). Using the layer bounds from `load_depth_levels`, only levels whose **full** bounds lie inside 0–50 m were retained, i.e. the 0, 10, 20 and 30 m levels (the 50 m level, whose bounds extend to 75 m, was excluded). No partial-layer weighting was applied.

**Models.** CESM2, GFDL-ESM4 and MPI-ESM1-2-LR `dissic`, `talk`, `ph`, `thetao`, `so` were loaded from `cmip6_historical` for 2002-01…2002-12 (`talk` and `ph` with `depth_range=[0,100]`), and `dissic`, `thetao`, `so` for 1971-01…2000-12 with `depth_range=[0,100]`. Annual/period means were formed with `temporal_aggregate(freq="mean", op="mean")`. Model DIC and TAlk (mol m⁻³) were converted to µmol kg⁻¹ by dividing by a density from `compute_potential_density_anomaly(thetao, so, ref_pressure=0)` (the anomaly as returned) via `transform_unit(from_unit="mol m-3", to_unit="umol kg-1", density=…)`; `ph` was used as provided. Each model's own levels wholly inside 0–50 m were selected with its `load_depth_levels` bounds, and the scalar means were reduced with `compute_basic_stat(op="volume_weighted_mean", dim=None)` directly on the selected fields (no separate `areacello`/`volcello`/`load_areas` operators were passed in the executed code). The 2-D 0–50 m DIC field is each cell's `op="mean", dim="depth"` over its own selected levels.

**Comparisons.** Maps: `regrid_horizontal(target="glodap", method="conservative")`, common mask from `build_comparison_mask(ocean_only=True)`, metrics from `compute_difference_metrics(dim=["lat","lon"], weighting="area")` and `analyze_correlation(method="pearson", weighting="area")`. Vertical: `regrid_horizontal` then `regrid_vertical(target=[0,100,500,1000,2000], method="linear")`, masked per depth. Period change: 2002 mean minus 1971–2000 mean, same level selection and weighting.

## 2. Key numerical results

**0–50 m global means (model − GLODAP)**

| Source | DIC (µmol kg⁻¹) | Δ | TAlk (µmol kg⁻¹) | Δ | pH (total) | Δ |
|---|---|---|---|---|---|---|
| GLODAP | 2029 | — | 2302 | — | 8.085 | — |
| CESM2 | 1997 | −32 | 2264 | −38 | 8.080 | −0.005 |
| GFDL-ESM4 | 2059 | +31 | 2336 | +34 | 8.088 | +0.004 |
| MPI-ESM1-2-LR | 1971 | −58 | 2216 | −85 | 8.061 | −0.024 |

**0–50 m DIC map on the GLODAP 1° grid**

| Model | Bias | RMSE | Spatial r |
|---|---|---|---|
| CESM2 | −28 | 53 | 0.836 |
| GFDL-ESM4 | +33 | 64 | 0.741 |
| MPI-ESM1-2-LR | −53 | 75 | 0.773 |

**Vertical DIC bias (µmol kg⁻¹)**

| Model | 0 m | 100 m | 500 m | 1000 m | 2000 m |
|---|---|---|---|---|---|
| CESM2 | −28 | −45 | −8 | +12 | +18 |
| GFDL-ESM4 | +31 | +20 | −6 | +9 | +9 |
| MPI-ESM1-2-LR | −57 | −58 | −73 | −57 | −41 |

**Period change (2002 − 1971–2000) and ratio to the 2002 bias**

| Model | ΔDIC (µmol kg⁻¹) | Ratio |
|---|---|---|
| CESM2 | +9 | −0.29 |
| GFDL-ESM4 | +12 | +0.41 |
| MPI-ESM1-2-LR | +11 | −0.19 |

All three models take up carbon over the ~30 years, at 9–12 µmol kg⁻¹. Because all |ratios| are below 1, the recent change is smaller than the bias, so the bias cannot be attributed *solely* to the temporal/climatological mismatch; but ratios of 0.19–0.41 mean that a non-negligible fraction (roughly a fifth to two-fifths) of the bias magnitude is of the same order as the model's own drift between the 1971–2000 average and 2002. The ratio is therefore a warning against reading the model−GLODAP offset as pure structural model error, especially for GFDL-ESM4 (largest ratio, and the only positive bias, so change and bias share sign).

## 3. Figures

*Figure 1* — three panels in one row, one per model, showing the 0–50 m DIC bias (model on the GLODAP 1° grid minus GLODAP, restricted to the common ocean mask) as pcolormesh maps, all on one shared symmetric diverging scale (`RdBu_r`) clipped at the 98th percentile of |bias| pooled across the three panels, with a single horizontal colourbar labelled µmol kg⁻¹. This construction makes the sign and amplitude of the bias directly comparable between models rather than rescaling each panel independently.
*Figure 2* — one panel with depth on the vertical axis increasing downward (0→2000 m), DIC bias on the horizontal axis, one marker-and-line series per model, a dashed vertical zero reference line, grid and legend, with ticks at the five sampled depths (0, 100, 500, 1000, 2000 m).

## 4. Interpretation

**Best model.** By the combined rank (|ΔDIC|, |ΔTAlk|, |ΔpH|, map RMSE, vertical RMS bias) GFDL-ESM4 is best (total rank 6), CESM2 second (9), MPI-ESM1-2-LR clearly worst (15). GFDL-ESM4 wins on the scalar 0–50 m differences and on the vertical profile (RMS 18 vs 26 and 58 µmol kg⁻¹), while CESM2 has the better spatial pattern skill (highest correlation 0.836, lowest RMSE 53, but a larger, sign-reversing vertical offset: −45 µmol kg⁻¹ at 100 m, +18 at 2000 m). MPI-ESM1-2-LR is biased low at every depth (−41 to −73 µmol kg⁻¹) and in both scalars.

**Sign consistency.** Across models the signs are *not* consistent for any variable: DIC, TAlk and pH are all negative for CESM2 and MPI-ESM1-2-LR but all positive for GFDL-ESM4. *Within* each model, however, the three variables always carry the same sign (CESM2 −/−/−, GFDL +/+/+, MPI −/−/−), as expected from carbonate-system coupling: a model high in DIC is also high in TAlk and pH. This means the model spread is mainly an offset in overall carbonate content, not a redistribution between DIC, alkalinity and pH.

**Caveats visible in the output.** The unweighted/regridded map means differ from the volume-weighted scalars (e.g. MPI-ESM1-2-LR: map mean 1999 vs scalar 1971 µmol kg⁻¹), so map bias and scalar bias are not interchangeable. The raw bias fields contain extreme outliers (minima near −1800 µmol kg⁻¹ for GFDL-ESM4 and MPI-ESM1-2-LR, with 0–50 m DIC values as low as ~320–350 µmol kg⁻¹), arising in marginal/coastal cells of the conservative regridding; this is why the final figure uses a 98th-percentile symmetric clip. Finally, GLODAP is a 1972–2013 climatology compared against a 2002 annual mean and a 1971–2000 model average, so part of every difference is temporal sampling mismatch, and the pH comparison presupposes that the CMIP6 `ph` is on the total scale, which the code took as given.
