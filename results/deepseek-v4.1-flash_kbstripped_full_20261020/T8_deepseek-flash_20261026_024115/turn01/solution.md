## Methodology (as actually run)

CMIP6 `historical` monthly fields were loaded on each model’s native grid with `load_regular_data(dataset="cmip6_historical", field="total")` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. For 2002 (Jan–Dec), `dissic`, `talk`, `ph`, `thetao` and `so` were loaded; `talk` and `ph` were requested for 0–50 m. For 1971–2000, `dissic`, `thetao` and `so` were loaded with depth range 0–50 m.

DIC and total alkalinity were converted from `mol m-3` to `umol kg-1` using `compute_potential_density_anomaly(thetao, so, ref_pressure=None)` and then `transform_unit(..., from_unit="mol m-3", to_unit="umol kg-1", density=<that density>)`. This conversion was done before the 2002 and 1971–2000 temporal aggregation. pH was not converted; it was averaged as supplied. Annual/period means were computed with `temporal_aggregate(freq="mean", op="mean")`.

GLODAPv2.2016b mapped fields were loaded with `load_climatology_data(dataset="glodap", ...)`: `TCO2` (then subset to 0–50 m), `TAlk` (0–50 m) and `pHtsinsitutp` (0–50 m, in-situ pH on the total scale). The GLODAP TCO2 profile was produced with `regrid_vertical(..., target=[0,100,500,1000,2000], method="linear")`.

The 0–50 m means were computed with `compute_basic_stat(..., op="volume_weighted_mean", dim=None)` on each product’s own grid. The shown code did not separately load or pass `areacello`/`volcello` arrays; it requested the `volume_weighted_mean` option. No unweighted cell-mean alternative was computed.

For DIC maps, each cell’s 0–50 m mean was computed with `compute_basic_stat(op="mean", dim="depth")`; model maps were regridded conservatively to the GLODAP 1° grid with `regrid_horizontal(target="glodap", method="conservative")`. A common comparison mask was built with `build_comparison_mask([...], ocean_only=False)` and applied before `compute_difference_metrics(..., weighting="area")` and `analyze_correlation(..., weighting="area")`.

For the vertical profile, annual-mean 2002 DIC was interpolated to 0, 100, 500, 1000 and 2000 m with `regrid_vertical(method="linear")`, regridded horizontally conservatively to the GLODAP 1° grid, masked per depth with `build_comparison_mask(..., ocean_only=False)`, and averaged with area weighting via `compute_difference_metrics`.

The 2002-minus-1971–2000 change was computed on each model’s own grid from the 2002 mean minus the model’s 1971–2000 mean, using `compute_basic_stat(op="volume_weighted_mean")`.

## Key results

### 0–50 m means and differences from GLODAP

| Quantity | GLODAP | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---|---:|---:|---:|---:|
| DIC mean (µmol kg⁻¹) | 2036 | 1997 | 2059 | 1974 |
| DIC difference | — | −39 | +24 | −61 |
| TAlk mean (µmol kg⁻¹) | 2304 | 2264 | 2336 | 2218 |
| TAlk difference | — | −40 | +31 | −87 |
| pH mean | 8.080 | 8.080 | 8.088 | 8.059 |
| pH difference | — | −0.0005 (≈−0.001 to 3 d.p.) | +0.008 | −0.022 |

### 0–50 m DIC map comparison on the GLODAP 1° grid

Model minus GLODAP, area-weighted on the common mask:

| Model | Mean bias (µmol kg⁻¹) | RMSE (µmol kg⁻¹) | Spatial correlation |
|---|---:|---:|---:|
| CESM2 | −33 | 60 | 0.817 |
| GFDL-ESM4 | +27 | 75 | 0.662 |
| MPI-ESM1-2-LR | −55 | 81 | 0.730 |

### Vertical profile of global mean DIC bias

Model minus GLODAP, area-weighted, with valid comparison cell counts:

| Depth (m) | CESM2 bias | CESM2 n | GFDL bias | GFDL n | MPI bias | MPI n |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | −29 | 40921 | +29 | 41301 | −58 | 40828 |
| 100 | −46 | 38318 | +19 | 38534 | −59 | 38355 |
| 500 | −12 | 36370 | −11 | 36561 | −78 | 36294 |
| 1000 | +2 | 35727 | −1 | 35975 | −67 | 35688 |
| 2000 | −2 | 33041 | −11 | 33795 | −61 | 33109 |

Coverage decreases with depth for all models, from about 41k cells near the surface to about 33k at 2000 m.

### 2002 minus 1971–2000 0–50 m DIC change

| Model | DIC change (µmol kg⁻¹) | Ratio to 0–50 m DIC difference from GLODAP |
|---|---:|---:|
| CESM2 | +9 | −0.24 |
| GFDL-ESM4 | +12 | +0.53 |
| MPI-ESM1-2-LR | +11 | −0.18 |

These changes are computed with each model’s own 1971–2000 baseline, so the model’s mean-state DIC bias cancels. The ratios therefore should not be read as model errors. The 2002-minus-1971–2000 change also comes from a single coupled `historical` member; internal variability and the fact that a model year need not align with the observed year limit interpretation of the period change as forced change.

## Ranking and sign consistency

By the specified ranking key—map RMSE first, then depth-profile MAE, then absolute bias, then correlation—the order is:

1. **CESM2**: RMSE 60, profile MAE 18, bias −33, r 0.817  
2. **GFDL-ESM4**: RMSE 75, profile MAE 14, bias +27, r 0.662  
3. **MPI-ESM1-2-LR**: RMSE 81, profile MAE 65, bias −55, r 0.730  

**CESM2 reproduces the observed carbon system best overall** under this ranking, mainly because it has the lowest DIC map RMSE and highest spatial correlation. GFDL-ESM4 has a smaller depth-profile MAE and a smaller absolute 0–50 m DIC mean difference, but its map RMSE and correlation are poorer. MPI-ESM1-2-LR has the largest biases and RMSE.

Sign consistency of model-minus-GLODAP errors:

- **CESM2**: DIC −39, TAlk −40, pH −0.001 → all negative, consistent.
- **GFDL-ESM4**: DIC +24, TAlk +31, pH +0.008 → all positive, consistent.
- **MPI-ESM1-2-LR**: DIC −61, TAlk −87, pH −0.022 → all negative, consistent.

Within each model, the DIC, TAlk and pH errors agree in sign. Across models, GFDL-ESM4 is positive whereas CESM2 and MPI-ESM1-2-LR are negative. MPI-ESM1-2-LR is persistently DIC-depleted over the full 0–2000 m profile, CESM2 has its largest negative DIC bias near 100 m and near-zero deep bias, and GFDL-ESM4 changes sign with depth, with positive near-surface bias and modest negative bias at 500 m and 2000 m.

## Figures generated

**Figure 1** shows three global panels, one per model, of the 0–50 m DIC bias (model minus GLODAP) on the masked GLODAP 1° comparison cells. All panels share one diverging `RdBu_r` colour scale centred on zero, labelled in µmol kg⁻¹. The shared scale is clipped at the 98th percentile of the pooled absolute bias to keep the maps legible.

**Figure 2** shows the vertical profile of global mean DIC bias. Depth is on the vertical axis increasing downward, DIC bias in µmol kg⁻¹ is on the horizontal axis, one line is drawn per model, a vertical zero reference line is included, and markers are placed at the five sampled depths: 0, 100, 500, 1000 and 2000 m.

## Caveats

- The model means are a single 2002 annual mean from `historical`, not an ensemble mean.
- The 0–50 m means are produced by the `volume_weighted_mean` option; the shown code did not compute an unweighted cell-mean alternative.
- GLODAP is a 1972–2013 mapped climatology; the model 2002 annual mean and 1971–2000 baseline do not represent the same observational sampling.
- The DIC, TAlk, pH and depth-profile differences from GLODAP carry the model’s mean-state bias. The 2002-minus-1971–2000 DIC change does not, because it uses the model’s own baseline.
