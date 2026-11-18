## 1. Direct answer

In the 200–1000 m layer, WOA23 has a volume-weighted mean N:P of 14.15. That is below the Redfield value of 16, the usual intermediate-water signature of denitrification. The models differ in how well they reproduce it.

- **GFDL-ESM4** matches the observed mean (14.13) and the global low-N:P volume (7.3 % vs 6.2 %). It puts that volume mostly in the Pacific and almost none in the Indian Ocean (0.5 % vs 10.3 %).
- **CESM2 and MPI-ESM1-2-LR** are biased low everywhere, Atlantic included (headline N:P 11.83 and 12.28). They have far too much water with N:P < 12 (36.7 % and 22.5 % vs 6.2 %). Their low N:P is a basin-wide offset, not a localized denitrification signal.

The calculation does not test mechanism. It cannot show whether a model's low N:P comes from denitrification or from other biogeochemical biases.

## 2. Methods as executed

**Data**
- Model `no3` and `po4` came from `load_regular_data` (cmip6_historical, monthly, 1971-01 to 2000-12, `depth_range=[200,1000]`).
- Each model was averaged with `temporal_aggregate(freq='mean')` into a 1971–2000 climatology.
- WOA23 `nitrate` and `phosphate` came from `load_climatology_data`. This is the all-period annual composite, so the periods are not identical to the models'.
- Each model's own `volcello` came from `load_volumes_for_cmip6`. `attrs['volume_source']` was "volcello published by cmip6_historical" for all three models.
- WOA23 has no volcello. The code built layer thickness from the WOA23 depth-level bounds, and the cell areas come from the tool.
- Native grids:
  - CESM2: `gn`, 384×320, 20 levels.
  - GFDL-ESM4: `gr`, the 1° archive grid, 180×360, 10 levels.
  - MPI-ESM1-2-LR: `gn`, 220×254, 11 levels.
  - WOA23: 1°, 23 levels.
- Original levels were used with no interpolation. The level bounds are not clipped exactly to 200–1000 m, and they differ by source:

| Source | Layer top (m) | Layer bottom (m) |
|---|---|---|
| CESM2 | 203.7 | 1045.5 |
| GFDL-ESM4 | 175 | 1050 |
| MPI-ESM1-2-LR | 200 | 1020 |
| WOA23 | 187.5 | 1025 |

- In every source, PO4 was set to NaN where it was ≤ 0, and NO3 and PO4 were forced to share one valid sample.

**Native-grid statistics** (each source on its own grid)
- N:P came from `compute_ratio(NO3, PO4)` cell by cell.
- The headline value is `compute_basic_stat(..., 'volume_weighted_mean', dim=['depth','lat','lon'])` applied to that ratio field.
- The secondary value is the ratio of the volume-weighted mean NO3 to the volume-weighted mean PO4.
- The volume weights were not passed explicitly in the code. `compute_basic_stat` and `compute_exceedance_fraction` resolved them internally from the datasets. The total volumes returned (about 2.7–2.9×10¹⁷ m³) are consistent with that.
- Volume shares used `compute_exceedance_fraction`:
  - "14 < N:P < 18" was `comparison='abs_less'`, `threshold=2`, `reference=16`, i.e. |N:P − 16| < 2 strictly.
  - "N:P < 12" was `comparison='less'`.
- Basins came from `select_basin` with the WOA23 basin mask, `include_marginal_seas=True` and `depth_dependent=False`. The Pacific, Atlantic and Indian Oceans were used. The Southern and Arctic Oceans were excluded.

**Regridding, mask and regression**
- Models were regridded to the WOA23 1° grid with `regrid_horizontal(method='conservative')`, applied to the 3-D NO3 and PO4 fields and to a thickness-weighted column-mean N:P field.
- `build_comparison_mask(ocean_only=True)` was built from the horizontal validity (any depth) of all four sources.
  - Cells valid per source: WOA23 38,759; CESM2 39,475; GFDL-ESM4 41,268; MPI 39,586.
  - Shared mask: **38,167 columns**. The ocean/land step removed 0 cells.
- Model concentrations were converted from mol m⁻³ to µmol kg⁻¹ with `transform_unit`, using ρ = 1025 kg m⁻³. WOA23 was already in µmol kg⁻¹.
- The regression NO3 = a·PO4 + b used `fit_multivariate_relation(method='linear')` on the regridded 3-D fields under the shared 2-D mask.
  - Every valid (column, level) cell counts equally.
  - Because the sources have different numbers of levels, the sample sizes differ: WOA23 856,031; CESM2 749,081; GFDL-ESM4 376,663; MPI 409,898.
  - A numpy `polyfit` cross-check gave the same coefficients.
- This regression is equal-weight and therefore not directly comparable with the volume-weighted statistics.

**Comparison with WOA23**
- `compute_difference_metrics` compared the regridded column-mean N:P with WOA23's on the shared mask (38,167 cells, equal weights).
- The Redfield judgement used a tolerance of ±0.5 around 16 on the headline value.
- The "reproduces" flag was a criterion I chose: the model's low-N:P share within 0.5–2× the WOA23 share.

## 3. Results

### Mean N:P, 200–1000 m (native grids, volume-weighted)

| Source | Headline: mean of cell-by-cell N:P | Secondary: mean NO3 / mean PO4 | Vol-wtd mean NO3 | Vol-wtd mean PO4 |
|---|---|---|---|---|
| WOA23 | **14.15** | 14.00 | 27.49 µmol/kg | 1.96 µmol/kg |
| CESM2 | **11.83** | 10.93 | 21.84 mmol/m³ | 2.00 mmol/m³ |
| GFDL-ESM4 | **14.13** | 13.65 | 25.08 mmol/m³ | 1.84 mmol/m³ |
| MPI-ESM1-2-LR | **12.28** | 11.77 | 20.89 mmol/m³ | 1.77 mmol/m³ |

(The model means are in native per-volume units, mol m⁻³ × 1000.)

### Regression NO3 = a·PO4 + b (shared mask, equal-weight cells, µmol kg⁻¹, models at ρ = 1025 kg m⁻³)

| Source | a | b (µmol/kg) | R² |
|---|---|---|---|
| WOA23 | 13.05 | 1.67 | 0.910 |
| CESM2 | 6.18 | 9.01 | 0.369 |
| GFDL-ESM4 | 11.38 | 4.05 | 0.820 |
| MPI-ESM1-2-LR | 8.09 | 6.22 | 0.643 |

### Volume with 14 < N:P < 18 (native, by volume)

| Source | Share | Volume (m³) | Total layer volume (m³) |
|---|---|---|---|
| WOA23 | 55.0 % | 1.50×10¹⁷ | 2.72×10¹⁷ |
| CESM2 | 17.2 % | 4.75×10¹⁶ | 2.77×10¹⁷ |
| GFDL-ESM4 | 74.6 % | 2.14×10¹⁷ | 2.87×10¹⁷ |
| MPI-ESM1-2-LR | 0.5 % | 1.49×10¹⁵ | 2.72×10¹⁷ |

### Basin headline N:P (native, volume-weighted; marginal seas included, depth_dependent=False)

| Source | Pacific | Atlantic | Indian |
|---|---|---|---|
| WOA23 | 13.57 | 15.55 | 13.73 |
| CESM2 | 10.58 | 13.55 | 11.18 |
| GFDL-ESM4 | 13.39 | 15.07 | 14.72 |
| MPI-ESM1-2-LR | 11.29 | 13.15 | 12.60 |

### Volume with N:P < 12 (percentage of the region's volume, with volume in m³)

| Source | Global | Pacific | Indian |
|---|---|---|---|
| WOA23 | 6.2 % (1.69×10¹⁶) | 8.1 % (9.97×10¹⁵) | 10.3 % (4.88×10¹⁵) |
| CESM2 | 36.7 % (1.02×10¹⁷) | 59.4 % (7.38×10¹⁶) | 35.8 % (1.72×10¹⁶) |
| GFDL-ESM4 | 7.3 % (2.08×10¹⁶) | 14.9 % (1.93×10¹⁶) | 0.5 % (2.6×10¹⁴) |
| MPI-ESM1-2-LR | 22.5 % (6.13×10¹⁶) | 40.0 % (4.88×10¹⁶) | 15.0 % (7.03×10¹⁵) |

### Difference from WOA23 (column-mean N:P, shared mask, 38,167 cells)

| Model | Bias | RMSE |
|---|---|---|
| CESM2 | −2.29 | 3.59 |
| GFDL-ESM4 | −0.02 | 1.69 |
| MPI-ESM1-2-LR | −1.84 | 2.72 |

## 4. Interpretation

**Redfield.**
- All four sources have headline N:P below 16, outside the ±0.5 tolerance. WOA23 is at 14.15 (about 2 below), GFDL-ESM4 at 14.13, MPI at 12.28 and CESM2 at 11.83.
- CESM2 and MPI are therefore more strongly sub-Redfield than observed. GFDL-ESM4 matches the observed mean, with a shared-mask bias of about 0.

**Regression.**
- WOA23 has the tightest and steepest NO3–PO4 relation (a = 13.05, R² = 0.91).
- GFDL-ESM4 is closest to it (a = 11.38, R² = 0.82).
- MPI and especially CESM2 have flatter slopes, larger positive intercepts and much weaker correlation (R² = 0.64 and 0.37). Their NO3 is less tightly tied to PO4.
- Part of CESM2's weak fit may come from slightly negative NO3 values in the model field (minimum about −0.0013 mol m⁻³).

**Low-N:P (denitrification) signal.**
- WOA23 has lower N:P in the Pacific (13.57) and Indian (13.73) than in the Atlantic (15.55). About 8 % of Pacific and 10 % of Indian volume has N:P < 12.
- All three models reproduce the sign of the Pacific–Atlantic contrast. Pacific minus Atlantic is −1.99 in WOA23, −2.97 in CESM2, −1.68 in GFDL-ESM4 and −1.86 in MPI.
- **CESM2** overestimates the low-N:P volume by about 5.9× globally, 7.4× in the Pacific and 3.5× in the Indian Ocean.
- **MPI** overestimates it by about 3.6× globally, 5.0× in the Pacific and 1.5× in the Indian Ocean.
- **GFDL-ESM4** reproduces the global volume (1.17× WOA23) but overestimates the Pacific share (1.85×) and nearly misses the Indian low-N:P water (0.05×, N:P mean 14.72 vs 13.73).
- By my 0.5–2× criterion, only GFDL-ESM4 reproduces the global volume, and no model reproduces the basin distribution.

## 5. Figures

**Figure 1** is a 2×2 grid of scatter plots in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR.
- Each panel plots NO3 (y) against PO4 (x) in µmol kg⁻¹ for the 200–1000 m cells on the shared 1° mask.
- Only a random subsample of up to 8,000 cells per source is drawn. The fits and statistics use all cells.
- All panels share the same axes and ranges, with limits set from the pooled 99.5th percentiles.
- Each panel shows its own OLS line (solid red) and the Redfield 16:1 line through the origin (dashed black), with a text box giving a, b and R².

**Figure 2** is a 2×2 grid of maps of the thickness-weighted 200–1000 m column-mean N:P minus 16.
- The fields are the conservatively regridded 1° versions on the shared mask, used for display only.
- They are drawn with `pcolormesh` on a regular longitude–latitude frame. This is plain lon/lat axes with equal aspect, not a cartographic projection.
- All panels share one RdBu_r diverging scale centred on zero. It is clipped at the 98th percentile of the pooled absolute values and has extend arrows.
- Panel titles give each source's headline N:P.

## 6. Caveats

- Native-grid numbers (means, volume shares, basin means) use each source's own grid and weights. GFDL-ESM4's "native" grid is its 1° `gr` archive grid with 10 levels, which is coarser vertically than the others. CESM2 and MPI are `gn`.
- Shared-mask numbers (regression, scatter, maps, bias/RMSE) use an equal-weight sample of 38,167 columns and differ from each source's own-ocean volume-weighted values. For example, the column-map mean is 12.29 for CESM2, against 11.83 for the headline value.
- The levels and layer bounds differ between sources, so the 200–1000 m layers are not strictly identical.
- The models cover 1971–2000, whereas WOA23 is an all-period composite.
- The WOA23 field contains some cells with zero nitrate and a few extreme ratios (up to about 139). These cells enter the cell-by-cell mean and the regression.
- The intercepts depend on the constant-density conversion (ρ = 1025 kg m⁻³). N:P, slope and R² do not.
