**Direct answer:** For the 200–1000 m layer, all four sources sit **below Redfield 16:1** in their headline mean cell-by-cell N:P. WOA23 has a mean N:P of 14.15 and a modest low-N:P volume: 6.2% of global water volume has N:P < 12, rising to 8.1% in the Pacific and 10.3% in the Indian Ocean. The three CMIP6 models do **not** reproduce both the observed low-N:P volume and its basin distribution. CESM2 and MPI-ESM1-2-LR strongly overproduce low-N:P water; GFDL-ESM4 is close to WOA23 in global mean and global N:P < 12 fraction, but places the low-N:P signal incorrectly among basins and has too much volume in the 14–18 range.

**Methodology actually executed**
- WOA23 annual objectively analysed nitrate and phosphate were loaded with `load_climatology_data(dataset="woa23", variable=..., time_slice="annual", depth_range=[200, 1000])`; WOA23 values were already in micromoles per kilogram. No temporal aggregation was applied to WOA23.
- CMIP6 historical monthly `no3` and `po4` were loaded with `load_regular_data(dataset="cmip6_historical", variable=..., model=..., time_range=["1971-01", "2000-12"], depth_range=[200, 1000])` for CESM2, GFDL-ESM4, and MPI-ESM1-2-LR. Each was temporally aggregated with `temporal_aggregate(..., freq="mean", op="mean")`. The executed code did not pass an additional day-weighting argument to this call.
- CMIP6 concentrations were converted with `transform_unit(..., to_unit="umol kg-1", density=1.025)`; the code comment equates this density argument to 1025 kg m−3. WOA23 was not converted.
- The valid-cell mask was defined by non-null nitrate and phosphate after loading/conversion: `no3.notnull() & po4.notnull()`. The executed code did not explicitly call `build_comparison_mask`, `load_depth_levels`, `load_basin_mask`, `load_areas`, or `load_volumes_for_cmip6`.
- Cell-by-cell N:P was computed with `compute_ratio`. The headline N:P is `compute_basic_stat(ratio, op="volume_weighted_mean", dim=None)`. The second N:P is `compute_ratio(mean_NO3, mean_PO4)`, where the means were also volume-weighted with `compute_basic_stat(..., op="volume_weighted_mean", dim=None)`.
- OLS was performed over all valid cells with equal cell weighting using `np.polyfit(PO4, NO3, 1)` for slope and intercept, and direct OLS for R². `fit_multivariate_relation` was also called with `target=no3`, `predictors=[po4]`, `method="linear"`, `validation="random"`, `test_size=0.2`, `standardize=False` for diagnostics; the reported R² was taken from that result when available.
- Volume fractions were computed with `compute_exceedance_fraction(..., weighting="volume")`: 14 < N:P < 18 via `threshold=2.0, comparison="abs_less", reference=16.0`; N:P < 12 via `threshold=12.0, comparison="less", reference=0.0`.
- Basin statistics used `select_basin(..., basin="Pacific"/"Atlantic"/"Indian", include_marginal_seas=True)`. Southern and Arctic basins were not selected for basin means. N:P < 12 was also computed for the selected Pacific and Indian basins.
- Maps were made with `compute_basic_stat(ratio, op="volume_weighted_mean", dim="depth")`, followed by subtracting 16. CMIP6 maps were displayed after `regrid_horizontal(..., target="woa23", method="bilinear")` for display only.

**Key numerical results (all for 200–1000 m; percentages are volume percentages)**

| Source | Headline mean N:P (mean of cell-by-cell N:P) | Ratio of mean NO3 to mean PO4 | OLS slope a | OLS intercept b (µmol kg⁻¹) | R² | 14 < N:P < 18 (%) | N:P < 12 global (%) | N:P < 12 Pacific (%) | N:P < 12 Indian (%) | Pacific mean N:P | Atlantic mean N:P | Indian mean N:P | vs Redfield 16 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| WOA23 | 14.15 | 14.00 | 13.04 | 1.69 | 0.908 | 55.0 | 6.2 | 8.1 | 10.3 | 13.57 | 15.55 | 13.73 | below |
| CESM2 | 11.83 | 10.93 | 4.89 | 10.84 | 0.238 | 17.2 | 36.7 | 59.4 | 35.8 | 10.58 | 13.55 | 11.18 | below |
| GFDL-ESM4 | 14.13 | 13.65 | 11.37 | 4.02 | 0.817 | 74.6 | 7.3 | 14.9 | 0.5 | 13.39 | 15.07 | 14.72 | below |
| MPI-ESM1-2-LR | 12.28 | 11.77 | 8.59 | 5.64 | 0.678 | 0.5 | 22.5 | 40.0 | 15.0 | 11.29 | 13.15 | 12.60 | below |

All three CMIP6 models published both nitrate and phosphate; no model result is missing.

**Interpretation**
- WOA23 shows a below-Redfield mean N:P of 14.15. Its Atlantic mean is highest at 15.55, close to Redfield, while the Pacific and Indian means are lower (13.57 and 13.73). The low-N:P (<12) volume is modest globally (6.2%) but larger in the Pacific (8.1%) and Indian Ocean (10.3%), consistent with a denitrification-related low-N:P imprint concentrated in those basins.
- CESM2 is far too low in mean N:P (11.83), has a very poor nitrate–phosphate relationship (R² = 0.238, slope 4.89, intercept 10.84 µmol kg⁻¹), and massively overproduces low-N:P water: 36.7% global, 59.4% Pacific, and 35.8% Indian. Its basin means are all much lower than WOA23. It does not reproduce the observed low-N:P volume or basin distribution.
- GFDL-ESM4 has a global mean N:P very close to WOA23 (14.13 vs 14.15) and a similar global N:P < 12 fraction (7.3% vs 6.2%), but it fails on distribution: the 14–18 fraction is much too large (74.6% vs 55.0%), Pacific N:P < 12 is too high (14.9% vs 8.1%), and Indian N:P < 12 is far too low (0.5% vs 10.3%). Its basin means are Pacific 13.39 (close), Atlantic 15.07 (too low), and Indian 14.72 (too high).
- MPI-ESM1-2-LR also has a low mean N:P (12.28), overproduces N:P < 12 globally and in both Pacific and Indian basins (22.5%, 40.0%, and 15.0%), and almost no water falls in 14 < N:P < 18 (0.5%). Its basin means are all lower than WOA23. It does not reproduce the observed pattern.
- All models have lower OLS slopes and higher intercepts than WOA23’s slope 13.04 and intercept 1.69 µmol kg⁻¹. Only GFDL-ESM4 has a moderately high R² (0.817) approaching the observed 0.908; CESM2 (0.238) and MPI-ESM1-2-LR (0.678) are poorer.

**Visualizations generated**
- **Figure 1:** A 2×2 scatter of nitrate (y) versus phosphate (x) for the four sources in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with identical axis ranges and ticks. Points were subsampled to 20,000 per source for display only. Each panel overlays its own OLS line and the Redfield 16:1 line through the origin, with distinguishable line styles, and is annotated with slope a, intercept b, and R².
- **Figure 2:** A 2×2 global map panel of the 200–1000 m volume-weighted mean N:P minus 16 for the four sources on a common longitude–latitude frame. A single shared diverging colour scale centred on zero is used, with the limit set from the 98th percentile of the pooled absolute values. CMIP6 fields were bilinearly regridded to the WOA23 regular grid for display only. The figure is intended to show where intermediate-water N:P deviates spatially from Redfield 16.
