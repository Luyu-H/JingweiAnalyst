## Assessment

The nitracline was assessed as the depth of the maximum vertical nitrate gradient within the top 1000 m, using WOA23 annual nitrate climatology and 1971–2000 CMIP6 historical `no3` means from CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. All four nitracline fields were compared on the common WOA23 1° ocean mask.

Implementation details from the executed code: WOA23 nitrate and model nitrate were loaded over 0–1500 m; the nitracline was computed with `compute_nutricline_depth(..., min_depth=0, max_depth=1000)`. Monthly model nitrate was aggregated to a 1971–2000 long-term mean using `temporal_aggregate(..., freq="mean", op="mean")`; no explicit weighting argument was supplied to that call. Model nitrate and nitracline fields were regridded to the WOA23 grid with `regrid_horizontal(..., target="woa23", method="conservative")`; no explicit native `areacello` argument was passed to the regridding call. Global and basin statistics used `compute_basic_stat(..., op="area_weighted_mean", dim=["lat","lon"])` on the common mask. Basin means used `select_basin(..., include_marginal_seas=True, depth_dependent=False)` for the Pacific, Atlantic and Indian only, thereby excluding Southern and Arctic basins. Correlations used `analyze_correlation(..., method="pearson", weighting="area")` and RMSE used `compute_difference_metrics(..., weighting="area")`.

## Key numerical results

### Global mean nitracline depth

| Source | Mean depth (m) | Difference from WOA23 (m) |
|---|---:|---:|
| WOA23 | 167 | — |
| CESM2 | 150 | -17 |
| GFDL-ESM4 | 135 | -32 |
| MPI-ESM1-2-LR | 135 | -32 |

All three models simulate a globally shallower nitracline than WOA23. CESM2 is closest globally, while GFDL-ESM4 and MPI-ESM1-2-LR are both about 32 m too shallow.

### Basin mean nitracline depth

| Basin | WOA23 | CESM2 | CESM2 diff | GFDL-ESM4 | GFDL diff | MPI-ESM1-2-LR | MPI diff |
|---|---:|---:|---:|---:|---:|---:|---:|
| Pacific | 203 | 164 | -39 | 160 | -42 | 185 | -18 |
| Atlantic | 166 | 110 | -56 | 91 | -75 | 111 | -55 |
| Indian | 181 | 166 | -15 | 171 | -10 | 122 | -59 |

All model basin means are shallower than WOA23. The largest underestimations are in the Atlantic for CESM2 and GFDL-ESM4, and in the Indian Ocean for MPI-ESM1-2-LR. MPI-ESM1-2-LR is closest in the Pacific but much too shallow in the Indian Ocean.

### Spatial correlation and RMSE against WOA23

| Model | Pearson correlation | RMSE (m) |
|---|---:|---:|
| CESM2 | 0.404 | 162.9 |
| GFDL-ESM4 | 0.499 | 161.3 |
| MPI-ESM1-2-LR | 0.598 | 151.7 |

MPI-ESM1-2-LR has both the highest spatial correlation and the lowest RMSE, so it reproduces the observed nitracline pattern best among the three models. However, even its correlation is only moderate (0.598) and its RMSE remains large (151.7 m).

### Basin ordering

WOA23 orders the basins as:

**Atlantic (166 m) < Indian (181 m) < Pacific (203 m)**

- CESM2 order: Atlantic < Pacific < Indian — does **not** match WOA23.
- GFDL-ESM4 order: Atlantic < Pacific < Indian — does **not** match WOA23.
- MPI-ESM1-2-LR order: Atlantic < Indian < Pacific — **matches** WOA23.

Thus MPI-ESM1-2-LR is the only model that reproduces the observed ordering of nitracline depth across the three basins.

## Visualization description

**Figure 1** is a 2-by-2 panel of global nitracline-depth maps in the order WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. All panels use the common WOA23 1° grid and common comparison mask, with one shared sequential `YlGnBu` colour scale from 0 m to the 98th percentile of all plotted values. Deeper nitraclines therefore appear darker/heavier. A single colourbar is shared across the four panels.

**Figure 2** is a 2-by-2 longitude–depth section of nitrate concentration through the equatorial Pacific, averaged over 2°S–2°N and covering 120°E–80°W and the top 1000 m. Panels are ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with depth increasing downward and one shared `YlGnBu` colour scale in mol m⁻³, clipped at the 98th percentile. WOA23 nitrate was converted from micromoles per kilogram to mol m⁻³ using `transform_unit(..., to_unit="mol m-3", density=1.025)`. A red nitracline-depth line is overlaid in each panel; in the executed code this line was taken from the precomputed nitracline map on the WOA23 grid averaged over 2°S–2°N and subset to 120°E–80°W, rather than being recomputed directly from each plotted nitrate section.

## Interpretation

All three CMIP6 models produce a nitracline that is too shallow relative to WOA23, both globally and in each of the three major basins. The bias is particularly strong in the Atlantic for CESM2 and GFDL-ESM4, and in the Indian Ocean for MPI-ESM1-2-LR. MPI-ESM1-2-LR is the best of the three by RMSE and correlation and is the only model that preserves the observed Atlantic–Indian–Pacific ordering. Nevertheless, the moderate correlations and large RMSE values indicate substantial regional mismatches in nitracline depth in all models.
