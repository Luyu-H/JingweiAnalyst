The executed analysis compares the depth of the maximum vertical nitrate gradient within the top 1000 m (the nitracline) for CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23. CMIP6 historical nitrate was used for January 1971–December 2000, aggregated with `temporal_aggregate(..., freq="mean", op="mean")`, loaded from 0–1300 m so the gradient at 1000 m could use the next deeper level, converted from mol m⁻³ to µmol kg⁻¹, and passed to `compute_nutricline_depth(min_depth=0, max_depth=1000)`. WOA23 used the all-time-periods annual nitrate climatology. Global means were computed with `compute_basic_stat(..., op="area_weighted_mean", dim=["lat","lon"])` on each source’s own grid. For basin statistics, models were conservatively regridded to the WOA23 1° grid. Basins were selected from the WOA23 basin mask for Pacific, Atlantic and Indian with `select_basin(..., include_marginal_seas=True)`; Southern and Arctic basins were not included. Correlation and RMSE were computed over the common ocean-valid cell mask obtained from `build_comparison_mask(..., ocean_only=True)` after regridding all model nitracline fields to the WOA23 grid. Differences are model minus WOA23 and were computed before rounding.

### Global mean nitracline depth

| Source | Global mean depth (m) | Model − WOA23 (m) |
|---|---:|---:|
| WOA23 | 166 | — |
| CESM2 | 147 | −18 |
| GFDL-ESM4 | 134 | −32 |
| MPI-ESM1-2-LR | 132 | −34 |

All three models simulate a globally shallower nitracline than WOA23.

### Basin mean nitracline depth
Values are mean depth in metres; parentheses give model minus WOA23, computed before rounding. Marginal seas are folded into the parent basin.

| Basin | WOA23 | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---|---:|---:|---:|---:|
| Pacific | 201 | 163 (−38) | 160 (−41) | 184 (−17) |
| Atlantic | 165 | 109 (−56) | 91 (−74) | 111 (−54) |
| Indian | 180 | 166 (−14) | 170 (−10) | 121 (−58) |

WOA23 basin ordering from deepest to shallowest is **Pacific > Indian > Atlantic**.  
CESM2 gives **Indian > Pacific > Atlantic**.  
GFDL-ESM4 gives **Indian > Pacific > Atlantic**.  
MPI-ESM1-2-LR gives **Pacific > Indian > Atlantic**.

Thus only MPI-ESM1-2-LR reproduces the WOA23 ordering of the three basins; CESM2 and GFDL-ESM4 both place the Indian Ocean slightly deeper than the Pacific, contrary to WOA23.

### Spatial pattern statistics on the WOA23 1° common mask

| Model | Pearson r | RMSE (m) | Bias (m) |
|---|---:|---:|---:|
| CESM2 | 0.398 | 162.9 | −16.7 |
| GFDL-ESM4 | 0.503 | 161.3 | −31.5 |
| MPI-ESM1-2-LR | 0.618 | 151.7 | −32.5 |

MPI-ESM1-2-LR has the highest spatial correlation and the lowest RMSE, so it is the best model by the RMSE criterion used in the code. However, all correlations are low to moderate and the RMSE values are very large relative to the mean nitracline depth, indicating substantial spatial mismatches. CESM2 has the smallest global mean bias but a weaker spatial correlation and larger RMSE than MPI-ESM1-2-LR.

### Equatorial Pacific section diagnostics
For the 120°E–280°E, 2°S–2°N section, after horizontal/vertical regridding to the WOA23 section grid and averaging over latitude, the mean section nitraclines are:

| Source | Mean section nitracline (m) |
|---|---:|
| WOA23 | 114 |
| CESM2 | 121 |
| GFDL-ESM4 | 119 |
| MPI-ESM1-2-LR | 39 |

MPI-ESM1-2-LR is much shallower than the other sources in this equatorial Pacific section.

### Figures
**Figure 1** is a 2×2 panel of global nitracline-depth maps in the order CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23. All panels use one shared sequential `YlGnBu` colour scale, with deeper values visually heavier/darker. The shared scale is clipped at the pooled 1st and 99th percentiles, and missing values are shown in grey. Panel titles include the global mean depth and, for models, the difference from WOA23, the spatial correlation and RMSE. The maps are shown over longitude −180° to 180° and latitude −70° to 75°.

**Figure 2** shows four longitude–depth nitrate sections for the equatorial Pacific from 120°E to 280°E, averaged over 2°S–2°N, in the order WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. All panels cover the top 1000 m with depth increasing downward and use a shared `YlGnBu` concentration scale in µmol kg⁻¹, clipped at the pooled 1st/99th percentiles. Each panel overlays its own nitracline depth as a red line. Model sections were regridded horizontally to the WOA23 section template and vertically to WOA23 depth levels before display.

### Overall assessment
The models capture only a modest amount of the observed nitracline spatial structure. All models are biased shallow globally, especially in the Atlantic. MPI-ESM1-2-LR performs best by spatial correlation and RMSE on the WOA23 grid and is the only model that preserves the observed WOA23 basin ordering, but it still has a large negative global bias, a very shallow equatorial Pacific section nitracline, and the Indian and Atlantic basins are much shallower than observed. CESM2 has the smallest global mean bias but weaker spatial agreement; GFDL-ESM4 has intermediate correlation and a strong shallow bias, particularly in the Atlantic and Pacific.
