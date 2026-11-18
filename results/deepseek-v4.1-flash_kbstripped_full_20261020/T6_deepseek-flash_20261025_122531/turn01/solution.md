**Method summary.** The nitracline was defined as the depth of the maximum vertical nitrate gradient searched only over levels in the top 1000 m, with the gradient at the deepest searched level using the level just below 1000 m; profiles were loaded to 0–2000 m for the nitracline calculation. The WOA23 reference was the annual objectively analysed nitrate climatology on the 1° grid. CMIP6 historical monthly `no3`, `thetao` and `so` were loaded on each model’s native grid for 1971-01 to 2000-12 and averaged with `temporal_aggregate(freq="mean", op="mean")`; the nitracline was then computed from the averaged nitrate field, not from monthly nitraclines. Global and basin means used `compute_basic_stat(op="area_weighted_mean")` over lat/lon. Basin means used `select_basin(..., include_marginal_seas=True, depth_dependent=False)` for the Pacific, Atlantic and Indian Oceans, so marginal seas were folded into their parent basins and the Southern and Arctic Oceans were excluded. For pattern metrics, each model nitracline was conservatively regridded to the WOA 1° grid, a single shared valid-ocean mask was built with `build_comparison_mask(ocean_only=True)` (40,623 cells), correlation used `analyze_correlation(..., weighting="area")`, and RMSE used `compute_difference_metrics(..., weighting="auto")`. Model sections for Figure 2 were converted from mol m−3 to µmol kg−1 using in-situ density from `compute_potential_density_anomaly(..., ref_pressure=None)`; the section latitude reduction used `compute_basic_stat(op="mean", dim="lat")` without an explicit spatial-weighting argument, and vertical interpolation used `regrid_vertical(method="linear")` to WOA levels including 1050 m.

**Global mean nitracline depth (1971–2000 climatology for models; WOA23 annual climatology), rounded to nearest metre.**

| Source | Global mean (m) | Difference from WOA23 (m) |
|---|---:|---:|
| WOA23 | 166 | — |
| CESM2 | 147 | -18 |
| GFDL-ESM4 | 134 | -32 |
| MPI-ESM1-2-LR | 132 | -34 |

All three models place the global-mean nitracline shallower than WOA23, by 18–34 m.

**Basin-mean nitracline depth, rounded to nearest metre; differences in parentheses are model minus WOA23.** Marginal seas are folded into each parent ocean using the WOA23 basin definitions; Southern and Arctic Oceans are excluded.

| Source | Pacific (m) | Atlantic (m) | Indian (m) |
|---|---:|---:|---:|
| WOA23 | 201 | 165 | 180 |
| CESM2 | 163 (-38) | 110 (-56) | 166 (-14) |
| GFDL-ESM4 | 160 (-41) | 91 (-74) | 170 (-10) |
| MPI-ESM1-2-LR | 184 (-17) | 112 (-54) | 121 (-58) |

WOA23’s shallow-to-deep basin ordering is **Atlantic (165 m) < Indian (180 m) < Pacific (201 m)**. Only **MPI-ESM1-2-LR** reproduces that ordering: Atlantic (112 m) < Indian (121 m) < Pacific (184 m). CESM2 and GFDL-ESM4 both produce Atlantic < Pacific < Indian, so they do not reproduce the observed ordering.

**Pattern statistics on the shared WOA 1° comparison mask (40,623 cells).**

| Model | Pearson correlation with WOA23 | RMSE (m) | Shared-mask bias (m, model − WOA23) |
|---|---:|---:|---:|
| CESM2 | 0.404 | 162.9 | -16.7 |
| GFDL-ESM4 | 0.499 | 161.3 | -31.5 |
| MPI-ESM1-2-LR | 0.598 | 151.7 | -32.5 |

The shared-mask biases differ from the global-mean differences above because the metrics use the common regridded valid-data mask rather than each model’s full native valid field.

**Best model.** **MPI-ESM1-2-LR** reproduces the observed nitracline best overall: it has the highest area-weighted spatial correlation (0.598) and the lowest RMSE (151.7 m) against WOA23, and it is the only model that reproduces WOA23’s ordering of Atlantic < Indian < Pacific basin-mean nitracline depth. However, it still underestimates the global mean by 34 m and is substantially too shallow in the Atlantic and Indian basins. CESM2 has the smallest global-mean bias but the lowest pattern correlation and highest RMSE; GFDL-ESM4 is intermediate in correlation and RMSE but has a strong shallow bias, especially in the Atlantic. The large RMSE values indicate substantial local depth mismatches, consistent with the nitracline being an extremum depth quantised by each model’s native vertical grid.

**Visualizations produced.** Figure 1 is a 2×2 global map panel of nitracline depth for WOA23, CESM2, GFDL-ESM4 and MPI-ESM1-2-LR on the WOA 1° grid, using one shared sequential `YlGnBu` scale spanning 0 m to the pooled 98th percentile so that deeper nitraclines are visually heavier; WOA23 land–sea masking (ocean > 1) is applied on the maps. Figure 2 is a four-panel longitude–depth section of nitrate through the equatorial Pacific (2°S–2°N, 120°E–80°W, top 1000 m, depth increasing downward), in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with one shared `viridis` scale in µmol kg−1 up to the pooled 98th percentile. Each panel overlays its own nitracline depth computed under the top-1000 m rule; model nitrate was converted to µmol kg−1 using in-situ density as described above.

**Key insights.** WOA23 shows the deepest basin-mean nitracline in the Pacific and the shallowest in the Atlantic. All models are biased shallow globally, but the basin-level biases differ: CESM2 and GFDL-ESM4 are especially shallow in the Atlantic, while MPI-ESM1-2-LR is especially shallow in the Indian and Atlantic basins. MPI-ESM1-2-LR’s better spatial pattern and correct basin ordering make it the best of the three for this diagnostic, even though its global-mean depth is too shallow.
