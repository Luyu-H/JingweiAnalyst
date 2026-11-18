**Direct answer**

The nitracline is defined here as the depth of the maximum vertical gradient of nitrate concentration within the top 1000 m; the profile data were loaded below 1000 m so that the deepest searched gradient could use the level just below 1000 m.

Against WOA23, all three CMIP6 models simulate a globally shallower nitracline. The WOA23 global mean is 166 m; CESM2 is 147 m (−18 m), GFDL-ESM4 is 134 m (−32 m), and MPI-ESM1-2-LR is 132 m (−34 m). According to the executed model-selection rule (lowest RMSE, tie-broken by absolute global bias), **MPI-ESM1-2-LR performs best**, with the highest pattern correlation and lowest RMSE, although it has a large shallow global bias. Only **MPI-ESM1-2-LR reproduces the observed basin ordering** Pacific > Indian > Atlantic.

**Methodology actually used in the executed code**

- WOA23 nitrate was loaded as the all-time-period annual objectively analysed climatology (`time_slice="annual"`), with attributes indicating the 1965–2022 composite period. The WOA23 land–sea mask was applied as ocean where `landsea > 1`.
- CMIP6 historical `no3` was loaded monthly for 1971-01 to 2000-12 for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. For the main nitracline calculation, the monthly nitrate was aggregated with `temporal_aggregate(no3_raw, freq="mean")` before computing the nitracline. The executed code did **not** pass an explicit day-weight or calendar argument in that call.
- Nitracline depth was calculated with `compute_nutricline_depth(min_depth=0, max_depth=1000)` on WOA23 and on each model after temporal aggregation. No density conversion was applied for the nitracline depth itself.
- Global and basin means were requested with `compute_basic_stat(op="area_weighted_mean")` over the horizontal dimensions. The executed code did not explicitly call `load_areas`/load `areacello`; the means are those produced by the `area_weighted_mean` call.
- Basin selections used `select_basin(..., include_marginal_seas=True, depth_dependent=False)` for Pacific, Atlantic and Indian; Southern and Arctic were excluded. This is the WOA basin labelling, not an OMIP 10-region definition.
- Model nitracline maps were conservatively regridded to the WOA23 1-degree grid with `regrid_horizontal(..., target="woa23", method="conservative")` only for map comparison and Figure 1. Global and basin means were computed on each model’s native grid before regridding.
- A common comparison mask was built with `build_comparison_mask(..., ocean_only=True)` over WOA23 and the three regridded model maps. It contains **40,623** common cells.
- Pattern correlations used `analyze_correlation(..., method="pearson", weighting="area")`. RMSE values used `compute_difference_metrics(..., weighting="auto")`.
- For the equatorial Pacific section, WOA23 nitrate was used directly in µmol kg⁻¹. CMIP6 `no3` was converted from mol m⁻³ to µmol kg⁻¹ using `transform_unit(..., from_unit="mol m-3", to_unit="umol/kg", density=...)`, where density came from `compute_potential_density_anomaly(temp_data=thetao_mean, sal_data=so_mean, ref_pressure=None)`. The density assumption reported by the tool was in-situ density anomaly (sigma + 1000), applied cell-by-cell at each cell’s own in-situ pressure: CESM2 1018.6–1032.4 kg m⁻³, GFDL-ESM4 1019.8–1033.0 kg m⁻³, MPI-ESM1-2-LR 1020.5–1032.1 kg m⁻³.
- The section averaging over 2°S–2°N used `compute_basic_stat(op="mean", dim="lat")` in the executed code, then linear vertical regridding to WOA23 depth levels for plotting. The overlaid nitracline lines were computed from the native-vertical-grid section, not from the vertically regridded plotting field.

**Global mean nitracline depth**

| Source | Global mean (m, nearest m) | Model − WOA23 (m) |
|---|---:|---:|
| WOA23 | 166 | — |
| CESM2 | 147 | −18 |
| GFDL-ESM4 | 134 | −32 |
| MPI-ESM1-2-LR | 132 | −34 |

All three models are globally too shallow, by 18–34 m.

**Basin mean nitracline depth**

| Basin | WOA23 (m) | CESM2 (m) | CESM2 diff (m) | GFDL-ESM4 (m) | GFDL diff (m) | MPI-ESM1-2-LR (m) | MPI diff (m) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Pacific | 201 | 163 | −38 | 160 | −41 | 184 | −17 |
| Atlantic | 165 | 110 | −56 | 91 | −74 | 112 | −54 |
| Indian | 180 | 166 | −14 | 170 | −10 | 121 | −58 |

All model basin means are shallower than WOA23. The Atlantic is strongly too shallow in all models, especially GFDL-ESM4 (−74 m). MPI-ESM1-2-LR is closest in the Pacific but much too shallow in the Indian Ocean.

**Spatial comparison on the common WOA23 1-degree mask**

Common-mask cells: **40,623**.

| Model | Pearson correlation | RMSE (m) | Bias on common mask (m) |
|---|---:|---:|---:|
| CESM2 | 0.404 | 162.9 | −16.7 |
| GFDL-ESM4 | 0.499 | 161.3 | −31.5 |
| MPI-ESM1-2-LR | 0.598 | 151.7 | −32.5 |

The correlations are modest and the RMSE values are large relative to the observed mean depth, so none of the models reproduces the spatial structure of the observed nitracline strongly. MPI-ESM1-2-LR is best by correlation and RMSE, but it is not the closest by global-mean bias; CESM2 has the closest global mean but weaker spatial agreement.

**Basin ordering**

WOA23 order, deepest first: **Pacific (201 m) > Indian (180 m) > Atlantic (165 m)**.

Model orders, deepest first:

- CESM2: Indian > Pacific > Atlantic — does **not** reproduce the observed order.
- GFDL-ESM4: Indian > Pacific > Atlantic — does **not** reproduce the observed order.
- MPI-ESM1-2-LR: Pacific > Indian > Atlantic — **does** reproduce the observed order.

Thus, **MPI-ESM1-2-LR is the best of the three models overall** under the executed RMSE/correlation assessment, and it is also the only model that reproduces the observed ordering of mean nitracline depth across the three basins.

**Figure descriptions**

- **Figure 1** is a 2-by-2 global map panel on the common WOA23 1-degree mask, ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. All four panels use one shared sequential `YlGnBu` colour scale for nitracline depth, clipped at the pooled 98th percentile, so deeper nitraclines are visually darker/heavier. Panel titles include each source’s global mean depth, and a shared horizontal colourbar is labelled in metres.
- **Figure 2** is a 2-by-2 equatorial Pacific nitrate section for 2°S–2°N and 120°E–80°W, plotted over 0–1000 m with depth increasing downward. Panels are ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. All panels share one `YlGnBu` colour scale in µmol kg⁻¹, clipped at the pooled 98th percentile. The nitrate section is vertically regridded to WOA23 levels for plotting, while each panel’s nitracline is overlaid as a red line and is based on the native vertical-grid diagnostic.

**Important caveats and insights**

- WOA23 nutrients are an all-period composite (1965–2022 in the loaded attributes), whereas the models are 1971–2000 historical means. This period mismatch is unavoidable with these files.
- Nitracline depths are quantized by each model’s native vertical grid; horizontal regridding to 1° was used only for maps/metrics, not for the global and basin means.
- The broad pattern is that all models put the nitracline too shallow, especially in the Atlantic. MPI-ESM1-2-LR has the best pattern correlation and lowest RMSE but also a large negative global-mean bias. CESM2 is closest in global mean depth but has the weakest pattern correlation and the highest RMSE. Only MPI-ESM1-2-LR preserves the observed Pacific > Indian > Atlantic ordering.
