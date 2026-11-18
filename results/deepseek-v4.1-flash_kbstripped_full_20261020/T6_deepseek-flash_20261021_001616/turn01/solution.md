**Assessment:** All three CMIP6 models place the global nitracline shallower than the WOA23 reference. **MPI-ESM1-2-LR performs best overall**: it has the highest spatial correlation and lowest RMSE against WOA23 on the shared 1° comparison grid, and it is the only model that preserves the WOA23 basin ordering of nitracline depth.

## Methodology used by the executed code

- **WOA23 reference:** annual nitrate climatology loaded with `load_climatology_data(dataset="woa23", variable="nitrate", time_slice="annual", depth_range=[0, 1200])`; masked with the WOA23 land–sea mask where `landsea > 1`. WOA23 areas, depth levels and the basin mask at `depth_levels="0m"` were loaded.
- **CMIP6 data:** historical monthly `no3`, `thetao` and `so` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR over **1971-01 to 2000-12**, loaded on each model’s native grid over 0–1200 m. The CMIP6 missing sentinel was removed with `da.where(da < 1e10)` before any arithmetic.
- **Nitrate unit conversion:** model nitrate was converted from mol m⁻³ to µmol kg⁻¹ using in-situ density computed from monthly `thetao` and `so` with `compute_potential_density_anomaly(ref_pressure=None)`, then passed to `transform_unit(..., to_unit="umol kg-1", density=sigma)`.
- **Temporal aggregation:** model nitrate was averaged over 1971–2000 with `temporal_aggregate(freq="mean", op="mean")`; the nitracline was then computed from the time-mean nitrate field.
- **Nitracline definition:** `compute_nutricline_depth(..., min_depth=0, max_depth=1000)`, i.e. the depth of the maximum vertical nitrate gradient searched within the top 1000 m, with the profile loaded below 1000 m so the deepest gradient can use the level just below 1000 m.
- **Global and basin means:** computed on each source’s native grid with `compute_basic_stat(op="area_weighted_mean")`. Basin selection used the WOA23 `basinmask_01` definition via `select_basin(..., basin=b, include_marginal_seas=True, depth_dependent=False)`, keeping only Pacific, Atlantic and Indian; Southern and Arctic were excluded.
- **Comparison grid and metrics:** nitracline fields were regridded to the WOA23 1° grid with `regrid_horizontal(..., method="conservative")`. A common mask was built with `build_comparison_mask(..., ocean_only=True)` and applied to every field. The shared comparison mask contains **40,621 cells**. Pearson correlation used `analyze_correlation(..., method="pearson", weighting="area")`; RMSE used `compute_difference_metrics(..., weighting="auto")`.

## Global mean nitracline depth

| Source | Global mean (m) | Difference from WOA23 (m) |
|---|---:|---:|
| WOA23 | **166** | — |
| CESM2 | **147** | **−19** |
| GFDL-ESM4 | **133** | **−33** |
| MPI-ESM1-2-LR | **132** | **−34** |

All three models simulate a global nitracline that is too shallow relative to WOA23, by roughly 19–34 m.

## Basin mean nitracline depths

Values are metres to the nearest whole metre; differences are model minus WOA23, computed from unrounded fields.

| Source | Pacific | Atlantic | Indian |
|---|---:|---:|---:|
| **WOA23** | **201** | **165** | **180** |
| **CESM2** | **162** (−38) | **110** (−56) | **165** (−14) |
| **GFDL-ESM4** | **159** (−41) | **91** (−74) | **169** (−11) |
| **MPI-ESM1-2-LR** | **184** (−17) | **111** (−54) | **121** (−59) |

The Atlantic nitracline is strongly underestimated by all three models, by 54–74 m. The Pacific is closest in MPI-ESM1-2-LR (−17 m), whereas CESM2 and GFDL-ESM4 are too shallow by about 38–41 m. The Indian basin is too shallow in MPI-ESM1-2-LR by about 59 m but relatively close in CESM2 and GFDL-ESM4.

## Spatial metrics on the shared WOA23 1° grid

Common comparison mask: **40,621 cells**. The shared-sample WOA23 mean is **167.0 m**.

| Model | Pearson r | RMSE (m) | Bias (m) | Mean model (m) | Mean WOA23 (m) |
|---|---:|---:|---:|---:|---:|
| CESM2 | **0.401** | **163.1** | −17.0 | 150.0 | 167.0 |
| GFDL-ESM4 | **0.499** | **161.1** | −32.0 | 135.0 | 167.0 |
| MPI-ESM1-2-LR | **0.599** | **151.6** | −32.7 | 134.3 | 167.0 |

**Best model:** **MPI-ESM1-2-LR**, by both the highest shared-mask correlation (0.599) and the lowest RMSE (151.6 m). The correlations are nevertheless only low-to-moderate and the RMSEs are large, indicating substantial spatial mismatch in all models.

## Basin ordering check

WOA23 orders the nitracline from deepest to shallowest as:

**Pacific (201 m) > Indian (180 m) > Atlantic (165 m).**

- **CESM2:** Indian > Pacific > Atlantic — does **not** match.
- **GFDL-ESM4:** Indian > Pacific > Atlantic — does **not** match.
- **MPI-ESM1-2-LR:** Pacific > Indian > Atlantic — **matches WOA23**.

Thus MPI-ESM1-2-LR is the only model that reproduces the observed three-basin ordering.

## Figures

- **Figure 1:** a 2×2 global map panel of regridded nitracline depth for CESM2, GFDL-ESM4, MPI-ESM1-2-LR and WOA23. All panels use one shared sequential colour scale, with the upper limit set from the pooled 98th percentile and deeper values shown as visually heavier/darker. The code uses the `viridis_r` colormap over a global longitude–latitude domain.
- **Figure 2:** a four-panel longitude–depth section of nitrate concentration through the equatorial Pacific, averaged over 2°S–2°N and covering 120°E–80°W (plotted as 120–280° longitude) and the top 1000 m, with depth increasing downward. Panels are ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, share one nitrate colour scale in µmol kg⁻¹, and each panel overlays its own nitracline depth as a black line. The latitude average was performed with `compute_basic_stat(op="mean", dim=lat)` for WOA23 and the analogous native latitude dimension for the models; the section nitracline was computed with the same 0–1000 m maximum-gradient rule before trimming the plotted section to 1000 m.

## Important insights

- The models are systematically biased toward a **shallower global nitracline** than WOA23.
- The **Atlantic** is the most poorly reproduced basin in all models, with especially large shallowing in GFDL-ESM4 (−74 m).
- **MPI-ESM1-2-LR** has the strongest spatial pattern agreement and the correct basin ordering, but it still underestimates the global mean depth by about 34 m and strongly underestimates the Indian and Atlantic nitraclines.
- The native-grid means and the shared-mask metrics differ slightly because the latter are computed after conservative regridding to the WOA23 1° grid and after applying the common comparison mask; for example, WOA23’s native global mean is 166 m while the shared-sample mean is 167.0 m.
- The nitracline is quantised by each product’s vertical grid, and the models have different vertical resolutions, so the maps and metrics should be interpreted alongside each model’s native level structure.
