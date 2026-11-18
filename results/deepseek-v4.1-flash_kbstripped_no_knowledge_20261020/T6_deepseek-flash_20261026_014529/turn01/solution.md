## 1. Direct answer

The nitracline depth (depth of the maximum vertical nitrate gradient in the top 1000 m) was derived for WOA23 and for the 1971–2000 climatologies of CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, and assessed globally, by WOA basin and spatially on the WOA23 1° grid.

All three CMIP6 models put the nitracline **too shallow** relative to WOA23, and all three are much too shallow in the Atlantic. **MPI-ESM1-2-LR reproduces the observed nitracline best overall**: it has the highest spatial correlation (r = 0.608), the lowest RMSE (156.3 m), a global mean within ~24 m of WOA23, and it is the only model that reproduces the observed basin ordering (Atlantic < Indian < Pacific). CESM2 has the smallest global-mean bias (−10 m) but a much weaker spatial pattern (r = 0.370), while GFDL-ESM4 has the largest Atlantic and global biases.

## 2. Methodology (as actually executed)

**WOA23 reference.** Nitrate was loaded as the WOA23 all-time-periods annual climatology (`load_climatology_data(dataset="woa23", variable="nitrate", time_slice="annual", depth_range=[0, 1200])`), together with its depth levels/bounds, the WOA23 land–sea mask (`load_landsea_mask`, ocean defined as value > 1) and the 1° cell areas (`load_areas`). Concentrations were converted to µmol kg⁻¹ with `transform_unit`. The profile was loaded over 0–1200 m so that the deepest retained top-1000 m gradient has a level below it.

**CMIP6 models.** Monthly `no3` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR was loaded from `cmip6_historical` over 1971-01 to 2000-12 on each native grid (native areacello provenance retained), over 0–1200 m. A day-weighted long-term mean was formed with `temporal_aggregate(freq="mean", op="mean")` and converted to µmol kg⁻¹.

**Nitracline rule.** For every source, `compute_vertical_gradient` was applied to the full loaded column, then gradients were restricted to levels with depth ≤ 1000 m and the depth of the maximum retained gradient (`idxmax`) was stored. Because the gradient is computed before the depth selection, the 1000 m gradient uses the next deeper level. Ocean masks were WOA `landsea > 1` and, for models, cells with a complete non-missing top-1000 m nitrate column.

**Averaging and regridding.** Global means are area-weighted on each native grid (model areacello; WOA 1° cell areas) via `compute_basic_stat(op="area_weighted_mean")`. Model nitraclines were conservatively regridded to the WOA 1° grid with `regrid_horizontal(method="conservative")`. Basin means use `select_basin(..., include_marginal_seas=True)` so marginal seas are folded into their parent ocean, with Southern and Arctic Oceans excluded, and are area-weighted with WOA 1° areas.

**Skill metrics.** `build_comparison_mask(fields=[WOA23, all three regridded models], ocean_only=True, extra_mask=(landsea > 1))` defined common valid ocean cells. Pearson correlation (`analyze_correlation`) and `compute_difference_metrics` (RMSE, bias, centred RMSE) were computed over those cells on the 1° grid.

## 3. Key numerical results

**Global mean nitracline depth (m; nearest metre) and difference from WOA23**

| Source | Global mean | Δ vs WOA23 |
|---|---|---|
| WOA23 | 166 | — |
| CESM2 | 156 | −10 |
| GFDL-ESM4 | 140 | −26 |
| MPI-ESM1-2-LR | 141 | −24 |

**Basin means (m; nearest metre), marginal seas folded in, Southern/Arctic excluded**

| Source | Pacific | Atlantic | Indian |
|---|---|---|---|
| WOA23 | 201 | 165 | 180 |
| CESM2 | 168 (−33) | 114 (−51) | 169 (−11) |
| GFDL-ESM4 | 163 (−38) | 94 (−71) | 172 (−7) |
| MPI-ESM1-2-LR | 190 (−11) | 117 (−48) | 123 (−56) |

Numbers in parentheses are differences from WOA23.

**Spatial skill on the WOA23 1° grid (common valid ocean cells)**

| Model | Pearson r | RMSE (m) | Bias (m) |
|---|---|---|---|
| CESM2 | 0.370 | 167.9 | −19.0 |
| GFDL-ESM4 | 0.493 | 165.1 | −34.8 |
| MPI-ESM1-2-LR | 0.608 | 156.3 | −34.4 |

**Basin ordering (shallowest → deepest).** WOA23: Atlantic < Indian < Pacific. CESM2 and GFDL-ESM4: Atlantic < Pacific < Indian (do not match). MPI-ESM1-2-LR: Atlantic < Indian < Pacific (matches WOA23). Thus only MPI reproduces the observed ordering.

## 4. Visualizations produced

- **Figure 1** — 2×2 global maps of nitracline depth (CESM2, GFDL-ESM4, MPI-ESM1-2-LR, WOA23), all on one shared sequential YlGnBu scale set from the pooled 2nd–98th percentiles so deeper values are visually heavier; land is grey. Titles include each source's native-grid global mean. This lets the spatial pattern and bias be compared directly across models and observations.
- **Figure 2** — four stacked longitude–depth sections of nitrate concentration along the equatorial Pacific (2°S–2°N, 120°E–80°W, crossing the dateline), in the order WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, shown for the top 1000 m with depth increasing downward and one shared colour scale in µmol kg⁻¹. A red line in each panel is that panel's nitracline computed by the same top-1000 m maximum-gradient rule. Section nitracline means are ~114 m (WOA23), ~120 m (CESM2), ~117 m (GFDL-ESM4) and ~34 m (MPI-ESM1-2-LR); the MPI section is notably coarse (only 16 longitude columns survive averaging on its native grid), so this apparent shallow section value should be treated with caution.
- **Figure 3** (additional summary generated by the code) — three bar panels: (a) global mean nitracline by source, (b) basin means by source, and (c) spatial Pearson r and RMSE per model.

## 5. Interpretation and caveats

- **Systematic shallow bias.** Every model places the nitracline shallower than observed, globally and in every basin. The bias is modest in the Indian Ocean for CESM2 (−11 m) and GFDL-ESM4 (−7 m) and for MPI in the Pacific (−11 m), but is large in the Atlantic for all models (−48 to −71 m).
- **Best model.** MPI-ESM1-2-LR wins on the two spatial metrics (r = 0.608, RMSE = 156.3 m) and is the only model matching the observed inter-basin ordering; CESM2 has the smallest global-mean bias but a weak spatial pattern and the second-largest RMSE.
- **Large RMSEs.** RMSE (~156–168 m) is comparable to the nitracline depth spread itself (0–900 m), so none of the models captures the full observed spatial structure; the low correlations (0.37–0.61) confirm only moderate pattern agreement.
- **Data caveats.** CESM2's native climatology contains small negative nitrate values (down to about −1.2 µmol kg⁻¹). The model nitracline fields were masked to cells with a complete top-1000 m column, and the skill statistics are computed over the intersection of valid cells across all sources, so regional coverage differs from the native-grid global means.
