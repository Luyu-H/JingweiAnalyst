### Overall assessment
All three CMIP6 models reproduce the broad spatial pattern of WOA23 sea-surface salinity (SSS) well, with spatial correlations of 0.893–0.926 on the common WOA23 1° comparison mask. Their 1971–2000 global-mean SSS is slightly fresher than WOA23, by 0.07–0.32 psu on their native grids. Historical 1950–2014 global-mean SSS trends are small and not statistically significant at the 95 % level in any model. In contrast, the 2071–2100 ssp245 segments show negative global-mean SSS trends in all models, significant for CESM2 and GFDL-ESM4 but not for MPI-ESM1-2-LR. End-of-century ΔSSS (ssp245 2071–2100 minus historical 1995–2014) is negative globally in all three models, with basin-dependent signs.

### Methodology actually applied
WOA23 salinity was loaded as the annual 1971–2000 climatology. SSS was computed as the full-cell-volume-weighted mean over the salinity levels whose depths lie in 0–10 m; WOA23 used the 0, 5 and 10 m levels with full layer thicknesses from the layer boundaries and cell areas, then restricted to `landsea > 1`. CMIP6 `so` was loaded monthly for historical 1950-01 to 2014-12 and ssp245 2071-01 to 2100-12 for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR. Model SSS was computed on each native grid as the full-`volcello`-weighted mean over model levels with depths in 0–10 m. Salinity remained practical salinity in psu; no conversion was applied.

The 1971–2000 model climatologies were formed with `temporal_aggregate` over 1971-01 to 2000-12 and global means with `compute_basic_stat(..., op="area_weighted_mean", dim=["lat","lon"])`. For comparison with WOA23, each model climatology was conservatively regridded with `regrid_horizontal` to the WOA23 1° grid. A single common comparison mask was built with `build_comparison_mask` over all three regridded model fields plus WOA23 (`ocean_only=True`), retaining **40,829 cells**. Bias and RMSE were computed with `compute_difference_metrics` using its `weighting="auto"` (WOA23 lat/lon areas), and spatial Pearson correlation with `analyze_correlation(..., weighting="area")`.

Trends of the annual global-mean SSS were computed from annual means (`temporal_aggregate(..., freq="year")`) followed by `fit_1D_trend(..., method="linear")`; slopes are reported per decade by multiplying the fitted per-year slope by 10. Confidence intervals are the serial-correlation-corrected 95 % intervals returned by the tool. Cell-by-cell historical trend maps used `fit_trend_map(..., alpha=0.05, per_decade=True)` on annual means; for display only the slopes were conservatively regridded to the WOA23 grid and the significance field was regridded by nearest neighbour and thresholded at 0.5. ΔSSS was computed on native grids as the ssp245 2071–2100 mean minus the historical 1995–2014 mean, then area-weighted globally and over the WOA23 basin mask using `select_basin(..., include_marginal_seas=True, depth_dependent=False)` for the Atlantic, Indian and Pacific only; the Southern and Arctic basins were excluded.

### Key numerical results

**1971–2000 global-mean SSS and native-grid difference from WOA23**

| Model | Global mean SSS (psu) | Model − WOA23 (psu) |
|---|---:|---:|
| WOA23 | 34.63 | — |
| CESM2 | 34.32 | −0.32 |
| GFDL-ESM4 | 34.56 | −0.07 |
| MPI-ESM1-2-LR | 34.36 | −0.27 |

On the common WOA23 comparison mask, the area-weighted biases were CESM2 **−0.31 psu**, GFDL-ESM4 **−0.03 psu**, and MPI-ESM1-2-LR **−0.27 psu**. These differ slightly from the native-grid global-mean differences because of conservative regridding and the common ocean mask.

**Bias, RMSE and spatial correlation against WOA23 on the common 40,829-cell WOA23 1° mask**

| Model | Bias (psu) | RMSE (psu) | Pearson r |
|---|---:|---:|---:|
| CESM2 | −0.31 | 0.77 | 0.903 |
| GFDL-ESM4 | −0.03 | 0.60 | 0.926 |
| MPI-ESM1-2-LR | −0.27 | 0.76 | 0.893 |

GFDL-ESM4 has the smallest bias and RMSE and the highest spatial correlation. CESM2 and MPI-ESM1-2-LR are both fresher than WOA23 on average and have slightly larger RMSE, while still reproducing the large-scale salinity pattern well.

**Linear trends of annual global-mean SSS**

Historical 1950–2014:

| Model | Slope (psu per decade) | 95 % CI | R² | Interval excludes zero? |
|---|---:|---:|---:|---|
| CESM2 | −2.68×10⁻³ | −6.93×10⁻³ to +1.57×10⁻³ | 0.245 | No |
| GFDL-ESM4 | −1.10×10⁻³ | −3.61×10⁻³ to +1.41×10⁻³ | 0.0639 | No |
| MPI-ESM1-2-LR | −5.51×10⁻⁴ | −5.04×10⁻³ to +3.93×10⁻³ | 0.00718 | No |

ssp245 2071–2100:

| Model | Slope (psu per decade) | 95 % CI | R² | Interval excludes zero? |
|---|---:|---:|---:|---|
| CESM2 | −2.09×10⁻² | −2.55×10⁻² to −1.64×10⁻² | 0.910 | Yes |
| GFDL-ESM4 | −5.34×10⁻³ | −1.06×10⁻² to −3.86×10⁻⁵ | 0.324 | Yes |
| MPI-ESM1-2-LR | −2.31×10⁻² | −5.13×10⁻² to +5.04×10⁻³ | 0.701 | No |

Thus, the historical global-mean SSS trends are not distinguishable from zero after accounting for serial correlation. The late-century ssp245 trends are negative in all three models, but only CESM2 and GFDL-ESM4 have 95 % intervals excluding zero; MPI-ESM1-2-LR has a steeper fitted slope but a broad interval because of low effective sample size.

**Ratio of ssp245 slope to historical slope**

| Model | Ratio (1 d.p.) | Meaning |
|---|---:|---|
| CESM2 | 7.8 | Descriptive only |
| GFDL-ESM4 | 4.9 | Descriptive only |
| MPI-ESM1-2-LR | 42.0 | Descriptive only; poorly constrained |

These ratios divide two independently fitted slopes from disjoint time periods and different forcing. They are **not physical climate sensitivities**. The MPI-ESM1-2-LR ratio is especially unstable because its historical slope is very small and its 95 % interval includes zero. For CESM2, the ratio is additionally affected by the historical/ssp245 ensemble-member mismatch noted below.

**ΔSSS: ssp245 2071–2100 minus historical 1995–2014**

| Model | Global | Atlantic | Indian | Pacific |
|---|---:|---:|---:|---:|
| CESM2 | −0.12 | +0.11 | +0.06 | −0.24 |
| GFDL-ESM4 | −0.01 | +0.19 | −0.00 | −0.04 |
| MPI-ESM1-2-LR | −0.09 | −0.04 | −0.05 | −0.11 |

All values are in psu. GFDL-ESM4’s Indian basin value rounds to −0.00 psu but is −0.004 psu before rounding. The global ΔSSS is negative in all three models, indicating net global surface freshening between the late historical period and the late twenty-first century under ssp245. Basin changes are not uniform: the Atlantic is saltier in CESM2 and GFDL-ESM4 but fresher in MPI-ESM1-2-LR; the Pacific freshens in all three models.

### Figures generated
**Figure 1** contains three side-by-side global maps, one per model, of the SSS bias defined as model 1971–2000 climatology minus WOA23 on the common WOA23 comparison cells. A single diverging colour scale centred on zero is shared across panels, with colour limits set from the 98th percentile of absolute bias across the models, and the colour bar is labelled in psu. The map isolates where each model is saltier or fresher than the WOA23 reference.

**Figure 2** shows three side-by-side longitude–latitude maps of the cell-by-cell linear trend of annual-mean SSS over 1950–2014, in psu per decade. The trend field is computed before display regridding; for display it is conservatively placed on the WOA23 regular grid. A shared diverging colour scale centred on zero is used, with limits from the 98th percentile of absolute slope across models. Cells significant at p < 0.05 after the serial-correlation-aware fit are stippled on top of the trend field.

**Figure 3** is one time-series panel of the annual global-mean SSS for all three models. Historical 1950–2014 and ssp245 2071–2100 are plotted as separate segments, leaving the 2015–2070 gap empty. Each segment has its own fitted linear line overlaid, and the y axis is labelled in psu. The same model colours are used consistently across Figures 1–3.

### Important patterns and caveats
- All three models have high spatial correlation with WOA23 and modest global-mean biases, with GFDL-ESM4 performing best on the common-mask RMSE and correlation metrics.
- Historical global-mean SSS trends over 1950–2014 are small and not significant in any model. The late-century ssp245 global-mean trends are negative, but significance depends on the model.
- The negative global ΔSSS in all models indicates net surface freshening, but the basin response is strongly model-dependent, and the global mean hides compensating regional changes.
- **CESM2 ensemble-member caveat:** CESM2’s historical data here are from r1i1p1f1, while its ssp245 data are from r10i1p1f1. Its ssp245 trend, slope ratio and ΔSSS therefore combine forced change with internal-variability differences between two members. GFDL-ESM4 and MPI-ESM1-2-LR use r1 in both experiments.
- The ssp245/historical slope ratio is a descriptive comparison of two disjoint, differently forced segments, not a physical sensitivity or a single continuous trend.
