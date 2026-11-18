**Direct answer.** In this analysis, all three CMIP6 models have a slightly lower area-weighted 1971–2000 global-mean SSS than WOA23: CESM2 by 0.32 psu, GFDL-ESM4 by 0.07 psu, and MPI-ESM1-2-LR by 0.27 psu. Their spatial patterns are strongly correlated with WOA23 (r = 0.893–0.926), while RMSE is 0.60–0.77 psu. The fitted 1950–2014 global-mean SSS trends are small and not significant once serial autocorrelation is accounted for. Late-century ssp245 trends over 2071–2100 are negative for all models, significant for CESM2 and GFDL-ESM4 but not for MPI-ESM1-2-LR. The global end-of-century ΔSSS is small and slightly negative for each model; basin changes are model-dependent.

### Methodology used by the executed code
- SSS was derived at each grid cell as the mean salinity of levels whose depths lie in 0–10 m, weighted by the **full volume** of each retained level. WOA23 used levels 0, 5, 10 m and layer thicknesses from `load_depth_levels`; for WOA23 the area factor cancels in the within-cell weighted mean. CMIP6 used each model’s `volcello`; if the volume field was monthly/time-varying it was collapsed to a static field before weighting and aligned to the salinity levels. No partial overlap weighting was used.
- Salinity was treated as practical salinity in psu directly; values of order 35 were compared with WOA23 without multiplication by 0.001 or 1000. No TEOS-10/Absolute Salinity conversion was performed.
- Historical data were monthly `so` from CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, 1950–2014. ssp245 data were monthly `so`, 2071–2100. The 1971–2000 climatology was formed with `temporal_aggregate(..., freq="mean")`. Annual fields used `freq="year"`.
- Model SSS climatologies were regridded conservatively to the WOA23 1° grid with `regrid_horizontal`. A single common comparison mask was built with `build_comparison_mask(..., ocean_only=True)` over all regridded model fields and WOA23, giving **40,829** comparison cells. All model and WOA fields were masked before `compute_difference_metrics` and `analyze_correlation(..., weighting="area")`.
- Global means were area-weighted with each product’s own `areacello` via `compute_basic_stat(..., op="area_weighted_mean")`, not simple latitude or grid-cell averages. The recorded source strings were `areacello published by cmip6_historical` and `volcello published by cmip6_historical` for all three models. Ocean areas behind the global means were approximately: WOA23 \(3.493\times10^{14}\) m²; CESM2 \(3.605\times10^{14}\) m²; GFDL-ESM4 \(3.611\times10^{14}\) m²; MPI-ESM1-2-LR \(3.613\times10^{14}\) m².
- Trends used `fit_1D_trend(method="linear")` with the effective-sample-size correction for serial autocorrelation; significance means the 95% confidence interval excludes zero. No trend was fitted across the 2014/2015 join.

### Mean state and WOA23 comparison
Area-weighted global means for 1971–2000 and native-grid differences from WOA23:

| Product | Global-mean SSS (psu) | Difference from WOA23 (psu) |
|---|---:|---:|
| WOA23 | 34.63 | — |
| CESM2 | 34.32 | -0.32 |
| GFDL-ESM4 | 34.56 | -0.07 |
| MPI-ESM1-2-LR | 34.36 | -0.27 |

On the common WOA23 1° comparison grid (40,829 cells), the WOA23 masked mean is 34.675 psu. Model masked means and metrics are:

| Model | Common-mask mean model (psu) | Common-mask mean WOA (psu) | Common-mask bias (psu) | RMSE (psu) | Spatial Pearson r |
|---|---:|---:|---:|---:|---:|
| CESM2 | 34.37 | 34.68 | -0.31 | 0.77 | 0.903 |
| GFDL-ESM4 | 34.64 | 34.68 | -0.03 | 0.60 | 0.926 |
| MPI-ESM1-2-LR | 34.41 | 34.68 | -0.27 | 0.76 | 0.893 |

The global-mean difference and the common-mask bias differ slightly because they are computed on different grids/masks. GFDL-ESM4 has the smallest mean-state difference, lowest RMSE and highest pattern correlation.

### Historical and ssp245 global-mean SSS trends
Slopes are in psu per decade, with 95% confidence intervals. “Significant” means the interval excludes zero.

| Model | Historical 1950–2014 slope [95% CI] | R² | Significant? | ssp245 2071–2100 slope [95% CI] | R² | Significant? | Ratio ssp/hist |
|---|---:|---:|---|---:|---:|---|---:|
| CESM2 | -2.68×10⁻³ [-6.93×10⁻³, +1.57×10⁻³] | 0.245 | No | -2.09×10⁻² [-2.55×10⁻², -1.64×10⁻²] | 0.910 | Yes | 7.8 |
| GFDL-ESM4 | -1.10×10⁻³ [-3.61×10⁻³, +1.41×10⁻³] | 0.064 | No | -5.34×10⁻³ [-1.06×10⁻², -3.85×10⁻⁵] | 0.324 | Yes | 4.9 |
| MPI-ESM1-2-LR | -5.51×10⁻⁴ [-5.04×10⁻³, +3.93×10⁻³] | 0.007 | No | -2.31×10⁻² [-5.13×10⁻², +5.04×10⁻³] | 0.701 | No | 42.0 |

The historical global-mean trends are not statistically significant for any model after accounting for autocorrelation. In ssp245, CESM2 and GFDL-ESM4 show significant negative trends; MPI-ESM1-2-LR does not, despite a large fitted slope, because its effective sample size is small.

The ratios are **not robust acceleration estimates**. The two trend windows are separated by a 56-year gap, so the ratio is only a scale comparison. For CESM2 it also mixes different realisations: historical r1i1p1f1 versus ssp245 r10i1p1f1 (parent r10i1p1f1). GFDL-ESM4 and MPI-ESM1-2-LR use continuous r1→r1 trajectories, making the comparison cleaner, but the MPI ratio is especially ill-conditioned because its historical slope is indistinguishable from zero.

### End-of-century ΔSSS: ssp245 2071–2100 minus historical 1995–2014
Values are area-weighted mean ΔSSS in psu, rounded to 2 decimal places.

| Model | Global | Atlantic | Indian | Pacific |
|---|---:|---:|---:|---:|
| CESM2 | -0.12 | +0.11 | +0.06 | -0.24 |
| GFDL-ESM4 | -0.01 | +0.19 | 0.00 | -0.04 |
| MPI-ESM1-2-LR | -0.09 | -0.04 | -0.05 | -0.11 |

Basins were selected with `select_basin(..., include_marginal_seas=True, depth_dependent=False)` using the WOA23 `basinmask_01`. Marginal seas were folded into their parent ocean; the Southern and Arctic Oceans were excluded as separate basins. The WOA23 basin mask was used, not the OMIP Ofx `basin` field. Some marginal seas in the WOA mask may have no model counterpart.

**Member caveat for ΔSSS:** For CESM2, the historical and ssp245 experiments are different realisations (r1i1p1f1 versus r10i1p1f1; no CESM2 historical r10 in the local archive), so its ΔSSS includes forced change plus an internal-variability difference. For GFDL-ESM4 and MPI-ESM1-2-LR the two experiments are the same realisation, so their ΔSSS is a cleaner forced-change estimate.

### Figures generated
- **Figure 1:** three global maps in one row of SSS bias (model 1971–2000 climatology minus WOA23) on the common 40,829 comparison cells. All panels share one diverging `RdBu_r` colour scale centred on zero, clipped at the 98th percentile of pooled |bias|, with colour bar labelled in psu. Panel titles include the spatial correlation and RMSE.
- **Figure 2:** three global maps in one row of the cell-by-cell linear trend of annual-mean SSS over 1950–2014 in psu per decade. The trend and significance fields were regridded bilinearly to a regular longitude–latitude frame for display only. Panels share one diverging colour scale centred on zero, clipped at the 98th percentile of pooled |trend|. Cells with p < 0.05 are stippled over the trend field.
- **Figure 3:** a single time-series panel of the global-mean SSS of all three models. Historical 1950–2014 and ssp245 2071–2100 segments are drawn separately, leaving the 2015–2070 gap empty/shaded. Each segment has its own dashed fitted linear trend line. The y-axis is labelled in psu, and each model uses one colour consistently across all three figures.

### Key insights
- Global-mean SSS is slightly fresh relative to WOA23 in all three models, by 0.07–0.32 psu, but spatial pattern agreement is high.
- GFDL-ESM4 performs best against WOA23 in this comparison (smallest global-mean difference, smallest RMSE, highest correlation).
- Historical global-mean SSS trends over 1950–2014 are weak and not significant after autocorrelation correction.
- Late-century ssp245 trends are more clearly negative for CESM2 and GFDL-ESM4, but MPI-ESM1-2-LR is not significant under the effective-sample-size test.
- Global ΔSSS is near zero to slightly negative, while basin responses differ: the Atlantic tends positive in CESM2/GFDL-ESM4, and the Pacific is negative in all three models, most strongly in CESM2.
