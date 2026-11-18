## Assessment

Using the common WOA23 comparison mask, the three CMIP6 models are all slightly fresh relative to WOA23 in 1971–2000, but by very different amounts: GFDL‑ESM4 is closest to WOA23 (global mean difference −0.01 psu), while CESM2 and MPI‑ESM1‑2‑LR are fresher by −0.45 and −0.28 psu. RMSE against WOA23 is about 0.92–1.06 psu and spatial pattern correlations are 0.798–0.841.

The historical 1950–2014 global-mean SSS trends are small and not significantly different from zero in all three models. The ssp245 2071–2100 trends are negative and significantly different from zero in all three models. The end-of-century delta SSS is negative globally for CESM2 and MPI‑ESM1‑2‑LR and near zero for GFDL‑ESM4, with substantial basin compensation.

## Methodology Summary

- **WOA23 reference:** The objectively analysed salinity field `s_an` was read from the 1971–2000 climate normal file, using time index 0. Levels at 0, 5 and 10 m were selected. Level weights were full cell volumes computed as WOA analytic 1° cell areas multiplied by full level thicknesses from `depth_bnds`, with no clipping to the 0–10 m interval. SSS was the volume-weighted mean over those levels. The land–sea mask used `landsea_01 > 1`. The basin mask used the surface WOA basin labels; higher-number marginal labels were folded to the nearest major basin. Salinity is practical salinity in psu.
- **CMIP6 processing:** Monthly `so` was taken for historical 1950–2014 and ssp245 2071–2100 for CESM2, GFDL‑ESM4 and MPI‑ESM1‑2‑LR. Model levels within 0–10 m were selected using the model level coordinate, with CESM2 levels converted from cm to m (`depth_scale=0.01`). Each selected level was weighted by its full `volcello` cell volume, without clipping to the 0–10 m overlap. Monthly SSS was formed and then averaged to annual means. Each model’s own `areacello` and `volcello` were opened once and used for both experiments; non-positive areas/volumes were masked, and CESM2 additionally used `sftof > 0`.
- **Regridding and comparison mask:** Model fields were mapped to the WOA23 1° longitude–latitude grid with a nearest-neighbour index built on Cartesian coordinates on the sphere. The common comparison mask required WOA `landsea_01 > 1`, valid WOA SSS, and valid regridded 1971–2000 SSS for all three models. It contains **40,778 cells**. Bias, RMSE, pattern correlation and comparison global means were computed only on this common mask, with area weighting: regridded model `areacello` for models and WOA analytic cell areas for WOA.
- **Trends:** Annual global-mean SSS was computed on each model’s native grid using native `areacello`. Separate linear trends were fitted for historical 1950–2014 and ssp245 2071–2100. The 95% confidence intervals use an effective-sample-size correction based on the lag-1 autocorrelation of residuals. Cell-by-cell historical trends were also computed with the same correction and p-values for Figure 2. No trend was fitted across the 2015–2070 gap.
- **Delta SSS:** Delta was ssp245 2071–2100 mean minus historical 1995–2014 mean, on cells valid in both periods. Global and basin means were area-weighted with the model’s regridded `areacello`. Basins were defined from the WOA23 basin mask after folding marginal seas into their parent ocean; reported basins are Atlantic (code 1), Pacific (code 2) and Indian (code 3). Southern and Arctic Oceans were excluded.

## Mean-State SSS, 1971–2000

All values below are on the common comparison mask with area weighting.

| Field | Global mean SSS (psu) | Difference vs WOA23 (psu) | RMSE (psu) | Spatial Pearson r | Mask cells |
|---|---:|---:|---:|---:|---:|
| WOA23 | **34.66** | — | — | — | 40,778 |
| CESM2 | **34.21** | **−0.45** | **1.06** | **0.836** | 40,778 |
| GFDL‑ESM4 | **34.65** | **−0.01** | **0.92** | **0.841** | 40,778 |
| MPI‑ESM1‑2‑LR | **34.38** | **−0.28** | **1.02** | **0.798** | 40,778 |

GFDL‑ESM4 has the smallest global-mean bias and the lowest RMSE, with CESM2 and MPI‑ESM1‑2‑LR both fresher than WOA23. The spatial correlations are all high but the residual RMSE near 1 psu indicates regional SSS errors remain important.

## Global-Mean SSS Trends

Slopes are in psu per decade. Confidence intervals are corrected for serial correlation using the effective sample size. “Excludes zero” refers to whether the 95% CI excludes zero.

| Model | Period | Slope (psu/decade) | 95% CI | R² | Excludes zero? |
|---|---|---:|---:|---:|---|
| CESM2 | historical 1950–2014 | **−0.00268** | [−0.00687, 0.00152] | 0.245 | No |
| GFDL‑ESM4 | historical 1950–2014 | **−0.00110** | [−0.00347, 0.00128] | 0.064 | No |
| MPI‑ESM1‑2‑LR | historical 1950–2014 | **−0.000687** | [−0.00511, 0.00373] | 0.011 | No |
| CESM2 | ssp245 2071–2100 | **−0.0209** | [−0.0252, −0.0166] | 0.910 | Yes |
| GFDL‑ESM4 | ssp245 2071–2100 | **−0.00534** | [−0.0104, −0.000253] | 0.325 | Yes |
| MPI‑ESM1‑2‑LR | ssp245 2071–2100 | **−0.0229** | [−0.0452, −0.000525] | 0.684 | Yes |

Slope ratios, ssp245/historical, to one decimal place:

| Model | ssp245/historical slope ratio | Meaningfulness |
|---|---:|---|
| CESM2 | **7.8** | Not meaningful as an acceleration estimate; the periods are separate and CESM2 uses different members. |
| GFDL‑ESM4 | **4.9** | Not meaningful as an acceleration estimate; the periods are separate and the 2015–2070 gap is not bridged. |
| MPI‑ESM1‑2‑LR | **33.3** | Not meaningful as an acceleration estimate; the periods are separate and the 2015–2070 gap is not bridged. |

The historical trends are weak and not significantly nonzero. The ssp245 trends are negative and significant for all three models. The numerical ratios are dominated by the near-zero historical slopes and should not be interpreted as continuous 1950–2100 trends, projections, or direct acceleration estimates. For CESM2, historical is `r1i1p1f1` and ssp245 is `r10i1p1f1`, so both its delta SSS and slope ratio include internal-variability differences between members. GFDL‑ESM4 and MPI‑ESM1‑2‑LR use `r1i1p1f1` for both experiments, but the two periods remain separate and non-contiguous.

## End-of-Century Delta SSS

Delta is ssp245 2071–2100 mean minus historical 1995–2014 mean. Global and basin means are area-weighted on each model’s delta-valid cells, not the common comparison mask.

| Model | Global ΔSSS (psu) | Atlantic ΔSSS (psu) | Indian ΔSSS (psu) | Pacific ΔSSS (psu) | Delta-valid cells |
|---|---:|---:|---:|---:|---:|
| CESM2 | **−0.12** | **+0.10** | **+0.08** | **−0.23** | 42,523 |
| GFDL‑ESM4 | **−0.01** | **+0.17** | **−0.00** | **−0.08** | 44,564 |
| MPI‑ESM1‑2‑LR | **−0.10** | **−0.00** | **−0.06** | **−0.13** | 42,679 |

CESM2 and MPI‑ESM1‑2‑LR show global freshening by the end of the century, while GFDL‑ESM4 has a near-zero global change. Basin changes are not uniform: CESM2 and GFDL‑ESM4 have Atlantic salinification, and CESM2 and MPI‑ESM1‑2‑LR have Pacific freshening. GFDL‑ESM4 illustrates strong basin compensation, with an Atlantic increase of +0.17 psu and a Pacific decrease of −0.08 psu nearly canceling globally. The CESM2 delta and slope ratio carry the additional caveat of different historical and ssp245 members; GFDL‑ESM4 and MPI‑ESM1‑2‑LR use matched members but still compare separate periods.

## Figures Produced

- **Figure 1:** Three global maps in one row of 1971–2000 SSS bias, model minus WOA23, on the common comparison cells. A single diverging `RdBu_r` colour scale centred on zero is shared across the panels, clipped at the 98th percentile of absolute bias, with the colour bar labelled in psu. The common-mask cell count is 40,778.
- **Figure 2:** Three global maps in one row of cell-by-cell linear SSS trend over 1950–2014 in psu per decade. Model trend and p-value fields are placed on the regular WOA23 1° longitude–latitude frame for display. A single diverging colour scale centred on zero is shared, clipped at the 98th percentile of absolute trend. Cells with p < 0.05 are stippled over the trend field.
- **Figure 3:** One time-series panel of annual global-mean SSS for all three models. Historical 1950–2014 and ssp245 2071–2100 are drawn as separate segments, leaving the 2015–2070 gap empty. Each segment has its own fitted line overlaid, the y-axis is labelled in psu, and each model uses a consistent colour across all three figures.

## Key Insights

- **Mean-state accuracy:** GFDL‑ESM4 has the smallest global-mean SSS bias and the lowest RMSE, while CESM2 and MPI‑ESM1‑2‑LR are substantially fresher than WOA23. All models capture the large-scale spatial pattern reasonably well, with correlations near 0.8.
- **Historical trends are not significant:** Over 1950–2014, none of the three models has a global-mean SSS trend whose 95% CI excludes zero. The best estimates are small and negative.
- **End-century ssp245 trends are significant:** All three models show significant negative global-mean SSS trends over 2071–2100, but with very different magnitudes: CESM2 and MPI‑ESM1‑2‑LR are much steeper than GFDL‑ESM4.
- **Delta SSS is basin-dependent:** Global near-zero change in GFDL‑ESM4 masks compensating Atlantic salinification and Pacific freshening. CESM2 and MPI‑ESM1‑2‑LR show global freshening, dominated in CESM2 by the Pacific.
- **Caveats:** Metrics and Figure 1 are conditional on the 40,778-cell common comparison mask. Delta SSS uses each model’s delta-valid cells and WOA23 surface basin labels after folding marginal seas into parent basins; Southern and Arctic Oceans are excluded from the reported basin means. CESM2’s historical/ssp245 member mismatch must be considered when interpreting its delta and slope ratio.
