## Direct answer

The analysis produced a complete comparison of sea-surface salinity (SSS) in CESM2, GFDL-ESM4 and MPI-ESM1-2-LR against WOA23 and between historical and ssp245 experiments. All SSS quantities are practical salinity in psu, with no factor-of-1000 rescaling. The WOA23 reference is the 1971–2000 annual climatology. Model SSS was formed from levels whose depths lie within 0–10 m, using full `volcello` weights for the models and full cell-volume weights for WOA23; global and basin means are area-weighted with `compute_basic_stat(op='area_weighted_mean')`. Bias/RMSE/pattern-correlation diagnostics use conservative regridding to the WOA23 1° grid and a common comparison mask.

## Mean-state SSS, 1971–2000

WOA23 global mean SSS: **34.63 psu**.

| Model | 1971–2000 global mean SSS (psu) | Difference from WOA23 (psu) | Common-mask mean bias (psu) | RMSE vs WOA23 (psu) | Spatial Pearson correlation |
|---|---:|---:|---:|---:|---:|
| CESM2 | 34.32 | −0.32 | −0.31 | 0.77 | 0.903 |
| GFDL-ESM4 | 34.56 | −0.07 | −0.03 | 0.60 | 0.926 |
| MPI-ESM1-2-LR | 34.36 | −0.27 | −0.27 | 0.76 | 0.893 |

The common comparison mask contains **40,829 cells**. GFDL-ESM4 has the smallest global-mean offset, the lowest RMSE and the highest spatial correlation. CESM2 and MPI-ESM1-2-LR are both slightly fresher than WOA23 in the global mean, with correlations around 0.89–0.90. The pattern correlation is area-weighted; no serial-correlation correction is applicable to a spatial correlation coefficient.

## Trends in annual global-mean SSS

Trends were fitted separately to the annual global-mean SSS series for historical 1950–2014 and ssp245 2071–2100. The 95% confidence intervals use the effective-sample-size correction for serial correlation. A trend is called significant only when its interval excludes zero.

| Model | Historical 1950–2014 slope (psu/decade) | 95% CI | R² | CI excludes zero? | ssp245 2071–2100 slope (psu/decade) | 95% CI | R² | CI excludes zero? | Ratio ssp/hist |
|---|---:|---:|---:|:---:|---:|---:|---:|:---:|---:|
| CESM2 | −2.68×10⁻³ | [−6.93×10⁻³, 1.57×10⁻³] | 0.245 | no | −2.09×10⁻² | [−2.55×10⁻², −1.64×10⁻²] | 0.910 | yes | 7.8 |
| GFDL-ESM4 | −1.10×10⁻³ | [−3.61×10⁻³, 1.41×10⁻³] | 0.0639 | no | −5.34×10⁻³ | [−1.06×10⁻², −3.86×10⁻⁵] | 0.324 | yes | 4.9 |
| MPI-ESM1-2-LR | −5.51×10⁻⁴ | [−5.04×10⁻³, 3.93×10⁻³] | 0.00718 | no | −2.31×10⁻² | [−5.13×10⁻², 5.04×10⁻³] | 0.701 | no | 42.0 |

All historical global-mean trends are statistically indistinguishable from zero at 95% after serial-correlation correction. Under ssp245, CESM2 and GFDL-ESM4 have significant negative global-mean trends; MPI-ESM1-2-LR has a large best-estimate negative slope but its corrected confidence interval includes zero, so it is not significant in this single-realization analysis.

The reported ssp/historical slope ratios are **not meaningful acceleration factors**. They compare a 65-year historical period with a 30-year late-century scenario period under different forcings and, for CESM2, different ensemble members. The two segments were never fitted as one continuous series across the 2014/2015 join.

## End-of-century SSS change: ssp245 2071–2100 minus historical 1995–2014

Basin means use the WOA23 basin labels via `select_basin` with `include_marginal_seas=True` and `depth_dependent=False`. Marginal seas are folded into their parent ocean; the Southern and Arctic Oceans are separate basins and are excluded from the Atlantic/Indian/Pacific values.

| Model | Global ΔSSS (psu) | Atlantic ΔSSS (psu) | Indian ΔSSS (psu) | Pacific ΔSSS (psu) |
|---|---:|---:|---:|---:|
| CESM2 | −0.12 | +0.11 | +0.06 | −0.24 |
| GFDL-ESM4 | −0.01 | +0.19 | 0.00 | −0.04 |
| MPI-ESM1-2-LR | −0.09 | −0.04 | −0.05 | −0.11 |

The GFDL Indian-basin value is **−0.004 psu before rounding**. In these single-realization calculations, the global-mean ΔSSS is slightly negative for all three models, but the basin responses differ: CESM2 and GFDL-ESM4 are positive in the Atlantic, while MPI-ESM1-2-LR is negative in all three reported basins. These changes should not be interpreted as robust forced signals without accounting for internal variability and ensemble spread.

**Ensemble-member caveat:** CESM2 ssp245 is r10i1p1f1, while only the historical r1 member is available locally. Thus the CESM2 ΔSSS contains forced change plus an internal-variability difference between two different members. GFDL-ESM4 and MPI-ESM1-2-LR use r1 for both experiments and are continuous.

## Visualizations

- **Figure 1:** Three global maps in one row, one per model, showing SSS bias (model 1971–2000 climatology minus WOA23) on the common comparison cells. A single diverging colour scale centred on zero is shared by all panels; the colour bar is labelled in psu. The plotted scale is clipped at the 98th percentile of absolute bias for display.
- **Figure 2:** Three global longitude–latitude maps of the cell-by-cell linear trend of annual-mean SSS over 1950–2014. The native model trend fields were conservatively regridded to the WOA23 1° grid for display only. A shared diverging colour scale centred on zero is used, with units psu/decade, and cells whose corrected trend is significant at p < 0.05 are stippled. The significance decision comes from `fit_trend_map` on the native grid and is carried to the display grid by nearest-neighbour mapping of the significance field.
- **Figure 3:** A single time-series panel of global-mean SSS for all three models. Historical 1950–2014 and ssp245 2071–2100 are drawn as separate segments with the 2015–2070 gap left empty. Each segment has its own fitted linear trend overlaid, the y axis is labelled in psu, and model colours are C0/C1/C2 within the panel.

## Important patterns and caveats

- GFDL-ESM4 agrees best with WOA23 in global mean, RMSE and spatial pattern correlation.
- The historical global-mean SSS trends over 1950–2014 are weak and not significant after serial-correlation correction for all three models.
- Under ssp245, the global-mean SSS trend is significant and negative for CESM2 and GFDL-ESM4; MPI-ESM1-2-LR is not significant despite a larger best-estimate negative slope.
- All reported global and basin means and ΔSSS values were computed on each model’s native grid with area-weighted statistics. Regridding was used only for the bias/RMSE/correlation comparison and for the Figure 2 display.
- The single-member design and the CESM2 historical/ssp245 member mismatch are the main limitations for interpreting ΔSSS and slope ratios.
