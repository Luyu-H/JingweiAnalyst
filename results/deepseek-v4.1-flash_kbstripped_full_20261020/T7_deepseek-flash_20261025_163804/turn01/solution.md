## Direct answer

All four 200–1000 m sources have a headline volume-weighted mean N:P below Redfield 16, but for different reasons and with very different fidelity. WOA23 is 14.15, GFDL-ESM4 is closest at 14.13, while CESM2 (11.83) and MPI-ESM1-2-LR (12.28) are substantially too low. The equal-weight NO3 = a·PO4 + b regressions also have slopes well below 16, with positive intercepts; GFDL-ESM4 best reproduces the observed NO3–PO4 scatter (R²=0.820 versus WOA23 R²=0.908), whereas CESM2 has a much weaker relationship (R²=0.369). The observed low-N:P volume (N:P<12) is modest globally and largest in the Indian Ocean (10.3%) versus the Pacific (8.1%); none of the models reproduces that basin distribution. GFDL-ESM4 gets the global N:P<12 percentage approximately right (7.3% versus 6.2%), but puts too little low-N:P water in the Indian Ocean (0.5%) and too much in the Pacific (14.9%). CESM2 and MPI-ESM1-2-LR strongly overestimate low-N:P water, especially in the Pacific.

## Methodology used by the executed code

- **WOA23:** loaded the annual objectively analysed nitrate and phosphate climatology, not a 1971–2000 normal, and restricted to 200–1000 m. WOA23 depth levels and layer boundaries came from `load_depth_levels`; WOA23 cell volumes were formed from the 1° cell areas (`load_areas`) times layer thickness. WOA23 nitrate and phosphate were already in µmol kg⁻¹.
- **CMIP6:** loaded historical `no3` and `po4` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, monthly 1971-01 to 2000-12, and restricted to 200–1000 m. Monthly fields were averaged with `temporal_aggregate` using day-weighted long-term means, not a simple time mean.
- **Units:** model concentrations were converted from mol m⁻³ to µmol kg⁻¹ using `transform_unit` with density `1.025` kg/L, i.e. the requested constant density of 1025 kg m⁻³. N:P was then formed from the 1971–2000 mean NO3 and PO4 fields.
- **Native-grid layer statistics:** volume-weighted means used each model’s own `volcello` from `load_volumes_for_cmip6`; WOA23 used area × thickness. The headline mean N:P is the volume-weighted mean of the cell-by-cell `compute_ratio` field. The “ratio of means” is the separately volume-weighted mean NO3 divided by the separately volume-weighted mean PO4.
- **Common-mask regression sample:** model NO3 and PO4 were regridded horizontally with bilinear interpolation to the WOA23 1° grid. A single common comparison mask was built with WOA23 as reference using `build_comparison_mask(..., ocean_only=True, reduce_other_dims="any")`, then applied to all four sources. The OLS fit was a plain equal-weight-per-cell fit of NO3 = a·PO4 + b; `fit_multivariate_relation` gave the same R² as the direct least-squares calculation.
- **Basin selection:** used the WOA23 `basinmask_01` via `load_basin_mask` and `select_basin` with `include_marginal_seas=True`, `depth_dependent=False`. Reported basins are Pacific, Atlantic and Indian; Southern and Arctic are separate basins and were excluded.
- **Threshold fractions:** computed with `compute_exceedance_fraction` using volume weights. The 14–18% was operationalized as |N:P − 16| < 2 with threshold 2.0, comparison `abs_less`, reference 16.0. The below-12 percentage used threshold 12.0, comparison `less`.

## 200–1000 m layer statistics

| Source | Headline mean N:P (volume-weighted mean of cell N:P) | Ratio of volume-weighted mean NO3 to mean PO4 | Mean NO3 (µmol kg⁻¹) | Mean PO4 (µmol kg⁻¹) | Volume with 14 < N:P < 18 |
|---|---:|---:|---:|---:|---:|
| WOA23 | **14.15** | **14.00** | 27.49 | 1.96 | **55.0%** (1.496×10¹⁷ m³) |
| CESM2 | **11.83** | **10.93** | 21.31 | 1.95 | **17.2%** (4.749×10¹⁶ m³) |
| GFDL-ESM4 | **14.13** | **13.65** | 24.47 | 1.79 | **74.6%** (2.137×10¹⁷ m³) |
| MPI-ESM1-2-LR | **12.28** | **11.77** | 20.38 | 1.73 | **0.5%** (1.487×10¹⁵ m³) |

All headline means are below Redfield 16. WOA23 is 14.15, GFDL-ESM4 is nearly identical at 14.13, while CESM2 and MPI-ESM1-2-LR are much lower. The 14–18 band is dominant in WOA23 and GFDL-ESM4, much smaller in CESM2, and almost absent in MPI-ESM1-2-LR.

## Equal-weight OLS regression on the common WOA23-grid sample

| Source | a (NO3 = a·PO4 + b) | b (µmol kg⁻¹) | R² | Common-mask cells |
|---|---:|---:|---:|---:|
| WOA23 | **13.03** | **1.71** | **0.908** | 855,988 |
| CESM2 | **6.17** | **9.02** | **0.369** | 748,981 |
| GFDL-ESM4 | **11.38** | **4.05** | **0.820** | 376,641 |
| MPI-ESM1-2-LR | **8.09** | **6.22** | **0.643** | 409,935 |

The observed WOA23 slope (13.03) is below Redfield 16 but the scatter is tight (R²=0.908). GFDL-ESM4 has the closest slope and highest model R², but still with a larger positive intercept. CESM2 and MPI-ESM1-2-LR have much lower slopes and larger intercepts, indicating a weaker and offset NO3–PO4 relationship.

## Basin mean N:P from the WOA23 basin mask

| Source | Pacific | Atlantic | Indian |
|---|---:|---:|---:|
| WOA23 | **13.57** | **15.55** | **13.73** |
| CESM2 | **10.58** | **13.55** | **11.18** |
| GFDL-ESM4 | **13.39** | **15.07** | **14.72** |
| MPI-ESM1-2-LR | **11.29** | **13.15** | **12.60** |

The observed order is Atlantic > Indian > Pacific, with the Atlantic clearly higher and the Pacific/Indian lower. All three models preserve that qualitative order, but CESM2 and MPI-ESM1-2-LR are too low in every basin. GFDL-ESM4 is close in the Pacific and Atlantic but too high in the Indian Ocean.

## Volume with N:P < 12

| Source | Global | Pacific | Indian |
|---|---:|---:|---:|
| WOA23 | **6.2%** (1.685×10¹⁶ m³) | **8.1%** (9.966×10¹⁵ m³) | **10.3%** (4.880×10¹⁵ m³) |
| CESM2 | **36.7%** (1.016×10¹⁷ m³) | **59.4%** (7.385×10¹⁶ m³) | **35.8%** (1.724×10¹⁶ m³) |
| GFDL-ESM4 | **7.3%** (2.083×10¹⁶ m³) | **14.9%** (1.929×10¹⁶ m³) | **0.5%** (2.593×10¹⁴ m³) |
| MPI-ESM1-2-LR | **22.5%** (6.133×10¹⁶ m³) | **40.0%** (4.884×10¹⁶ m³) | **15.0%** (7.031×10¹⁵ m³) |

WOA23 has a modest global low-N:P volume, with the Indian Ocean having more N:P<12 water than the Pacific. All models fail to reproduce that basin ordering: in every model the Pacific N:P<12 percentage exceeds the Indian percentage. GFDL-ESM4 matches the global percentage reasonably well but not the basin distribution; CESM2 and MPI-ESM1-2-LR overproduce low-N:P water substantially.

## Interpretation against Redfield and the observed low-N:P signal

- **WOA23:** below Redfield overall (14.15), with a relatively tight NO3–PO4 relationship and a modest N:P<12 volume. The low-N:P signal is stronger in the Indian and Pacific than in the Atlantic by basin mean.
- **CESM2:** clearly below Redfield (11.83), far below WOA23. It produces too much low-N:P water globally and in both Pacific and Indian basins, and its NO3–PO4 relationship is much weaker (R²=0.369). It does not reproduce the observed magnitude of the low-N:P volume.
- **GFDL-ESM4:** closest to WOA23 in headline mean N:P (14.13 versus 14.15), ratio of means (13.65 versus 14.00), global N:P<12 percentage (7.3% versus 6.2%), and Pacific/Atlantic basin means. However, it misses the Indian Ocean low-N:P volume strongly (0.5% versus 10.3%) and overproduces N:P<12 in the Pacific (14.9% versus 8.1%). It best captures the overall stoichiometry but not the full basin distribution of the low-N:P signal.
- **MPI-ESM1-2-LR:** below Redfield (12.28) and much lower than WOA23. It overproduces N:P<12 water globally and in the Pacific, and its basin means are all too low. It does not reproduce the observed low-N:P volume or its basin distribution.

Overall, the models reasonably capture that intermediate-water N:P is below Redfield, but only GFDL-ESM4 is close to the observed global mean and global low-N:P volume. None of the three models reproduces the observed basin pattern of low-N:P water, especially the relatively high Indian Ocean N:P<12 volume and the observed Indian > Pacific ordering for that threshold.

## Visualization description

- **Figure 1:** four panels of nitrate versus phosphate from the common WOA23-grid sample, ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR. Panels share identical x and y limits and ticks. Each panel shows a blue scatter of cells, a red fitted OLS line from the equal-weight regression, and a black dashed Redfield 16:1 line through the origin. Each panel is annotated with the fitted slope a, intercept b, R² and common-mask cell count.
- **Figure 2:** four longitude–latitude pcolormesh panels of the 200–1000 m vertical volume-weighted mean N:P minus 16, using a common longitude–latitude frame and one shared diverging colour scale centred on zero. WOA23 is shown on its 1° grid; model fields were bilinearly regridded to the WOA23 1° grid for display only. The shared symmetric colour limit was set from the 98th percentile of absolute N:P−16 values across the maps. The figure uses a plain longitude–latitude plotting frame, not a projected cartographic map with coastlines.
