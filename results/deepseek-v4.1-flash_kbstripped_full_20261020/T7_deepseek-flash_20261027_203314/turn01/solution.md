**Direct answer.** Over 200–1000 m, all four sources have mean N:P below Redfield 16:1. The WOA23 observed headline N:P (volume-weighted mean of cell-by-cell N:P) is **14.18**; the models give **14.17** (GFDL-ESM4), **12.28** (MPI-ESM1-2-LR) and **11.83** (CESM2). GFDL-ESM4 reproduces the observed global mean and global low-N:P fraction best (bias −0.01; N:P<12 fraction 7.0% vs observed 5.7%), but **no model reproduces the observed low-N:P volume and its basin distribution**. CESM2 and MPI-ESM1-2-LR strongly overproduce low-N:P water, while GFDL-ESM4 nearly removes the Indian Ocean low-N:P signal.

**Methodology actually used.** CMIP6 historical `no3` and `po4` for CESM2, GFDL-ESM4 and MPI-ESM1-2-LR were loaded monthly for **1971-01 to 2000-12** with `load_regular_data` and depth range 200–1000 m. They were converted with `transform_unit` from mol m−3 to µmol kg−1 using the constant-density conversion argument `density=1.025`; WOA23 nitrate/phosphate were already µmol kg−1 and were not converted. Annual means were formed with `temporal_aggregate(freq='mean')` on `no3` and `po4`, and N:P was computed afterwards with `compute_ratio`, not from monthly ratios. Levels whose depth coordinate lies strictly inside 200–1000 m were retained with `sub_layer`; the 200 m and 1000 m interfaces themselves were excluded. This left 21 WOA23 levels, 20 CESM2 levels, 8 GFDL-ESM4 levels and 11 MPI-ESM1-2-LR levels. The WOA23 input was the annual objectively analyzed climatology on the 1° grid. Volume-weighted diagnostics used `compute_basic_stat(op='volume_weighted_mean')` and `compute_exceedance_fraction(weighting='volume')`; recorded volume sources were published `volcello` for each CMIP6 model and WOA23 area × layer thickness. The OLS fit used `fit_multivariate_relation(method='linear')` over all valid cells with equal cell weight, in µmol kg−1. Basin means used `select_basin(..., include_marginal_seas=True, depth_dependent=False)` for Pacific, Atlantic and Indian; Southern and Arctic basins were excluded from these three basin means but remain in the global numbers. The headline mean is the **volume-weighted mean of the cell-by-cell N:P ratio**; the ratio of the volume-weighted mean nitrate to the volume-weighted mean phosphate is reported separately because N:P is nonlinear.

### Key stoichiometry results

| Source | Headline mean N:P: mean of cell-by-cell N:P | Ratio of VW mean NO3 to VW mean PO4 | OLS slope a | OLS intercept b (µmol kg−1) | R² | N cells | % volume 14<N:P<18 |
|---|---:|---:|---:|---:|---:|---:|---:|
| WOA23 | **14.18** | 14.01 | 12.96 | 1.87 | 0.907 | 784,057 | 55.2 |
| CESM2 | **11.83** | 10.93 | 4.89 | 10.84 | 0.238 | 1,558,515 | 17.2 |
| GFDL-ESM4 | **14.17** | 13.69 | 11.29 | 4.23 | 0.818 | 315,011 | 75.3 |
| MPI-ESM1-2-LR | **12.28** | 11.77 | 8.59 | 5.64 | 0.678 | 330,410 | 0.5 |

All sources sit **below Redfield 16:1** in the 200–1000 m mean. The WOA23 observed mean is 14.18; GFDL-ESM4 is almost identical at 14.17, while CESM2 and MPI-ESM1-2-LR are substantially lower. The fitted slopes are also below 16 in every source, with positive intercepts; CESM2 has the weakest NO3–PO4 relationship (R² = 0.238), while WOA23 has the strongest (R² = 0.907).

### Basin means and low-N:P volume fractions

| Source | Pacific mean N:P | Atlantic mean N:P | Indian mean N:P | % volume N:P<12 global | % volume N:P<12 Pacific | % volume N:P<12 Indian |
|---|---:|---:|---:|---:|---:|---:|
| WOA23 | 13.61 | 15.59 | 13.75 | 5.7 | 7.3 | 10.0 |
| CESM2 | 10.58 | 13.55 | 11.18 | 36.7 | 59.4 | 35.8 |
| GFDL-ESM4 | 13.46 | 15.06 | 14.76 | 7.0 | 14.5 | 0.3 |
| MPI-ESM1-2-LR | 11.29 | 13.15 | 12.60 | 22.5 | 40.0 | 15.0 |

The observed low-N:P signal is modest globally (5.7%) and is larger in the Indian Ocean (10.0%) than in the Pacific (7.3%). CESM2 and MPI-ESM1-2-LR greatly overproduce N:P<12 water globally and especially in the Pacific. GFDL-ESM4 roughly matches the global low-N:P fraction (7.0%) and the global mean, but has too much low-N:P water in the Pacific (14.5%) and almost none in the Indian Ocean (0.3%). Thus the models do not reproduce the observed basin distribution of the low-N:P denitrification signal. The mean N:P ordering with Atlantic highest and Pacific/Indian lower is broadly present, but the N:P<12 volume distribution is wrong.

Volumes behind the shares: layer volume is 2.478×10^17 m³ (WOA23), 2.767×10^17 m³ (CESM2), 2.376×10^17 m³ (GFDL-ESM4) and 2.724×10^17 m³ (MPI-ESM1-2-LR). Global N:P<12 volume is 1.408×10^16, 1.016×10^17, 1.670×10^16 and 6.133×10^16 m³ respectively. Pacific N:P<12 volume is 8.204×10^15, 7.385×10^16, 1.560×10^16 and 4.884×10^16 m³; Indian N:P<12 volume is 4.315×10^15, 1.724×10^16, 1.276×10^14 and 7.031×10^15 m³.

### Model–observation agreement

Each model’s 200–1000 m column N:P was conservatively regridded to the WOA23 1° grid for comparison and display; the metrics use a shared comparison mask built with `build_comparison_mask(ocean_only=True)`. The mask retained **38,084 of 64,800** cells. The WOA23 mean on that shared sample is 14.179.

| Model | Bias (model − WOA23) | RMSE | Centred RMSE | Max absolute error and location | Compared mean N:P | n valid |
|---|---:|---:|---:|---|---:|---:|
| CESM2 | −2.32 | 3.60 | 2.75 | 24.21 at 23.5°N, 79.5°W | 11.86 | 38,084 |
| GFDL-ESM4 | −0.01 | 1.66 | 1.66 | 26.26 at 31.5°N, 17.5°E | 14.17 | 38,084 |
| MPI-ESM1-2-LR | −1.87 | 2.73 | 1.98 | 29.33 at 31.5°N, 17.5°E | 12.31 | 38,084 |

GFDL-ESM4 has by far the smallest global bias and a mean N:P almost identical to WOA23, but it still has large local errors and the wrong spatial low-N:P distribution. CESM2 is strongly biased low and has the largest RMSE. MPI-ESM1-2-LR is also biased low and overproduces low-N:P volume.

### Figures

- **Figure 1** is a four-panel scatter of nitrate (y) versus phosphate (x) in µmol kg−1 over 200–1000 m, ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR, with identical axis ranges and ticks. Each panel shows its own OLS fit, the Redfield 16:1 line through the origin, and annotations of slope a, intercept b and R². For drawing only, the scatter uses up to 20,000 randomly selected valid cells per source; the fitted statistics use all valid cells.
- **Figure 2** is a 2×2 map of the 200–1000 m mean N:P minus 16 for the four sources on a common longitude–latitude frame, using one shared diverging colour scale centred on zero. The model fields were put on the WOA23 1° grid for display only; all means, volume shares, basin statistics and comparison masks were computed on each model’s archive/native grid with its own area/volume information.

### Important insights and caveats

The observed low-N:P signal is small in total volume but basin-dependent. None of the three CMIP6 models reproduces both the magnitude and the Pacific/Indian distribution of N:P<12 water: CESM2 and MPI-ESM1-2-LR generate far too much low-N:P water, whereas GFDL-ESM4 matches the global mean and global low-N:P fraction but has an almost absent Indian Ocean low-N:P signal. CESM2 also contains slightly negative nitrate (minimum about −1.24 µmol kg−1) and negative N:P in some cells, which weakens its fit and inflates low-N:P fractions. The CMIP6 concentrations were placed on a constant-density basis for the regression; using in-situ density below the surface would change the µmol kg−1 concentrations and therefore the fitted a and b, while ratio-based quantities would be unchanged because the same density factor cancels.
