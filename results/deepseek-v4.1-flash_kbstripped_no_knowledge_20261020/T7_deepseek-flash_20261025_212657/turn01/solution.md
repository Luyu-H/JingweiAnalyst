## Direct answer

For the 200–1000 m layer, the WOA23 observed climatology and all three CMIP6 models have mean N:P below the Redfield value of 16. The observed headline mean cell-by-cell N:P is **14.15** (ratio of means **14.00**). Among the models, **GFDL-ESM4** is closest to WOA23 in global headline mean N:P (14.13), global low-N:P volume (7.3% vs 6.2% observed), and Pacific mean N:P (13.39 vs 13.57 observed). However, GFDL-ESM4 almost completely misses the observed Indian Ocean low-N:P signal (0.5% of volume with N:P < 12 versus 10.3% in WOA23; Indian mean 14.72 versus 13.73). CESM2 and MPI-ESM1-2-LR are both too low in mean N:P and substantially overproduce low-N:P volume, especially in the Pacific. Thus **no model fully reproduces the observed low-N:P volume and its basin distribution**; GFDL-ESM4 is closest globally and in the Pacific, but it fails to reproduce the Indian Ocean denitrification-related low-N:P signature, while CESM2 and MPI-ESM1-2-LR overdo the low-N:P signal.

## Methodology as executed

- **WOA23** nitrate and phosphate were loaded with `load_climatology_data(dataset='woa23', variable='nitrate'/'phosphate', time_slice='annual', depth_range=[200,1000])`, aligned, and used to compute cell-by-cell N:P. The executed code did not explicitly load WOA23 areas or depth-layer boundaries; the 200–1000 m restriction was imposed through `depth_range`.
- **CMIP6 historical** nitrate and phosphate were loaded with `load_regular_data(dataset='cmip6_historical', variable='no3'/'po4', model=..., time_range=['1971-01','2000-12'], depth_range=[200,1000])`. Monthly fields were aggregated with `temporal_aggregate(freq='mean', op='mean')`; the code did not pass explicit day weights. Concentrations were converted with `transform_unit(from_unit='mol m-3', to_unit='umol/kg', density=1.025)` as executed.
- Cell-by-cell N:P was computed only where both nutrients were valid and phosphate was nonzero.
- Headline and basin means were requested with `compute_basic_stat(..., op='volume_weighted_mean', dim=None)`. The executed code did not pass an explicit `weights` argument to those calls; the volume-weighted-mean operation was requested through the `op` argument. Exceedance fractions explicitly used `weighting='volume'`.
- The 14–18% volume statistic was computed as `compute_exceedance_fraction(..., threshold=2, comparison='abs_less', reference=16)`, which is equivalent to strict \(14 < \mathrm{N:P} < 18\). The low-N:P fraction used `threshold=12, comparison='less'`, i.e. strictly below 12.
- Basin statistics used `select_basin(..., basin=..., include_marginal_seas=True, depth_dependent=False)` on the N:P fields. The executed code did not explicitly call `load_basin_mask`, `load_volumes_for_cmip6`, `load_areas`, or `load_depth_levels`; it relied on `load_regular_data`/`load_climatology_data` with `depth_range` and on `select_basin`.
- OLS regressions were performed with `fit_multivariate_relation(target=NO3, predictors=[PO4], method='linear', validation='random')`, with each valid cell weighted equally. Units for both NO3 and PO4 are micromoles per kilogram.
- Figure 2 model depth-mean fields were regridded only for display with `regrid_horizontal(depth_maps[model], target='woa23', method='bilinear')`.

## Key numerical results

### Global 200–1000 m statistics

| Source | Mean N:P, cell-by-cell ratio (headline) | Mean N:P, ratio of means | OLS slope \(a\) | OLS intercept \(b\) (\(\mu\)mol kg\(^{-1}\)) | \(R^2\) | % volume \(14<\mathrm{N:P}<18\) | % volume N:P < 12 | Relative to Redfield (16) |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| WOA23 | **14.15** | **14.00** | **13.04** | **1.69** | **0.908** | **55.0%** | **6.2%** | below (−1.85) |
| CESM2 | **11.83** | **10.93** | **4.89** | **10.84** | **0.238** | **17.2%** | **36.7%** | below (−4.17) |
| GFDL-ESM4 | **14.13** | **13.65** | **11.37** | **4.02** | **0.817** | **74.6%** | **7.3%** | below (−1.87) |
| MPI-ESM1-2-LR | **12.28** | **11.77** | **8.59** | **5.64** | **0.678** | **0.5%** | **22.5%** | below (−3.72) |

The regression equation is \(\mathrm{NO3} = a\,\mathrm{PO4} + b\). The fitted slopes are all below the Redfield slope of 16, and the intercepts are positive. WOA23 has the strongest linear NO3–PO4 relationship (\(R^2=0.908\)), GFDL-ESM4 is intermediate (\(R^2=0.817\)), MPI-ESM1-2-LR is weaker (\(R^2=0.678\)), and CESM2 is much weaker (\(R^2=0.238\)) with a large positive intercept. This indicates that CESM2’s NO3–PO4 relationship is poorly described by a single linear slope over this layer.

For every source, the headline cell-by-cell mean N:P is slightly higher than the ratio of the mean nitrate to mean phosphate.

### Basin means and low-N:P volume fractions

| Source | Pacific mean N:P | Atlantic mean N:P | Indian mean N:P | % volume N:P < 12 in Pacific | % volume N:P < 12 in Indian |
|---|---:|---:|---:|---:|---:|
| WOA23 | **13.57** | **15.55** | **13.73** | **8.1%** | **10.3%** |
| CESM2 | **10.58** | **13.55** | **11.18** | **59.4%** | **35.8%** |
| GFDL-ESM4 | **13.39** | **15.07** | **14.72** | **14.9%** | **0.5%** |
| MPI-ESM1-2-LR | **11.29** | **13.15** | **12.60** | **40.0%** | **15.0%** |

WOA23 shows the lowest basin mean N:P in the Pacific (13.57) and Indian Ocean (13.73), with a higher Atlantic mean (15.55), and nonzero low-N:P volume in both the Pacific (8.1%) and Indian Ocean (10.3%). This is the observed basin pattern expected from denitrification-related low N:P.

## Interpretation against observations

- **WOA23**: Headline mean N:P **14.15**, ratio of means **14.00**, below Redfield by about 1.85. The low-N:P signal is concentrated in the Pacific and Indian basins, consistent with a denitrification influence.
- **GFDL-ESM4**: Global mean N:P **14.13**, almost identical to WOA23, and global low-N:P volume **7.3%** is close to the observed **6.2%**. Its Pacific mean (**13.39**) is also close to WOA23 (**13.57**). However, the Pacific low-N:P volume is higher than observed (**14.9%** versus **8.1%**), while the Indian low-N:P volume is nearly absent (**0.5%** versus **10.3%**) and the Indian mean is too high (**14.72** versus **13.73**). Thus GFDL-ESM4 reproduces the global mean and part of the Pacific signal but does not reproduce the observed Indian Ocean low-N:P distribution.
- **CESM2**: Mean N:P is far too low (**11.83** versus **14.15**), global low-N:P volume is much too large (**36.7%** versus **6.2%**), and Pacific/Indian low-N:P fractions are strongly overproduced (**59.4%** and **35.8%** versus **8.1%** and **10.3%**). Basin means are also too low in all three basins. CESM2 does not reproduce the observed low-N:P volume or basin distribution; it overdoes the low-N:P signal.
- **MPI-ESM1-2-LR**: Mean N:P is too low (**12.28**), global low-N:P volume is too large (**22.5%** versus **6.2%**), and Pacific/Indian low-N:P fractions are too high (**40.0%** and **15.0%** versus **8.1%** and **10.3%**). Basin means are lower than observed. MPI-ESM1-2-LR also overproduces low-N:P water and does not match the observed basin distribution.

Overall, **all sources are below Redfield (N:P = 16)**. The closest match to the observed global mean is GFDL-ESM4, but none of the models fully reproduces the observed low-N:P volume and its basin-specific distribution. GFDL-ESM4 captures the global low-N:P fraction and Pacific mean reasonably well but misses the Indian Ocean low-N:P signal; CESM2 and MPI-ESM1-2-LR overproduce low-N:P water, especially in the Pacific.

## Visualization descriptions

**Figure 1** is a 2×2 scatter plot of nitrate (y-axis) versus phosphate (x-axis) in \(\mu\)mol kg\(^{-1}\) for the 200–1000 m layer, with panels ordered WOA23, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR. All panels use identical axis limits and ticks, set from pooled data approximately to phosphate 0–3.7 and nitrate 0–65. Each panel shows the data as gray points, the panel-specific OLS regression line, and the Redfield 16:1 line through the origin with distinguishable styles. Each panel is annotated with the fitted slope \(a\), intercept \(b\), and \(R^2\). For readability, the plotted scatter uses at most 5000 subsampled points per panel, while the reported regressions use all valid cells.

**Figure 2** is a 2×2 global map of the 200–1000 m mean N:P minus 16 for WOA23, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR. Each panel is on a common longitude–latitude geographic frame. The WOA23 field is already regular; the CMIP6 depth-mean fields were regridded to the WOA23 regular grid with bilinear interpolation for display only, while all statistics remained on native grids. All panels share one diverging `RdBu_r` color scale centered at zero, with the color limit clipped at about ±9.98 based on the 98th percentile of absolute values. The maps are intended to show where the 200–1000 m mean N:P lies above or below Redfield and how those spatial patterns differ among the observations and models.
