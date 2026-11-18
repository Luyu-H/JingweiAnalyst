# Tool verification report

Run directory: `results/tool_verification`

**591 metrics**: 498 pass, 16 fail, 0 error, 77 info (reported against a stated range, no threshold).

## Per tool (pass/fail/error, info in brackets)

| tool | 1 analytic | 2 reference impl. | 3 published products | 4 magnitudes | 5 consistency | 6 Monte Carlo | 7 safeguards | verdict |
|---|---|---|---|---|---|---|---|---|
| `analyze_correlation` |  | 12/0/0 |  |  |  | 2/0/0 (2) | 3/0/0 | pass |
| `build_comparison_mask` |  |  |  |  |  |  | 5/0/0 | pass |
| `compute_anomaly_vs_own_climatology` | 4/0/0 | 2/0/0 |  |  |  |  | 4/0/0 | pass |
| `compute_anomaly_vs_woa23` | 9/0/0 |  |  |  |  |  | 4/0/0 | pass |
| `compute_basic_stat` | 16/0/0 | 3/3/0 |  |  | 4/0/0 |  | 8/0/0 | FAIL |
| `compute_difference_metrics` | 8/0/0 | 7/0/0 |  |  | 3/0/0 |  | 2/0/0 | pass |
| `compute_enso_nino_index` | 4/0/0 | 4/0/0 | 0/12/0 (8) |  | 2/0/0 |  | 4/0/0 | FAIL |
| `compute_exceedance_fraction` | 8/0/0 |  |  |  |  |  | 2/0/0 | pass |
| `compute_mixed_layer_depth` | 16/0/0 |  | 1/0/0 (4) |  |  |  | 3/0/0 | pass |
| `compute_nutricline_depth` | 4/0/0 |  |  |  |  |  | 1/0/0 | pass |
| `compute_pdo_index` | 4/0/0 (1) | 5/0/0 | 3/0/0 (9) |  |  |  | 3/0/0 | pass |
| `compute_potential_density_anomaly` |  | 20/0/0 (1) |  |  |  |  | 3/0/0 | pass |
| `compute_ratio` | 3/0/0 |  |  |  |  |  | 1/0/0 | pass |
| `compute_vertical_gradient` | 6/0/0 (1) |  |  |  |  |  | 2/0/0 | pass |
| `convert_concentration_to_content` | 2/0/0 |  |  |  |  |  | 2/0/0 | pass |
| `detect_outliers` | 4/0/0 | 4/0/0 |  |  |  | 4/0/0 | 3/0/0 | pass |
| `detect_threshold_events` | 5/0/0 |  | 4/0/0 |  |  |  | 2/0/0 | pass |
| `fit_1D_trend` | 7/0/0 | 12/0/0 |  |  |  | 3/1/0 (3) | 3/0/0 | FAIL |
| `fit_multivariate_relation` |  | 9/0/0 |  |  |  |  | 2/0/0 | pass |
| `fit_trend_map` | 4/0/0 | 9/0/0 |  |  |  | 6/0/0 (10) | 2/0/0 | pass |
| `identify_oxygen_minimum_zone` | 3/0/0 |  |  | 0/0/0 (4) | 1/0/0 |  | 1/0/0 | pass |
| `load_areas` | 9/0/0 |  |  | 5/0/0 (3) |  |  |  | pass |
| `load_basin_mask` |  | 2/0/0 |  |  | 3/0/0 |  | 1/0/0 | pass |
| `load_climatology_data` |  | 9/0/0 |  |  |  |  | 2/0/0 | pass |
| `load_data_with_geolocator` |  | 10/0/0 (10) |  |  | 0/0/0 (1) |  | 3/0/0 | pass |
| `load_depth_levels` |  | 7/0/0 |  |  | 12/0/0 |  | 2/0/0 | pass |
| `load_landsea_mask` |  | 1/0/0 |  |  | 1/0/0 |  |  | pass |
| `load_regular_data` |  | 20/0/0 |  |  | 2/0/0 |  | 8/0/0 | pass |
| `load_volumes_for_cmip6` |  |  |  | 6/0/0 (3) |  |  | 2/0/0 | pass |
| `locate_extreme` | 3/0/0 |  |  |  |  |  | 2/0/0 | pass |
| `regrid_horizontal` | 15/0/0 (4) | 2/0/0 (1) |  |  |  |  | 6/0/0 | pass |
| `regrid_vertical` | 4/0/0 | 1/0/0 (1) |  |  |  |  | 4/0/0 | pass |
| `select_basin` |  |  |  |  | 7/0/0 (6) |  | 4/0/0 | pass |
| `temporal_aggregate` | 35/0/0 | 8/0/0 | 0/0/0 (4) |  |  |  | 5/0/0 | pass |
| `transform_unit` | 7/0/0 |  |  |  |  |  | 3/0/0 (1) | pass |

## Failures and errors (16)

- **compute_basic_stat** [cat2_reference] area- and volume-weighted means on model grids vs numpy with the tools' own weights — CESM2 tos: global mean vs float64 numpy with load_areas: value 3.902e-05, threshold <= 1e-06. float32 inputs: the tool accumulates in float32
- **compute_basic_stat** [cat2_reference] area- and volume-weighted means on model grids vs numpy with the tools' own weights — GFDL-ESM4 tos: global mean vs float64 numpy with load_areas: value 1.328e-05, threshold <= 1e-06. float32 inputs: the tool accumulates in float32
- **compute_basic_stat** [cat2_reference] area- and volume-weighted means on model grids vs numpy with the tools' own weights — MPI-ESM1-2-LR tos: global mean vs float64 numpy with load_areas: value 3.703e-06, threshold <= 1e-06. float32 inputs: the tool accumulates in float32
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO1+2 TOTAL (box mean SST, edge cells included): share of months within 0.01 degC: value 0.05263, threshold >= 0.95. 912 months; RMSE 0.1341, max |d| 0.320 (2017-09), mean d +0.0921
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO1+2 ANOM, 1991-2020 base: share of months within 0.01 degC: value 0.1557, threshold >= 0.95. 912 months; RMSE 0.0517, max |d| 0.203 (1997-11), mean d +0.0033
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO1+2 ANOM, 1981-2010 base: share of months within 0.01 degC: value 0.162, threshold >= 0.95. 852 months; RMSE 0.0527, max |d| 0.195 (1997-11), mean d +0.0084
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO3 TOTAL (box mean SST, edge cells included): share of months within 0.01 degC: value 0.4243, threshold >= 0.95. 912 months; RMSE 0.0150, max |d| 0.037 (2024-12), mean d +0.0100
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO3 ANOM, 1991-2020 base: share of months within 0.01 degC: value 0.8059, threshold >= 0.95. 912 months; RMSE 0.0076, max |d| 0.029 (1983-01), mean d -0.0004
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO3 ANOM, 1981-2010 base: share of months within 0.01 degC: value 0.7969, threshold >= 0.95. 852 months; RMSE 0.0082, max |d| 0.076 (2020-12), mean d +0.0015
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO3.4 TOTAL (box mean SST, edge cells included): share of months within 0.01 degC: value 0.7675, threshold >= 0.95. 912 months; RMSE 0.0084, max |d| 0.027 (2018-08), mean d +0.0017
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO3.4 ANOM, 1991-2020 base: share of months within 0.01 degC: value 0.807, threshold >= 0.95. 912 months; RMSE 0.0075, max |d| 0.023 (1992-01), mean d -0.0028
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO3.4 ANOM, 1981-2010 base: share of months within 0.01 degC: value 0.8627, threshold >= 0.95. 852 months; RMSE 0.0069, max |d| 0.029 (2020-12), mean d -0.0011
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO4 TOTAL (box mean SST, edge cells included): share of months within 0.01 degC: value 0.5768, threshold >= 0.95. 912 months; RMSE 0.0119, max |d| 0.033 (1991-10), mean d -0.0058
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO4 ANOM, 1991-2020 base: share of months within 0.01 degC: value 0.7204, threshold >= 0.95. 912 months; RMSE 0.0091, max |d| 0.029 (2021-02), mean d +0.0008
- **compute_enso_nino_index** [cat3_products] CPC ERSSTv5 Nino indices: area-mean SST (TOTAL) and anomalies on two baselines — NINO4 ANOM, 1981-2010 base: share of months within 0.01 degC: value 0.7171, threshold >= 0.95. 852 months; RMSE 0.0090, max |d| 0.033 (2018-07), mean d +0.0010
- **fit_1D_trend** [cat6_montecarlo] AICc model selection on a positive series generated with a linear trend — linear truth, positive series: share selected as linear: value 0.56, threshold >= 0.7. {'quadratic': 15, 'linear': 112, 'sinusoidal': 5, 'exponential': 68}

## Reported values (77)

- **compute_vertical_gradient** [cat1_analytic] smooth profile, levels below 500 m: max relative error vs dc/dz: 0.06803 — report only: coarse 50-100 m levels
- **compute_pdo_index** [cat1_analytic] explained variance: 97.38
- **regrid_horizontal** [cat1_analytic] CESM2: source cells' share of the sphere (sum of areas / 4 pi - 1): -0.009497 — the grid stops short of the South Pole
- **regrid_horizontal** [cat1_analytic] CESM2: ocean cells, error vs centre values (conservative / bilinear p99): 2.22e-03 / 1.85e-04 — first-order conservative onto a finer grid is O(grad f * dx_source), not O(dx^2); bilinear max over all cells 1.35e-03
- **regrid_horizontal** [cat1_analytic] MPI-ESM1-2-LR: source cells' share of the sphere (sum of areas / 4 pi - 1): -0.01261 — the grid stops short of the South Pole
- **regrid_horizontal** [cat1_analytic] MPI-ESM1-2-LR: ocean cells, error vs centre values (conservative / bilinear p99): 8.71e-03 / 4.70e-04 — first-order conservative onto a finer grid is O(grad f * dx_source), not O(dx^2); bilinear max over all cells 1.82e-01
- **compute_potential_density_anomaly** [cat2_reference] |sigma0 error| from treating thetao as in-situ, at 2000 / 4000 m: 0.0116 / 0.0215 kg m-3 — docstring states about 0.01 at 2000 m and 0.026 at 4000 m
- **regrid_horizontal** [cat2_reference] real WOA23 SST with land: coastal cells filled only by the tool: 272 — a great-circle edge reaches slivers of a neighbouring wet cell that parallel edges do not
- **regrid_vertical** [cat2_reference] cells defined by numpy.interp but missing in the tool: 0 — nonzero only where a column has interior gaps
- **load_data_with_geolocator** [cat2_reference] Mediterranean Sea: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 0.9796 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Gulf of Mexico: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 0.9752 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Caribbean Sea: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 0.968 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Bay of Bengal: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 0.9765 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Red Sea: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 1 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Persian Gulf: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 1 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Barents Sea: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 0.9858 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Bering Sea: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 0.9795 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Gulf of Guinea: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 0.9538 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **load_data_with_geolocator** [cat2_reference] Arabian Sea: diagnostic - Jaccard against cells displaced +0.5 deg in longitude: 0.952 — 1.0 here means the mask was built on [c, c+1] instead of [c-0.5, c+0.5]
- **compute_enso_nino_index** [cat3_products] NINO1+2: box-edge convention, RMSE with edges included / excluded: 0.1341 / 0.2023 — 6x6 vs 4x4 cells
- **compute_enso_nino_index** [cat3_products] NINO1+2: diagnosis - box mean weighting each cell by its fraction inside the box: RMSE 0.0141, within 0.01: 0.85 — the tool gives a cell centred on the box edge full weight; CPC's values match partial (fractional-overlap) weighting far better
- **compute_enso_nino_index** [cat3_products] NINO3: box-edge convention, RMSE with edges included / excluded: 0.0150 / 0.0158 — 5x31 vs 5x29 cells
- **compute_enso_nino_index** [cat3_products] NINO3: diagnosis - box mean weighting each cell by its fraction inside the box: RMSE 0.0037, within 0.01: 0.98 — the tool gives a cell centred on the box edge full weight; CPC's values match partial (fractional-overlap) weighting far better
- **compute_enso_nino_index** [cat3_products] NINO3.4: box-edge convention, RMSE with edges included / excluded: 0.0084 / 0.0089 — 5x26 vs 5x24 cells
- **compute_enso_nino_index** [cat3_products] NINO3.4: diagnosis - box mean weighting each cell by its fraction inside the box: RMSE 0.0045, within 0.01: 0.95 — the tool gives a cell centred on the box edge full weight; CPC's values match partial (fractional-overlap) weighting far better
- **compute_enso_nino_index** [cat3_products] NINO4: box-edge convention, RMSE with edges included / excluded: 0.0119 / 0.0137 — 5x26 vs 5x24 cells
- **compute_enso_nino_index** [cat3_products] NINO4: diagnosis - box mean weighting each cell by its fraction inside the box: RMSE 0.0052, within 0.01: 0.94 — the tool gives a cell centred on the box edge full weight; CPC's values match partial (fractional-overlap) weighting far better
- **temporal_aggregate** [cat3_products] ONI, fixed 1991-2020 base: share of months within 0.01 degC: 0.03183 — 911 months; RMSE 0.2711, max |d| 0.925 (1954-11), mean d -0.1942; report only: early decades sit on a cooler base
- **temporal_aggregate** [cat3_products] ONI, CPC rolling 30-year bases: share of months within 0.01 degC: 0.08123 — 911 months; RMSE 0.0952, max |d| 0.441 (1952-06), mean d +0.0120; report only; bases 1936-65 for 1950-55, then moved every 5 years, latest 1996-2025
- **temporal_aggregate** [cat3_products] diagnosis: published ONI rebuilt from CPC's own monthly ERSSTv5 Nino3.4 (no tool involved): RMSE 0.0853, max 0.360 — 1966-2025, documented base rule; a large residual means the ONI file rests on another SST input or vintage, so the tool cannot be judged against it
- **temporal_aggregate** [cat3_products] ONI rolling base: max |d| by decade: {"1950": 0.441, "1960": 0.303, "1970": 0.351, "1980": 0.237, "1990": 0.291, "2000": 0.234, "2010": 0.194, "2020": 0.121}
- **compute_pdo_index** [cat3_products] 5-year running-mean correlation: 0.9428
- **compute_pdo_index** [cat3_products] share of months with the same sign: 0.8915
- **compute_pdo_index** [cat3_products] RMS difference of the standardised series: 0.3938
- **compute_pdo_index** [cat3_products] phase 1925-1946: mean index, tool / NCEI: +0.38 / +0.73
- **compute_pdo_index** [cat3_products] phase 1947-1976: mean index, tool / NCEI: -0.43 / -0.27
- **compute_pdo_index** [cat3_products] phase 1977-1998: mean index, tool / NCEI: +0.62 / +0.66
- **compute_pdo_index** [cat3_products] phase 1999-2013: mean index, tool / NCEI: -0.36 / -0.55
- **compute_pdo_index** [cat3_products] phase 2014-2016: mean index, tool / NCEI: +1.06 / +1.01
- **compute_pdo_index** [cat3_products] phase 2020-2025: mean index, tool / NCEI: -1.17 / -1.78
- **compute_mixed_layer_depth** [cat3_products] area-weighted geometric-mean ratio tool / SEANOE: 0.6531 — below 1 expected: monthly-mean profiles are smoother than individual ones
- **compute_mixed_layer_depth** [cat3_products] Northern winter (JFM, 40-70N): median MLD tool / SEANOE (m): 67.8 / 100.8
- **compute_mixed_layer_depth** [cat3_products] Southern winter (JAS, 40-70S): median MLD tool / SEANOE (m): 98.5 / 129.8
- **compute_mixed_layer_depth** [cat3_products] tropics (20S-20N, all months): median MLD tool / SEANOE (m): 19.6 / 31.0
- **load_areas** [cat4_magnitude] WOA23: area of cells the land-sea mask calls ocean / 3.61e14: 0.9677 — 3.4933e+14 m2; a binary 1-degree mask has no partial coastal cells, so it is not held to the +-1% meant for a model's fractional areacello
- **load_areas** [cat4_magnitude] 1-degree grid total / WGS84 surface area: 1 — short only by the two 0.01-degree polar caps the area solver clips
- **load_areas** [cat4_magnitude] MPI-ESM1-2-LR: the same with its repeated columns [] counted once: 1.001
- **load_volumes_for_cmip6** [cat4_magnitude] CESM2: rebuilt / published global volume: 1 — docstring: GFDL-ESM4 1.53e18 vs 1.33e18; MPI partial bottom cells; CESM2 agree
- **load_volumes_for_cmip6** [cat4_magnitude] GFDL-ESM4: rebuilt / published global volume: 1 — docstring: GFDL-ESM4 1.53e18 vs 1.33e18; MPI partial bottom cells; CESM2 agree
- **load_volumes_for_cmip6** [cat4_magnitude] MPI-ESM1-2-LR: rebuilt / published global volume: 1 — docstring: GFDL-ESM4 1.53e18 vs 1.33e18; MPI partial bottom cells; CESM2 agree
- **identify_oxygen_minimum_zone** [cat4_magnitude] share of the WOA23 ocean volume with O2 < 20 umol/kg: 1.033 — 1.417e+16 m3; literature puts the < 20 umol/kg core at the order of 1% of the ocean
- **identify_oxygen_minimum_zone** [cat4_magnitude] share of the WOA23 ocean volume with O2 < 60 umol/kg: 5.854 — 8.033e+16 m3; literature puts the < 20 umol/kg core at the order of 1% of the ocean
- **identify_oxygen_minimum_zone** [cat4_magnitude] eastern tropical North Pacific: maximum OMZ (<20) thickness (m): 1088 — these are the two thickest OMZs of the world ocean
- **identify_oxygen_minimum_zone** [cat4_magnitude] Arabian Sea: maximum OMZ (<20) thickness (m): 1212 — these are the two thickest OMZs of the world ocean
- **select_basin** [cat5_consistency] CESM2: ocean area left out of every basin (%): 3.5 — cells whose nearest 1-degree cell is land in the WOA23 mask
- **select_basin** [cat5_consistency] CESM2: basin area / WOA23 basin area: {"Atlantic": 0.999, "Pacific": 1.0003, "Indian": 1.0015, "Southern": 0.9832, "Arctic": 0.981}
- **select_basin** [cat5_consistency] GFDL-ESM4: ocean area left out of every basin (%): 3.464 — cells whose nearest 1-degree cell is land in the WOA23 mask
- **select_basin** [cat5_consistency] GFDL-ESM4: basin area / WOA23 basin area: {"Atlantic": 1.0001, "Pacific": 1.0021, "Indian": 1.0004, "Southern": 0.9917, "Arctic": 0.9842}
- **select_basin** [cat5_consistency] MPI-ESM1-2-LR: ocean area left out of every basin (%): 3.842 — cells whose nearest 1-degree cell is land in the WOA23 mask
- **select_basin** [cat5_consistency] MPI-ESM1-2-LR: basin area / WOA23 basin area: {"Atlantic": 0.9944, "Pacific": 0.9985, "Indian": 0.9962, "Southern": 0.9952, "Arctic": 0.9769}
- **load_data_with_geolocator** [cat5_consistency] Mediterranean Sea: wet area on each grid (10^6 km2): {"woa23": 1.447, "ersst": 1.765, "CESM2": 1.788, "GFDL-ESM4": 1.77, "MPI-ESM1-2-LR": 1.828} — Natural Earth polygon 2.5e6 km2 approx.; coarse grids gain coastal cells
- **fit_trend_map** [cat6_montecarlo] n=480, phi=0.9: false-positive rate (nominal 0.05): 0.06925 — without the correction 0.671
- **fit_trend_map** [cat6_montecarlo] n=480, phi=0.9: 95% interval coverage: 0.9337
- **fit_trend_map** [cat6_montecarlo] n=120, phi=0.0: false-positive rate (nominal 0.05): 0.04225 — without the correction 0.048
- **fit_trend_map** [cat6_montecarlo] n=120, phi=0.0: 95% interval coverage: 0.953
- **fit_trend_map** [cat6_montecarlo] n=120, phi=0.3: false-positive rate (nominal 0.05): 0.05875 — without the correction 0.148
- **fit_trend_map** [cat6_montecarlo] n=120, phi=0.3: 95% interval coverage: 0.9387
- **fit_trend_map** [cat6_montecarlo] n=120, phi=0.6: false-positive rate (nominal 0.05): 0.062 — without the correction 0.331
- **fit_trend_map** [cat6_montecarlo] n=120, phi=0.6: 95% interval coverage: 0.9335
- **fit_trend_map** [cat6_montecarlo] n=120, phi=0.9: false-positive rate (nominal 0.05): 0.09725 — without the correction 0.656
- **fit_trend_map** [cat6_montecarlo] n=120, phi=0.9: 95% interval coverage: 0.898
- **fit_1D_trend** [cat6_montecarlo] theil_sen, n=240, phi=0.5: false-positive rate: 0.2617 — 600 replicates (s.e. ~0.009); documented: Theil-Sen and Mann-Kendall assume independence, so expected high
- **fit_1D_trend** [cat6_montecarlo] theil_sen, n=240, phi=0.5: 95% interval coverage: 0.7433
- **fit_1D_trend** [cat6_montecarlo] linear truth crossing zero (no exponential candidate): share selected as linear: 0.865 — {'quadratic': 22, 'linear': 173, 'sinusoidal': 5}
- **analyze_correlation** [cat6_montecarlo] n=240, phi=0.8: rejection rate with n_eff (nominal 0.05): 0.05 — without the correction 0.349
- **analyze_correlation** [cat6_montecarlo] n=120, phi=0.9: rejection rate with n_eff (nominal 0.05): 0.0455 — without the correction 0.493
- **transform_unit** [cat7_safeguards] documented limitation: sigma + 1000 still labelled 'anomaly' -> converted value ratio: 1.976 — the docstring warns of this; ~2 means 1000 was added twice
