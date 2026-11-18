I assessed nitracline depth in WOA23 and in three CMIP6 historical simulations: CESM2, GFDL-ESM4, and MPI-ESM1-2-LR.

Method used
- WOA23: I loaded the annual nitrate climatology from WOA23 and masked it with the WOA23 land–sea mask, keeping ocean where `landsea_01 > 1`.
- WOA23 nitracline: I diagnosed nitracline depth with `compute_nutricline_depth(…, min_depth=0, max_depth=1000)`, so the search was within the top 1000 m.
- CMIP6 models: I loaded monthly `no3`, `thetao`, and `so` from `cmip6_historical` for 1971-01 to 2000-12, then formed a long-term mean with `temporal_aggregate(..., freq="mean", op="mean")`.
- Unit harmonization: Model nitrate was converted to `umol/kg` using `transform_unit(..., to_unit="umol/kg", density=sigma_insitu)`, where `sigma_insitu` came from `compute_potential_density_anomaly(..., ref_pressure=None, temperature_type="auto")`.
- Model nitracline: I then diagnosed nitracline depth for each model with `compute_nutricline_depth(..., min_depth=0, max_depth=1000)`.
- Means: Global and basin means were computed as area-weighted means with `compute_basic_stat(..., op="area_weighted_mean", dim=["lat","lon"])`.
- Basins: Pacific, Atlantic, and Indian basin masks were applied with `select_basin(..., include_marginal_seas=True, depth_dependent=False)`. Southern and Arctic Oceans were not included in those basin summaries.
- Map skill metrics: Each model nitracline map was regridded to the WOA23 1° grid using `regrid_horizontal(..., target="woa23", method="bilinear", skipna=True, na_thres=0.5)`. Metrics were then computed against WOA23 on a single shared ocean mask common to all four fields, using 40,614 grid cells. Spatial Pearson correlation was area-weighted; RMSE came from `compute_difference_metrics(..., weighting="area")`.

Direct answer

Global mean nitracline depth
- WOA23: 166 m
- CESM2: 147 m, bias relative to WOA23 = -19 m
- GFDL-ESM4: 133 m, bias = -33 m
- MPI-ESM1-2-LR: 132 m, bias = -34 m

So all three models simulate a shallower global mean nitracline than observed.

Basin mean nitracline depth

WOA23
- Pacific: 201 m
- Atlantic: 165 m
- Indian: 180 m

CESM2
- Pacific: 162 m, bias = -38 m
- Atlantic: 110 m, bias = -56 m
- Indian: 165 m, bias = -14 m

GFDL-ESM4
- Pacific: 159 m, bias = -41 m
- Atlantic: 91 m, bias = -74 m
- Indian: 169 m, bias = -11 m

MPI-ESM1-2-LR
- Pacific: 184 m, bias = -17 m
- Atlantic: 111 m, bias = -54 m
- Indian: 121 m, bias = -59 m

This shows that all three models are too shallow in all three major basins, but the size of the bias differs by basin and model.

Spatial skill against WOA23 on the WOA23 1° grid
- CESM2: Pearson r = 0.398, RMSE = 163.7 m
- GFDL-ESM4: Pearson r = 0.499, RMSE = 161.2 m
- MPI-ESM1-2-LR: Pearson r = 0.602, RMSE = 150.5 m

Best-performing model
- MPI-ESM1-2-LR reproduces the observed nitracline best overall.
- It has both the highest spatial correlation and the lowest RMSE of the three models.

Observed basin ordering and whether models reproduce it
- WOA23 ordering from shallowest to deepest: Atlantic < Indian < Pacific
- CESM2: Atlantic < Pacific < Indian → does not match
- GFDL-ESM4: Atlantic < Pacific < Indian → does not match
- MPI-ESM1-2-LR: Atlantic < Indian < Pacific → matches

Key interpretation
- The observed nitracline is deepest in the Pacific, intermediate in the Indian, and shallowest in the Atlantic.
- Every model underestimates nitracline depth globally.
- CESM2 is closest to WOA23 in the global mean, but not in map pattern skill.
- GFDL-ESM4 improves somewhat on CESM2 in spatial pattern metrics, but still has a large shallow bias, especially in the Atlantic.
- MPI-ESM1-2-LR has the best spatial agreement with WOA23 and is the only model to reproduce the observed ordering across the three basins, although it still has strong shallow biases in the Atlantic and Indian.

Important numerical patterns
- Atlantic biases are large and negative in all three models: -56 m, -74 m, and -54 m.
- Pacific bias is smallest in MPI-ESM1-2-LR (-17 m), much larger in CESM2 and GFDL-ESM4 (-38 m and -41 m).
- Indian bias is relatively modest in CESM2 and GFDL-ESM4 (-14 m and -11 m), but large in MPI-ESM1-2-LR (-59 m).
- Thus, MPI-ESM1-2-LR wins overall on full-map skill, but not because it is uniformly best in every basin mean.

Figures generated

Figure 1
- A 2×2 panel of global nitracline-depth maps was generated for WOA23, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR.
- The three model maps were first regridded to the WOA23 grid.
- All four panels use one shared sequential color scale (`viridis`), with the upper end clipped at the 98th percentile of pooled values.
- The figure is designed to support direct visual comparison of the geographic structure of nitracline depth across observation and models.

Figure 2
- A 4-panel longitude–depth nitrate section was generated for the equatorial Pacific, ordered WOA23, CESM2, GFDL-ESM4, MPI-ESM1-2-LR.
- The longitude range is 120°E to 280°E (80°W), and the latitude average is over 2°S–2°N.
- WOA23 used the annual climatology; each CMIP6 model used the 1971-01 to 2000-12 mean.
- All sections are in `umol kg^-1`, plotted over the top 1000 m with depth increasing downward.
- Each panel includes that source’s own nitracline depth line overlaid on the section.
- One shared color scale was used across all four panels, clipped at the pooled 2nd and 98th percentiles.

Overall conclusion
- All three CMIP6 models simulate a nitracline that is too shallow relative to WOA23, both globally and in the Pacific, Atlantic, and Indian basin means.
- Among them, MPI-ESM1-2-LR gives the best overall reproduction of observed nitracline depth based on the spatial correlation and RMSE on the common WOA23-grid mask.
- MPI-ESM1-2-LR is also the only model that reproduces the observed inter-basin ordering: Atlantic shallowest, Indian intermediate, Pacific deepest.
