I assessed intermediate-water nitrate-to-phosphate stoichiometry over 200–1000 m for WOA23 and for three CMIP6 historical models: CESM2, GFDL-ESM4, and MPI-ESM1-2-LR.

Method actually used
- WOA23 nitrate and phosphate were loaded from the annual climatology over the 200–1000 m depth range.
- The WOA23 basin mask was loaded.
- CMIP6 historical monthly no3 and po4 were loaded for 1971-01 to 2000-12 over 200–1000 m, then aggregated to a climatology and averaged over month.
- Model cell volumes were taken from each model’s volcello. For GFDL-ESM4, time-varying volcello was loaded over 1971-01 to 2000-12; for CESM2 and MPI-ESM1-2-LR, volcello was loaded without time selection.
- Cellwise N:P was computed as nitrate/phosphate on each source’s native grid.
- Global and basin means of N:P were computed with volume-weighted means using the source data as passed to compute_basic_stat.
- For regressions and scatter plots, CMIP6 nitrate and phosphate were converted from native units to µmol kg-1 using density=1.025; WOA23 was used in matching µmol kg-1 units.
- Linear fits were obtained from fit_multivariate_relation with target = nitrate, predictor = phosphate, method="linear", validation="blocked". Equal-cell weighting was therefore implicit in the fit inputs rather than volume weighting.
- Fractions of water volume satisfying 14 < N:P < 18 and N:P < 12 were computed with compute_exceedance_fraction using weighting="auto".
- Basin summaries used select_basin for Pacific, Atlantic, and Indian with include_marginal_seas=True and depth_dependent=False.
- For the map figure, a 200–1000 m mean N:P−16 field was formed on each native grid; model maps were then bilinearly regridded horizontally to the WOA23 grid for display only.

Direct answer
WOA23 shows that intermediate waters are globally below Redfield in N:P, with a clear low-N:P signal consistent with denitrification, especially in the Pacific and Indian Oceans. Among the three models, GFDL-ESM4 comes closest to the observed mean stoichiometry and reproduces the global amount of low-N:P water reasonably well, but it misses the observed basin distribution. CESM2 and MPI-ESM1-2-LR are too N-poor overall and strongly overproduce low-N:P water, especially in the Pacific, so they do not reproduce the observed denitrification imprint well.

Key numerical results

1) Global mean N:P over 200–1000 m
Headline value = volume-weighted mean of cellwise N:P
Ratio-of-means = volume-weighted mean nitrate divided by volume-weighted mean phosphate

- WOA23:
  - Headline mean N:P: 14.15
  - Ratio of mean NO3 to mean PO4: 14.00
- CESM2:
  - Headline mean N:P: 11.83
  - Ratio of mean NO3 to mean PO4: 10.93
- GFDL-ESM4:
  - Headline mean N:P: 14.13
  - Ratio of mean NO3 to mean PO4: 13.65
- MPI-ESM1-2-LR:
  - Headline mean N:P: 12.28
  - Ratio of mean NO3 to mean PO4: 11.77

Interpretation:
- All four sources sit below the Redfield ratio of 16.
- WOA23 indicates moderately sub-Redfield intermediate waters.
- GFDL-ESM4 matches the observed global mean N:P closely.
- CESM2 and MPI-ESM1-2-LR are much more sub-Redfield than observed.

2) OLS regression of NO3 = a·PO4 + b over all valid 200–1000 m cells
Units are µmol kg-1 for both variables in the regression.

The executed fit object reported the slope in its equation and coefficients, but the post-processing extracted the wrong slope key, so the summary table showed NaN for slope. The actual fitted slopes are recoverable from the regression output:

- WOA23:
  - a = 13.04
  - b = 1.69 µmol kg-1
  - R² = 0.908
- CESM2:
  - a = 4.89
  - b = 10.84 µmol kg-1
  - R² = 0.238
- GFDL-ESM4:
  - a = 11.37
  - b = 4.02 µmol kg-1
  - R² = 0.817
- MPI-ESM1-2-LR:
  - a = 8.59
  - b = 5.64 µmol kg-1
  - R² = 0.678

Interpretation:
- None of the fitted slopes reaches 16, again consistent with sub-Redfield stoichiometry.
- WOA23 has the tightest NO3–PO4 relationship.
- GFDL-ESM4 is the closest model to observations in both slope and goodness of fit.
- CESM2 has a very weak phosphate control on nitrate in this layer, with a very low slope and low R².
- MPI-ESM1-2-LR is intermediate but still too low in slope.

3) Volume fraction with 14 < N:P < 18
- WOA23: 55.0%
- CESM2: 17.2%
- GFDL-ESM4: 74.6%
- MPI-ESM1-2-LR: 0.5%

Interpretation:
- Observed intermediate waters are often near, but generally a bit below, Redfield.
- GFDL-ESM4 places too much volume in the near-Redfield range.
- CESM2 and especially MPI-ESM1-2-LR place far too little volume there.

4) Basin headline mean N:P over 200–1000 m
Pacific / Atlantic / Indian

- WOA23:
  - Pacific: 13.57
  - Atlantic: 15.55
  - Indian: 13.73
- CESM2:
  - Pacific: 10.58
  - Atlantic: 13.55
  - Indian: 11.18
- GFDL-ESM4:
  - Pacific: 13.39
  - Atlantic: 15.07
  - Indian: 14.72
- MPI-ESM1-2-LR:
  - Pacific: 11.29
  - Atlantic: 13.15
  - Indian: 12.60

Interpretation:
- In WOA23, the Atlantic is closest to Redfield, while the Pacific and Indian are lower-N:P, consistent with a denitrification imprint being strongest outside the Atlantic.
- GFDL-ESM4 reproduces that broad ordering best.
- CESM2 and MPI-ESM1-2-LR are too low in all three basins, especially the Pacific.

5) Volume fraction with N:P < 12
Global / Pacific / Indian

- WOA23:
  - Global: 6.2%
  - Pacific: 8.1%
  - Indian: 10.3%
- CESM2:
  - Global: 36.7%
  - Pacific: 59.4%
  - Indian: 35.8%
- GFDL-ESM4:
  - Global: 7.3%
  - Pacific: 14.9%
  - Indian: 0.5%
- MPI-ESM1-2-LR:
  - Global: 22.5%
  - Pacific: 40.0%
  - Indian: 15.0%

Interpretation:
- WOA23 shows a real but limited low-N:P volume, with the strongest occurrence in the Indian and Pacific.
- CESM2 greatly overestimates low-N:P water globally and in both target basins.
- MPI-ESM1-2-LR also overestimates low-N:P water strongly, especially in the Pacific.
- GFDL-ESM4 gets the global low-N:P volume close to observed, but shifts the basin distribution: too much in the Pacific and far too little in the Indian.

Assessment against Redfield and denitrification signal

- WOA23:
  - Below Redfield.
  - Shows a clear low-N:P denitrification-like signal, particularly in the Pacific and Indian Oceans.

- CESM2:
  - Well below Redfield.
  - Does not reproduce the observed low-N:P signal realistically; it produces far too much low-N:P water globally and in both the Pacific and Indian Oceans.

- GFDL-ESM4:
  - Below Redfield, close to WOA23 in global mean.
  - Best of the three models overall.
  - Reproduces the global low-N:P volume reasonably well, but not its basin distribution.

- MPI-ESM1-2-LR:
  - Below Redfield.
  - Does not reproduce the observed low-N:P signal well; low-N:P waters are much too extensive, especially in the Pacific.

Model-by-model verdict on reproducing observed low-N:P water
Using the executed code’s comparison to WOA23:
- CESM2: does not reproduce observed low-N:P volume and basin distribution
- GFDL-ESM4: reproduces global low-N:P volume but not basin distribution
- MPI-ESM1-2-LR: does not reproduce observed low-N:P volume and basin distribution

What the visualizations show

Figure 1
- A four-panel scatterplot of nitrate versus phosphate over 200–1000 m, one panel each for WOA23, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR.
- All panels use common axis limits derived from pooled data percentiles.
- Each panel includes:
  - scatter points from valid cells,
  - a fitted linear regression line from the regression step,
  - a dashed Redfield 16:1 line through the origin,
  - an annotation box with fitted slope, intercept, and R².
- The purpose is to compare the NO3–PO4 relationship structure across observations and models on the same scale.

Figure 2
- A 2×2 set of global longitude–latitude maps of mean N:P minus 16 over 200–1000 m for WOA23, CESM2, GFDL-ESM4, and MPI-ESM1-2-LR.
- WOA23 is shown on its native grid; model fields are regridded horizontally to the WOA23 grid for display only.
- All panels use one shared diverging colour scale centered on zero, clipped using the pooled 98th percentile of absolute values.
- The purpose is to compare the spatial pattern of departures from Redfield.

Important patterns and insights
- Observed intermediate waters are substantially sub-Redfield on average, but not extremely so: the global mean N:P is about 14.15.
- The observed basin structure matters: Atlantic intermediate waters are much closer to Redfield than Pacific and Indian waters.
- GFDL-ESM4 is the only model that gets close to the observed global mean stoichiometry and global low-N:P volume.
- CESM2 and MPI-ESM1-2-LR both appear too nitrate-depleted relative to phosphate through much of the intermediate ocean.
- Model spread is large, and differences likely reflect both circulation and biogeochemical model structure; the mismatches should not be attributed to transport alone.

One technical note
- The regression slopes printed in the summary table were NaN because the extraction code looked for the wrong coefficient key. However, the regression output itself clearly reported the fitted equations, allowing the actual slopes above to be read directly and reported.
