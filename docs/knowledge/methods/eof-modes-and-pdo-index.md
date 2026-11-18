---
id: eof-modes-and-pdo-index
title: "EOF modes and the PDO index: global-mean removal, area weighting, sign, explained variance and pattern amplitude"
type: method
summary: "The PDO is the leading EOF of North Pacific SST anomalies with the global-mean anomaly removed; EOFs are area-weighted, their sign is a convention, explained variance depends on grid, domain and record length, and a mode's amplitude is read from its regression pattern on the standardised index, compared between sources with an area-weighted correlation."
keywords: [PDO, Pacific Decadal Oscillation, 太平洋年代际振荡, EOF, empirical orthogonal function, 经验正交函数, principal component, 主成分, leading mode, 主模态, explained variance, 解释方差, regression pattern, 回归型, North Pacific, 北太平洋, global mean removal, 去除全球平均, standardized index, 标准化指数]
datasets: [ersst, cmip6_historical, cpc_indices]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [tos, thetao]
tools: [compute_pdo_index, analyze_correlation, regrid_horizontal]
severity: should
related: [anomaly-baselines-and-climate-indices, serial-correlation-and-significance, omip-weighted-means-and-integrals, cmip6-tos-bulk-top-cell]
sources:
  - citation: "National Center for Atmospheric Research Staff (Eds.). The Climate Data Guide: Pacific Decadal Oscillation (PDO): Definition and Indices."
    locator: "Definition; expert guidance"
    url: https://climatedataguide.ucar.edu/climate-data/pacific-decadal-oscillation-pdo-definition-and-indices
  - citation: "National Center for Atmospheric Research Staff (Eds.). The Climate Data Guide: Empirical Orthogonal Function (EOF) Analysis and Rotated EOF Analysis."
    locator: "Spatial weighting; domain dependence of the patterns"
    url: https://climatedataguide.ucar.edu/climate-data-tools-and-analysis/empirical-orthogonal-function-eof-analysis-and-rotated-eof-analysis
---

## Key points

- The PDO is "defined by the leading pattern (EOF) of sea surface temperature
  (SST) anomalies in the North Pacific basin". The anomalies are formed by
  removing both the climatological annual cycle and the global-mean SST anomaly
  at each grid point, "as a way of removing any long-term trends that might be
  associated with changes in external forcing". [CDG PDO]
- "Positive values of the PDO index correspond with negative SST anomalies in
  central and western North Pacific", with positive anomalies in the eastern
  North Pacific. [CDG PDO]
- EOFs are computed on area-weighted data; "Most commonly, these spatial
  weights are sqrt(cos(latitude))". The derived patterns "may be domain
  dependent". [CDG EOF]

## Implications for this project

- `compute_pdo_index` removes the global-mean SST anomaly so that the index
  measures the North Pacific pattern rather than global warming, and weights
  each cell by the square root of its area so that the EOF is taken on a field
  whose variance is area-proportional.
- The sign of an EOF is arbitrary, so the tool fixes it by the convention above
  (positive phase cool in the central and western North Pacific) for both the
  pattern and the index. A negative correlation between two such patterns is
  then a real disagreement, not an arbitrary flip.
- The explained variance of a mode depends on the grid and on the length of the
  record, as well as on the domain: compare it only between fields on the same
  grid over the same period.
- The EOF itself is dimensionless and carries no amplitude. The regression
  pattern (the anomalies regressed on the standardised index, in units per
  standard deviation) is the map with a physical amplitude and the form a PDO
  pattern is normally drawn in; its area-weighted RMS is the mode's strength in
  degrees. The index is normally published standardised (zero mean, unit
  variance).
- Compare patterns from different sources on one grid with an area-weighted
  correlation (`analyze_correlation` with `weighting='area'`): the domain
  reaches high latitudes, where cells are much smaller.

## Common mistakes

- Comparing explained variances computed on different grids or over records of
  different length.
- Reading a mode's amplitude off the EOF pattern instead of the regression
  pattern.

## Source excerpts

> "Positive values of the PDO index correspond with negative SST anomalies in
> central and western North Pacific" (Climate Data Guide, PDO)

> "Most commonly, these spatial weights are sqrt(cos(latitude))." (Climate Data
> Guide, EOF analysis)
