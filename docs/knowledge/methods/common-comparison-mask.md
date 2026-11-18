---
id: common-comparison-mask
title: "One comparison sample: a shared mask for model–observation and multi-model metrics"
type: method
summary: "A bias, RMSE or pattern correlation between two fields exists only where both have data, and several models' cell-by-cell scores can be ranked only on the cells all of them cover; for a cell-by-cell comparison of one field put the fields on one grid, build one mask over all of them and apply it to every field (a pattern that an analysis produced on each source's own domain, such as a variability mode, is correlated with the reference over the cells the two share), and keep WOA23's land-sea mask out of comparisons WOA23 is not part of (GLODAP, for instance)."
keywords: [comparison mask, 比较掩膜, common mask, 公共掩膜, shared mask, valid cells, 有效格点, ranking, 排名, model-observation comparison, 模式与观测对比, bias, 偏差, RMSE, 均方根误差, pattern correlation, 空间相关, land-sea mask, 海陆掩膜, landsea_01, ocean_only, coastline, 海岸线]
datasets: [woa23, glodap, ersst, cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [tos, thetao, so, o2, no3, po4, dissic, talk, ph]
tools: [build_comparison_mask, compute_difference_metrics, analyze_correlation, regrid_horizontal]
severity: must
related: [omip-native-grid-vs-regridded, omip-weighted-means-and-integrals, missing-values-and-land-cells]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix A3 (p. 3239)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "Lauvset, S. K., et al. (2016). A new global interior ocean mapped climatology: the 1° × 1° GLODAP version 2. Earth Syst. Sci. Data, 8, 325–340."
    locator: "Abstract"
    url: https://doi.org/10.5194/essd-8-325-2016
---

## Key points

- Regridding model output to a spherical grid "generally involves a land–sea
  mask distinct from the model native mask". [OMIP-phys A3] After remapping,
  each source covers its own set of cells, and coastal cells differ most.
- The GLODAPv2 mapped climatology is a global 1° × 1° product on 33 standard
  depth surfaces, mapped from the GLODAPv2 bottle data. [Lauvset et al. 2016]
  It shares the horizontal geometry of the WOA 1° grid, but not WOA's data
  coverage or its land-sea mask.

## Implications for this project

- A model is comparable with an observation, or with another model, only where
  both have data, and several models can be ranked against each other on a
  cell-by-cell score (a bias, an RMSE or the pattern correlation of one field)
  only where all of them do. Left to itself each pair covers a different set of
  cells, so such a score from one pair describes a different piece of ocean
  from the next, and a ranking between them means nothing. Put the fields on
  one grid, build one mask over all of them with `build_comparison_mask`, and
  apply it with `.where(mask)` to every field before any metric is taken.
- Where the shared mask does not apply: a pattern that is itself the result of
  an analysis on each source's own domain, such as the regression pattern of a
  variability mode (see `eof-modes-and-pdo-index`). Each mode was computed on
  its own source's cells, so masking the patterns afterwards to the other
  sources' coverage does not make them comparable. Correlate each with the
  reference pattern over the cells the two share (`analyze_correlation` leaves
  out the cells either one lacks), and compare the models on those
  correlations. Nor does a mask built on one level carry to another: coverage
  shrinks with depth, so a statistic at several depths takes its cells from the
  fields at each depth, not from the surface mask.
- The numbers of one cell-by-cell comparison come from one sample. Report the
  means that `compute_difference_metrics` returns over the shared sample
  (`mean_data`, `mean_reference`) together with the bias, not means computed
  separately over different cells. A global or regional mean of each source
  over its own ocean, with its own weights, is a different quantity (see
  `omip-native-grid-vs-regridded`): when the difference of such means is asked
  for, take each mean on its own grid and say that the difference is not the
  bias on the shared mask.
- `ocean_only=True`, the default of `build_comparison_mask`, is right only when
  WOA23 is part of the comparison. The check is made on the grid, so a
  GLODAP-against-models comparison on the shared 1° geometry silently inherits
  WOA23's coastline and drops cells where GLODAP has data and WOA23 says land.
  Those are among the hardest cells for a model to match, so dropping them
  flatters every model. Pass `ocean_only=False` when the reference is not WOA23.
- Report the cell counts the mask records alongside the metrics. A mask much
  smaller than the individual fields means the products disagree about where
  the ocean is, not that the models are poor.

## Checked in the local files

- GLODAPv2.2016b surface `TCO2` has 41320 valid cells on its 1° grid, and
  WOA23 `landsea_01` calls 1459 of them land: 442 north of 60°N, 179 south of
  60°S, the rest coastal cells at all latitudes. `build_comparison_mask` on the
  GLODAP field keeps 39861 cells with `ocean_only=True` and all 41320 with
  `ocean_only=False` (checked 2026-10-20).

## Common mistakes

- Scoring each model on its own pairwise mask in a cell-by-cell comparison of
  one field, and then ranking the models on those scores.
- Applying the WOA23 land-sea mask to a comparison against GLODAP.

## Source excerpts

> "Unfortunately, regridding from model native grids to spherical grids
> generally involves a land–sea mask distinct from the model native mask."
> (OMIP-phys, A3)

> "Climatologies were created for 33 standard depth surfaces." (Lauvset et al.
> 2016, abstract)
