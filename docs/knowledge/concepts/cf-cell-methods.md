---
id: cf-cell-methods
title: "Reading CF cell_methods on CMIP6 ocean variables (mean vs sum, \"where sea\", climatologies)"
type: concept
summary: "cell_methods says how each cell value was formed; \"mean\" fields are intensive and are averaged with a weight, \"sum\" fields are extensive and are added, \"where sea\" restricts the statistic to the ocean part of the cell, and climatologies use \"within years / over years\"."
keywords: [cell_methods, cell_measures, CF conventions, CF 约定, area mean, where sea, extensive, intensive, 广延量, 强度量, climatology bounds, 气候态, within years, over years, time mean, 时间平均]
datasets: [cmip6_historical, cmip6_ssp245, woa23, glodap, seanoe_mld, ersst]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, volcello, areacello, wmo, intdic, mlotstmax]
tools: [compute_basic_stat, temporal_aggregate]
severity: info
related: [omip-weighted-means-and-integrals, cmip6-static-grid-fields, cmip6-integrated-carbon-content]
sources:
  - citation: "Eaton, B., et al. (2017). NetCDF Climate and Forecast (CF) Metadata Conventions, Version 1.7."
    locator: "Sections 7.2 Cell Measures, 7.3 Cell Methods (7.3.1, 7.3.3), 7.4 Climatological Statistics; Appendix E"
    url: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json and CMIP6_Ofx.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "cell_methods / cell_measures fields of ocean variables"
    url: https://github.com/PCMDI/cmip6-cmor-tables
---

## Key points

- `cell_methods` is a list of `name: method` pairs, applied **left to right**.
  `area` stands for the horizontal dimensions together. [CF 7.3, 7.3.1]
- `point` is the default for intensive quantities and `sum` for extensive ones.
  Methods include mean, sum, maximum, minimum, standard_deviation and variance.
  variance has units u². [CF 7.3, Appendix E]
- `area: mean where sea` means the statistic was evaluated **only over the sea
  portion of the cell**. `mean where type1 over type2` distinguishes the
  averaging domain. [CF 7.3.3]
- `cell_measures` (for example `area: areacello volume: volcello`) names the
  variables that hold the cell area and volume to weight with. These may live
  in another file listed in `external_variables`. [CF 7.2]
- Climatologies use `time: mean within years time: mean over years`, with the
  intervals in a `climatology` bounds variable. A DJF interval runs from
  1 Dec to 1 Mar across the year boundary. [CF 7.4]
- The ocean variables in the CMIP6 tables show these patterns [CMOR Omon/Ofx]:

  | pattern | examples | how to aggregate |
  |---|---|---|
  | `area: mean where sea time: mean` | thetao, so, tos, o2, no3, dissic, mlotst | weighted mean (area or volume) |
  | `area: sum` / `area: sum where sea` | areacello, volcello, wmo | add up; never average without weights |
  | `... depth: sum where sea ...` | intdic, intdoc, intpoc | already integrated over depth |
  | `time: maximum` | mlotstmax | max of instantaneous values within the month |
  | `time: mean` only | uo, vo, wo, umo, vmo | vector components; see `velocity-components-are-grid-relative` |

## Implications for this project

- Before aggregating a variable, read its cell_methods. A field marked `sum`
  (such as `wmo` or `volcello`) is extensive. Its regional total is a plain sum,
  and averaging it is meaningless unless the result is divided by a measure.
- The monthly `time: mean` fields are averages of all model time steps, not
  snapshots (see `monthly-means-are-time-averages`).

## Source excerpts

> "If more than one cell method is to be indicated, they should be arranged in
> the order they were applied. The left-most operation is assumed to have been
> applied first." (CF-1.7, 7.3.1)
