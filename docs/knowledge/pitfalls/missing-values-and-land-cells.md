---
id: missing-values-and-land-cells
title: "CMIP6 missing value 1e20, land cells, and fill values under scaling"
type: pitfall
summary: "CMIP6 files mark land and below-sea-floor cells with the missing value 1e20; any read path that skips CF masking (raw netCDF4 arrays, scale/offset before masking, float32 overflow) turns those into huge numbers that silently dominate sums and means."
keywords: [missing value, 缺测值, _FillValue, 1e20, land mask, 陆地掩膜, NaN, fill value, 填充值, scale_factor, add_offset, 掩膜, below sea floor, 海底以下, weighted sum, 加权求和, skipna, threshold fraction, 超阈值比例]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, so, o2, uo, vo, deptho, areacello]
tools: [load_regular_data, temporal_aggregate, compute_exceedance_fraction]
severity: should
related: [cmip6-static-grid-fields, cf-cell-methods]
sources:
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "Header: missing_value 1e20, int_missing_value -999"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix A7 (p. 3242); Appendix G7 (p. 3253)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "Eaton, B., et al. (2017). NetCDF Climate and Forecast (CF) Metadata Conventions, Version 1.7."
    locator: "Section 2.5.1 Missing data, valid and actual range of data"
    url: https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/cf-conventions.html
---

## Key points

- The CMIP6 Omon table header sets `missing_value: 1e20` (integers: −999).
  [CMOR Omon header]
- "All fields are reported as 'missing' over grid cells that are entirely
  land." `areacello` "should be zero or missing over land", and `deptho` is
  missing over land. [OMIP-phys A7, G7; CMOR Ofx]
- CF: missing data is any value outside `valid_range` or equal to
  `_FillValue`/`missing_value`. For packed data, test for missing values
  **before** applying `scale_factor`/`add_offset`. Scaling the default float
  fill value can overflow to infinity. Missing data is not allowed in
  coordinate variables. [CF 2.5.1]

## Implications for this project

- xarray decodes 1e20 to NaN by default. The local files trigger a harmless
  "multiple fill values" warning, and both declared fill values decode to NaN.
  Code that reads with `netCDF4` directly, or opens with
  `mask_and_scale=False`, sees 1e20 at land and sub-bottom points. One such
  value makes a sum or mean meaningless.
- NaN-aware reductions are required, and the weights must be masked to the
  same cells (see `omip-weighted-means-and-integrals`).
- A value near 1e20 (or 9.97e36, the netCDF float default) in any output is
  an unmasked fill value, not a result.
- Masks can differ between variables of the same model. In the local
  historical files at the top level, first month, the valid cells of `uo` and
  `thetao` differ in 3019 cells for CESM2 (89115 vs 86096 valid) and in 968
  for MPI-ESM1-2-LR (35823 vs 36791). Build each field's valid mask from the
  field itself.
- A hand-rolled weighted reduction loses the land mask. Written as
  `(data * weights).sum('time') / weights.sum()`, an all-land column returns
  0.0 rather than NaN, because xarray's sum skips missing values; that zero
  then spreads into neighbouring ocean cells the moment the field is regridded.
  `temporal_aggregate` masks each cell's weights to the steps it has data, so a
  land cell stays missing.
- A hand-written threshold comparison counts land as "not exceeding": NaN
  compares as False, so land enters the denominator of any fraction.
  `compute_exceedance_fraction` leaves missing cells out of the numerator and
  the denominator.

## Source excerpts

> "All fields are reported as 'missing' over grid cells that are entirely
> land." (OMIP-phys, A7)
