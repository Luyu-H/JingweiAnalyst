---
id: derive-unavailable-diagnostics
title: "Standard CMIP6 diagnostics missing from the local archive, and how to derive them from what exists"
type: method
summary: "A lookup from commonly requested CMIP6 ocean variables that are not stored locally (sos, tob, sob, chlos, o2min, zo2min, msftmz, obvfsq, thkcello, zos, spco2, fgco2, intpp, epc100, hfds) to their CMOR definition and to a derivation from local variables, or to \"not derivable\"."
keywords: [missing variable, 缺失变量, derive, 推导, sos, sea surface salinity, 海表盐度, tob, bottom temperature, 底层温度, sob, chlos, o2min, zo2min, msftmz, AMOC, obvfsq, N2, buoyancy frequency, 浮力频率, thkcello, zos, sea level, 海平面, spco2, fgco2, air-sea CO2 flux, intpp, primary production, 初级生产力, epc100, export, hfds, heat flux]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [sos, tob, sob, chlos, o2min, zo2min, msftmz, obvfsq, thkcello, zos, spco2, fgco2, intpp, epc100, hfds]
tools: [load_regular_data, load_volumes_for_cmip6, load_areas, compute_vertical_gradient, locate_extreme, compute_potential_density_anomaly]
severity: should
related: [cmip6-tos-bulk-top-cell, cmip6-static-grid-fields, overturning-from-mass-transport, cmip6-o2sat-and-derived-oxygen]
sources:
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry sos, tob, sob, chlos, o2min, zo2min, msftmz, obvfsq, thkcello, zos, spco2, fgco2, intpp, epc100, hfds"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix H16, H18, H21, H22 (pp. 3261–3262); Table H1 (p. 3254)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
---

## Key points (CMOR definitions)

| requested | CMOR definition | derivable locally? |
|---|---|---|
| `sos` | Sea surface salinity, units 0.001 | yes: top level of `so` |
| `tob` | "Potential temperature at the ocean bottom-most grid cell" | yes: `thetao` at each column's deepest valid level |
| `sob` | "Model prognostic salinity at bottom-most model grid cell" | yes: `so` at the deepest valid level |
| `chlos` | surface chlorophyll, sum of all groups, kg m-3 | yes: top level of `chl` |
| `o2min` / `zo2min` | O2 minimum concentration / its depth ("if two, then the shallower") | yes, from `o2`. Take the first minimum from the top on ties. |
| `msftmz` / `msftyz` | overturning mass streamfunction by basin, kg s-1 | partly: MPI from `vmo`, CESM2 from `vo` (see `overturning-from-mass-transport`), GFDL-ESM4 not at all |
| `obvfsq` | squared buoyancy frequency N², s-2. OMIP recommends **locally referenced** potential density | yes: from `thetao` and `so` with TEOS-10 |
| `thkcello` | cell thickness, m | yes: `volcello / areacello` |
| `zos` | dynamic sea level, "should have zero global area mean" | no |
| `spco2` | surface aqueous pCO2, Pa | no (it needs full carbonate chemistry from dissic, talk, T, S and nutrients) |
| `fgco2` | air–sea CO2 flux, kg C m-2 s-1, positive into the ocean | no |
| `intpp` | vertically integrated primary production, mol m-2 s-1 | no |
| `epc100` | POC flux at 100 m, mol m-2 s-1 | no |
| `hfds` | net downward heat flux at the sea-water surface, W m-2 | no |

[CMOR Omon; OMIP-phys H16, H18, H21, H22]

## Implications for this project

- If a question asks for a "not derivable" variable, say that the local
  archive does not hold it. Do not substitute a loosely related field without
  saying so.
- Derived versions of archived diagnostics differ from what the model would
  have archived, because they are built from monthly means (see
  `monthly-means-are-time-averages`). Label them as derived, for example
  "sos (derived from so level 1)".
- Bottom values: the "deepest valid level" differs per column. Find it from
  the field's own NaN pattern or from `deptho` with the level bounds. Do not use
  a fixed level index.
- For N², OMIP recommends locally referenced potential density rather than σ0,
  which earlier CMIPs used. Compute N² from adjacent levels referenced to their
  mid-point pressure, and state the method if σ0 is used instead.

## Source excerpts

> "We recommend use of locally referenced potential density for computing
> this measure of vertical gravitational stability." (OMIP-phys, H22)
