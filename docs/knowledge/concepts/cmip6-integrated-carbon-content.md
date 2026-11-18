---
id: cmip6-integrated-carbon-content
title: "intdic, intdoc and intpoc are full-depth carbon mass per area (kg C m-2)"
type: concept
summary: "intdic/intdoc/intpoc are vertical integrals over the whole water column in kg of carbon per m2 (cell_methods \"depth: sum where sea\"), not mol m-2 and not partial-depth integrals."
keywords: [intdic, intdoc, intpoc, carbon inventory, 碳储量, 碳含量, column integral, 垂直积分, kg m-2, mol m-2, dissolved inorganic carbon content, particulate organic carbon]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [intdic, intdoc, intpoc, dissic, dissoc]
tools: [load_volumes_for_cmip6, compute_basic_stat]
severity: should
related: [cmip6-bgc-units-per-volume, cmip6-bgc-variable-scope, cf-cell-methods]
sources:
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry intdic, intdoc, intpoc, dissic"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Orr, J. C., et al. (2017). Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model Intercomparison Project (OMIP). Geosci. Model Dev., 10, 2169–2199."
    locator: "Table 11 and footnotes c–e (p. 2189)"
    url: https://doi.org/10.5194/gmd-10-2169-2017
---

## Key points

- `intdic` = "Dissolved Inorganic Carbon Content" (standard name
  `ocean_mass_content_of_dissolved_inorganic_carbon`), "Vertically integrated
  DIC", units **kg m-2**. [CMOR Omon]
- `intdoc` = vertically integrated DOC, "explicit pools only". `intpoc` =
  vertically integrated POC (standard name
  `ocean_mass_content_of_particulate_organic_matter_expressed_as_carbon`).
  Both are kg m-2. [CMOR Omon; OMIP-BGC Table 11]
- cell_methods for all three: `area: mean where sea depth: sum where sea time:
  mean`. The vertical operation is a **sum over the whole water column**.
  [CMOR Omon]
- The 3-D counterparts (`dissic`, `dissoc`) are in mol m-3. [CMOR Omon]

## Implications for this project

- Converting a per-volume field to the content variable:
  intdic ≈ Σ_k dissic_k · Δz_k · M_C, with M_C = 12.011e-3 kg mol-1 (the
  molar mass of carbon, a standard constant not given in the sources). The
  layer thickness can be taken as Δz_k = volcello / areacello. Minor
  differences from the model's own `intdic` are expected from time averaging
  (a product of monthly means versus a monthly mean of the product).
- To report inventories: global inventory (kg C) = Σ intdic · areacello over
  ocean cells. Divide by 1e12 for Pg C.
- Availability: CESM2 lacks `intdic` and `intpoc` in historical, and GFDL-ESM4
  and MPI have them (see `configs/datasets.yaml`). For CESM2 historical, derive
  the content from `dissic` as above and say so.

## Common mistakes

- Treating `intdic` as mol m-2 and multiplying by 12 again, or dividing a
  derived mol m-2 value by 1000 instead of converting by molar mass.
- Integrating `dissic` only over the upper ocean and comparing the result with
  `intdic`, which covers the full depth.

## Source excerpts

> intdic: units "kg m-2", cell_methods "area: mean where sea depth: sum where
> sea time: mean", comment "Vertically integrated DIC" (CMOR Omon)
