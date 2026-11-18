---
id: cmip6-o2sat-and-derived-oxygen
title: "CMIP6 o2 and o2sat, and deriving AOU and percent saturation from them"
type: concept
summary: "o2sat is the model's oxygen saturation concentration (mol m-3, per volume) on the full 3-D grid; AOU = o2sat − o2 and percent saturation = 100·o2/o2sat can be derived per model, but observational AOU and saturation use their own solubility formula and mass units."
keywords: [o2, oxygen, 溶解氧, o2sat, oxygen saturation, 氧饱和度, AOU, apparent oxygen utilization, 表观耗氧量, percent saturation, Garcia and Gordon, solubility, 溶解度, deoxygenation, 脱氧, hypoxia, 低氧, suboxia, oxygen minimum zone, OMZ, 最小含氧带]
datasets: [cmip6_historical, cmip6_ssp245, woa23]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [o2, o2sat, o2min, zo2min]
tools: [compute_ratio, identify_oxygen_minimum_zone, transform_unit]
severity: should
related: [cmip6-bgc-units-per-volume, derive-unavailable-diagnostics]
sources:
  - citation: "Orr, J. C., et al. (2017). Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model Intercomparison Project (OMIP). Geosci. Model Dev., 10, 2169–2199."
    locator: "Section 2.1.2 (p. 2174); Section 2.5.3, Eq. 16 (pp. 2179–2180); Tables 7, 10, 12, 15"
    url: https://doi.org/10.5194/gmd-10-2169-2017
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry o2, o2sat, o2min, zo2min"
    url: https://github.com/PCMDI/cmip6-cmor-tables
---

## Key points

- `o2` is "Dissolved Oxygen Concentration" and `o2sat` is "Dissolved Oxygen
  Concentration at Saturation". Both are mol m-3 on (lon, lat, olevel, time).
  `o2sat` means "the mole concentration in a saturated solution". [CMOR Omon]
- For oxygen, OMIP takes the saturation concentration [O2]⁰_sat from **Garcia
  and Gordon (1992)** (Eq. 8, Table 1). The atmospheric pressure correction is
  [A]_sat = (P_a / P⁰_a) · [A]⁰_sat, with atmospheric O2 fixed at a mole
  fraction of 0.20946. [OMIP-BGC §2.1.2, §2.5.3]
- Surface pressure variations "must not be neglected". Between 60°S and 30°S the
  mean surface pressure is about 3 % below the global mean, which lowers
  [O2]_sat there by about 10 µmol kg-1. [OMIP-BGC §2.5.3]
- The data request also defines `o2min` (oxygen minimum concentration) and
  `zo2min` (its depth, the shallower one if there are two). [CMOR Omon;
  OMIP-BGC Table 15]

## Implications for this project

- Per model, AOU = `o2sat` − `o2` and percent saturation = 100 · `o2` /
  `o2sat` are consistent derivations (same grid, same member, same per-volume
  basis). The volume basis cancels in the percent saturation but **not** in
  AOU.
- Against WOA23 AOU or percent saturation, the model and observed saturation
  values may rest on different solubility conventions and different pressure
  assumptions. The model value follows the OMIP recipe above, while the WOA
  convention has to be taken from WOA documentation. Treat the differences as
  including a definitional component. For AOU, convert to µmol kg-1 first.
- Low-oxygen thresholds are usually stated in µmol kg-1, commonly 60 for
  hypoxia and 5 for suboxia. `identify_oxygen_minimum_zone` compares the field
  with the threshold in the field's own units, so model `o2` (mol m-3) has to
  be converted first (see `cmip6-bgc-units-per-volume`).
- Availability is uneven: CESM2 has no `o2`/`o2sat` in historical, only in
  ssp245 (`configs/datasets.yaml`). `o2min`/`zo2min` are not in the local
  archive and must be derived from `o2` (see `derive-unavailable-diagnostics`).

## Common mistakes

- Computing AOU from a model `o2` field and an observational saturation field
  (or the reverse). The two fields are on different bases and conventions.
- Reporting model AOU in mol m-3 next to WOA AOU in µmol kg-1.

## Source excerpts

> "For O2, it is not φ⁰_A that is available, but rather [O2]⁰_sat (Garcia and
> Gordon, 1992)" (OMIP-BGC, §2.5.3)
