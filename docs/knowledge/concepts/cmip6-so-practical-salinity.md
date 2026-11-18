---
id: cmip6-so-practical-salinity
title: "CMIP6 so is (an approximation of) Practical Salinity, units \"0.001\""
type: concept
summary: "CMIP6 salinity is compared as Practical Salinity whatever the model's thermodynamics; \"0.001\" means parts per thousand; TEOS-10 routines need Absolute Salinity derived from it."
keywords: [so, sos, salinity, 盐度, practical salinity, 实用盐度, absolute salinity, 绝对盐度, PSS-78, psu, 0.001, units, 单位, TEOS-10, gsw]
datasets: [cmip6_historical, cmip6_ssp245, woa23, glodap]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [so, sos, sob]
tools: [compute_potential_density_anomaly, compute_mixed_layer_depth, transform_unit]
severity: should
related: [cmip6-thetao-potential-temperature, ocean-heat-and-salt-content]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix D4 (pp. 3247–3248); Appendix D5 (p. 3248); Table A1 (p. 3241); Table H1 (p. 3254)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry so, sos, sob"
    url: https://github.com/PCMDI/cmip6-cmor-tables
---

## Key points

- Most CMIP6 ocean models use pre-TEOS-10 thermodynamics. Their prognostic
  salinity approximates **Practical Salinity S_P** (PSS-78), the quantity
  stored in observational databases such as WOA. [OMIP-phys D4]
- The protocol states that CMIP6 salinity comparisons are made against Practical
  Salinity whatever the model thermodynamics. This matches earlier CMIPs.
  [OMIP-phys D4]
- TEOS-10 recommends **Absolute Salinity S_A** for thermodynamic
  calculations. S_A reflects the mass of dissolved constituents, while S_P
  depends only on conductivity. Adding non-ionic silicate raises S_A and density
  but leaves S_P unchanged. [OMIP-phys D4]
- CMOR units for `so`, `sos` and `sob` are the string **`0.001`**, meaning parts per
  thousand (dimensionless). The CMOR comment warns that the unqualified term
  "salinity" "does not necessarily imply any particular method of
  calculation". [CMOR Omon]
- Salt content in a cell is formally cell mass × S_A. Pre-TEOS-10 models,
  whose conserved salinity is S_P, report salt content as cell mass × S_P.
  [OMIP-phys D5]

## Implications for this project

- Model `so` and WOA23 salinity (practical salinity, units "1") can be
  differenced directly. No unit conversion is needed: "0.001", "1e-3", "1" and
  "psu" all denote the same practical-salinity number.
- For density or other TEOS-10 routines (gsw), convert S_P to S_A with the
  location-dependent function before calling routines that expect S_A.
  Passing `so` straight in as S_A biases density slightly, and the error varies
  by region.
- A unit parser that rejects "0.001" is not a sign of bad data. Map it to
  practical salinity explicitly.

## Common mistakes

- Multiplying `so` by 0.001 or 1000 because of the units string. The stored
  numbers are already of order 35.
- Calling a salinity difference "psu" in a report while treating the field as
  Absolute Salinity in a density calculation.

## Source excerpts

> "salinity comparisons in CMIP6 will be made against Practical Salinity
> regardless of the model thermodynamics." (OMIP-phys, D4)

> "The units of salinity are dimensionless and the units attribute should
> normally be given as 1e-3 or 0.001 i.e. parts per thousand." (CMOR Omon, so)
