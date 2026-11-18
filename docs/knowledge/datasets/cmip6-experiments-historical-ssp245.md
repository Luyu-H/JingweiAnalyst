---
id: cmip6-experiments-historical-ssp245
title: "What the CMIP6 historical and ssp245 experiments are"
type: dataset
summary: "historical is the all-forcing 1850–2014 run branched from piControl; ssp245 is the concentration-driven SSP2-4.5 scenario (an update of RCP4.5, medium forcing) for 2015–2100 branched from historical; both are coupled runs with free internal variability, so their years are not synchronised with the observed years."
keywords: [historical, 历史试验, ssp245, SSP2-4.5, scenario, 情景, ScenarioMIP, concentration-driven, 浓度驱动, RCP4.5, piControl, esm-hist, experiment, 试验, internal variability, 内部变率]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: []
tools: [load_regular_data]
severity: should
related: [historical-ssp245-member-continuity, model-drift-in-deep-trends]
sources:
  - citation: "WCRP-CMIP CMIP6_CVs, CMIP6_experiment_id.json, CV collection version 6.2.60.0."
    locator: "experiment_id historical, ssp245, piControl, esm-hist"
    url: https://github.com/WCRP-CMIP/CMIP6_CVs
  - citation: "Orr, J. C., et al. (2017). Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model Intercomparison Project (OMIP). Geosci. Model Dev., 10, 2169–2199."
    locator: "Section 2 (p. 2173): forced OMIP runs vs coupled historical runs"
    url: https://doi.org/10.5194/gmd-10-2169-2017
---

## Key points

| experiment | description | years | parent | activity |
|---|---|---|---|---|
| historical | "all-forcing simulation of the recent past" | 1850–2014 (min 165 yr) | piControl (or past1000/past2k) | CMIP |
| ssp245 | "Future scenario with medium radiative forcing by the end of century. Following approximately RCP4.5 global forcing pathway but with new forcing based on SSP2. **Concentration-driven**" | 2015–2100 (min 86 yr) | historical | ScenarioMIP |
| piControl | "pre-industrial control" | ≥ 500 yr | piControl-spinup | CMIP |
| esm-hist | historical with **CO2 emission-driven**, atmospheric CO2 calculated | 1850–2014 | esm-piControl | CMIP |

[CVs experiment_id]

- The coupled historical run differs from the forced ocean-only OMIP runs
  (omip1) "only" in initialisation and forcing. omip1 is initialised from
  observations and forced by reanalysis, while historical is coupled after a
  spin-up. [OMIP-BGC §2]

## Implications for this project

- Both local experiments are **concentration-driven** (not esm-hist). Atmospheric
  CO2 is prescribed, so ocean carbon uptake differences between models do not
  feed back on atmospheric CO2.
- Coupled runs have their own internal variability. A given model year (e.g.
  the 1997 El Niño year) is **not** expected to match the observed year.
  Compare statistics (means, variances, spectra, trends over long windows),
  not individual years or event timing.
- A single member per model (see `historical-ssp245-member-continuity`)
  cannot separate forced change from internal variability. Short-window trends
  should carry that caveat.

## Source excerpts

> ssp245: "Future scenario with medium radiative forcing by the end of century.
> Following approximately RCP4.5 global forcing pathway but with new forcing
> based on SSP2. Concentration-driven" (CMIP6_experiment_id.json)
