---
id: cmip6-thetao-potential-temperature
title: "CMIP6 thetao is potential temperature referenced to the sea surface"
type: concept
summary: "thetao is potential temperature (reference pressure 0 dbar) in every CMIP6 model, even TEOS-10 models; at the surface it equals in-situ temperature, below it does not."
keywords: [thetao, potential temperature, 位温, conservative temperature, 保守温度, bigthetao, in-situ temperature, 现场温度, TEOS-10, temperature comparison, 温度对比, density, 密度, sigma0, potential density, 位密, in-situ density, 现场密度]
datasets: [cmip6_historical, cmip6_ssp245, woa23, glodap]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, tos, bigthetao]
tools: [compute_potential_density_anomaly, compute_mixed_layer_depth, compute_anomaly_vs_woa23]
severity: should
related: [cmip6-so-practical-salinity, cmip6-tos-bulk-top-cell, ocean-heat-and-salt-content]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix D2 (p. 3246); Appendix H10–H11 (p. 3259); Appendix H13 (p. 3260)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry thetao, bigthetao"
    url: https://github.com/PCMDI/cmip6-cmor-tables
---

## Key points

- The CMIP6 diagnostic protocol requires every model to archive **potential
  temperature θ**, `thetao`, referenced to the ocean surface. This holds whatever
  the model's prognostic temperature is. [OMIP-phys D2, H10]
- Models whose prognostic variable is TEOS-10 **Conservative Temperature Θ**
  still convert it online to θ for `thetao`. They archive Θ separately as
  `bigthetao`. The CMOR comment on `bigthetao` reads: "this should be contributed only
  for models using conservative temperature as prognostic field". [OMIP-phys
  D2, H11; CMOR Omon]
- CMOR standard_name of `thetao`: `sea_water_potential_temperature`, units
  `degC`, cell_methods `area: mean where sea time: mean`. The CMOR comment
  says: "Diagnostic should be contributed even for models using conservative
  temperature as prognostic field."
- Potential temperature is referenced to 0 dbar, so **surface θ equals surface
  in-situ temperature**. SST (`tos`) is also potential temperature, not
  Conservative Temperature. [OMIP-phys H13]

## Implications for this project

- None of the three local models provides `bigthetao`. Any TEOS-10 calculation
  (density, N², mixed layer) must treat `thetao` as potential temperature and
  convert it with the potential-temperature route, not the in-situ route.
  `temperature_kind()` in `src/tools/scientific_calculator.py` already reads this
  from the metadata. Taking the wrong route biases the density increasingly
  with depth.
- Which density for which purpose: a sigma referenced to a fixed pressure
  (sigma0 at the surface, sigma2 or sigma4 for the deep ocean) is what labels a
  water mass, because it removes the compression that makes any water denser
  with depth. The in-situ density, at each cell's own pressure
  (`compute_potential_density_anomaly` with `ref_pressure=None`), is what
  mass-per-volume arithmetic needs: converting a concentration between molar
  bases, or a volume into a mass (see `cmip6-bgc-units-per-volume`).
- Comparing `thetao` with an **in-situ** temperature climatology (WOA23 and
  GLODAP report in-situ temperature) is exact at the surface. The two
  quantities diverge with depth. For surface-only comparisons (SST maps,
  top-level bias) no conversion is needed. For sub-surface comparisons one side
  must be converted.

## Common mistakes

- Treating `thetao` as Conservative Temperature because the model is "modern".
  None of the three local models archived Θ.
- Converting SST from potential to in-situ temperature. The two are identical
  at 0 dbar.

## Source excerpts

> "Regardless of the model thermodynamics, modellers should archive potential
> temperature, θ." (OMIP-phys, D2)

> "Since the potential temperature is referenced to 0 gauge pressure, the
> surface potential temperature is the same as surface in situ temperature."
> (OMIP-phys, H13)
