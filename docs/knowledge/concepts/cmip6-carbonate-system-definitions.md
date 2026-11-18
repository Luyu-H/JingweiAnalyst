---
id: cmip6-carbonate-system-definitions
title: "CMIP6 dissic, talk and ph: what they contain and on which pH scale"
type: concept
summary: "dissic is total dissolved inorganic carbon (CO2*+HCO3-+CO3--), talk includes carbonate, borate, water, phosphate and silicate (and nitrogen) terms, and ph is in-situ pH on the TOTAL scale; all three models follow OMIP best-practice constants."
keywords: [dissic, DIC, TCO2, 溶解无机碳, talk, alkalinity, 总碱度, TAlk, ph, pH, total scale, 总标度, seawater scale, free scale, carbonate chemistry, 碳酸盐体系, mocsy, GLODAP pH]
datasets: [cmip6_historical, cmip6_ssp245, glodap]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [dissic, talk, ph]
tools: [transform_unit, compute_difference_metrics]
severity: must
related: [cmip6-bgc-units-per-volume, cmip6-integrated-carbon-content, model-gfdl-esm4-configuration, model-mpi-esm1-2-lr-configuration, model-cesm2-configuration]
sources:
  - citation: "Orr, J. C., et al. (2017). Biogeochemical protocols and diagnostics for the CMIP6 Ocean Model Intercomparison Project (OMIP). Geosci. Model Dev., 10, 2169–2199."
    locator: "Section 2 intro (p. 2173); Section 2.6, Eqs. 32–38 (pp. 2182–2183); Table 10 (p. 2188)"
    url: https://doi.org/10.5194/gmd-10-2169-2017
  - citation: "PCMDI cmip6-cmor-tables, Tables/CMIP6_Omon.json, data_specs_version 01.00.33 (18 Nov 2020)."
    locator: "variable_entry dissic, talk, ph"
    url: https://github.com/PCMDI/cmip6-cmor-tables
  - citation: "Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1 (GFDL-ESM 4.1). J. Adv. Model. Earth Syst., 12, e2019MS002015."
    locator: "Section 2.6 (p. 7)"
    url: https://doi.org/10.1029/2019MS002015
  - citation: "Mauritsen, T., et al. (2019). Developments in the MPI-M Earth System Model version 1.2 (MPI-ESM1.2). J. Adv. Model. Earth Syst., 11, 998–1038."
    locator: "Section 5.2 (p. 1014)"
    url: https://doi.org/10.1029/2018MS001400
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 2.2 (pp. 4–5)"
    url: https://doi.org/10.1029/2019MS001916
---

## Key points

- **dissic**: "Dissolved inorganic carbon (CO3+HCO3+H2CO3) concentration",
  mol m-3. It is C_T = CO2* + HCO3- + CO3--. [CMOR Omon; OMIP-BGC §2.6]
- **talk**: "total alkalinity equivalent concentration (including carbonate,
  nitrogen, silicate, and borate components)", mol m-3. OMIP requires
  A_T = A_C + A_B + A_W + A_P + A_Si (+ others). It notes that neglecting the
  phosphate and silicate terms "has been common among model groups but leads to
  systematic errors in computed pCO2". [CMOR Omon; OMIP-BGC Eq. 32]
- **ph**: standard name `sea_water_ph_reported_on_total_scale`. The comment
  is "negative log of hydrogen ion concentration with the concentration
  expressed as mol H kg-1". OMIP asks for the **total pH scale** and the
  Dickson et al. (2007) best-practice constants. Derived quantities are
  computed through the whole water column (in situ), not only at the surface.
  [CMOR Omon; OMIP-BGC §2, §2.6]
- Model implementations:
  - GFDL-ESM4 (COBALTv2) uses **mocsy 2.0** with Lee et al. (2010) total boron,
    Millero (2010) K1/K2 and Dickson & Riley (1979) Kf. [GFDL §2.6]
  - MPI-ESM1.2 (HAMOCC6) moved to the total scale and the Dickson constants, and
    added phosphate and silicate alkalinity. [MPI §5.2]
  - CESM2 (MARBL) has "fully prognostic carbonate chemistry". The paper does
    not name the constants. [CESM2 §2.2]

## Implications for this project

- Compare model `ph` only with an observational pH on the **total scale at
  in-situ temperature and pressure**. A pH at 25 °C / 0 dbar is a different
  quantity. Among GLODAP's pH fields, check which one is in-situ
  (`configs/datasets.yaml` lists `pHtsinsitutp` and `pHts25p0`).
- `dissic` and `talk` are per volume. See `cmip6-bgc-units-per-volume` for the
  conversion to µmol kg-1.
- The CESM2 paper does not document its carbonate constants. Mention this as
  an unverified assumption if the comparison depends on them.

## Source excerpts

> "For carbonate chemistry, groups should strive to use the constants
> recommended for best practices (Dickson et al., 2007) on the total pH scale
> and to avoid common modeling assumptions that lead to significant biases,
> notably an oversimplified alkalinity equation." (OMIP-BGC, §2)
