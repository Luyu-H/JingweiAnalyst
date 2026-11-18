---
id: model-drift-in-deep-trends
title: "Control-run drift contaminates deep-ocean trends in historical and scenario runs"
type: pitfall
summary: "All three models still drift in their pre-industrial controls (deep warming, slow approach to equilibrium); without a piControl in the local archive, deep-ocean or whole-ocean trends from historical/ssp245 contain an unremoved drift component that must be stated."
keywords: [drift, 漂移, model drift, 模式漂移, piControl, spin-up, 起转, deep ocean trend, 深海趋势, dedrift, 去漂移, global mean temperature, thetaoga, heat uptake, 热吸收]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, so, o2, dissic]
tools: [fit_1D_trend, fit_trend_map, compute_basic_stat]
severity: should
related: [ocean-heat-and-salt-content, historical-ssp245-member-continuity, model-cesm2-configuration, model-gfdl-esm4-configuration, model-mpi-esm1-2-lr-configuration]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix H12 (pp. 3259–3260): global means as a measure of drift"
    url: https://doi.org/10.5194/gmd-9-3231-2016
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 6 (p. 14): control-run energy imbalance and ocean warming"
    url: https://doi.org/10.1029/2019MS001916
  - citation: "Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1 (GFDL-ESM 4.1). J. Adv. Model. Earth Syst., 12, e2019MS002015."
    locator: "Section 4.1 (p. 14); Section 4.10 (pp. 29–30); Section 4.11 (p. 32)"
    url: https://doi.org/10.1029/2019MS002015
  - citation: "Mauritsen, T., et al. (2019). Developments in the MPI-M Earth System Model version 1.2 (MPI-ESM1.2). J. Adv. Model. Earth Syst., 11, 998–1038."
    locator: "Section 8.1 (p. 1021)"
    url: https://doi.org/10.1029/2018MS001400
---

## Key points

- OMIP: global-mean temperature and salinity time series "provide a measure of
  the model drift". [OMIP-phys H12]
- **CESM2**: in the pre-industrial control the ocean gains about 0.07 W m-2
  (last 500 years). Global-mean ocean potential temperature rises from 3.92 to
  about 4.23 °C by year 1200. Depths below about 2 km warm in all major basins.
  Global-mean salinity falls by about 5 × 10^-5 psu per century. [CESM2 §6]
- **GFDL-ESM4.1**: the volume-averaged ocean in the final configuration warms
  0.017 °C over the first 300 years and 0.013 °C over the last 300 years of
  piControl, "illustrating the slow approach to equilibrium". Interior drift
  concentrates "at the base of the thermocline and deep Atlantic and Indian
  basins". The historical run began after 500 years of spin-up, "so ESM4.1 has
  had more opportunity to accumulate bias on its path to equilibrium".
  [GFDL §4.1, §4.10, §4.11]
- **MPI-ESM1.2**: spin-up continues from predecessor ocean states. The
  cumulative spin-up since the ocean was initialised from observations is
  "well in excess of 20,000 years", aimed at a quasi-stationary control.
  [MPI §8.1]

## Implications for this project

- The local archive has no piControl, so drift cannot be removed by
  subtracting a parallel control segment. Deep (> ~2 km) and whole-ocean trends
  from historical/ssp245 therefore include an unknown drift term. Its size
  differs by model. Going by the papers, CESM2's control warms the global-mean
  ocean by about 0.03 °C per century averaged over its 1200 years (3.92 → 4.23 °C).
  GFDL-ESM4's warms by about 0.004 °C per century over its last 300 years. The
  MPI paper gives no rate. State the limitation,
  and prefer upper-ocean or anomaly-based statements.
- Upper-ocean forced trends over recent decades are usually much larger than
  the drift rates above. Deep-ocean historical trends can be of the same order.
- Mean-state biases in a historical run reflect both model error and
  accumulated drift from the spin-up path. Do not attribute them only to the
  physics of the historical period.

## Source excerpts

> "the global‐mean ocean potential temperature … increases from its initial
> value of 3.92 to about 4.23 °C by Year 1,200 in CESM2(CAM6)" (CESM2, §6)
