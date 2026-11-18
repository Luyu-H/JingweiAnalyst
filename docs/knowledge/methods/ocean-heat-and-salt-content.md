---
id: ocean-heat-and-salt-content
title: "Computing ocean heat content and salt content from CMIP6 thetao and so"
type: method
summary: "Model heat content is c_p × prognostic temperature × ρ dV with the model's own constant heat capacity (TEOS-10 value 3991.87 J kg-1 K-1) and, for Boussinesq models, a constant reference density; changes are scale-independent, absolute values are not; 1 W m-2 over a decade warms the global ocean by ~0.02 °C."
keywords: [ocean heat content, 海洋热含量, OHC, heat capacity, 比热容, cp, 3991.87, reference density, 参考密度, rhozero, 1035, salt content, 盐含量, heat uptake, 热吸收, warming trend, 增温趋势]
datasets: [cmip6_historical, cmip6_ssp245, woa23]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, so, volcello, areacello]
tools: [load_volumes_for_cmip6, compute_basic_stat, fit_1D_trend, fit_trend_map]
severity: should
related: [cmip6-thetao-potential-temperature, cmip6-so-practical-salinity, model-drift-in-deep-trends, omip-weighted-means-and-integrals]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix D3–D3.3 (pp. 3246–3247); D5 (p. 3248); Appendix E (pp. 3248–3249); G3–G4 (p. 3251); H12.2, Eqs. H30–H31 (p. 3260)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
---

## Key points

- TEOS-10 heat content of a parcel is h° ρ dV with h° = c_p° Θ.
  c_p° = 3991.86795711963 J kg-1 K-1, and "there is no physics in c_p° beyond
  roughly four significant digits". [OMIP-phys D3, Eqs. D1–D3]
- Models measure heat content as **c_p × prognostic temperature × ρ dV** with a
  **constant** c_p. That c_p "is not always equal to the TEOS-10 recommended value".
  Pre-TEOS-10 models use their prognostic potential temperature. They should
  not diagnose Conservative Temperature for this. [OMIP-phys D3.2–D3.3, G4]
- Boussinesq models replace ρ by a constant ρ_o. Choices differ between
  models: "ρ_o = 1035 kg m−3 is a sensible choice", but "early versions of the GFDL
  ocean model … set ρ_o = 1000 kg m−3". [OMIP-phys D3.1, G3]
- Ocean heat content depends on the temperature scale (°C versus K). Its **time
  evolution does not, as long as ocean mass is constant**. [OMIP-phys E]
- Scale estimate: V ρ_o c_p ≈ 5.4 × 10^24 J °C-1 and A = 3.6 × 10^14 m2, so
  Δ⟨θ⟩ per decade ≈ 0.02 · Q_H. A 1 W m-2 ocean-mean heating over a decade gives
  ≈ 0.02 °C per decade, or 0.2 °C per century. [OMIP-phys H12.2]
- Salt content = cell mass × S_A formally. Pre-TEOS-10 models use cell mass ×
  S_P. [OMIP-phys D5]

## Implications for this project

- `cpocean` and `rhozero` are not in the local archive. Use c_p = 3991.87 and
  ρ_o = 1035 kg m-3 **and state both**. Report heat content **changes or
  anomalies** relative to a stated baseline, not absolute values. Changes are
  robust to the scale offset and much less sensitive to the constants.
- OHC(t) = c_p ρ_o Σ thetao · volcello over the layer of interest (J).
  For a 0–700 m or 0–2000 m layer, cut at the layer boundary using the level
  bounds. Include only the part of the straddling cell above the cut.
- Convert a heat-content trend to W m-2 by dividing by the seconds per year and
  by the **Earth's** or the **ocean's** surface area. Say which one; the
  literature uses both.
- Deep-layer heat content trends in a historical run include model drift
  (see `model-drift-in-deep-trends`).

## Source excerpts

> "For example, with a 1 W m−2 ocean area average heating of the ocean over the
> course of a decade, we expect a global mean temperature trend of roughly
> 0.02 °C per decade, or 0.2 °C per century." (OMIP-phys, H12.2)
