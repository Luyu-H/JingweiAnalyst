---
id: model-cesm2-reported-biases
title: "CESM2 ocean biases and variability as reported in its description paper"
type: dataset
summary: "CESM2 is ~0.4 °C too warm in global-mean SST vs ERSSTv5 (1985–2014) with strong eastern-boundary upwelling warm biases and a Gulf Stream dipole, a ~24 Sv AMOC (above RAPID's 17 Sv), a very weak AABW cell and weak deep ventilation; ENSO patterns are realistic."
keywords: [CESM2, bias, 偏差, SST bias, 海温偏差, warm bias, upwelling, 上升流, Gulf Stream, AMOC, AABW, ventilation, 通风, ENSO, Niño3.4, model evaluation, 模式评估]
datasets: [cmip6_historical, ersst]
models: [CESM2]
variables: [tos, thetao, vo]
tools: [compute_difference_metrics, compute_enso_nino_index]
severity: info
related: [model-cesm2-configuration, overturning-from-mass-transport, model-drift-in-deep-trends]
sources:
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 7 (pp. 14–15); Section 8.4 and Fig. 14 caption (p. 21); Section 9 (pp. 22–23)"
    url: https://doi.org/10.1029/2019MS001916
---

The paper evaluates the CESM2(CAM6) ensemble. Each number below comes with the
period and reference data set the paper used. It is context for interpreting
results, not a target value.

## Key points

**SST (1985–2014 mean, top 0–10 m, against ERSSTv5)** [CESM2 §7]
- Global-mean SST is **+0.39 °C** warmer than observed (the CESM2(WACCM6)
  variant is +0.33 °C). The RMS difference is 0.77 °C. The warm offset largely
  reflects a PI control that is also warm (+0.33 °C versus ERSSTv5 pre-industrial).
- Persistent patterns:
  - a warm–cold dipole in the North Atlantic from the Gulf Stream / North
    Atlantic Current separation and path;
  - **warm biases in the eastern boundary upwelling regions** off California,
    South America and South Africa, larger than in CESM1 and hard to reduce at
    1° resolution.

**ENSO (piControl composites)** [CESM2 §8.4]
- Spatial patterns are "generally well simulated". The DJF Niño3.4 El Niño
  minus La Niña composite is 1.38 °C (CAM6) against 1.26 °C (WACCM6). The two
  differ significantly, "for reasons that remain to be understood".
- The composite uses normalised December Niño3.4 with a ±1 standard deviation
  threshold. Observations are ERSSTv5 1920–2017.

**AMOC (1995–2014, Eulerian mean)** [CESM2 §9]
- The NADW maximum is **about 24 Sv**. At 26.5°N the model exceeds RAPID
  (17.0 ± 3.3 Sv, 2004–2017) by about 2–3 Sv. The NADW cell is shallower than
  in CESM1.
- The **AABW cell is about 0.5 Sv**, against an observed 1–5 Sv at 24.5°N. The
  result is "weak ventilation at depth evident in tracer distributions
  particularly in the North Pacific".

## Implications for this project

- A warm CESM2 SST bias and warm eastern-boundary upwelling regions are
  expected features of the model. They do not by themselves indicate a
  processing error.
- Weak deep ventilation implies that deep-water tracer ages and
  remineralisation signals (lower O2, higher nutrients and DIC at depth,
  notably in the North Pacific) may be biased. This is an inference from the
  paper's statement. Check it against the data before asserting it.

## Source excerpts

> "The warm biases in the Eastern boundary upwelling regions, that is, off the
> coasts of California, South America, and South Africa, are more pronounced
> in CESM2(CAM6) and CESM2(WACCM6) than in CESM1(LENS)." (CESM2, §7)
