---
id: model-gfdl-esm4-reported-biases
title: "GFDL-ESM4 ocean biases and variability as reported in its description paper"
type: dataset
summary: "GFDL-ESM4.1 has a small global SST bias (−0.07 °C vs WOA13) but cold subtropical gyres and subpolar North Atlantic, warm eastern boundaries and Southern Ocean, shallow tropical/subtropical mixed layers, a diffuse thermocline, a shallow AMOC, too much abyssal remineralisation, and a slightly strong, ~4-year ENSO."
keywords: [GFDL-ESM4, ESM4.1, bias, 偏差, SST bias, SSS bias, 盐度偏差, mixed layer depth bias, 混合层偏差, thermocline, 温跃层, AMOC, ENSO, oxygen, nutrient, remineralization, 再矿化, model evaluation, 模式评估]
datasets: [cmip6_historical, woa23, glodap]
models: [GFDL-ESM4]
variables: [tos, thetao, so, mlotst, o2, po4, dissic, talk]
tools: [compute_difference_metrics, compute_enso_nino_index]
severity: info
related: [model-gfdl-esm4-configuration, overturning-from-mass-transport, model-drift-in-deep-trends, cmip6-mlotst-definition]
sources:
  - citation: "Dunne, J. P., et al. (2020). The GFDL Earth System Model Version 4.1 (GFDL-ESM 4.1): Overall coupled model description and simulation characteristics. J. Adv. Model. Earth Syst., 12, e2019MS002015."
    locator: "Section 4.9 (pp. 24–28); 4.10 (pp. 29–32); 4.11 (pp. 32–33); 4.12 (pp. 34–35)"
    url: https://doi.org/10.1029/2019MS002015
---

Each number below comes with the period and reference data set the paper used.
It is context for interpreting results, not a target value.

## Key points

**Surface (against WOA13)** [GFDL §4.9]
- SST: RMS 0.68 °C and mean bias −0.07 °C. Cold biases in the **subtropical
  gyres** and an "**extreme cold bias** in ESM4.1 subpolar North Atlantic". Warm
  biases in **eastern boundary current** regions, the **Southern Ocean** and the
  Labrador Sea.
- SSS: RMSE 0.47 ppt. High bias in the **Arctic**. Low bias in the equatorial and
  south-east Atlantic (Benguela). Low bias in Oceania. High bias off the north-east US
  coast and in the Bay of Bengal (weak freshwater input).
- MLD (0.03 kg m-3 criterion relative to the upper 10 m, monthly min and max
  against WOA13): pattern r² = 0.68 (min) and 0.51 (max). There is a **shallow bias
  in the tropics and subtropics** against the 15–25 m bands in WOA13. Deep winter
  mixing occurs along east Siberia, where observations show none.

**Interior** [GFDL §4.10]
- "Diffuse thermocline" bias: warm at the thermocline base and in the deep
  Atlantic and Indian basins, cold in the upper thermocline (notably the
  Pacific).
- NADW is too warm and salty, leaving a salt deficit in the Pacific
  thermocline. Red Sea (200–400 m) and Mediterranean (600–1000 m) overflows are
  too strong and too shallow.
- AMOC at 24°N is shallower than observed and ends at 3000 m. AABW transport
  into the southern Indo-Pacific is weak.

**Biogeochemistry (tropical 30°S–30°N, 1990–2009)** [GFDL §4.11]
- Good agreement for PO4, O2, DIC and Alk−DIC in the upper
  nutricline/oxycline. In the tropical Atlantic, PO4 and O2 are underestimated
  in the main nutricline/oxycline, and there are mid-depth biases linked to the
  shallow AMOC.
- "the signature of **excess remineralization in the abyss** in all four
  variables" in both basins.

**ENSO** [GFDL §4.12]
- Slightly **stronger** than observed, with a **longer period (~4 years**
  against ~3.3 in observations and CM4.0). The warmest eastern-equatorial
  anomalies sit slightly farther west. Tropical mixed layers are shallower than
  observed, and there is an equatorial Pacific subsurface cold bias.

## Implications for this project

- Tropical and subtropical MLD biases and the diffuse thermocline affect
  nutricline and oxycline depths. A GFDL-ESM4 nutricline or oxygen bias is
  consistent with these documented physical biases.
- Abyssal DIC and nutrients that are too high and O2 that is too low are
  expected in GFDL-ESM4.

## Source excerpts

> "In both the Atlantic and Indo‐Pacific, ESM4.1 exhibits the signature of
> excess remineralization in the abyss in all four variables." (GFDL, §4.11)
