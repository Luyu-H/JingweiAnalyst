---
id: anomaly-baselines-and-climate-indices
title: "Anomalies and climate indices: which climatology to remove, which baseline, which season convention"
type: method
summary: "An anomaly from a field's own climatology removes its bias (what an index, a trend or a comparison between two products needs), one from an observed climatology keeps it; an anomaly is defined only relative to its baseline period, so series on different baselines differ by an offset; DJF and JFM seasons must not be mixed; events and outliers are judged on anomalies, with persistence."
keywords: [anomaly, 距平, climatology, 气候态, baseline, 基准期, base period, 参考期, climate normal, index, 指数, ENSO, El Niño, La Niña, Niño 3.4, Nino 3.4, Nino3.4, ONI, CPC, event, 事件, persistence, 持续性, season, 季节, DJF, JFM, seasonal cycle, 季节循环, outlier, 异常值, marine heatwave, 海洋热浪]
datasets: [cpc_indices, ersst, woa23, cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [tos, thetao, so]
tools: [compute_anomaly_vs_own_climatology, compute_anomaly_vs_woa23, compute_enso_nino_index, temporal_aggregate, detect_threshold_events, detect_outliers, load_regular_data]
severity: must
related: [cmip6-tos-bulk-top-cell, serial-correlation-and-significance, eof-modes-and-pdo-index]
sources:
  - citation: "NOAA Climate Prediction Center. Changes to the Oceanic Niño Index (ONI)."
    locator: "Explanation of the base-period strategy (first three paragraphs)"
    url: https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_change.shtml
  - citation: "Reagan, J. R., et al. (2024). World Ocean Atlas 2023: Product Documentation. NOAA National Centers for Environmental Information, Silver Spring, MD."
    locator: "File naming conventions, averaging period codes 13–16 (p. 13)"
    url: https://www.ncei.noaa.gov/data/oceans/woa/WOA23/DOCUMENTATION/WOA23_Product_Documentation.pdf
---

## Key points

- CPC changed the base period of its Oceanic Niño Index because, with a warming
  trend in the Niño-3.4 region since 1950, episodes defined against a fixed base
  period "are increasingly incorporating longer-term trends that do not reflect
  interannual ENSO variability". The ONI now uses centred 30-year base periods,
  each applied to successive 5-year spans of the record. [CPC ONI change]
  The same index therefore exists on several baselines; `configs/datasets.yaml`
  lists the local CPC products and the baseline of each.
- The WOA23 seasonal climatologies are calendar quarters: winter January–March,
  spring April–June, summer July–September, autumn October–December. [WOA23
  documentation]

## Implications for this project

- Which climatology is removed decides what an anomaly means.
  - Against the field's **own** climatology (`compute_anomaly_vs_own_climatology`,
    and `compute_enso_nino_index` on model output) whatever bias the field
    carries cancels. This is what an index, a trend, a composite or a comparison
    between two products needs.
  - Against the **observed** climatology (`compute_anomaly_vs_woa23`) the
    field's bias against observations stays inside the anomaly. This answers
    "how far is this field from the observed ocean". An observed climatology
    forced on a model index puts the model's mean-state bias into the index
    (see `cmip6-tos-bulk-top-cell`).
- The baseline is part of the answer. An anomaly is defined only relative to
  the period its climatology came from, and two series built on different
  baselines differ by an offset (one constant per calendar month) even when
  their variability is identical. Compute an index on the baseline of any
  published series it is compared with (`base_period`), and state the
  climatology and baseline used; the tools record both in the result's attrs.
  The CPC products in the registry differ in baseline and are not
  interchangeable.
- Niño 3.4 is the standard region for ENSO monitoring. A monthly index is
  matched against the climatology of its own calendar month, which removes the
  seasonal cycle.
- Seasons: DJF/MAM/JJA/SON is the usual (meteorological) climate definition,
  while the WOA23 seasonal fields are the calendar quarters JFM/AMJ/JAS/OND.
  The two are offset by a month and must not be mixed in one comparison
  (`season_convention` of `temporal_aggregate`; `freq='season'` of
  `compute_anomaly_vs_woa23` uses the WOA23 quarters).
- Events are defined operationally by persistence, not by eye: a warm or cold
  phase of a climate index, a marine heatwave or a hypoxic spell counts only
  when the series stays beyond the threshold for a minimum number of
  consecutive steps. A single step beyond the threshold is noise;
  `min_duration` of `detect_threshold_events` is what separates an event from
  a wobble.
- Scan a time series for outliers on its anomalies, not on its raw values. A
  raw ocean series is dominated by its seasonal cycle, so every summer maximum
  sits in the tail of the distribution and the "outliers" returned are simply
  the warm season.

## Common mistakes

- Reading the offset between a reconstructed index and a published one built on
  another baseline as an error of the reconstruction.
- Comparing seasonal means built on the DJF and the JFM conventions.

## Source excerpts

> "In order to remove this warming trend, CPC is adopting a new strategy to
> update the base period." (CPC, Changes to the ONI)

> "Season 13 – North Hemisphere winter (January - March); Season 14 – North
> Hemisphere spring (April - June); Season 15 – North Hemisphere summer (July -
> September); Season 16 – North Hemisphere autumn (October - December)" (WOA23
> Product Documentation, p. 13)
