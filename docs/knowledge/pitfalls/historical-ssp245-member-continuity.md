---
id: historical-ssp245-member-continuity
title: "Joining historical and ssp245 only makes one continuous run when the member matches the scenario's parent"
type: pitfall
summary: "ssp245 branches from historical (parent_experiment_id) of the member named in parent_variant_label; locally CESM2 ssp245 is r10i1p1f1 continuing historical r10i1p1f1, while the local CESM2 historical is r1i1p1f1, so historical→ssp245 changes for CESM2 mix forced change with internal variability between members."
keywords: [historical, ssp245, scenario, 情景, concatenate, 拼接, ensemble member, 集合成员, variant label, r1i1p1f1, r10i1p1f1, parent_variant_label, branch_time, future change, 未来变化, 2015, continuity, 连续性]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, so, o2, dissic, ph]
tools: [load_regular_data, fit_1D_trend, fit_trend_map, compute_difference_metrics]
severity: must
related: [cmip6-experiments-historical-ssp245, model-drift-in-deep-trends]
sources:
  - citation: "WCRP-CMIP CMIP6_CVs, CMIP6_experiment_id.json, CV collection version 6.2.60.0."
    locator: "experiment_id historical (parent piControl/past1000/past2k), ssp245 (parent historical)"
    url: https://github.com/WCRP-CMIP/CMIP6_CVs
  - citation: "Danabasoglu, G., et al. (2020). The Community Earth System Model Version 2 (CESM2). J. Adv. Model. Earth Syst., 12, e2019MS001916."
    locator: "Section 4 (p. 10): historical members initialised from different piControl years"
    url: https://doi.org/10.1029/2019MS001916
---

## Key points

- `ssp245` has `parent_experiment_id: historical` and runs 2015–2100.
  `historical` runs 1850–2014 with parent `piControl` (or past1000/past2k).
  [CVs experiment_id]
- The 11 CESM2(CAM6) historical members start from **different years of the
  piControl run** (years 501, 601, 631, …, 871). They are independent
  realisations of internal variability, not copies. [CESM2 §4]

## Checked in the local files

These are the global attributes of the first `thetao`/`o2` file in each
directory, read on 2026-09-24:

| model | historical variant (parent) | ssp245 variant → parent_variant_label |
|---|---|---|
| CESM2 | r1i1p1f1 (piControl r1i1p1f1) | **r10i1p1f1 → historical r10i1p1f1** |
| GFDL-ESM4 | r1i1p1f1 | r1i1p1f1 → historical r1i1p1f1 |
| MPI-ESM1-2-LR | r1i1p1f1 | r1i1p1f1 → historical r1i1p1f1 |

The local archive does **not** contain CESM2 historical r10i1p1f1.

## Implications for this project

- For GFDL-ESM4 and MPI-ESM1-2-LR, historical r1 followed by ssp245 r1 is one
  continuous trajectory. A time series across 2014/2015 is valid.
- For CESM2, historical (r1) and ssp245 (r10) are **different realisations**.
  - A future-minus-historical difference (a late-21st-century window minus a
    late-20th-century window)
    contains forced change **plus** the difference in internal variability
    between two members. For long averaging windows and large forced signals
    the extra term is small, but it is not zero. State it.
  - A trend fitted across the join, or a plot that joins the two, shows a
    spurious discontinuity at 2015. Do not fit across the join for CESM2.
- The same mismatch occurs within one experiment: CESM2 historical `phycos` is
  r6 and MPI historical `ph` is r10 (see `configs/datasets.yaml`). Joint
  point-by-point analysis with other variables mixes members.

## Common mistakes

- Concatenating CESM2 historical and ssp245 and fitting one trend over
  1850–2100.
- Reporting a CESM2 scenario change without noting that the baseline comes
  from another member.

## Source excerpts

> "The initial conditions for the 11 members of the CESM2(CAM6) historical
> simulations come from Years 501, 601, 631, 661, 691, 721, 751, 781, 811, 841,
> and 871 of the corresponding PI control simulation" (CESM2, §4)
