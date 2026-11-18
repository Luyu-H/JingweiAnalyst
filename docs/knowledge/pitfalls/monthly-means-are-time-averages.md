---
id: monthly-means-are-time-averages
title: "CMIP6 monthly fields are time averages of every model time step, so nonlinear diagnostics from them are not the model's own"
type: pitfall
summary: "Monthly CMIP6 values are averages over all time steps (products averaged as products, max/min from snapshots); any nonlinear diagnostic (MLD, density, N², variance, extremes) computed offline from monthly means differs from the online one and from the mean of the diagnostic."
keywords: [monthly mean, 月平均, time average, 时间平均, aliasing, 混叠, nonlinear, 非线性, offline diagnostic, 离线诊断, online, snapshot, 瞬时, variance, 方差, climatology first, 先平均后计算]
datasets: [cmip6_historical, cmip6_ssp245]
models: [CESM2, GFDL-ESM4, MPI-ESM1-2-LR]
variables: [thetao, so, mlotst, mlotstmax, mlotstmin, o2]
tools: [temporal_aggregate, compute_mixed_layer_depth, compute_potential_density_anomaly, compute_nutricline_depth]
severity: must
related: [cmip6-mlotst-definition, cf-cell-methods, derive-unavailable-diagnostics]
sources:
  - citation: "Griffies, S. M., et al. (2016). OMIP contribution to CMIP6: experimental and diagnostic protocol for the physical component of the Ocean Model Intercomparison Project. Geosci. Model Dev., 9, 3231–3296."
    locator: "Appendix A1 (p. 3239); Appendix C2 (p. 3244); Appendix F3 (p. 3250); Appendix H24.2 (p. 3263); Table H1 caption (p. 3254)"
    url: https://doi.org/10.5194/gmd-9-3231-2016
---

## Key points

- Each requested period (day, month, year) "is sampled as a time average. Time
  averages include all model time steps over the given period; there is no
  sub-sampling. Products of time-dependent fields are time-averaged as a
  product". [OMIP-phys A1]
- `mlotstmax`/`mlotstmin` are the extremes over snapshots within the month.
  Adjacent points "need not have extrema registered at the same time step".
  [OMIP-phys A1, H24.1]
- Squared fields (e.g. `tossq`, `mlotstsq`) are accumulated online so that
  variance can be computed. Squaring a monthly mean does not give them.
  [OMIP-phys Table H1 caption, H25]
- MLD is recommended online "in order to avoid aliasing". [OMIP-phys H24.2]
- For extensive quantities, "experience at GFDL with Boussinesq z* models
  suggests that errors are minor when working offline using time averaged
  thickness and time averaged tracer concentration". Errors are "generally
  unacceptable" for isopycnal models. [OMIP-phys F3]

## Implications for this project

- f(monthly mean) ≠ monthly mean of f for any nonlinear f. Examples are density
  from T and S, MLD by a threshold crossing, depth of a maximum gradient, a
  threshold exceedance (hypoxia) or a ratio of two tracers. The order of
  operations therefore changes the answer and must be chosen on purpose:
  - **diagnostic first, then average** reproduces the mean of the diagnostic
    and matches how archived online diagnostics are defined;
  - **average first, then diagnostic** gives the diagnostic of the mean state.
    This is the like-for-like choice when the reference exists only as a
    climatology.
  State which order was used and why.
- Monthly variance within a month cannot be recovered from monthly means. Only
  interannual or inter-month variance can.
- Averaging monthly fields over a period needs day weights in the model's own
  calendar. `temporal_aggregate` does this, and equal weights do not: calendar
  months are 28 to 31 days long, so an equal-weight average of monthly means is
  not a time average. `compute_basic_stat(op='mean', dim='time')` gives every
  month the same weight. `temporal_aggregate` is the only tool that weights by
  month length and so the only correct way to average monthly data over a
  period that spans different calendar months; its `freq='mean'` gives the
  mean state of a period, the field a bias or an RMSE is measured on. Day
  weights do not enter a monthly climatology, each calendar month averaged
  across the years of a baseline (what an anomaly or a climate index is taken
  against): that is the plain mean of that month's values over the years, as
  climate normals are defined, and `compute_anomaly_vs_own_climatology` and the
  index tools build it so.

## Common mistakes

- Calling an MLD from climatological T/S profiles "the model's climatological
  MLD". It is a different quantity from the time mean of the monthly MLDs.
- Reporting a variance of monthly means as if it were sub-monthly variability.

## Source excerpts

> "Time averages include all model time steps over the given period; there is
> no sub-sampling. Products of time-dependent fields are time-averaged as a
> product, again using all model time steps to build the average."
> (OMIP-phys, A1)
