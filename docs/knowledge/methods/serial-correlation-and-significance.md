---
id: serial-correlation-and-significance
title: "Autocorrelated ocean data: significance of trends and correlations, and validation of fitted relations"
type: method
summary: "Monthly ocean series and neighbouring ocean cells are strongly autocorrelated, so naive sample sizes overstate the information: trend intervals and correlation p-values need an effective sample size, a trend question needs an explicit straight-line fit quoted with its interval, shared trends inflate correlations, and fitted relations need spatially blocked validation."
keywords: [trend, 趋势, linear trend, 线性趋势, confidence interval, 置信区间, significance, 显著性, p-value, serial correlation, autocorrelation, 自相关, effective sample size, 有效样本量, degrees of freedom, 自由度, correlation, 相关, lag correlation, 滞后相关, detrend, 去趋势, regression, 回归, Theil-Sen, per decade, 每十年, stippling, cross-validation, 交叉验证, collinearity, 共线性]
tools: [fit_1D_trend, fit_trend_map, analyze_correlation, fit_multivariate_relation]
severity: should
related: [anomaly-baselines-and-climate-indices, eof-modes-and-pdo-index]
sources:
  - citation: "Bretherton, C. S., Widmann, M., Dymnikov, V. P., Wallace, J. M., and Bladé, I. (1999). The effective number of spatial degrees of freedom of a time-varying field. J. Climate, 12, 1990–2009."
    locator: "Section 5, Relation to effective sample size for serial correlation, Eqs. (30)–(31) (pp. 2003–2004)"
    url: https://doi.org/10.1175/1520-0442(1999)012<1990:TENOSD>2.0.CO;2
---

## Key points

- Serially correlated observations carry fewer independent samples than they
  have time steps. For red noise with lag-one autocorrelation r, the effective
  sample size for the sample mean (order 1) is T(1 − r)/(1 + r). For the
  correlation between two series with lag-one autocorrelations r_X and r_Y it
  is T(1 − r_X r_Y)/(1 + r_X r_Y) (Eq. 31). [Bretherton et al. 1999, §5]
- In practice the autocorrelation has to be estimated from the one series at
  hand, so the effective sample size is itself uncertain, the more so when the
  series has substantial power at periods near its length. [Bretherton et al.
  1999, §5]

## Implications for this project

- A monthly ocean series is strongly autocorrelated, so its number of steps
  overstates the independent information by a large factor. Uncorrected, a test
  uses far too many degrees of freedom: weak correlations come out highly
  significant, and a monthly trend map is significant almost everywhere, which
  says more about the persistence of the ocean than about any trend.
  `fit_1D_trend`, `fit_trend_map` and `analyze_correlation` apply an effective
  sample size along time by default and return the unadjusted values alongside;
  quote the corrected ones.
- For the question "is this quantity changing, and how fast", fit a straight
  line explicitly: `method='linear'` (the standard choice) or
  `method='theil_sen'` (the same question, resistant to outliers and non-normal
  residuals). Automatic model selection answers "which curve describes this
  series", not "what is the trend": a seasonal series will legitimately select
  the sinusoid. Quote the slope with its confidence interval, and call it
  significant only if the interval excludes zero.
- A trend map is a straight-line fit at every cell, so every cell is estimated
  the same way and read as a rate. Over whole years the raw field can be
  fitted, but its seasonal cycle inflates the lag-1 autocorrelation and makes
  the test conservative; removing the seasonal cycle first gives a sharper test
  and the same trend. Draw the slope everywhere and stipple the significant
  cells, rather than blanking the insignificant ones: a coherent pattern that
  misses significance is still worth seeing. A low R² with a significant slope
  is normal for a noisy field, a real drift buried in variability.
- Correlations along time: when both series contain a long-term trend, detrend
  them first (`detrend=True`), because two independently rising series
  correlate through their trends alone. That is for asking whether two
  different quantities are related. When two estimates of the same quantity are
  compared, such as an index rebuilt from a gridded field against the
  published index, or one series against its counterpart in another product,
  correlate them as they are: their shared trend is part of the agreement being
  checked, not a spurious link. When the effective sample size is much
  smaller than the number of pairs, the correlation rests on that many
  independent observations; a correlation from a handful of points is unstable.
- Across space no such correction exists for a pattern correlation: read its
  coefficient, not its p-value.
- Fitted relations between fields (`fit_multivariate_relation`): neighbouring
  ocean cells are near-duplicates, so a random train/test split scores a model
  against its own neighbours and reports a skill it does not have. Hold out
  whole spatial blocks larger than the field's correlation length scale; a
  random split suits only data with no spatial structure. Ocean tracers are
  often strongly collinear, which makes individual coefficients unreliable even
  when the fit is good; the variance inflation factors show when that happens.

## Common mistakes

- Quoting a trend's p-value or interval from the naive sample size of a monthly
  series.
- Letting automatic model selection answer a trend question.

## Source excerpts

> "For instance, for red noise with lag-one autocorrelation r, the ESS of order
> 1 is T* = T(1 − r)/(1 + r) [...], while (31) shows that the ESS of order 2 is
> T* = T(1 − r²)/(1 + r²), which is up to twice as large as the Jones–Leith ESS
> when r approaches 1, that is, if the autocorrelation is strong."
> (Bretherton et al. 1999, §5)
