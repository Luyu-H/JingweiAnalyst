# Reference answer - T9: PDO pattern and decadal statistics

Window 1950-01 to 2014-12. Leading EOF of North Pacific SST anomalies, 20N-70N / 110E-100W.
Every source, the observations included, regridded to a common 1-degree grid before the EOF.

## Steps 2-3. Variance, amplitude and persistence

| Source | EOF1 variance | amplitude (degC per sd) | lag-1-year autocorrelation |
|---|---|---|---|
| ERSSTv5 | 23.0% | 0.330 | +0.342 |
| CESM2 | 26.1% | 0.364 | +0.223 |
| GFDL-ESM4 | 14.2% | 0.291 | +0.305 |
| MPI-ESM1-2-LR | 22.9% | 0.409 | +0.425 |

The standard deviation of the standardized index is 1 for every source by construction and
is deliberately not reported as a result. The amplitude comes from the regression pattern,
the only output carrying physical units.

## Step 4. Reconstruction against the published index

Pearson r = **0.962** over 780 months.

The published PDO is built from ERSST **v6** while the gridded field used here is **v5**, so
this is not a pure operator validation: the residual contains the difference between two SST
product versions as well as any difference in method.

## Step 5. Pattern agreement

| Model | area-weighted spatial correlation with ERSSTv5 |
|---|---|
| CESM2 | 0.925 |
| GFDL-ESM4 | 0.891 |
| MPI-ESM1-2-LR | 0.750 |

Closest to observations: **CESM2**.

## Step 6. Observed phase transitions

Sign changes of the 9-year centred running mean: [1961, 1962, 1976, 2002].
After merging crossings less than 3 years apart: **[1961, 1976, 2002]**.

## What the comparison shows

- The ERSSTv5 reconstruction correlates with the published series at r = 0.962. That is high, but it is NOT a clean operator validation in the way case T2 is: the published index comes from ERSST v6 and the gridded field here is v5, so the residual mixes a product-version difference with any difference in method. It shows the reconstruction is sound; it cannot separate the two sources of disagreement.
- The models place the mode's pattern well: correlations with observations are CESM2 0.925, GFDL-ESM4 0.891, MPI-ESM1-2-LR 0.750, with CESM2 closest.
- Reporting the standard deviation of the standardized index would have been vacuous - it is 1 for every source by construction. The amplitude has to come from the regression pattern, which is the only output in physical units: observations show 0.330 degC per standard deviation against CESM2 0.364, GFDL-ESM4 0.291, MPI-ESM1-2-LR 0.409.
- The observed phase transitions are [1961, 1976, 2002]. The 1976 one is the documented Pacific climate shift; a 25-year window such as 1990-2014 would not contain it at all, which is why the full record is used.
- The models' index series are NOT expected to track the observed one in time, and figure 2 shows that they do not. These are free-running coupled simulations, so their internal variability has its own phasing and is not synchronised with the real world; only a run nudged to observations could be. What is comparable is the PATTERN and the STATISTICS - the variance share, the amplitude, the persistence - which is what steps 2, 3 and 5 measure. A low time correlation between a model index and the observed one is the expected result and not a model failure; the published series is the only line in figure 2 that should track the observed reconstruction, because it is built from nearly the same data.

## Figures

- `fig1_pdo_patterns.png` - the regression pattern for the four sources on one shared scale.
- `fig2_pdo_timeseries.png` - the standardized index for four sources plus the published
  series, with the smoothed observed curve and the transition years marked.
