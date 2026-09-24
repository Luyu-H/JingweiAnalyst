# Reference answer - T8: carbon system against GLODAPv2

Models sampled for calendar year 2002 (day-weighted annual mean). GLODAPv2.2016b mapped
climatology, input data spanning 1972-2013. Model dissic and talk converted from mol m-3 to
umol/kg with an in-situ density field from each model's own thetao and so.

## Steps 1-2. Surface (0-50 m) volume-weighted global means

| Variable | GLODAP | CESM2 | GFDL-ESM4 | MPI-ESM1-2-LR |
|---|---|---|---|---|
| DIC (umol/kg) | 2036 | 1997 (-39) | 2059 (+23) | 1974 (-62) |
| TAlk (umol/kg) | 2304 | 2264 (-40) | 2336 (+31) | 2217 (-87) |
| pH | 8.080 | 8.080 (-0.001) | 8.088 (+0.008) | 8.059 (-0.022) |

Model columns show the value with its difference from GLODAP in brackets.

## Step 3. Surface DIC pattern agreement (41229 shared cells)

| Model | bias (umol/kg) | RMSE (umol/kg) | spatial correlation |
|---|---|---|---|
| CESM2 | -37 | 60 | 0.828 |
| GFDL-ESM4 | +25 | 65 | 0.709 |
| MPI-ESM1-2-LR | -59 | 83 | 0.736 |

## Step 4. Global mean DIC bias by depth (umol/kg)

| Model | 0 m | 100 m | 500 m | 1000 m | 2000 m |
|---|---|---|---|---|---|
| CESM2 | -29 | -46 | -12 | +2 | -2 |
| GFDL-ESM4 | +29 | +19 | -11 | -1 | -11 |
| MPI-ESM1-2-LR | -58 | -59 | -78 | -67 | -61 |

## Step 5. How much the sampled year matters

| Model | 1971-2000 mean | 2002 | change | bias vs GLODAP | change / bias |
|---|---|---|---|---|---|
| CESM2 | 1987 | 1997 | +9 | -39 | -0.24 |
| GFDL-ESM4 | 2047 | 2059 | +12 | +23 | +0.53 |
| MPI-ESM1-2-LR | 1963 | 1974 | +11 | -62 | -0.18 |

## What the comparison shows

- The model fields had to be converted with an in-situ density field rather than a constant. DIC and alkalinity are absolute concentrations, so unlike a ratio or a depth the conversion moves every number, and a constant near-surface density would err by about the size of the biases reported here.
- On surface DIC the smallest bias against GLODAP is GFDL-ESM4 (+23 micromoles per kilogram).
- Step 5 is what decides how much of the DIC bias can be read as model error. The period change each model shows between its 1971-2000 mean and 2002 is CESM2 +9, GFDL-ESM4 +12, MPI-ESM1-2-LR +11 micromoles per kilogram, against biases of CESM2 -39, GFDL-ESM4 +23, MPI-ESM1-2-LR -62. Where those are comparable in size the bias is sensitive to the year sampled.
- That period change is mostly the forced rise of surface DIC as the ocean takes up anthropogenic CO2 between the two periods, plus internal variability. It is not model drift in the CMIP6 sense - the residual trend of the pre-industrial control run - which cannot be measured here without a piControl run.
- Alkalinity is far less sensitive to anthropogenic CO2 than DIC, so its bias cannot be explained by the period mismatch and is the cleaner indicator of a structural model error.
- The DIC and alkalinity biases carry the SAME SIGN in every model - CESM2 -39 and -40, GFDL-ESM4 +23 and +31, MPI-ESM1-2-LR -62 and -87 micromoles per kilogram - which is why pH comes out far better than either variable alone. pH responds to the balance between the two, so errors that move DIC and alkalinity together largely cancel in it: CESM2 is off by -39 in DIC and -40 in alkalinity yet reproduces pH to 0.001. Agreement in pH alone is therefore NOT evidence that a model's carbon system is right, and the three variables have to be read together - which is what this case is built to show.

## Figures

- `fig1_dic_bias_maps.png` - 0-50 m DIC bias against GLODAP, three panels on one shared scale.
- `fig2_dic_bias_profile.png` - global mean DIC bias against depth, one line per model.
