# Reference answer - T10: global oxygen minimum zones and the ssp245 change

Thresholds 20 and 60 umol/kg, strictly less than. Model oxygen converted from mol m-3 with an
in-situ density field from each model's own thetao and so. Observations are the WOA23
1971-2000 climate normal. Every volume and inventory on the source's own grid and levels.

CESM2 publishes no historical o2, so it appears only under ssp245.

## Step 1. OMZ volume (1e15 m3)

| Source | < 20 umol/kg | < 60 umol/kg |
|---|---|---|
| WOA23 | 14.2 | 80.3 |
| GFDL-ESM4 / historical | 43.1 | 105 |
| MPI-ESM1-2-LR / historical | 92 | 214 |
| CESM2 / ssp245 | 60.9 | 165 |
| GFDL-ESM4 / ssp245 | 41.3 | 110 |
| MPI-ESM1-2-LR / ssp245 | 89.2 | 212 |

## Step 2. Thickest point at the 60 threshold

| Source | max thickness (m) | latitude | longitude |
|---|---|---|---|
| WOA23 | 3125.0 | +7.50 | +121.50 |
| GFDL-ESM4 / historical | 4343.2 | +13.50 | -93.50 |
| MPI-ESM1-2-LR / historical | 4960.0 | -9.39 | -81.22 |
| CESM2 / ssp245 | 3631.9 | +7.45 | +121.44 |
| GFDL-ESM4 / ssp245 | 4318.2 | +13.50 | -93.50 |
| MPI-ESM1-2-LR / ssp245 | 4950.0 | -9.39 | -81.22 |

## Step 3. Scenario change in OMZ volume, ssp245 minus historical

| Model | < 20 umol/kg | < 60 umol/kg |
|---|---|---|
| GFDL-ESM4 | -4.06% | +4.41% |
| MPI-ESM1-2-LR | -3.03% | -0.72% |

## Step 4. Global dissolved oxygen inventory

| Model and experiment | inventory (mol) |
|---|---|
| GFDL-ESM4 / historical | 2.2e+17 |
| MPI-ESM1-2-LR / historical | 2.25e+17 |
| CESM2 / ssp245 | 2.02e+17 |
| GFDL-ESM4 / ssp245 | 2.15e+17 |
| MPI-ESM1-2-LR / ssp245 | 2.2e+17 |

## Step 5. OMZ volume by basin at the 60 threshold (1e15 m3)

| Source | Pacific | Indian | Atlantic |
|---|---|---|---|
| WOA23 | 67.4 | 11.6 | 1.29 |
| GFDL-ESM4 / historical | 92.8 | 5.43 | 6.39 |
| MPI-ESM1-2-LR / historical | 182 | 19.2 | 11.3 |

Marginal seas are folded into their parent ocean. The Southern and Arctic Oceans are separate basins in the WOA23 mask and are not part of these numbers.

## Step 6. What kind of OMZ, from the observations

Inside the observed OMZ: mean AOU 271.3 umol/kg, mean saturation 11.3%.

## What the comparison shows

- Both historical models overestimate the OMZ. At the 60 threshold the observed volume is 80.3e15 m3 against GFDL-ESM4 105e15, MPI-ESM1-2-LR 214e15; GFDL-ESM4 is the closer of the two. Overestimating the OMZ is a documented common bias.
- Under ssp245 the two models do not agree on the sign of the change: GFDL-ESM4 +4.4% at 60 and -4.1% at 20; MPI-ESM1-2-LR -0.7% at 60 and -3.0% at 20. The threshold matters as much as the model: the suboxic core and the hypoxic envelope need not move together, which is why both thresholds are asked for. Each change compares two periods of one continuous run (r1i1p1f1 in both experiments), but it also contains whatever drift the model's pre-industrial control still carries, which cannot be removed without a piControl run.
- Inside the observed OMZ the mean apparent oxygen utilization is 271.3 micromoles per kilogram and the mean saturation is 11.3 per cent. A high AOU with a low saturation says the zone is set by the respiration of sinking organic matter in water that has been out of contact with the atmosphere for a long time, rather than by cold, poorly ventilated water that was already low in oxygen when it left the surface.
- CESM2 could not be assessed against observations or given a scenario change at all, because the historical experiment publishes no o2 for it. Its ssp245 OMZ is reported and drawn, but it cannot be told whether that state is a change or a pre-existing bias.

## Figures

- `fig1_omz_thickness_historical.png` - OMZ thickness for the two historical models and
  WOA23 on one shared scale.
- `fig2_omz_thickness_change.png` - the ssp245 minus historical change for the two models
  that have both experiments, plus CESM2's ssp245 state on its own scale.
