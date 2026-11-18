**Assessment of sea-surface salinity (SSS) in CESM2, GFDL-ESM4 and MPI-ESM1-2-LR**

SSS was defined as the volume-weighted mean salinity of model/WOA levels whose depths lie within 0–10 m, using each model’s own `volcello` and `areacello`, and compared on the WOA23 1° grid.

### 1. Mean state, 1971–2000

WOA23 global mean SSS: **34.63 psu**.

| Model | Global mean SSS (psu) | Bias vs WOA23 (psu) | RMSE vs WOA23 (psu) | Spatial correlation |
|---|---:|---:|---:|---:|
| CESM2 | 24.50 | −10.13 | 1.89 | 0.751 |
| GFDL-ESM4 | 34.56 | −0.07 | 0.79 | 0.920 |
| MPI-ESM1-2-LR | 34.36 | −0.27 | 0.89 | 0.906 |

GFDL-ESM4 and MPI-ESM1-2-LR have physically plausible global means and small global biases. GFDL-ESM4 has the lowest RMSE and highest spatial correlation, followed closely by MPI-ESM1-2-LR. The CESM2 result is not physically credible: a global SSS of 24.50 psu and a bias of −10.13 psu indicate that land or non-ocean grid cells are very likely entering the area-weighted global mean or regridding. Therefore the CESM2 mean-state bias/RMSE/correlation should be treated as unreliable until the land/area masking is fixed.

### 2. Linear trends of annual global-mean SSS

Units: **psu decade⁻¹**. Confidence intervals are 95%.

| Model | Period | Slope | 95% CI | R² | CI excludes zero? |
|---|---|---:|---:|---:|---|
| CESM2 | Historical 1950–2014 | −1.91×10⁻³ | [−2.76×10⁻³, −1.07×10⁻³] | 0.245 | Yes |
| GFDL-ESM4 | Historical 1950–2014 | −1.10×10⁻³ | [−2.15×10⁻³, −3.76×10⁻⁵] | 0.064 | Yes |
| MPI-ESM1-2-LR | Historical 1950–2014 | −6.87×10⁻⁴ | [−2.35×10⁻³, +9.77×10⁻⁴] | 0.011 | No |
| CESM2 | ssp245 2071–2100 | −1.49×10⁻² | [−1.68×10⁻², −1.31×10⁻²] | 0.910 | Yes |
| GFDL-ESM4 | ssp245 2071–2100 | −5.34×10⁻³ | [−8.32×10⁻³, −2.36×10⁻³] | 0.325 | Yes |
| MPI-ESM1-2-LR | ssp245 2071–2100 | −2.29×10⁻² | [−2.89×10⁻², −1.68×10⁻²] | 0.684 | Yes |

Historical trends are weak. CESM2 and GFDL-ESM4 show a small but statistically significant freshening; MPI-ESM1-2-LR has no significant historical trend. Under ssp245, all three models show statistically significant freshening, with MPI-ESM1-2-LR strongest, then CESM2, then GFDL-ESM4. However, CESM2’s trend amplitude may be affected by the same masking/mean-state problem.

Ratio of ssp245 slope to historical slope:

| Model | Ratio | Meaningfulness |
|---|---:|---|
| CESM2 | 7.8 | Descriptive only: the two slopes are from disjoint periods and different forcings; not a continuation ratio. CESM2 also mixes ensemble members. |
| GFDL-ESM4 | 4.9 | Descriptive only for the same reason; both individual slopes exclude zero. |
| MPI-ESM1-2-LR | 33.3 | Not meaningful: the historical 95% CI includes zero, so the ratio is numerically unstable. |

### 3. End-of-century SSS change

ΔSSS = ssp245 2071–2100 mean minus historical 1995–2014 mean, in psu.

| Model | Global | Atlantic | Indian | Pacific |
|---|---:|---:|---:|---:|
| CESM2 | −0.12 | +0.09 | +0.06 | −0.24 |
| GFDL-ESM4 | −0.01 | +0.19 | −0.00 | −0.04 |
| MPI-ESM1-2-LR | −0.08 | −0.05 | −0.03 | −0.10 |

The global end-of-century changes are small. GFDL-ESM4 is nearly unchanged globally but has Atlantic salinification and slight Pacific freshening. MPI-ESM1-2-LR freshens in all three basins. CESM2 shows Pacific freshening and modest Atlantic/Indian salinification, but its absolute values are compromised by the apparent land-mask/area-weighting issue.

**Ensemble-member caveat:** CESM2 uses historical `r1i1p1f1` but ssp245 `r10i1p1f1`, so its ΔSSS mixes different ensemble members. GFDL-ESM4 and MPI-ESM1-2-LR use `r1i1p1f1` in both experiments, so their deltas do not have this member mismatch.

### Overall assessment

- **GFDL-ESM4** performs best in mean-state SSS: near-zero global bias, lowest RMSE, and highest spatial correlation with WOA23.
- **MPI-ESM1-2-LR** is also close in mean state but slightly fresher overall; its historical trend is not significant, while its ssp245 freshening is the strongest of the three.
- **CESM2** produces an implausible global mean SSS in this calculation, so its mean-state metrics and trend magnitudes should not be interpreted scientifically without correcting the masking/area-weighting problem. Its ssp245 delta additionally mixes ensemble members.
