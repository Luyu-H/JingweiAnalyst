The internal-consistency test was performed by computing an offline MLD from each model’s own monthly `thetao` and `so` with a sigma-theta threshold of 0.03 kg m⁻³ relative to 10 m, interpolating the first crossing between model levels, and comparing it with the same model’s published `mlotst`. All statistics below use only cells where both fields are defined for all 30 Marches or all 30 Septembers. Mean difference, RMSD and correlation are cell-wise over that common mask; the percentage-area statistics are area-weighted with each model’s `areacello`.

### Numerical results: operator MLD minus published `mlotst`

| Model | Month | Mean difference (m) | RMSD (m) | Pearson correlation | Area with \|difference\| > 20 m (%) | Common-valid area (% of model ocean) |
|---|---:|---:|---:|---:|---:|---:|
| CESM2 | March | 1.6 | 35.7 | 0.959 | 3.8 | 100.0 |
| CESM2 | September | 0.3 | 7.8 | 0.992 | 1.4 | 100.0 |
| GFDL-ESM4 | March | -5.5 | 33.7 | 0.969 | 6.9 | 100.0 |
| GFDL-ESM4 | September | -4.1 | 19.4 | 0.973 | 7.4 | 100.0 |
| MPI-ESM1-2-LR | March | 2.8 | 28.9 | 0.993 | 1.7 | 100.0 |
| MPI-ESM1-2-LR | September | -5.0 | 36.5 | 0.985 | 4.0 | 100.0 |

Positive mean differences mean the offline operator MLD is deeper than the published `mlotst`.

### Ranking and criterion

Using the mean of the March and September RMSD as the ranking criterion, lower RMSD being better, the models rank from most to least internally consistent as:

1. **CESM2** — mean RMSD = 21.8 m  
2. **GFDL-ESM4** — mean RMSD = 26.5 m  
3. **MPI-ESM1-2-LR** — mean RMSD = 32.7 m  

This ranking is metric-dependent. By mean correlation, MPI-ESM1-2-LR is highest (mean r ≈ 0.989, versus CESM2 ≈ 0.975 and GFDL-ESM4 ≈ 0.971). By the percentage of area with |difference| > 20 m, GFDL-ESM4 is clearly the worst, while CESM2 and MPI-ESM1-2-LR are closer. The robust conclusion is that CESM2 performs best or near-best by RMSD and has a very small September error; MPI-ESM1-2-LR has the highest spatial correlation and a small March >20 m area but a large September RMSD in the Southern Ocean; GFDL-ESM4 has the most widespread area exceeding 20 m.

### Where the disagreement concentrates

The coarse 10° area-weighted diagnostics show that disagreement is not spatially uniform:

- **CESM2:** largest March errors are localized in high northern latitudes (especially 60–70°N near 150–160°E and North Atlantic/Nordic boxes), the Arabian/Red Sea region around 20–30°N, 40–50°E, and a Southern Indian sector near 30–20°S, 140–150°E. September is much better, with residual boxes near the Gulf of Mexico/western Atlantic, the western Pacific, and high-latitude North Atlantic.
- **GFDL-ESM4:** the most widespread disagreement, concentrated in western boundary current regions (Gulf Stream and Kuroshio), North Atlantic/North Pacific midlatitudes, Arctic/Nordic seas in March, and the Southern Ocean/Antarctic plus South Atlantic/Indian sectors in September.
- **MPI-ESM1-2-LR:** March disagreement is concentrated in the high-latitude North Atlantic/Nordic/Arctic region and western North Pacific; September disagreement is strongly concentrated in the Southern Ocean/Antarctic, especially 60–80°S.

### Why the operator and `mlotst` are not expected to agree exactly

Even though both quantities come from the same simulation, they are not expected to be identical:

1. **Online versus offline calculation:** `mlotst` is the model’s online mixed-layer diagnostic, computed from the model state and its own algorithm at model time steps. The offline operator uses monthly-mean `thetao` and `so`, so it loses sub-monthly covariation and nonlinear averaging effects.
2. **Different density and threshold implementation:** `mlotst` is defined with a density threshold relative to 10 m, but models may use slightly different density formulas, reference treatment, and vertical interpolation. The offline operator here uses EOS-80 from potential temperature and practical salinity.
3. **Vertical discretization:** The offline MLD interpolates the first threshold crossing between model levels. The published `mlotst` may use a different interpolation, a different reference level, or the model’s native vertical-grid treatment.
4. **Masking, sea ice, and marginal seas:** Mixed-layer diagnostics are sensitive near coasts, marginal seas, high latitudes, and sea-ice edges, where monthly means, missing values, and ice/ocean masks can affect the computed crossing.
5. **Highly nonlinear threshold behavior:** Small differences in near-surface stratification can move the first density-threshold crossing by many metres, especially in regions with weak or shallow stratification and strong vertical gradients, such as western boundary currents, the Southern Ocean, and high-latitude convection regions.

### Figures

The requested figures were generated:

- **Figure 1:** 3 × 2 global maps of operator MLD minus published `mlotst`, with rows CESM2, GFDL-ESM4 and MPI-ESM1-2-LR, columns March and September, all panels on a common longitude–latitude display grid and using one shared diverging colour scale centred on zero (approximately ±52.7 m from the 98th percentile of the displayed differences).
- **Figure 2:** 3 × 2 log–log scatter plots of operator MLD versus published `mlotst`, with common axes (0.1–10,000 m), a dashed y = x line, and each panel annotated with its mean difference, RMSD and correlation. The points lie close to the one-to-one line overall, consistent with the high correlations, with the largest departures in regions and months with deep MLDs and strong localized biases.
