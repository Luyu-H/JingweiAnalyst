**Short answer.** The offline density-threshold MLD reproduces the published `mlotst` field **broadly but not exactly**. Spatial correlations are high (0.966–0.993), but there are systematic negative mean differences and model/month-dependent RMSEs, with the largest disagreements concentrated in high-latitude and Southern Ocean regions. This is consistent with an internal-consistency test of the diagnostic: both fields come from the same simulation, so the residual disagreement diagnoses algorithm, density-formulation, vertical-discretization and time-averaging differences rather than model physics.

**Area-weighted metrics for March and September, 1971–2000.**  
Mean difference is operator MLD minus published `mlotst`; negative means the operator MLD is, on average, shallower than `mlotst`.

| Model | Month | Mean diff (m) | RMSE (m) | Spatial corr. | Area with \|diff\|>20 m (%) | Area both defined (% ocean area) |
|---|---:|---:|---:|---:|---:|---:|
| CESM2 | March | −2.4 | 25.2 | 0.966 | 5.3 | 100.0 |
| CESM2 | September | −2.3 | 8.0 | 0.993 | 1.6 | 100.0 |
| GFDL-ESM4 | March | −8.4 | 28.3 | 0.975 | 7.3 | 100.0 |
| GFDL-ESM4 | September | −7.5 | 16.9 | 0.981 | 8.1 | 100.0 |
| MPI-ESM1-2-LR | March | −2.1 | 12.6 | 0.992 | 2.2 | 100.0 |
| MPI-ESM1-2-LR | September | −4.0 | 28.4 | 0.978 | 4.0 | 100.0 |

The common valid mask covers 100.0% of each model’s ocean area in these diagnostics, because requiring both fields to be defined in all 30 Marches or 30 Septembers did not remove any additional ocean cells beyond the model ocean mask.

**Ranking of internal consistency.**  
Using the stated criterion — the mean area-weighted RMSE of operator MLD minus `mlotst` over March and September, with lower RMSE meaning more internally consistent — the ranking is:

1. **CESM2** — mean RMSE ≈ 16.6 m  
2. **MPI-ESM1-2-LR** — mean RMSE ≈ 20.5 m  
3. **GFDL-ESM4** — mean RMSE ≈ 22.6 m  

This ranking is metric-dependent, but under the chosen RMSE criterion CESM2 is the most internally consistent overall. CESM2 has a much smaller September RMSE than the other two; GFDL-ESM4 has the largest March bias and RMSE, while MPI-ESM1-2-LR has a very large September RMSE concentrated in the Southern Ocean.

**Where the disagreements concentrate.**  
The latitude-band diagnostics show that the differences are not globally uniform:
- **CESM2:** March disagreement is concentrated in the Northern Hemisphere high and mid-latitudes (60–90°N mean absolute difference ≈ 26.0 m; 30–60°N ≈ 14.2 m). September disagreement is weaker and shifts to the Southern subtropics/mid-latitudes (30°S–0 ≈ 5.8 m; 60–30°S ≈ 5.5 m).
- **GFDL-ESM4:** March disagreement is largest in northern mid/high latitudes (30–60°N ≈ 28.2 m; 60–90°N ≈ 27.9 m) and also in the southern high latitudes. September disagreement is strongest in the Southern Ocean and southern mid/high latitudes (60–30°S ≈ 16.5 m; 90–60°S ≈ 14.8 m).
- **MPI-ESM1-2-LR:** March disagreement is moderate and concentrated in northern high/mid-latitudes. September is dominated by a strong Southern Ocean signal (90–60°S mean absolute difference ≈ 48.9 m), with secondary disagreement at 60–30°S.

**Why exact agreement is not expected.**  
The published `mlotst` is a model-diagnosed field computed online with each model’s own mixed-layer algorithm, reference-depth convention, potential-density formulation, vertical discretization and temporal sampling. The offline operator instead uses monthly-mean `thetao`/`so`, EOS-80 density at zero pressure, a surface reference density, the first crossing of the 0.03 kg m⁻³ threshold, and linear interpolation between level centres. It also cannot capture sub-monthly mixing events or nonlinear averaging of the mixed-layer diagnostic. Therefore nonzero differences are expected even for the same simulation; they diagnose the sensitivity of the MLD definition and implementation, not a failure of the model physics.

The requested figures were generated: Figure 1 shows common-colour-scale maps of operator MLD minus `mlotst` for the three models and two months, and Figure 2 shows the corresponding log–log scatter plots against the one-to-one line, with bias, RMSE and correlation annotated in each panel.
