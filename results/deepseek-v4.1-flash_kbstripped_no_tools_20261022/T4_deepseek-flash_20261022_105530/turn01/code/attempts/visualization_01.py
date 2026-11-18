import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Figure 1: maps of D = operator MLD - published mlotst
# -----------------------------------------------------------------------------
model_order = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

# Shared symmetric colour limits from pooled remapped fields.
pooled = []
for model in model_order:
    for m in (3, 9):
        arr = results[model][m]["D_display"]
        pooled.append(arr[np.isfinite(arr)].ravel())
pooled = np.concatenate(pooled) if pooled else np.array([0.0])
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if not np.isfinite(vmax) or vmax <= 0:
    vmax = 1.0

fig1, axes1 = plt.subplots(
    3, 2, figsize=(14, 11), constrained_layout=True
)

for i, model in enumerate(model_order):
    lon_c = results[model]["display"]["lon"]
    lat_c = results[model]["display"]["lat"]
    ocean_frac = results[model]["display"]["ocean_frac"]

    for j, m in enumerate((3, 9)):
        ax = axes1[i, j]
        D_disp = results[model][m]["D_display"]

        im = ax.pcolormesh(
            lon_c, lat_c, D_disp,
            cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto"
        )

        # Coastline-like contour from the remapped ocean fraction.
        ax.contour(
            lon_c, lat_c, ocean_frac,
            levels=[0.5], colors="k", linewidths=0.4
        )

        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_xticks(np.arange(-180, 181, 60))
        ax.set_yticks(np.arange(-90, 91, 30))
        ax.set_title(f"{model} — {month_names[m]}", fontsize=10)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

# One shared colorbar.
cbar = fig1.colorbar(
    im, ax=axes1.ravel().tolist(),
    orientation="horizontal", fraction=0.03, pad=0.04,
    extend="both"
)
cbar.set_label("Operator MLD − mlotst (m)")
fig1.suptitle(
    "Operator MLD minus published mlotst, 1971–2000 monthly climatologies\n"
    f"Display: conservative area-weighted remap to 2° grid; colour clipped at pooled 98th percentile |D| = {vmax:.1f} m",
    fontsize=11
)

# -----------------------------------------------------------------------------
# Figure 2: scatter operator MLD vs published mlotst on common-mask cells
# -----------------------------------------------------------------------------
# Pool all positive common-mask values to set identical log limits.
all_positive = []
for model in model_order:
    for m in (3, 9):
        x = results[model][m]["op_m"].values.ravel()
        y = results[model][m]["ml_m"].values.ravel()
        good = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        if good.any():
            all_positive.append(x[good])
            all_positive.append(y[good])

if all_positive:
    all_positive = np.concatenate(all_positive)
    lim_lo = max(0.1, float(np.nanpercentile(all_positive, 1.0)) * 0.8)
    lim_hi = float(np.nanpercentile(all_positive, 99.5)) * 1.2
    if not np.isfinite(lim_hi) or lim_hi <= lim_lo:
        lim_hi = 5000.0
else:
    lim_lo, lim_hi = 0.1, 5000.0

tick_candidates = np.array(
    [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
)
ticks = tick_candidates[(tick_candidates >= lim_lo) & (tick_candidates <= lim_hi)]
if len(ticks) < 2:
    ticks = np.array([lim_lo, lim_hi])

fig2, axes2 = plt.subplots(
    3, 2, figsize=(12, 12), constrained_layout=True
)

for i, model in enumerate(model_order):
    for j, m in enumerate((3, 9)):
        ax = axes2[i, j]
        r = results[model][m]["metrics"]

        x = results[model][m]["op_m"].values.ravel()
        y = results[model][m]["ml_m"].values.ravel()
        good = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        omitted = int(len(x) - good.sum())

        ax.scatter(
            x[good], y[good],
            s=1.0, alpha=0.25, edgecolors="none", rasterized=True
        )
        ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "k--", lw=0.9)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lim_lo, lim_hi)
        ax.set_ylim(lim_lo, lim_hi)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks], fontsize=7)
        ax.set_yticklabels([f"{t:g}" for t in ticks], fontsize=7)
        ax.set_xlabel("Operator MLD (m)")
        ax.set_ylabel("Published mlotst (m)")
        ax.set_title(
            f"{model} — {month_names[m]}\n"
            f"bias={r['mean']:.1f} m, RMSE={r['rmse']:.1f} m, r={r['corr']:.3f}\n"
            f"non-positive cells omitted: {omitted}",
            fontsize=9
        )

fig2.suptitle(
    "Operator MLD vs published mlotst on common-mask ocean cells\n"
    "Logarithmic axes with identical limits and ticks; dashed line is y = x",
    fontsize=11
)
