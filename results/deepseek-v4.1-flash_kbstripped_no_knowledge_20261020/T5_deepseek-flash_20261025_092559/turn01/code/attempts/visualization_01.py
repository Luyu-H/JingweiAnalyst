import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Figure 1: SSS bias maps
fig, axes = plt.subplots(1, 3, figsize=(16, 4), constrained_layout=True)
vals = []
for model in models:
    arr = bias_fields[model].compute().values.ravel()
    vals.append(arr[np.isfinite(arr)])
allv = np.concatenate(vals)
vmax1 = float(np.nanpercentile(np.abs(allv), 98))
if not np.isfinite(vmax1) or vmax1 == 0:
    vmax1 = 0.1

for ax, model in zip(axes, models):
    b = bias_fields[model]
    lon = b.lon.values
    lat = b.lat.values
    im = ax.pcolormesh(
        lon, lat, b.values, shading="auto",
        cmap="RdBu_r", vmin=-vmax1, vmax=vmax1
    )
    ax.set_title(model, color=model_colors[model], fontweight="bold")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cb = fig.colorbar(
    im, ax=axes, orientation="horizontal",
    fraction=0.05, pad=0.05, extend="both"
)
cb.set_label("SSS bias (psu)")
fig.suptitle(
    "Model minus WOA23 1971-2000 SSS bias (shared scale, clipped at 98th percentile)",
    fontsize=12
)

# Figure 2: cell-by-cell historical SSS trend maps
fig, axes = plt.subplots(1, 3, figsize=(16, 4), constrained_layout=True)
vals2 = []
for model in models:
    arr = trend_slope_regrid[model].compute().values.ravel()
    vals2.append(arr[np.isfinite(arr)])
allv2 = np.concatenate(vals2)
vmax2 = float(np.nanpercentile(np.abs(allv2), 98))
if not np.isfinite(vmax2) or vmax2 == 0:
    vmax2 = 0.1

for ax, model in zip(axes, models):
    sl = trend_slope_regrid[model]
    p = trend_p_regrid[model]
    lon = sl.lon.values
    lat = sl.lat.values

    im = ax.pcolormesh(
        lon, lat, sl.values, shading="auto",
        cmap="RdBu_r", vmin=-vmax2, vmax=vmax2
    )

    sig = xr.where(np.isfinite(p), (p < 0.05), 0.0).values
    ax.contourf(
        lon, lat, sig, levels=[0.5, 1.5],
        colors="none", hatches=["...."]
    )

    ax.set_title(model, color=model_colors[model], fontweight="bold")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cb = fig.colorbar(
    im, ax=axes, orientation="horizontal",
    fraction=0.05, pad=0.05, extend="both"
)
cb.set_label("SSS trend (psu decade$^{-1}$)")
fig.suptitle(
    "Historical 1950-2014 annual-mean SSS trend; stippling p < 0.05 "
    "(shared scale, clipped at 98th percentile)",
    fontsize=12
)

# Figure 3: global-mean SSS time series, historical and ssp245 segments
fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)

for model in models:
    color = model_colors[model]

    h = annual_hist_global[model]
    s = annual_ssp_global[model]

    ax.plot(
        h["year"], h, color=color, lw=1.2, alpha=0.6,
        label=f"{model} historical"
    )
    ax.plot(
        s["year"], s, color=color, lw=1.2, alpha=0.6, ls="--",
        label=f"{model} ssp245"
    )

    hf = hist_fit[model]["fitted_data"]
    sf = ssp_fit[model]["fitted_data"]
    ax.plot(hf["year"], hf, color=color, lw=2.2)
    ax.plot(sf["year"], sf, color=color, lw=2.2, ls="--")

ax.set_xlabel("Year")
ax.set_ylabel("Global mean SSS (psu)")
ax.set_xlim(1950, 2100)
ax.grid(True, alpha=0.3)
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=9)
ax.set_title("Global mean SSS: historical 1950-2014 and ssp245 2071-2100")
