# ============================================================
# Visualization
# ============================================================

# Figure 1: SSS bias maps (model 1971-2000 minus WOA23)
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
vmax1 = max(np.nanmax(np.abs(bias_maps[m])) for m in MODEL_NAMES)
if not np.isfinite(vmax1) or vmax1 == 0:
    vmax1 = 1.0
im1 = None
for ax, model in zip(axes1, MODEL_NAMES):
    im1 = ax.pcolormesh(
        woa_lon, woa_lat, bias_maps[model],
        cmap="RdBu_r", vmin=-vmax1, vmax=vmax1, shading="auto"
    )
    ax.set_title(model)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
fig1.colorbar(im1, ax=axes1, shrink=0.8, label="psu")

# Figure 2: cell-by-cell historical SSS trend, 1950-2014
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
vmax2 = max(np.nanmax(np.abs(trend_maps[m])) for m in MODEL_NAMES)
if not np.isfinite(vmax2) or vmax2 == 0:
    vmax2 = 0.1
im2 = None
for ax, model in zip(axes2, MODEL_NAMES):
    im2 = ax.pcolormesh(
        woa_lon, woa_lat, trend_maps[model],
        cmap="RdBu_r", vmin=-vmax2, vmax=vmax2, shading="auto"
    )
    sig = np.nan_to_num(p_maps[model], nan=1.0) < 0.05
    ax.contourf(
        woa_lon, woa_lat, sig.astype(float),
        levels=[0.5, 1.5], colors="none", hatches=["..."]
    )
    ax.set_title(model)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
fig2.colorbar(im2, ax=axes2, shrink=0.8, label="psu per decade")

# Figure 3: global mean SSS time series, historical and ssp245
fig3, ax3 = plt.subplots(figsize=(10, 5))
for model in MODEL_NAMES:
    color = MODEL_COLORS[model]
    hist = gm_hist_dict[model]
    ssp = gm_ssp_dict[model]

    ax3.plot(hist["year"], hist.values, color=color, linewidth=1.5,
             label=f"{model} historical")
    ax3.plot(ssp["year"], ssp.values, color=color, linewidth=1.5,
             label=f"{model} ssp245")

    for da, ls in [(hist, "--"), (ssp, ":")]:
        x = da["year"].values.astype(float)
        y = da.values.astype(float)
        coef = np.polyfit(x, y, 1)
        ax3.plot(x, np.polyval(coef, x), color=color, linestyle=ls,
                 linewidth=1.0)

ax3.set_xlabel("Year")
ax3.set_ylabel("SSS (psu)")
ax3.set_title("Global mean sea surface salinity")
ax3.grid(True, alpha=0.3)
ax3.legend(ncol=2, fontsize=8)
