# Visualization code only

model_colors = {
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}

# -----------------------------
# Figure 1: SSS bias maps
# -----------------------------
bias_arrays = [bias_maps[m].values.ravel() for m in models]
bias_pool = np.concatenate([a[np.isfinite(a)] for a in bias_arrays])
bias_vmax = float(np.nanpercentile(np.abs(bias_pool), 98))

fig1, axes1 = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
for ax, model in zip(axes1, models):
    da = bias_maps[model]
    im = ax.pcolormesh(
        da["lon"], da["lat"], da,
        cmap="RdBu_r", vmin=-bias_vmax, vmax=bias_vmax, shading="auto"
    )
    ax.set_title(f"{model}")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
cbar1 = fig1.colorbar(im, ax=axes1, orientation="vertical", shrink=0.95, extend="both")
cbar1.set_label("SSS bias (psu)")
fig1.suptitle("Figure 1. 1971–2000 SSS bias vs WOA23 (shared 98th-percentile scale)", y=1.02)

# -----------------------------
# Figure 2: Historical trend maps with stippling
# -----------------------------
trend_arrays = [trend_maps_rg[m].values.ravel() for m in models]
trend_pool = np.concatenate([a[np.isfinite(a)] for a in trend_arrays])
trend_vmax = float(np.nanpercentile(np.abs(trend_pool), 98))

fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
for ax, model in zip(axes2, models):
    slope = trend_maps_rg[model]
    sig = signif_maps_rg[model]

    im2 = ax.pcolormesh(
        slope["lon"], slope["lat"], slope,
        cmap="RdBu_r", vmin=-trend_vmax, vmax=trend_vmax, shading="auto"
    )

    lon2d, lat2d = np.meshgrid(slope["lon"].values, slope["lat"].values)
    sig_mask = sig.values
    ax.scatter(
        lon2d[sig_mask], lat2d[sig_mask],
        s=2, c="k", marker=".", linewidths=0, alpha=0.5
    )

    ax.set_title(f"{model}")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

cbar2 = fig2.colorbar(im2, ax=axes2, orientation="vertical", shrink=0.95, extend="both")
cbar2.set_label("SSS trend (psu/decade)")
fig2.suptitle("Figure 2. Historical 1950–2014 annual-mean SSS trend (shared 98th-percentile scale)", y=1.02)

# -----------------------------
# Figure 3: Global-mean SSS time series and fitted segments
# -----------------------------
fig3, ax3 = plt.subplots(figsize=(11, 5), constrained_layout=True)

for model in models:
    c = model_colors[model]

    hist = global_annual_hist[model]
    ssp = global_annual_ssp[model]

    hist_x = hist["year"].values if "year" in hist.coords else hist.coords[hist.dims[0]].values
    ssp_x = ssp["year"].values if "year" in ssp.coords else ssp.coords[ssp.dims[0]].values

    ax3.plot(hist_x, hist.values, color=c, lw=1.8, label=model)
    ax3.plot(ssp_x, ssp.values, color=c, lw=1.8)

    hist_fit = trend_hist[model]["fitted_data"]
    ssp_fit = trend_ssp[model]["fitted_data"]

    hist_fit_x = hist_fit["year"].values if "year" in hist_fit.coords else hist_fit.coords[hist_fit.dims[0]].values
    ssp_fit_x = ssp_fit["year"].values if "year" in ssp_fit.coords else ssp_fit.coords[ssp_fit.dims[0]].values

    ax3.plot(hist_fit_x, hist_fit.values, color=c, lw=2.5, ls="--")
    ax3.plot(ssp_fit_x, ssp_fit.values, color=c, lw=2.5, ls="--")

ax3.set_xlabel("Year")
ax3.set_ylabel("SSS (psu)")
ax3.set_title("Figure 3. Global mean SSS: historical and ssp245 with separate fitted segments")
ax3.legend(loc="best")
plt.show()
