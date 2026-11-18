# -----------------------------------------------------------------------------
# Shared plotting setup
# -----------------------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
model_colors = {"CESM2": "C0", "GFDL-ESM4": "C1", "MPI-ESM1-2-LR": "C2"}

lon = lon_woa
lat = lat_woa
lon2d, lat2d = np.meshgrid(lon, lat)

# -----------------------------------------------------------------------------
# Figure 1: 1971-2000 SSS bias, model minus WOA23
# -----------------------------------------------------------------------------
all_bias = np.concatenate([bias_maps[m].values.ravel() for m in models])
vmax_bias = float(np.nanpercentile(np.abs(all_bias), 98))

fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
for ax, model in zip(axes1, models):
    data = bias_maps[model].values
    im = ax.pcolormesh(
        lon2d, lat2d, data, cmap="RdBu_r",
        vmin=-vmax_bias, vmax=vmax_bias, shading="auto",
    )
    ax.set_title(f"{model}  1971-2000 SSS bias")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cbar1 = fig1.colorbar(
    im, ax=axes1, orientation="horizontal",
    fraction=0.05, pad=0.05, extend="both",
)
cbar1.set_label("SSS bias (psu)")
fig1.suptitle(
    f"Model minus WOA23 SSS on common mask (cells={common_mask_cells}); "
    "shared scale clipped at 98th percentile",
    fontsize=12,
)

# -----------------------------------------------------------------------------
# Figure 2: 1950-2014 annual SSS trend, p < 0.05 stippled
# -----------------------------------------------------------------------------
all_trend = np.concatenate([trend_maps[m].values.ravel() for m in models])
vmax_trend = float(np.nanpercentile(np.abs(all_trend), 98))

fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
for ax, model in zip(axes2, models):
    data = trend_maps[model].values
    im = ax.pcolormesh(
        lon2d, lat2d, data, cmap="RdBu_r",
        vmin=-vmax_trend, vmax=vmax_trend, shading="auto",
    )
    sig = (pval_maps[model].values < 0.05).astype(float)
    ax.contourf(
        lon2d, lat2d, sig,
        levels=[0.5, 1.5], colors="none", hatches=["...."],
    )
    ax.set_title(f"{model}  1950-2014 SSS trend")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cbar2 = fig2.colorbar(
    im, ax=axes2, orientation="horizontal",
    fraction=0.05, pad=0.05, extend="both",
)
cbar2.set_label("SSS trend (psu/decade)")
fig2.suptitle(
    "Annual SSS trend, 1950-2014; p<0.05 stippled; "
    "shared scale clipped at 98th percentile",
    fontsize=12,
)

# -----------------------------------------------------------------------------
# Figure 3: global mean SSS time series, separate segments and fitted lines
# -----------------------------------------------------------------------------
fig3, ax3 = plt.subplots(figsize=(10, 5.2), constrained_layout=True)

for model in models:
    color = model_colors[model]

    gm_h = global_mean_hist[model]
    years_h = gm_h["year"].values
    vals_h = gm_h.values
    ax3.plot(years_h, vals_h, color=color, lw=1.5, label=f"{model} historical")

    hstats = trend_stats_hist[model]
    fit_h = hstats["slope"] * years_h + hstats["intercept"]
    ax3.plot(years_h, fit_h, color=color, ls="--", lw=1.2)

    gm_s = global_mean_ssp[model]
    years_s = gm_s["year"].values
    vals_s = gm_s.values
    ax3.plot(years_s, vals_s, color=color, lw=1.5, label=f"{model} ssp245")

    sstats = trend_stats_ssp[model]
    fit_s = sstats["slope"] * years_s + sstats["intercept"]
    ax3.plot(years_s, fit_s, color=color, ls="--", lw=1.2)

ax3.set_xlabel("Year")
ax3.set_ylabel("Global mean SSS (psu)")
ax3.set_xlim(1945, 2105)
ax3.set_title(
    "Annual global mean SSS: historical 1950-2014 and ssp245 2071-2100 "
    "(gap not bridged; dashed = fitted trend)"
)
ax3.legend(
    loc="center left", bbox_to_anchor=(1.01, 0.5),
    borderaxespad=0.0, fontsize=8,
)
