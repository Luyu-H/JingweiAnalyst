# Visualization code only.
# Uses variables from the computation block.

colors = {
    "CESM2": "C0",
    "GFDL-ESM4": "C1",
    "MPI-ESM1-2-LR": "C2",
}

# -----------------------------------------------------------------------------
# Figure 1: SSS bias vs WOA23 on the common comparison cells
# -----------------------------------------------------------------------------
bias_abs_pool = np.concatenate(
    [bias_maps[m][np.isfinite(bias_maps[m])] for m in models]
)
vmax_bias = float(np.nanpercentile(np.abs(bias_abs_pool), 98))

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
for ax, model in zip(axes, models):
    im = ax.pcolormesh(
        woa_lon, woa_lat, bias_maps[model],
        shading="auto", cmap="RdBu_r", vmin=-vmax_bias, vmax=vmax_bias
    )
    ax.set_title(model, color=colors[model])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
fig.colorbar(
    im, ax=axes, orientation="horizontal", fraction=0.05, pad=0.05,
    extend="both", label="SSS bias (psu)"
)
fig.suptitle(
    f"SSS bias vs WOA23, 1971-2000, shared cells; shared diverging scale, "
    f"98th pct |bias|={vmax_bias:.3f} psu",
    fontsize=10
)

# -----------------------------------------------------------------------------
# Figure 2: 1950-2014 cell-by-cell SSS trend with p<0.05 stippling
# -----------------------------------------------------------------------------
trend_abs_pool = np.concatenate(
    [trend_maps_woa[m][np.isfinite(trend_maps_woa[m])] for m in models]
)
vmax_trend = float(np.nanpercentile(np.abs(trend_abs_pool), 98))
lon_grid, lat_grid = np.meshgrid(woa_lon, woa_lat)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), constrained_layout=True)
for ax, model in zip(axes, models):
    im = ax.pcolormesh(
        woa_lon, woa_lat, trend_maps_woa[model],
        shading="auto", cmap="RdBu_r", vmin=-vmax_trend, vmax=vmax_trend
    )
    p = p_maps_woa[model]
    sig = np.isfinite(p) & (p < 0.05) & np.isfinite(trend_maps_woa[model])
    ax.scatter(
        lon_grid[sig], lat_grid[sig], s=0.15, c="k", marker=".",
        linewidths=0, alpha=0.55
    )
    ax.set_title(model + "\nstippled: p<0.05", color=colors[model])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
fig.colorbar(
    im, ax=axes, orientation="horizontal", fraction=0.05, pad=0.05,
    extend="both", label="SSS trend (psu/decade)"
)
fig.suptitle(
    f"1950-2014 SSS trend, shared diverging scale; "
    f"98th pct |trend|={vmax_trend:.3g} psu/decade; stippled p<0.05",
    fontsize=10
)

# -----------------------------------------------------------------------------
# Figure 3: global mean SSS time series, historical and ssp245 segments
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9.2, 4.8), constrained_layout=True)

for model in models:
    color = colors[model]
    ts = time_series[model]

    ax.plot(
        ts["hist_years"], ts["hist_vals"],
        color=color, lw=1.8, label=model
    )
    ax.plot(
        ts["ssp_years"], ts["ssp_vals"],
        color=color, lw=1.8
    )

    fh = ts["fit_h"]
    fs = ts["fit_s"]
    yh = fh["intercept"] + fh["slope_per_year"] * ts["hist_years"]
    ys = fs["intercept"] + fs["slope_per_year"] * ts["ssp_years"]
    ax.plot(ts["hist_years"], yh, color=color, ls="--", lw=1.0)
    ax.plot(ts["ssp_years"], ys, color=color, ls="--", lw=1.0)

ax.set_xlim(1945, 2105)
ax.set_xlabel("Year")
ax.set_ylabel("Global mean SSS (psu)")
ax.set_title(
    "Global mean SSS: historical 1950-2014 and ssp245 2071-2100\n"
    "dashed lines are ESS-corrected linear fits"
)
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0)
