# -----------------------------------------------------------------------------
# Common plotting setup
# -----------------------------------------------------------------------------
dx = woa_lon[1] - woa_lon[0]
dy = woa_lat[1] - woa_lat[0]
lon_edges = np.r_[woa_lon - dx / 2.0, woa_lon[-1] + dx / 2.0]
lat_edges = np.r_[woa_lat - dy / 2.0, woa_lat[-1] + dy / 2.0]

model_colors = {
    "CESM2": "C0",
    "GFDL-ESM4": "C1",
    "MPI-ESM1-2-LR": "C2",
}

# -----------------------------------------------------------------------------
# Figure 1: SSS bias maps, 1971-2000 climatology minus WOA23
# -----------------------------------------------------------------------------
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
vmax1 = max(np.nanmax(np.abs(bias_maps[m])) for m in MODELS)
vmax1 = float(vmax1)

for ax, model in zip(axes1, MODELS):
    im1 = ax.pcolormesh(
        lon_edges, lat_edges, bias_maps[model],
        cmap="RdBu_r", vmin=-vmax1, vmax=vmax1, shading="auto",
    )
    ax.set_title(model)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect("equal")
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-90, -45, 0, 45, 90])

cbar1 = fig1.colorbar(
    im1, ax=axes1, orientation="horizontal", fraction=0.05, pad=0.05
)
cbar1.set_label("SSS bias (psu)")

# -----------------------------------------------------------------------------
# Figure 2: cell-by-cell historical SSS trend, 1950-2014
# -----------------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)
vmax2 = max(np.nanmax(np.abs(trend_maps[m])) for m in MODELS)
vmax2 = float(vmax2)

for ax, model in zip(axes2, MODELS):
    im2 = ax.pcolormesh(
        lon_edges, lat_edges, trend_maps[model],
        cmap="RdBu_r", vmin=-vmax2, vmax=vmax2, shading="auto",
    )

    p = trend_p_maps[model]
    sig = (p < 0.05) & np.isfinite(p)
    yy, xx = np.where(sig)
    if len(yy) > 0:
        step = max(1, len(yy) // 6000)
        ax.scatter(
            woa_lon[xx[::step]], woa_lat[yy[::step]],
            s=0.5, c="k", marker=".", linewidths=0,
        )

    ax.set_title(model)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect("equal")
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_yticks([-90, -45, 0, 45, 90])

cbar2 = fig2.colorbar(
    im2, ax=axes2, orientation="horizontal", fraction=0.05, pad=0.05
)
cbar2.set_label(r"SSS trend (psu decade$^{-1}$)")

# -----------------------------------------------------------------------------
# Figure 3: global mean SSS time series, historical and ssp245
# -----------------------------------------------------------------------------
fig3, ax3 = plt.subplots(figsize=(10, 5), constrained_layout=True)

for model in MODELS:
    c = model_colors[model]
    d = global_series[model]

    ax3.plot(
        d["hist_years"], d["hist_sss"],
        color=c, lw=1.5, label=f"{model} historical",
    )
    ax3.plot(
        d["ssp_years"], d["ssp_sss"],
        color=c, lw=1.5, ls="--", label=f"{model} ssp245",
    )
    ax3.plot(
        d["hist_years"], d["hist_fit"],
        color=c, lw=1.0, ls="-",
    )
    ax3.plot(
        d["ssp_years"], d["ssp_fit"],
        color=c, lw=1.0, ls="--",
    )

ax3.set_xlim(1945, 2105)
ax3.set_xlabel("Year")
ax3.set_ylabel("Global mean SSS (psu)")
ax3.grid(True, alpha=0.3)
ax3.legend(fontsize=8, ncol=2)
