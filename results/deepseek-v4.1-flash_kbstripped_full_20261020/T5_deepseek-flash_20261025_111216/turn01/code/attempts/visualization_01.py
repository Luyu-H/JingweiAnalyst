# ============================================================
# Figure 1: SSS bias (model 1971-2000 minus WOA23), shared scale
# ============================================================
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

all_bias_values = np.concatenate(
    [bias_maps[m].values.ravel() for m in models]
)
vmax_bias = float(np.nanpercentile(np.abs(all_bias_values), 98))

lon = woa_sss.lon.values
lat = woa_sss.lat.values

pcm1 = None
for ax, m in zip(axes1, models):
    pcm1 = ax.pcolormesh(
        lon,
        lat,
        bias_maps[m].values,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax_bias,
        vmax=vmax_bias,
    )
    ax.set_title(m)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig1.colorbar(
    pcm1,
    ax=axes1,
    label="SSS bias (psu)",
    extend="both",
)
fig1.suptitle(
    "SSS bias: model 1971-2000 minus WOA23 "
    f"(shared scale clipped at 98th percentile, {vmax_bias:.2f} psu)",
    fontsize=12,
)

# ============================================================
# Figure 2: cell-by-cell historical SSS trend, 1950-2014
# ============================================================
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

all_slope_values = np.concatenate(
    [model_results[m]["slope_regrid"].values.ravel() for m in models]
)
vmax_slope = float(np.nanpercentile(np.abs(all_slope_values), 98))

pcm2 = None
lon2d, lat2d = np.meshgrid(lon, lat)

for ax, m in zip(axes2, models):
    slope = model_results[m]["slope_regrid"].values
    sig = model_results[m]["sig_regrid"].values

    pcm2 = ax.pcolormesh(
        lon,
        lat,
        slope,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax_slope,
        vmax=vmax_slope,
    )
    ax.scatter(
        lon2d[sig],
        lat2d[sig],
        s=1,
        marker=".",
        color="k",
        alpha=0.6,
    )
    ax.set_title(m)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig2.colorbar(
    pcm2,
    ax=axes2,
    label="SSS trend (psu/decade)",
    extend="both",
)
fig2.suptitle(
    "Cell-by-cell SSS trend 1950-2014 "
    f"(shared scale clipped at 98th percentile, {vmax_slope:.3f} psu/decade; "
    "black dots: corrected p < 0.05)",
    fontsize=12,
)

# ============================================================
# Figure 3: global mean SSS time series, historical and ssp245
# ============================================================
fig3, ax3 = plt.subplots(figsize=(10, 5), constrained_layout=True)

model_colors = {
    "CESM2": "C0",
    "GFDL-ESM4": "C1",
    "MPI-ESM1-2-LR": "C2",
}

for m in models:
    c = model_colors[m]
    hist = model_results[m]["gm_hist_series"]
    ssp = model_results[m]["gm_ssp_series"]

    x_hist = hist.coords[hist.dims[0]].values
    x_ssp = ssp.coords[ssp.dims[0]].values

    ax3.plot(
        x_hist,
        hist.values,
        color=c,
        label=f"{m} historical",
    )
    ax3.plot(
        x_ssp,
        ssp.values,
        color=c,
        linestyle="--",
        label=f"{m} ssp245",
    )

    th_fit = model_results[m]["trend_hist"]["fitted_data"]
    ts_fit = model_results[m]["trend_ssp"]["fitted_data"]

    x_th = th_fit.coords[th_fit.dims[0]].values
    x_ts = ts_fit.coords[ts_fit.dims[0]].values

    ax3.plot(
        x_th,
        th_fit.values,
        color=c,
        linestyle=":",
        linewidth=1.0,
    )
    ax3.plot(
        x_ts,
        ts_fit.values,
        color=c,
        linestyle=":",
        linewidth=1.0,
    )

ax3.set_xlabel("Year")
ax3.set_ylabel("Global mean SSS (psu)")
ax3.set_title(
    "Global mean SSS: historical 1950-2014 and ssp245 2071-2100 "
    "(dotted: fitted linear trends; 2015-2070 gap left empty)"
)
ax3.legend(loc="best", fontsize=8)
ax3.grid(True, alpha=0.3)
