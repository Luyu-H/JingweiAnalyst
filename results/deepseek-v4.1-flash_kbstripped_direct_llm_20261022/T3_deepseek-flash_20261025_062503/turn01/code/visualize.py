# ----------------------------------------------------------------------
# Figure 1: 2x2 panel of model-minus-WOA bias fields and WOA reference
# ----------------------------------------------------------------------
bias_fields = {
    m: np.where(common_mask, model_on_woa[m] - woa_sst_vals, np.nan)
    for m in models
}
woa_ref_common = np.where(common_mask, woa_sst_vals, np.nan)

bias_abs_max = max(np.nanmax(np.abs(bias_fields[m])) for m in models)

fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)

panel_specs = [
    (axes[0, 0], "CESM2 SST bias (model - WOA23)", bias_fields["CESM2"], "RdBu_r", -bias_abs_max, bias_abs_max, "SST bias (°C)"),
    (axes[0, 1], "GFDL-ESM4 SST bias (model - WOA23)", bias_fields["GFDL-ESM4"], "RdBu_r", -bias_abs_max, bias_abs_max, "SST bias (°C)"),
    (axes[1, 0], "MPI-ESM1-2-LR SST bias (model - WOA23)", bias_fields["MPI-ESM1-2-LR"], "RdBu_r", -bias_abs_max, bias_abs_max, "SST bias (°C)"),
    (axes[1, 1], "WOA23 0 m annual SST reference", woa_ref_common, "viridis", np.nanmin(woa_ref_common), np.nanmax(woa_ref_common), "SST (°C)"),
]

for ax, title, data, cmap, vmin, vmax, cbar_label in panel_specs:
    im = ax.pcolormesh(
        woa_lon,
        woa_lat,
        data,
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(woa_lon.min(), woa_lon.max())
    ax.set_ylim(woa_lat.min(), woa_lat.max())
    cbar = fig.colorbar(im, ax=ax, shrink=0.82)
    cbar.set_label(cbar_label)

# ----------------------------------------------------------------------
# Figure 2: global mean SST over each product's native ocean grid
# ----------------------------------------------------------------------
bar_labels = models + ["WOA23"]
bar_values = [model_data[m]["global_mean"] for m in models] + [woa_global_mean]

fig2, ax2 = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
bars = ax2.bar(bar_labels, bar_values, color=["C0", "C1", "C2", "C3"])

for bar, val in zip(bars, bar_values):
    ax2.text(
        bar.get_x() + bar.get_width() / 2.0,
        bar.get_height(),
        f"{val:.2f}",
        ha="center",
        va="bottom",
    )

ax2.set_ylabel("SST (°C)")
ax2.set_title("Global mean SST over native ocean grids, 1971-2000")
ax2.set_ylim(0.0, max(bar_values) * 1.12)
