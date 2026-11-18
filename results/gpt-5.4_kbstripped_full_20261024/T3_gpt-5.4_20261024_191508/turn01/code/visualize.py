# Visualization code only

# --- Figure 1: 2x2 panel of bias maps and WOA23 SST map ---
fig1, axes = plt.subplots(2, 2, figsize=(15, 8), constrained_layout=True)

# Robust shared color scale for bias panels
all_bias_vals = np.concatenate([
    bias_fields["CESM2"].values.ravel(),
    bias_fields["GFDL-ESM4"].values.ravel(),
    bias_fields["MPI-ESM1-2-LR"].values.ravel(),
])
all_bias_vals = all_bias_vals[np.isfinite(all_bias_vals)]
bias_vmax = float(np.nanpercentile(np.abs(all_bias_vals), 98))

# Robust scale for WOA23 SST
woa_vals = woa_comp.values.ravel()
woa_vals = woa_vals[np.isfinite(woa_vals)]
woa_vmin = float(np.nanpercentile(woa_vals, 2))
woa_vmax = float(np.nanpercentile(woa_vals, 98))

panel_order = [
    ("CESM2", axes[0, 0]),
    ("GFDL-ESM4", axes[0, 1]),
    ("MPI-ESM1-2-LR", axes[1, 0]),
]

bias_mappables = []
for model, ax in panel_order:
    im = ax.pcolormesh(
        bias_fields[model]["lon"],
        bias_fields[model]["lat"],
        bias_fields[model],
        cmap="RdBu_r",
        shading="auto",
        vmin=-bias_vmax,
        vmax=bias_vmax,
    )
    bias_mappables.append(im)
    ax.set_title(f"{model} - WOA23 SST bias")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(float(bias_fields[model]["lon"].min()), float(bias_fields[model]["lon"].max()))
    ax.set_ylim(float(bias_fields[model]["lat"].min()), float(bias_fields[model]["lat"].max()))

im_woa = axes[1, 1].pcolormesh(
    woa_comp["lon"],
    woa_comp["lat"],
    woa_comp,
    cmap="viridis",
    shading="auto",
    vmin=woa_vmin,
    vmax=woa_vmax,
)
axes[1, 1].set_title("WOA23 0 m annual SST")
axes[1, 1].set_xlabel("Longitude")
axes[1, 1].set_ylabel("Latitude")
axes[1, 1].set_xlim(float(woa_comp["lon"].min()), float(woa_comp["lon"].max()))
axes[1, 1].set_ylim(float(woa_comp["lat"].min()), float(woa_comp["lat"].max()))

cbar_bias = fig1.colorbar(
    bias_mappables[0],
    ax=[axes[0, 0], axes[0, 1], axes[1, 0]],
    orientation="vertical",
    shrink=0.9,
    extend="both",
)
cbar_bias.set_label("SST bias (°C)")

cbar_woa = fig1.colorbar(
    im_woa,
    ax=axes[1, 1],
    orientation="vertical",
    shrink=0.9,
    extend="both",
)
cbar_woa.set_label("SST (°C)")

fig1.suptitle("Figure 1. 1971-2000 SST climatology biases on shared WOA23 1° comparison cells (98th percentile clipping)")

# --- Figure 2: bar chart of own-grid global means ---
fig2, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

bar_labels = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR", "WOA23"]
bar_values = [
    native_global_means["CESM2"],
    native_global_means["GFDL-ESM4"],
    native_global_means["MPI-ESM1-2-LR"],
    woa_native_global_mean,
]

bars = ax.bar(bar_labels, bar_values, color=["tab:blue", "tab:orange", "tab:green", "tab:gray"])
ax.set_ylabel("Global mean SST (°C)")
ax.set_title("Figure 2. Global mean 1971-2000 SST climatology on each native ocean grid")

for bar, val in zip(bars, bar_values):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height(),
        f"{val:.2f}",
        ha="center",
        va="bottom",
    )
