# ---------------------------------------------------------------------------
# Figure 1: bias maps on the shared WOA23 comparison cells
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)

model_list = list(models.keys())

bias_values = []
for m in model_list:
    v = bias_maps[m].values
    bias_values.append(v[np.isfinite(v)])
pooled_bias = np.concatenate(bias_values)

vmax = float(np.nanpercentile(np.abs(pooled_bias), 98))
if not np.isfinite(vmax) or vmax <= 0:
    vmax = 1.0

for ax, m in zip(axes.flat[:3], model_list):
    pcm = ax.pcolormesh(
        woa_temp.lon,
        woa_temp.lat,
        bias_maps[m].values,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_title(f"{m} minus WOA23 SST")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    cb = fig.colorbar(pcm, ax=ax, orientation="vertical", shrink=0.85, extend="both")
    cb.set_label("SST bias (°C)")

woa_plot = woa_temp.where(comparison_mask)
woa_vals = woa_plot.values[np.isfinite(woa_plot.values)]
vmin_w, vmax_w = np.nanpercentile(woa_vals, [2, 98])

ax = axes.flat[3]
pcm = ax.pcolormesh(
    woa_temp.lon,
    woa_temp.lat,
    woa_plot.values,
    cmap="viridis",
    vmin=vmin_w,
    vmax=vmax_w,
    shading="auto",
)
ax.set_title("WOA23 1971-2000 SST at 0 m")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
cb = fig.colorbar(pcm, ax=ax, orientation="vertical", shrink=0.85, extend="both")
cb.set_label("WOA23 SST (°C)")

fig.suptitle(
    "SST on shared WOA23 comparison cells; bias scale clipped at 98th percentile of |bias|",
    fontsize=12,
)

# ---------------------------------------------------------------------------
# Figure 2: global mean SST on each product's own native ocean grid
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

labels = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR", "WOA23"]
values = [
    model_native_means["CESM2"],
    model_native_means["GFDL-ESM4"],
    model_native_means["MPI-ESM1-2-LR"],
    woa_global_mean,
]
colors = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd"]

bars = ax.bar(labels, values, color=colors)

for bar, val in zip(bars, values):
    ax.annotate(
        f"{val:.2f}",
        (bar.get_x() + bar.get_width() / 2.0, bar.get_height()),
        textcoords="offset points",
        xytext=(0, 3),
        ha="center",
        va="bottom",
    )

ax.set_ylabel("Mean SST (°C)")
ax.set_title("1971-2000 mean SST on each native ocean grid")
ax.set_ylim(min(values) - 0.8, max(values) + 0.8)
