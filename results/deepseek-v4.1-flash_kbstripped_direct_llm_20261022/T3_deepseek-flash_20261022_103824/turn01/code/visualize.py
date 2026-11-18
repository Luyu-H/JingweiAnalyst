# ============================================================
# Visualization block
# ============================================================

# -----------------------------
# Figure 1: 2x2 maps
# -----------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)

# Shared diverging scale for the three bias panels.
vmax_bias = 0.0
for model in models:
    vals = np.abs(model_bias_on_woa[model][common_mask])
    if vals.size > 0:
        vmax_bias = max(vmax_bias, float(np.nanmax(vals)))
if vmax_bias == 0.0:
    vmax_bias = 1.0

for ax, model in zip(axes.flat[:3], models):
    data = np.ma.masked_where(~common_mask, model_bias_on_woa[model])
    im = ax.pcolormesh(
        lon_woa, lat_woa, data,
        shading="auto",
        cmap="RdBu_r",
        vmin=-vmax_bias,
        vmax=vmax_bias,
    )
    ax.set_title(f"{model} bias (model - WOA23)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(lon_woa[0], lon_woa[-1])
    ax.set_ylim(lat_woa[0], lat_woa[-1])
    cb = fig.colorbar(im, ax=ax, orientation="vertical", shrink=0.8)
    cb.set_label("Delta SST (degC)")

# WOA23 reference SST on the same common cells.
ax = axes.flat[3]
data = np.ma.masked_where(~common_mask, woa_sst)
im = ax.pcolormesh(
    lon_woa, lat_woa, data,
    shading="auto",
    cmap="viridis",
)
ax.set_title("WOA23 SST (reference, 1971-2000 annual, 0 m)")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_xlim(lon_woa[0], lon_woa[-1])
ax.set_ylim(lat_woa[0], lat_woa[-1])
cb = fig.colorbar(im, ax=ax, orientation="vertical", shrink=0.8)
cb.set_label("SST (degC)")

# -----------------------------
# Figure 2: native-grid global mean bar chart
# -----------------------------
fig2, ax2 = plt.subplots(figsize=(7, 5))

names = models + ["WOA23"]
values = [global_mean_models[m] for m in models] + [woa_global_mean]
colors = ["C0", "C1", "C2", "0.5"]

bars = ax2.bar(names, values, color=colors, edgecolor="black")
for bar, val in zip(bars, values):
    ax2.text(
        bar.get_x() + bar.get_width() / 2.0,
        val,
        f"{val:.2f}",
        ha="center",
        va="bottom",
        fontsize=10,
    )

ax2.set_ylabel("SST (degC)")
ax2.set_title("Global mean SST, 1971-2000, native ocean grids")
ax2.set_ylim(min(values) - 1.0, max(values) + 1.0)
ax2.tick_params(axis="x", rotation=20)
