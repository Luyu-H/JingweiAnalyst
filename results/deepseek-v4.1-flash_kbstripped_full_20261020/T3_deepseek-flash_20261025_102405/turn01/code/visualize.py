# Figure 1: 2x2 maps of SST biases and WOA23 reference SST
fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
axes = axes.ravel()

bias_fields = [bias_maps[m] for m in models]

# Robust shared diverging scale for the three bias panels.
all_bias = np.concatenate([np.asarray(f.values).ravel() for f in bias_fields])
all_bias = all_bias[np.isfinite(all_bias)]
vmax = float(np.nanpercentile(np.abs(all_bias), 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0
vmin = -vmax

lon = woa_shared_c["lon"].values
lat = woa_shared_c["lat"].values

for ax, model, bias in zip(axes[:3], models, bias_fields):
    im = ax.pcolormesh(
        lon,
        lat,
        bias.values,
        cmap="RdBu_r",
        vmin=vmin,
        vmax=vmax,
        shading="auto",
    )
    ax.set_title(f"{model} - WOA23 SST bias")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    cb = fig.colorbar(im, ax=ax, extend="both")
    cb.set_label("SST bias (degC)")

# WOA23 reference panel.
woa_vals = woa_shared_c.values
finite_woa = woa_vals[np.isfinite(woa_vals)]
vmin_woa = float(np.nanpercentile(finite_woa, 1))
vmax_woa = float(np.nanpercentile(finite_woa, 99))

im = axes[3].pcolormesh(
    lon,
    lat,
    woa_vals,
    cmap="viridis",
    vmin=vmin_woa,
    vmax=vmax_woa,
    shading="auto",
)
axes[3].set_title("WOA23 0 m SST (1971-2000)")
axes[3].set_xlabel("Longitude")
axes[3].set_ylabel("Latitude")
cb = fig.colorbar(im, ax=axes[3], extend="both")
cb.set_label("WOA23 SST (degC)")

fig.suptitle(
    "SST on shared WOA23 1° comparison cells; "
    "bias colour limits clipped at the 98th percentile of |bias|",
    fontsize=12,
)

# Figure 2: native-grid global mean SST bar chart
fig2, ax2 = plt.subplots(figsize=(8, 5), constrained_layout=True)

labels = models + ["WOA23"]
values = [model_means_native[m] for m in models] + [woa_mean_native]
areas = [model_areas_native[m] for m in models] + [woa_area_native]

bars = ax2.bar(labels, values, color=["C0", "C1", "C2", "C3"])
ax2.set_ylabel("SST (degC)")
ax2.set_title("1971-2000 global mean SST over each native ocean grid")

for bar, value, area in zip(bars, values, areas):
    ax2.text(
        bar.get_x() + bar.get_width() / 2.0,
        bar.get_height() + 0.15,
        f"{value:.2f} degC\n({area:.2e} m²)",
        ha="center",
        va="bottom",
        fontsize=9,
    )

ax2.set_ylim(min(values) - 1.5, max(values) + 2.5)
