# Figure 1: 2x2 maps of model biases and WOA23 reference
lon = woa_cmp.coords["lon"].values
lat = woa_cmp.coords["lat"].values

bias_arrays = {m: bias_fields[m].compute().values for m in models}
woa_vals = woa_cmp.compute().values

# Shared diverging scale for the three bias panels
pooled = np.concatenate([bias_arrays[m].ravel() for m in models])
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0

# Robust sequential scale for WOA23 reference
vmin_woa = float(np.nanpercentile(woa_vals, 2))
vmax_woa = float(np.nanpercentile(woa_vals, 98))
if not np.isfinite(vmin_woa) or not np.isfinite(vmax_woa) or vmin_woa == vmax_woa:
    vmin_woa = float(np.nanmin(woa_vals))
    vmax_woa = float(np.nanmax(woa_vals))

fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes = axes.ravel()

for ax, model in zip(axes[:3], models):
    im = ax.pcolormesh(
        lon, lat, bias_arrays[model],
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto",
    )
    ax.set_title(f"{model} SST bias (model - WOA23)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("SST bias (°C)")

ax = axes[3]
im = ax.pcolormesh(
    lon, lat, woa_vals,
    cmap="viridis", vmin=vmin_woa, vmax=vmax_woa, shading="auto",
)
ax.set_title("WOA23 0 m SST (reference)")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
cb = fig.colorbar(im, ax=ax, shrink=0.8)
cb.set_label("SST (°C)")

fig.suptitle(
    f"1971-2000 SST climatology on common comparison cells ({mask_ncells} cells); "
    f"bias panels share a scale clipped at |bias|={vmax:.2f} °C (98th percentile)",
    fontsize=11,
)

# Figure 2: native-grid global mean SST bar chart
labels = models + ["WOA23"]
values = [native_means[m] for m in models] + [woa_native_mean]

fig2, ax2 = plt.subplots(figsize=(8, 5), constrained_layout=True)
bars = ax2.bar(labels, values, color=["C0", "C1", "C2", "C3"])
ax2.set_ylabel("SST (°C)")
ax2.set_title("Native-grid global mean SST, 1971-2000")
ax2.grid(axis="y", linestyle="--", alpha=0.5)

for bar, val in zip(bars, values):
    ax2.annotate(
        f"{val:.2f} °C",
        xy=(bar.get_x() + bar.get_width() / 2, val),
        xytext=(0, 3),
        textcoords="offset points",
        ha="center", va="bottom",
    )

ax2.set_ylim(min(values) - 1.0, max(values) + 1.0)
