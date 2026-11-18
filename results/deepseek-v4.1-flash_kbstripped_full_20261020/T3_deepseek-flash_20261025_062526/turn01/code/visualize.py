lon = woa_ref_map["lon"].values
lat = woa_ref_map["lat"].values

bias_arrays = [bias_map[m].values.ravel() for m in models]
vmax_bias = float(np.nanpercentile(np.abs(np.concatenate(bias_arrays)), 98))
if not np.isfinite(vmax_bias) or vmax_bias == 0:
    vmax_bias = 1.0

woa_vals = woa_ref_map.values
vmin_woa = float(np.nanpercentile(woa_vals, 2))
vmax_woa = float(np.nanpercentile(woa_vals, 98))
if not np.isfinite(vmin_woa):
    vmin_woa = float(np.nanmin(woa_vals))
if not np.isfinite(vmax_woa):
    vmax_woa = float(np.nanmax(woa_vals))

fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
axes = axes.ravel()

for ax, model in zip(axes[:3], models):
    pcm = ax.pcolormesh(
        lon, lat, bias_map[model].values,
        cmap="RdBu_r", vmin=-vmax_bias, vmax=vmax_bias,
        shading="auto"
    )
    ax.set_title(f"{model} SST bias (model - WOA23)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    cb = fig.colorbar(pcm, ax=ax, orientation="vertical", shrink=0.9, extend="both")
    cb.set_label("SST bias (°C)")

ax = axes[3]
pcm = ax.pcolormesh(
    lon, lat, woa_ref_map.values,
    cmap="viridis", vmin=vmin_woa, vmax=vmax_woa,
    shading="auto"
)
ax.set_title("WOA23 0 m SST (shared cells)")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)
cb = fig.colorbar(pcm, ax=ax, orientation="vertical", shrink=0.9, extend="both")
cb.set_label("SST (°C)")

fig.suptitle(
    f"SST bias on shared comparison cells; bias scale clipped at 98th percentile "
    f"of |bias| = {vmax_bias:.2f} °C. WOA23 panel uses 2nd-98th percentile.",
    fontsize=11
)

fig2, ax2 = plt.subplots(figsize=(8, 5), constrained_layout=True)
labels = models + ["WOA23"]
vals = [native_gm[m] for m in models] + [woa_global_gm]
colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
bars = ax2.bar(labels, vals, color=colors)
for bar, v in zip(bars, vals):
    ax2.text(
        bar.get_x() + bar.get_width() / 2, bar.get_height(),
        f"{v:.2f}", ha="center", va="bottom", fontsize=10
    )

ax2.set_ylabel("Global mean SST (°C)")
ax2.set_title(
    "Native-grid global mean SST, 1971-2000\n"
    "(each bar over its own native ocean grid; not the shared comparison cells)"
)
ax2.set_ylim(min(vals) - 0.5, max(vals) + 0.5)
ax2.grid(axis="y", alpha=0.3)
ax2.tick_params(axis="x", rotation=15)
