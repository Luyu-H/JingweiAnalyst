# ---------------------------------------------------------------
# Figure 1: 2x2 bias maps (shared diverging scale) + WOA23 reference
# ---------------------------------------------------------------
lon = woa_sst["lon"].values
lat = woa_sst["lat"].values

pooled = np.concatenate([bias_field[m].values.ravel() for m in models])
vmax = float(np.nanpercentile(np.abs(pooled), 98))

ref_vals = woa_sst.where(mask).values
rvmin = float(np.nanpercentile(ref_vals, 1))
rvmax = float(np.nanpercentile(ref_vals, 99))

fig, axes = plt.subplots(2, 2, figsize=(14, 8.0), constrained_layout=True)
axes = axes.ravel()

for i, m in enumerate(models):
    ax = axes[i]
    pcm = ax.pcolormesh(lon, lat, bias_field[m].values, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    ax.set_title(f"{m} − WOA23")
    cb = fig.colorbar(pcm, ax=ax, extend="both")
    cb.set_label("SST bias (°C)")

ax = axes[3]
pcm = ax.pcolormesh(lon, lat, ref_vals, cmap="viridis",
                   vmin=rvmin, vmax=rvmax, shading="auto")
ax.set_title("WOA23 SST (0 m, 1971–2000)")
cb = fig.colorbar(pcm, ax=ax)
cb.set_label("SST (°C)")

for ax in axes:
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(float(np.nanmin(lon)), float(np.nanmax(lon)))
    ax.set_ylim(float(np.nanmin(lat)), float(np.nanmax(lat)))

fig.suptitle("SST bias vs WOA23 on the common comparison cells "
             f"(diverging scale shared across models, clipped at 98th pct = ±{vmax:.1f} °C)",
             fontsize=11)

# ---------------------------------------------------------------
# Figure 2: global mean SST of each model (native grid) and WOA23
# ---------------------------------------------------------------
labels = models + ["WOA23"]
values = [native_mean[m] for m in models] + [woa_native_mean]
colors = ["#1f77b4", "#2ca02c", "#9467bd", "#d62728"]

fig2, ax2 = plt.subplots(figsize=(7.5, 5.2), constrained_layout=True)
bars = ax2.bar(labels, values, color=colors)
for b, v in zip(bars, values):
    ax2.annotate(f"{v:.2f}",
                 xy=(b.get_x() + b.get_width() / 2.0, v),
                 xytext=(0, 3), textcoords="offset points",
                 ha="center", va="bottom", fontsize=10)

pad = 0.6
ax2.set_ylim(min(values) - pad, max(values) + pad)
ax2.set_ylabel("Global mean SST (°C)")
ax2.set_title("1971–2000 mean SST over each model's native ocean grid")
ax2.grid(axis="y", alpha=0.3)
