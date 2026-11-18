fig, axes = plt.subplots(2, 2, figsize=(15, 8.5), constrained_layout=True)
pooled = np.concatenate([np.abs(diff_maps[m].values).ravel() for m in models])
vmax = float(np.nanpercentile(pooled, 98))

panels = [(axes[0, 0], "CESM2"), (axes[0, 1], "GFDL-ESM4"), (axes[1, 0], "MPI-ESM1-2-LR")]
for ax, m in panels:
    d = diff_maps[m]
    im = ax.pcolormesh(d["lon"], d["lat"], d.values, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    ax.set_title(f"{m} SST bias (model - WOA23), mean {scores[m]['bias']:.2f} °C")
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
bias_im = im
cb = fig.colorbar(bias_im, ax=[axes[0, 0], axes[0, 1], axes[1, 0]], extend="both",
                  orientation="vertical", shrink=0.85, pad=0.02)
cb.set_label("SST bias (°C)")

ax = axes[1, 1]
im2 = ax.pcolormesh(woa_m["lon"], woa_m["lat"], woa_m.values, cmap="viridis", shading="auto")
ax.set_title("WOA23 0 m SST, 1971-2000 (common mask)")
ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
cb2 = fig.colorbar(im2, ax=ax, orientation="vertical", shrink=0.85, pad=0.02)
cb2.set_label("SST (°C)")
fig.suptitle(f"Biases shared scale ±{vmax:.2f} °C (98th percentile of |bias|, clipped)")

# Figure 2
fig2, ax2 = plt.subplots(figsize=(8, 5), constrained_layout=True)
names = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR", "WOA23"]
vals = [native_mean[n] for n in names]
bars = ax2.bar(names, vals, color=["tab:blue", "tab:orange", "tab:green", "tab:gray"])
for b, v in zip(bars, vals):
    ax2.annotate(f"{v:.2f}", (b.get_x() + b.get_width() / 2, v),
                 ha="center", va="bottom", xytext=(0, 3), textcoords="offset points")
ax2.set_ylabel("Global mean SST (°C)")
ax2.set_title("1971-2000 global mean SST on each product's native ocean grid")
ax2.set_ylim(0, max(vals) * 1.15)
plt.show()
