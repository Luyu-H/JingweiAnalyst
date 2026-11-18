# ===============================================================
# FIGURE 1 - 2x2 panel: three bias maps + WOA23 reference
# ===============================================================
fig, axes = plt.subplots(2, 2, figsize=(15, 9), constrained_layout=True)
panels = axes.ravel()

pooled = np.concatenate([np.abs(np.asarray(bias_fields[m].values)).ravel()
                         for m in models])
vmax = float(np.nanpercentile(pooled, 98))
vmax = max(vmax, 0.1)

for ax, m in zip(panels[:3], models):
    pcm = ax.pcolormesh(lon1, lat1, bias_fields[m].values,
                        cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                        shading="auto")
    ax.set_title(f"{m} − WOA23   (mean bias {bias[m]:+.2f} °C)", fontsize=11)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    cb = fig.colorbar(pcm, ax=ax, extend="both")
    cb.set_label("SST bias (°C)")

# fourth panel: WOA23 reference SST on the same cells
v1, v2 = np.nanpercentile(woa_ref_field.values, [1, 99])
pcm4 = panels[3].pcolormesh(lon1, lat1, woa_ref_field.values,
                            cmap="magma", vmin=float(v1), vmax=float(v2),
                            shading="auto")
panels[3].set_title(f"WOA23 0 m SST, 1971-2000  (mean {comp_mean_woa:.2f} °C)",
                    fontsize=11)
panels[3].set_xlabel("Longitude")
panels[3].set_ylabel("Latitude")
panels[3].set_xlim(-180, 180)
panels[3].set_ylim(-90, 90)
cb4 = fig.colorbar(pcm4, ax=panels[3])
cb4.set_label("SST (°C)")

fig.suptitle("Figure 1. SST climatology 1971–2000: CMIP6 historical vs WOA23 on the common "
             "1° comparison cells\n(Bias maps share a diverging scale clipped at the pooled "
             f"98th percentile, ±{vmax:.2f} °C)", fontsize=12)

# ===============================================================
# FIGURE 2 - bar chart of the four native-grid global means
# ===============================================================
labels = models + ["WOA23"]
values = [native_mean[m] for m in models] + [woa_global_mean]
colors = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]

fig2, ax2 = plt.subplots(figsize=(9, 6), constrained_layout=True)
bars = ax2.bar(labels, values, color=colors, edgecolor="black", linewidth=0.6)
ax2.set_ylabel("Sea surface temperature (°C)")
ax2.set_title("Figure 2. 1971–2000 global mean SST on each native ocean grid\n"
              "(CMIP6 day-weighted monthly climatologies; WOA23 annual 0 m normal, all valid 1° cells)")
ax2.set_ylim(min(values) - 1.0, max(values) + 0.8)
ax2.grid(axis="y", alpha=0.3)

for bar, val in zip(bars, values):
    ax2.annotate(f"{val:.2f} °C",
                 xy=(bar.get_x() + bar.get_width() / 2.0, bar.get_height()),
                 xytext=(0, 5), textcoords="offset points",
                 ha="center", va="bottom", fontsize=10, fontweight="bold")

for lbl in ax2.get_xticklabels():
    lbl.set_rotation(15)
    lbl.set_ha("right")
