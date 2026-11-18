# Visualization code only

# Figure 1: 3x2 maps of offline MLD minus mlotst on a common display grid
fig1, axes1 = plt.subplots(3, 2, figsize=(14, 11), constrained_layout=True)

pooled = np.concatenate([x[np.isfinite(x)] for x in display_maps if np.isfinite(x).any()])
vmax = float(np.nanpercentile(np.abs(pooled), 98))
vmin = -vmax

for i, model in enumerate(models):
    for j, month in enumerate(target_months):
        ax = axes1[i, j]
        da = diff_maps_native[model][f"{month}_display"]
        lon = da["lon"].values
        lat = da["lat"].values
        im = ax.pcolormesh(lon, lat, da.values, cmap="RdBu_r", vmin=vmin, vmax=vmax, shading="auto")
        ax.set_title(f"{model} - {month_names[month]}")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_xlim(lon.min(), lon.max())
        ax.set_ylim(lat.min(), lat.max())

cbar1 = fig1.colorbar(im, ax=axes1, orientation="vertical", shrink=0.9, extend="both")
cbar1.set_label("Offline MLD - mlotst (m)")
fig1.suptitle("Figure 1. CMIP6 internal MLD consistency: native-grid differences regridded for display only\nShared color scale clipped at pooled 98th percentile", y=1.02)

# Figure 2: 3x2 scatter plots with identical logarithmic axes
fig2, axes2 = plt.subplots(3, 2, figsize=(12, 12), constrained_layout=True)

all_xy = []
for model in models:
    for month in target_months:
        x = scatter_inputs[model][month]["x_offline"]
        y = scatter_inputs[model][month]["y_mlotst"]
        valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
        all_xy.append(np.concatenate([x[valid], y[valid]]))
all_xy = np.concatenate(all_xy)
xy_min = float(np.nanpercentile(all_xy, 1))
xy_max = float(np.nanpercentile(all_xy, 99.5))

ticks = np.array([1, 3, 10, 30, 100, 300, 1000], dtype=float)
ticks = ticks[(ticks >= xy_min) & (ticks <= xy_max)]

for i, model in enumerate(models):
    for j, month in enumerate(target_months):
        ax = axes2[i, j]
        x = scatter_inputs[model][month]["x_offline"]
        y = scatter_inputs[model][month]["y_mlotst"]
        valid = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)

        ax.scatter(x[valid], y[valid], s=3, alpha=0.25, linewidths=0)
        ax.plot([xy_min, xy_max], [xy_min, xy_max], "--", color="k", lw=1)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(xy_min, xy_max)
        ax.set_ylim(xy_min, xy_max)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks])
        ax.set_yticklabels([f"{t:g}" for t in ticks])

        rr = results[model][month]
        ax.set_title(f"{model} - {month_names[month]}")
        ax.set_xlabel("Offline MLD (m)")
        ax.set_ylabel("mlotst (m)")
        ax.text(
            0.03,
            0.97,
            f"Bias = {rr['bias_m']:.1f} m\nRMSD = {rr['rmsd_m']:.1f} m\nr = {rr['corr']:.3f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
        )

fig2.suptitle("Figure 2. Offline MLD versus published mlotst\nOne point per shared valid ocean cell; identical log axes in all panels", y=1.02)
