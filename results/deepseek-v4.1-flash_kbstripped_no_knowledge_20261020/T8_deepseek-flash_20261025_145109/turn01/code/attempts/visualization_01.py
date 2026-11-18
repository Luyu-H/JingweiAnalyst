# ------------------------- Figure 1: DIC bias maps -------------------------
all_diffs = np.concatenate(
    [dic_bias_maps[model].values.ravel() for model in models]
)
vmax = float(np.nanpercentile(np.abs(all_diffs), 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0

fig, axes = plt.subplots(
    1, 3, figsize=(15, 4.2), constrained_layout=True
)
for ax, model in zip(axes, models):
    d = dic_bias_maps[model]
    lon = d.lon.values
    lat = d.lat.values
    im = ax.pcolormesh(
        lon,
        lat,
        d.values,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_title(model)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim([lon.min(), lon.max()])
    ax.set_ylim([lat.min(), lat.max()])

fig.colorbar(im, ax=axes, shrink=0.82, label="DIC bias (umol kg-1)")
fig.suptitle(
    "0–50 m DIC bias vs GLODAP "
    f"(shared diverging scale, 98th-percentile |bias| = {vmax:.0f} umol kg-1)"
)

# --------------------- Figure 2: vertical bias profile ---------------------
fig2, ax2 = plt.subplots(figsize=(5.4, 6.0), constrained_layout=True)
for model in models:
    prof = profile_bias[model]
    depths = [p[0] for p in prof]
    biases = [p[1] for p in prof]
    ax2.plot(biases, depths, marker="o", linewidth=1.6, label=model)

ax2.axvline(0.0, color="k", linestyle="--", linewidth=0.9)
ax2.invert_yaxis()
ax2.set_yticks(depths_profile)
ax2.set_xlabel("DIC bias (umol kg-1)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Global mean DIC bias vs GLODAP")
ax2.legend(loc="best")
