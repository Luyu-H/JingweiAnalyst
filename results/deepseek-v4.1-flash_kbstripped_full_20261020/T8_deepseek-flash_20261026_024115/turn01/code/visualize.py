# ---------------------------------------------------------------------------
# Figure 1: 0-50 m DIC bias maps (model minus GLODAP), shared diverging scale
# ---------------------------------------------------------------------------
fig1, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

# Robust shared colour limit from the pooled bias fields
pooled = np.concatenate(
    [np.asarray(dic_diff_fields[m].values).ravel() for m in models]
)
pooled = pooled[np.isfinite(pooled)]
vmax = float(np.nanpercentile(np.abs(pooled), 98)) if pooled.size else 1.0
vmax = max(vmax, 1.0)

for ax, m in zip(axes, models):
    field = dic_diff_fields[m]
    pcm = ax.pcolormesh(field.lon,
                        field.lat,
                        field,
                        cmap="RdBu_r",
                        vmin=-vmax,
                        vmax=vmax,
                        shading="auto")
    ax.set_title(m)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cbar = fig1.colorbar(pcm, ax=axes, orientation="horizontal",
                     fraction=0.05, pad=0.05, extend="both")
cbar.set_label("DIC bias (µmol kg-1)")
fig1.suptitle(
    "0–50 m DIC bias (model - GLODAP); shared scale clipped at "
    "98th percentile of |bias|",
    fontsize=12,
)

# ---------------------------------------------------------------------------
# Figure 2: vertical profile of global mean DIC bias
# ---------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(7, 5), constrained_layout=True)

for m in models:
    ax2.plot(results[m]["prof_bias"],
             profile_depths,
             marker="o",
             label=m)

ax2.axvline(0, color="k", lw=1, ls="--")
ax2.invert_yaxis()
ax2.set_yticks(profile_depths)
ax2.set_xlabel("DIC bias (µmol kg-1)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Global mean DIC bias profile (model - GLODAP)")
ax2.legend(loc="best")

plt.show()
