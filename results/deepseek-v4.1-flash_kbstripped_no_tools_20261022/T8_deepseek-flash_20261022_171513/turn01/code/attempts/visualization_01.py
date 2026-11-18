# Visualization code only

# ---------------------------------------------------------------------------
# Figure 1: 0-50 m DIC bias maps against GLODAP
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(
    1, 3, figsize=(15, 4.2), constrained_layout=True
)

pooled_abs = np.concatenate(
    [np.abs(dic_bias_maps[m].values).ravel() for m in MODEL_LIST]
)
vmax = float(np.nanpercentile(pooled_abs, 98.0))

for ax, model in zip(axes, MODEL_LIST):
    im = ax.pcolormesh(
        g_lon_vals,
        g_lat_vals,
        dic_bias_maps[model].values,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_title(model)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(float(g_lon_vals.min()), float(g_lon_vals.max()))
    ax.set_ylim(float(g_lat_vals.min()), float(g_lat_vals.max()))

fig.colorbar(
    im,
    ax=axes.ravel().tolist(),
    label="DIC bias (µmol kg$^{-1}$)",
    extend="both",
)
fig.suptitle(
    "0–50 m mean DIC bias vs GLODAPv2.2016b; "
    "shared scale clipped at pooled 98th percentile of |bias|"
)

# ---------------------------------------------------------------------------
# Figure 2: vertical profile of global mean DIC bias
# ---------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(7.2, 5.0), constrained_layout=True)

for model in MODEL_LIST:
    ax2.plot(
        model_dic_profile_bias[model],
        DEPTHS,
        marker="o",
        linewidth=1.8,
        markersize=6,
        label=model,
    )

ax2.axvline(0.0, color="black", linewidth=0.9, linestyle="--")
ax2.invert_yaxis()
ax2.set_yticks(DEPTHS)
ax2.set_xlabel("DIC bias (µmol kg$^{-1}$)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Global mean DIC bias profile (model minus GLODAP)")
ax2.grid(True, alpha=0.25)
ax2.legend(loc="best")
