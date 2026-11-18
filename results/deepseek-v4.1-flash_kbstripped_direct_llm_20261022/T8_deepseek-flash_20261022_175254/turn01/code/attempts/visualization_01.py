import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Figure 1: 0-50 m DIC bias maps, shared diverging scale
# ---------------------------------------------------------------------------
all_abs = np.concatenate(
    [np.abs(dic_bias_maps[m][np.isfinite(dic_bias_maps[m])]).ravel() for m in MODELS]
)
vmax = np.nanpercentile(all_abs, 98) if all_abs.size else 50.0
vmax = max(vmax, 1.0)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
for ax, m in zip(axes, MODELS):
    pm = ax.pcolormesh(
        glodap_lon,
        glodap_lat,
        dic_bias_maps[m],
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_title(m)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(0, 360)
    ax.set_ylim(-90, 90)

cbar = fig.colorbar(
    pm, ax=axes, orientation="horizontal", fraction=0.05, pad=0.08, shrink=0.85
)
cbar.set_label(r"DIC bias (model $-$ GLODAP), $\mu$mol kg$^{-1}$")
fig.suptitle("0-50 m mean DIC bias against GLODAPv2.2016b", y=1.02)


# ---------------------------------------------------------------------------
# Figure 2: vertical profile of global mean DIC bias
# ---------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(5.5, 5.5))
for m in MODELS:
    ax2.plot(
        vertical_bias_df.loc[m].values,
        profile_depths,
        marker="o",
        label=m,
    )
ax2.axvline(0.0, color="k", linewidth=0.8)
ax2.invert_yaxis()
ax2.set_xlabel(r"DIC bias (model $-$ GLODAP), $\mu$mol kg$^{-1}$")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Global mean DIC bias profile")
for d in profile_depths:
    ax2.axhline(d, color="0.7", linewidth=0.4, linestyle=":")
ax2.legend(frameon=False)
