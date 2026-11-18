# Figure 1: 2002 0-50 m DIC bias maps on the GLODAP 1-degree grid
finite_abs = []
for m in models:
    b = bias_maps[m]
    finite_abs.append(b[np.isfinite(b)])
finite_abs = np.concatenate(finite_abs)

if finite_abs.size > 0:
    vmax = float(np.nanpercentile(np.abs(finite_abs), 98))
else:
    vmax = 1.0

fig1, axes1 = plt.subplots(1, 3, figsize=(18, 4.5), constrained_layout=True)
pcm = None
for ax, m in zip(axes1, models):
    pcm = ax.pcolormesh(
        glodap_lon, glodap_lat, bias_maps[m],
        shading="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax
    )
    ax.set_title(m)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_yticks(np.arange(-90, 91, 30))
    ax.set_xlim(0, 360)
    ax.set_ylim(-90, 90)

cb = fig1.colorbar(
    pcm, ax=axes1, orientation="horizontal", shrink=0.85,
    pad=0.08, extend="both"
)
cb.set_label("DIC bias (µmol kg$^{-1}$)")
fig1.suptitle(
    "2002 0–50 m DIC bias (model − GLODAPv2.2016b); shared diverging scale, "
    "clipped at pooled 98th percentile of |bias|",
    fontsize=12
)

# Figure 2: vertical profile of global mean DIC bias
fig2, ax2 = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
for m in models:
    ax2.plot(
        profile_bias[m], depths,
        marker="o", linewidth=1.6, label=m
    )
ax2.axvline(0.0, color="k", linewidth=0.8)
ax2.set_yticks(depths)
ax2.set_ylim(2100, -50)
ax2.set_xlabel("DIC bias (µmol kg$^{-1}$)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Global mean DIC bias profile (model − GLODAPv2.2016b)")
ax2.grid(True, axis="y", linestyle=":", alpha=0.6)
ax2.legend(loc="best", fontsize=9)
