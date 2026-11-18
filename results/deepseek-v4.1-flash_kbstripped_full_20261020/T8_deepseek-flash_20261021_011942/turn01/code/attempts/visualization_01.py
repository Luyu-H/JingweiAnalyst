# ---------------------------------------------------------------------
# Figure 1: 0-50 m DIC bias maps (model - GLODAP) on the shared mask
# ---------------------------------------------------------------------
for m in models:
    bias_maps[m] = bias_maps[m].compute()

pooled = np.concatenate([np.asarray(bias_maps[m].values).ravel() for m in models])
pooled = pooled[np.isfinite(pooled)]
vmax = float(np.nanpercentile(np.abs(pooled), 98))

lons = np.asarray(bias_maps[models[0]]["lon"].values)
lats = np.asarray(bias_maps[models[0]]["lat"].values)

fig1, axes = plt.subplots(1, 3, figsize=(17, 4.2), constrained_layout=True)
for ax, m in zip(axes, models):
    pm = ax.pcolormesh(lons, lats, np.asarray(bias_maps[m].values),
                       cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto")
    ax.set_title(f"{m}  (bias {summary['map_bias'][m]:+.0f}, "
                 f"RMSE {summary['map_rmse'][m]:.0f}, r {summary['map_corr'][m]:.3f})",
                 fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig1.colorbar(pm, ax=list(axes), shrink=0.85, extend="both",
              label="DIC bias (µmol kg$^{-1}$)")
fig1.suptitle("0–50 m DIC bias vs GLODAPv2 (2002 annual mean); shared scale clipped at "
              "the pooled 98th percentile of |bias|", fontsize=11)

# ---------------------------------------------------------------------
# Figure 2: vertical profile of the global mean DIC bias
# ---------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(6.2, 5.2), constrained_layout=True)
for m in models:
    vals = [profile_bias[m][d] for d in target_depths]
    ax2.plot(vals, target_depths, marker="o", label=m)
ax2.axvline(0.0, color="k", linestyle="--", linewidth=1)
ax2.invert_yaxis()
ax2.set_yticks(target_depths)
ax2.set_ylabel("Depth (m)")
ax2.set_xlabel("DIC bias (µmol kg$^{-1}$)")
ax2.set_title("Global mean DIC bias vs GLODAPv2 (area-weighted on the GLODAP 1° grid)")
ax2.legend(loc="best", frameon=True)
ax2.grid(alpha=0.3)
