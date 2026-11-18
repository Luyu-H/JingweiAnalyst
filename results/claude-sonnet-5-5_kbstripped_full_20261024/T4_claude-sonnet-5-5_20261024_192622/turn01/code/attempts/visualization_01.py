keys = [(m, mn) for m in models for mn in ["March", "September"]]

# ---------- Figure 1: difference maps ----------
pooled = np.concatenate([np.abs(diff_map[k].values.ravel()) for k in keys])
vmax = float(np.nanpercentile(pooled, 98))
fig1, axes1 = plt.subplots(3, 2, figsize=(14, 11), constrained_layout=True)
mesh = None
for i, m in enumerate(models):
    for j, mn in enumerate(["March", "September"]):
        ax = axes1[i, j]
        d = diff_map[(m, mn)]
        ax.set_facecolor("0.85")
        mesh = ax.pcolormesh(d["lon"].values, d["lat"].values, d.values, cmap="RdBu_r",
                             vmin=-vmax, vmax=vmax, shading="auto")
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        ax.set_aspect("equal")
        r = results[(m, mn)]
        ax.set_title(f"{m}, {mn}  (bias {r['mean_diff_m']:.1f} m, RMSE {r['rmse_m']:.1f} m)", fontsize=10)
        if j == 0:
            ax.set_ylabel("Latitude (°N)")
        if i == 2:
            ax.set_xlabel("Longitude (°E)")
cb = fig1.colorbar(mesh, ax=axes1, orientation="horizontal", fraction=0.04, pad=0.03, extend="both")
cb.set_label("Operator MLD - mlotst (m)")
fig1.suptitle(f"Operator MLD minus mlotst, 1971-2000 climatology, shared mask "
              f"(colour scale clipped at the 98th percentile of |difference| = {vmax:.0f} m)", fontsize=12)

# ---------- Figure 2: scatter plots ----------
lim = (1, 1e4)
ticks = [1, 10, 100, 1000, 10000]
fig2, axes2 = plt.subplots(3, 2, figsize=(10, 14), constrained_layout=True)
for i, m in enumerate(models):
    for j, mn in enumerate(["March", "September"]):
        ax = axes2[i, j]
        x, y = scatter[(m, mn)]
        ax.scatter(x, y, s=1, alpha=0.25, color="tab:blue", rasterized=True)
        ax.plot(lim, lim, "k--", lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xticklabels([str(t) for t in ticks])
        ax.set_yticklabels([str(t) for t in ticks])
        ax.set_aspect("equal")
        r = results[(m, mn)]
        ax.text(0.04, 0.96, f"bias = {r['mean_diff_m']:.1f} m\nRMSE = {r['rmse_m']:.1f} m\nr = {r['corr']:.3f}",
                transform=ax.transAxes, va="top", ha="left", fontsize=9,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.6"))
        ax.set_title(f"{m}, {mn}", fontsize=10)
        ax.set_xlabel("Operator MLD (m)")
        ax.set_ylabel("mlotst (m)")
        ax.grid(True, which="major", alpha=0.3)
fig2.suptitle("mlotst vs operator MLD, one point per shared-mask ocean cell", fontsize=12)
