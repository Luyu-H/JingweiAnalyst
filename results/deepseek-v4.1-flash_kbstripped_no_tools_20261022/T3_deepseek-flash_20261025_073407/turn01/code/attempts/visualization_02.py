# ----------------------------------------------------------------------
# Figure 1 -- 2x2 global map panel of SST bias (model - WOA23) on the
# shared comparison cells + the WOA23 0 m reference on the same cells
# ----------------------------------------------------------------------
fig1, axs = plt.subplots(2, 2, figsize=(15, 9), constrained_layout=True)
axs = axs.ravel()

_pool = np.concatenate([np.asarray(bias_maps[m].values, float).ravel() for m in models])
_pool = _pool[np.isfinite(_pool)]
vmax = float(np.nanpercentile(np.abs(_pool), 98))       # robust shared limit

for ax, m in zip(axs[:3], models):
    # NOTE: `extend` is not a valid QuadMesh keyword in current matplotlib;
    #       pass it to the colorbar instead so the triangular ends are drawn.
    pm = ax.pcolormesh(woa_lon, woa_lat, bias_maps[m].values,
                       cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       shading="auto")
    ax.set_title(f"{m} - WOA23   (delta_SST = {delta_sst[m]:+.2f} degC, "
                 f"RMSE = {rmse_sst[m]:.3f} degC, r = {corr_sst[m]:.3f})", fontsize=10)
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
    cb = fig1.colorbar(pm, ax=ax, orientation="vertical", shrink=0.9, pad=0.02,
                       extend="both")
    cb.set_label("SST bias (degC)")

_ref = woa_sst_shared.values
vo = np.nanpercentile(_ref, 2); ve = np.nanpercentile(_ref, 98)
pm4 = axs[3].pcolormesh(woa_lon, woa_lat, _ref, cmap="viridis",
                        vmin=vo, vmax=ve, shading="auto")
axs[3].set_title("WOA23 0 m reference SST (1971-2000 mean, shared cells)", fontsize=10)
axs[3].set_xlabel("Longitude"); axs[3].set_ylabel("Latitude")
axs[3].set_xlim(-180, 180); axs[3].set_ylim(-90, 90)
cb4 = fig1.colorbar(pm4, ax=axs[3], orientation="vertical", shrink=0.9, pad=0.02,
                    extend="both")
cb4.set_label("SST (degC)")

fig1.suptitle("Bias colour scale clipped at the pooled 98th percentile of |bias| "
              f"(+/-{vmax:.2f} degC); shared mask = {mask_count} cells", fontsize=11)

# ----------------------------------------------------------------------
# Figure 2 -- four global mean SST values, each over its own grid
# ----------------------------------------------------------------------
fig2, ax = plt.subplots(figsize=(8.5, 6), constrained_layout=True)
labels = models + ["WOA23"]
values = [native_mean[m] for m in models] + [woa_global_mean_unweighted]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#9467bd"]
bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.6)

for b, v in zip(bars, values):
    ax.annotate(f"{v:.2f}", (b.get_x() + b.get_width() / 2.0, b.get_height()),
                xytext=(0, 4), textcoords="offset points",
                ha="center", va="bottom", fontsize=11, fontweight="bold")

ax.set_ylabel("Global mean SST (degC)")
ax.set_title("Global mean SST 1971-2000\n"
             "CMIP6: native-grid area-weighted (own areacello); "
             "WOA23: whole-1deg-grid unweighted mean")
ax.set_ylim(0, max(values) * 1.18)
ax.grid(axis="y", alpha=0.3)
fig2.text(0.5, -0.02,
          "Note: model bars are area-weighted over each native ocean grid; the WOA23 bar "
          "is the unweighted whole-1deg-grid mean, so the four bars are not identical quantities.",
          ha="center", fontsize=8, wrap=True)
