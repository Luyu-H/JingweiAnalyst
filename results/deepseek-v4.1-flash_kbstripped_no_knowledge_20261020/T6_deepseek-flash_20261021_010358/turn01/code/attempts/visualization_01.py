# ============================================================================
# Figure 1 - global nitracline depth maps (2 x 2, shared sequential scale)
# ============================================================================
panels1 = [("CESM2",         model_nit_1deg["CESM2"]),
           ("GFDL-ESM4",     model_nit_1deg["GFDL-ESM4"]),
           ("MPI-ESM1-2-LR", model_nit_1deg["MPI-ESM1-2-LR"]),
           ("WOA23 (obs)",   woa_nit)]

pooled1 = np.concatenate([np.asarray(p.transpose("lat", "lon").values).ravel()
                          for _, p in panels1])
vmin1 = max(0.0, float(np.nanpercentile(pooled1, 1)))
vmax1 = float(np.nanpercentile(pooled1, 99))

cmap1 = plt.get_cmap("YlGnBu").copy()
cmap1.set_bad("0.85")

fig1, axes1 = plt.subplots(2, 2, figsize=(14, 8.0), constrained_layout=True)
for ax, (name, field) in zip(axes1.ravel(), panels1):
    da = field.transpose("lat", "lon")
    mesh1 = ax.pcolormesh(da["lon"].values, da["lat"].values, da.values,
                          cmap=cmap1, vmin=vmin1, vmax=vmax1, shading="auto")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-70, 75)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-60, -30, 0, 30, 60])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    if name == "WOA23 (obs)":
        ax.set_title(f"{name}\nglobal mean {woa_nit_global:.0f} m", fontweight="bold")
    else:
        ax.set_title(f"{name}\nglobal mean {model_nit_global[name]:.0f} m "
                     f"({model_nit_global[name] - woa_nit_global:+.0f} m), "
                     f"r={corr_val[name]:.3f}, RMSE={rmse_val[name]:.1f} m")

cb1 = fig1.colorbar(mesh1, ax=axes1.ravel().tolist(), shrink=0.9, extend="both")
cb1.set_label("Nitracline depth (m)")
fig1.suptitle("Figure 1 - Nitracline depth (depth of max vertical nitrate gradient, "
              "top 1000 m); shared scale clipped at 1st/99th percentile",
              fontsize=12, fontweight="bold")

# ============================================================================
# Figure 2 - equatorial Pacific nitrate sections (2S-2N mean), 120E-280E
# ============================================================================
panels2 = [("WOA23 (obs)",   woa_sec_mean,   woa_sec_nit),
           ("CESM2",         model_sec["CESM2"],         model_sec_nit["CESM2"]),
           ("GFDL-ESM4",     model_sec["GFDL-ESM4"],     model_sec_nit["GFDL-ESM4"]),
           ("MPI-ESM1-2-LR", model_sec["MPI-ESM1-2-LR"], model_sec_nit["MPI-ESM1-2-LR"])]

pooled2 = np.concatenate([np.asarray(da.sel(depth=slice(0, 1000))
                                     .transpose("depth", "lon").values).ravel()
                          for _, da, _ in panels2])
vmin2 = float(np.nanpercentile(pooled2, 1))
vmax2 = float(np.nanpercentile(pooled2, 99))

cmap2 = plt.get_cmap("YlGnBu").copy()
cmap2.set_bad("0.85")

xt = [120, 150, 180, 210, 240, 270]
xl = ["120°E", "150°E", "180°", "150°W", "120°W", "90°W"]

fig2, axes2 = plt.subplots(4, 1, figsize=(12, 13), constrained_layout=True)
for ax, (name, da, nit) in zip(axes2, panels2):
    d = da.sel(depth=slice(0, 1000)).transpose("depth", "lon")
    mesh2 = ax.pcolormesh(d["lon"].values, d["depth"].values, d.values,
                          cmap=cmap2, vmin=vmin2, vmax=vmax2, shading="auto")
    nl = nit["lon"].values
    nd = nit.values
    ax.plot(nl, nd, color="red", lw=1.8, label="nitracline")
    ax.set_ylim(0, 1000)
    ax.invert_yaxis()
    ax.set_xlim(float(sec_lon[0]), float(sec_lon[-1]))
    ax.set_xticks(xt)
    ax.set_xticklabels(xl)
    ax.set_ylabel("Depth (m)")
    ax.set_title(f"{name} - mean section nitracline {float(np.nanmean(nd)):.0f} m")
    ax.legend(loc="lower left", framealpha=0.85)

axes2[-1].set_xlabel("Longitude (equatorial Pacific, averaged over 2°S-2°N)")
cb2 = fig2.colorbar(mesh2, ax=axes2.tolist(), shrink=0.9, extend="both")
cb2.set_label(r"Nitrate ($\mu$mol kg$^{-1}$)")
fig2.suptitle("Figure 2 - Equatorial Pacific nitrate sections (1971-2000), "
              "top 1000 m; shared scale clipped at 1st/99th percentile",
              fontsize=12, fontweight="bold")
