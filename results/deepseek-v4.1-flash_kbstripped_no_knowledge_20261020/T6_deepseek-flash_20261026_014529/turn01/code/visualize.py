# ============================================================================
# Figure 1 - 2x2 global maps of nitracline depth (shared sequential scale)
# ============================================================================
fig1, axes1 = plt.subplots(2, 2, figsize=(15, 9), constrained_layout=True)
fig1.suptitle("Nitracline depth (depth of maximum vertical nitrate gradient, top 1000 m)\n"
              "CMIP6 1971-2000 climatologies vs WOA23 observed climatology - shared scale (2-98th pct)",
              fontsize=13)

panels1 = [("CESM2", reg["CESM2"]),
           ("GFDL-ESM4", reg["GFDL-ESM4"]),
           ("MPI-ESM1-2-LR", reg["MPI-ESM1-2-LR"]),
           ("WOA23 (observed)", nc_woa)]

pooled1 = np.concatenate([np.asarray(f.values, dtype=float).ravel() for _, f in panels1])
pooled1 = pooled1[np.isfinite(pooled1)]
vmin1, vmax1 = np.nanpercentile(pooled1, [2, 98])

cmap1 = plt.get_cmap("YlGnBu").copy()
cmap1.set_bad("0.85")   # land
mesh1 = None
for ax, (name, f) in zip(axes1.ravel(), panels1):
    lon = np.asarray(f["lon"].values, dtype=float)
    lat = np.asarray(f["lat"].values, dtype=float)
    mesh1 = ax.pcolormesh(lon, lat, np.asarray(f.values, dtype=float),
                          cmap=cmap1, vmin=vmin1, vmax=vmax1,
                          shading="auto", rasterized=True)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-75, 80)
    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-60, 61, 30))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(f"{name}   (global mean {glob[name]:.0f} m)" if name in glob
                 else f"{name}   (global mean {glob['WOA23']:.0f} m)", fontsize=11)
    ax.grid(True, linewidth=0.25, alpha=0.4)

cb1 = fig1.colorbar(mesh1, ax=axes1.ravel().tolist(), shrink=0.85, extend="both",
                    label="Nitracline depth (m)")

# ============================================================================
# Figure 2 - 4-panel equatorial Pacific longitude-depth nitrate sections
# ============================================================================
fig2, axes2 = plt.subplots(4, 1, figsize=(13, 13), constrained_layout=True,
                           sharex=True)
fig2.suptitle("Nitrate along the equatorial Pacific (2\u00b0S-2\u00b02N), 120\u00b0E-80\u00b0W, top 1000 m\n"
              "CMIP6 1971-2000 climatology vs WOA23 annual climatology - "
              "shared scale, red line = nitracline (max gradient in top 1000 m)",
              fontsize=13)

sec_plot = {name: sec_fields[name].sel(depth=slice(0, 1000)) for name in panel_order}
pooled2 = np.concatenate([np.asarray(s.values, dtype=float).ravel() for s in sec_plot.values()])
pooled2 = pooled2[np.isfinite(pooled2)]
vmin2, vmax2 = np.nanpercentile(pooled2, [2, 98])

cmap2 = plt.get_cmap("YlGnBu").copy()
cmap2.set_bad("0.85")
mesh2 = None
for ax, name in zip(axes2, panel_order):
    s = sec_plot[name]
    lname = lon_name_of(s)
    x = np.asarray(s[lname].values, dtype=float)
    z = np.asarray(s["depth"].values, dtype=float)
    mesh2 = ax.pcolormesh(x, z, np.asarray(s.values, dtype=float),
                          cmap=cmap2, vmin=vmin2, vmax=vmax2,
                          shading="auto", rasterized=True)

    nl = sec_nc[name]
    xl = np.asarray(nl[lon_name_of(nl)].values, dtype=float)
    yl = np.asarray(nl.values, dtype=float)
    ax.plot(xl, yl, color="crimson", linewidth=2.0, zorder=5,
            label="nitracline")

    ax.set_ylim(1000, 0)            # depth increasing downward
    ax.set_xlim(120, 280)
    ax.set_xticks(np.arange(120, 281, 20))
    ax.set_ylabel("Depth (m)")
    ax.set_title(name, fontsize=11, loc="left")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.8)

axes2[-1].set_xlabel("Longitude (\u00b0E)")
cb2 = fig2.colorbar(mesh2, ax=axes2.tolist(), shrink=0.85, extend="both",
                    label="Nitrate (\u00b5mol kg$^{-1}$)")

# ============================================================================
# Figure 3 - summary of global and basin differences and skill scores
# ============================================================================
fig3, axes3 = plt.subplots(1, 3, figsize=(16, 4.5), constrained_layout=True)
fig3.suptitle("CMIP6 nitracline depth assessed against WOA23", fontsize=13)

# (a) global means
ax = axes3[0]
names = ["WOA23"] + models
vals = [glob[n] for n in names]
cols = ["0.35", "#4C72B0", "#DD8452", "#55A868"]
ax.bar(names, vals, color=cols)
ax.axhline(glob["WOA23"], color="0.35", linestyle="--", linewidth=1)
for i, v in enumerate(vals):
    ax.text(i, v + 2, f"{v:.0f}", ha="center", fontsize=9)
ax.set_ylabel("Nitracline depth (m)")
ax.set_title("(a) Global mean")
ax.tick_params(axis="x", rotation=20)

# (b) basin means
ax = axes3[1]
w = 0.2
xb = np.arange(len(basin_names))
for k, src in enumerate(["WOA23"] + models):
    ax.bar(xb + (k - 1.5) * w, [basin[(src, b)] for b in basin_names], width=w,
           label=src, color=cols[k])
ax.set_xticks(xb)
ax.set_xticklabels(basin_names)
ax.set_ylabel("Nitracline depth (m)")
ax.set_title("(b) Basin means")
ax.legend(fontsize=8, loc="upper left")

# (c) skill
ax = axes3[2]
r = [skill[m]["correlation"] for m in models]
rm = [skill[m]["rmse"] for m in models]
xp = np.arange(len(models))
ax.bar(xp - 0.2, r, width=0.4, color="#4C72B0", label="Pearson r")
ax.set_ylim(0, 1.05)
ax.set_ylabel("Pearson r")
ax.set_xticks(xp)
ax.set_xticklabels(models, rotation=15)
for i, v in enumerate(r):
    ax.text(i - 0.2, v + 0.02, f"{v:.3f}", ha="center", fontsize=9)
ax2 = ax.twinx()
ax2.bar(xp + 0.2, rm, width=0.4, color="#DD8452", label="RMSE (m)")
ax2.set_ylabel("RMSE (m)")
for i, v in enumerate(rm):
    ax2.text(i + 0.2, v * 1.02, f"{v:.1f}", ha="center", fontsize=9)
ax.set_title("(c) Spatial skill vs WOA23")
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")

plt.show()
