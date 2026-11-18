# ----------------------------------------------------------------------
# Figure 1: global nitracline depth maps, shared sequential scale
# ----------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

map_items = [("WOA23", woa_nc.where(woa_landsea > 1))]
for m in models:
    map_items.append((m, model_nc_regrid[m].where(woa_landsea > 1)))

pool1_list = []
for _, da in map_items:
    arr = np.asarray(da.values)
    pool1_list.append(arr[np.isfinite(arr)])
pool1 = np.concatenate(pool1_list)

vmin1 = 0.0
vmax1 = float(np.nanpercentile(pool1, 98))
cmap1 = plt.get_cmap("YlGnBu")

p1_last = None
for ax, (title, da) in zip(axes1.flat, map_items):
    lon = da.lon.values
    lat = da.lat.values
    p1_last = ax.pcolormesh(
        lon,
        lat,
        da.values,
        cmap=cmap1,
        vmin=vmin1,
        vmax=vmax1,
        shading="auto",
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig1.colorbar(p1_last, ax=axes1, label="Nitracline depth (m)", shrink=0.9)
fig1.suptitle(
    "Nitracline depth (top 1000 m search; shared scale capped at 98th percentile)",
    fontsize=14,
)


# ----------------------------------------------------------------------
# Figure 2: equatorial Pacific nitrate sections and nitracline overlay
# ----------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, 4, figsize=(18, 4.8), constrained_layout=True)

sec_items = [("WOA23", woa_sec, woa_sec_nc)]
for m in models:
    sec_items.append((m, model_sec[m], model_sec_nc[m]))

pool2_list = []
for _, sec, _ in sec_items:
    sec_top = sec.where(sec.depth <= 1000.0, drop=True)
    arr = np.asarray(sec_top.values)
    pool2_list.append(arr[np.isfinite(arr)])
pool2 = np.concatenate(pool2_list)

vmin2 = 0.0
vmax2 = float(np.nanpercentile(pool2, 98))
cmap2 = plt.get_cmap("viridis")

p2_last = None
for ax, (title, sec, sec_nc) in zip(axes2, sec_items):
    sec_top = sec.where(sec.depth <= 1000.0, drop=True)
    lon = sec_top.lon.values
    dep = sec_top.depth.values
    p2_last = ax.pcolormesh(
        lon,
        dep,
        sec_top.values,
        cmap=cmap2,
        vmin=vmin2,
        vmax=vmax2,
        shading="auto",
    )
    ax.plot(
        sec_nc.lon.values,
        sec_nc.values,
        color="k",
        linewidth=1.4,
        label="Nitracline",
    )
    ax.set_ylim(1000, 0)
    ax.set_xlim(120, 280)
    ax.set_xticks([120, 160, 200, 240, 280])
    ax.set_xticklabels(["120°E", "160°E", "160°W", "120°W", "80°W"])
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Depth (m)")
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)

fig2.colorbar(
    p2_last,
    ax=axes2,
    label="Nitrate (µmol kg$^{-1}$)",
    shrink=0.9,
)
fig2.suptitle(
    "Equatorial Pacific (2°S–2°N) nitrate sections; shared scale capped at 98th percentile",
    fontsize=14,
)
