# ------------------------------------------------------------------
# Figure 1: global nitracline-depth maps on the common WOA 1-degree mask
# ------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes1 = axes1.ravel()

fields_fig1 = [("WOA23", woa_nit_depth.where(common_mask))]
for m in models:
    fields_fig1.append((m, model_regrid[m].where(common_mask)))

cmap1 = plt.get_cmap("YlGnBu")
pooled1 = np.concatenate(
    [da.values.ravel() for _, da in fields_fig1]
)
vmax1 = float(np.nanpercentile(pooled1, 98))
vmin1 = 0.0

for ax, (label, da) in zip(axes1, fields_fig1):
    lon = da["lon"].values
    lat = da["lat"].values
    im = ax.pcolormesh(
        lon,
        lat,
        da.values,
        shading="auto",
        cmap=cmap1,
        vmin=vmin1,
        vmax=vmax1,
    )
    if label == "WOA23":
        mean_val = woa_global
    else:
        mean_val = model_global[label]
    ax.set_title(f"{label}  (mean {mean_val:.0f} m)")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig1.colorbar(
    axes1[0].collections[0],
    ax=axes1,
    orientation="horizontal",
    shrink=0.7,
    pad=0.03,
    label="Nitracline depth (m)",
    extend="max",
)
fig1.suptitle(
    "Global nitracline depth on the common WOA23 1-degree mask; "
    "shared sequential scale clipped at pooled 98th percentile (deeper = darker)",
    fontsize=12,
)

# ------------------------------------------------------------------
# Figure 2: equatorial Pacific nitrate section and nitracline
# ------------------------------------------------------------------
fig2, axes2 = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
axes2 = axes2.ravel()

order_fig2 = ["WOA23"] + models
cmap2 = plt.get_cmap("YlGnBu")

pooled2 = []
for label in order_fig2:
    sec = sections[label]
    sec_top = sec.sel(depth=slice(0, 1000))
    pooled2.append(sec_top.values.ravel())

vmax2 = float(np.nanpercentile(np.concatenate(pooled2), 98))
vmin2 = 0.0

for ax, label in zip(axes2, order_fig2):
    sec = sections[label]
    nit = nit_sections[label]

    xdim = [d for d in sec.dims if d != "depth"][0]
    lon_name = sec.attrs.get("lon_coord", xdim)
    if lon_name not in sec.coords:
        lon_name = xdim

    sec_top = sec.sel(depth=slice(0, 1000))
    lon_vals = sec_top[lon_name].values
    depth_vals = sec_top["depth"].values
    data2d = sec_top.transpose("depth", xdim).values

    ax.pcolormesh(
        lon_vals,
        depth_vals,
        data2d,
        shading="auto",
        cmap=cmap2,
        vmin=vmin2,
        vmax=vmax2,
    )

    if lon_name in nit.coords:
        line_lon = nit[lon_name].values
    else:
        line_lon = lon_vals
    ax.plot(line_lon, nit.values, color="red", linewidth=1.3)

    ax.set_title(label)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Depth (m)")
    ax.set_xlim(eq_lon_range[0], eq_lon_range[1])
    ax.set_ylim(1000, 0)   # depth increases downward
    ax.set_xticks([120, 150, 180, 210, 240, 270])

fig2.colorbar(
    axes2[0].collections[0],
    ax=axes2,
    orientation="horizontal",
    shrink=0.7,
    pad=0.03,
    label="Nitrate (µmol kg$^{-1}$)",
    extend="max",
)
fig2.suptitle(
    "Equatorial Pacific nitrate (2°S–2°N, 120°E–80°W); "
    "shared scale clipped at pooled 98th percentile; red line = nitracline",
    fontsize=12,
)
