# ------------------------------------------------------------------
# Figure 1: global nitracline-depth maps on the common WOA 1-degree mask
# ------------------------------------------------------------------
fig1, axes1 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes1 = axes1.ravel()

fields_fig1 = [("WOA23", woa_nit_depth.where(common_mask))]
for m in models:
    fields_fig1.append((m, model_regrid[m].where(common_mask)))

cmap1 = plt.get_cmap("YlGnBu")


def _to_float(a):
    """Plain float array with masked / missing entries turned into NaN."""
    return np.ma.filled(np.ma.asarray(a, dtype=float), np.nan)


pooled1 = np.concatenate(
    [_to_float(da.values).ravel() for _, da in fields_fig1]
)
vmax1 = float(np.nanpercentile(pooled1, 98))
vmin1 = 0.0

for ax, (label, da) in zip(axes1, fields_fig1):
    lon = _to_float(da["lon"].values)
    lat = _to_float(da["lat"].values)
    data = _to_float(da.values)

    # keep only the finite coordinate rows / columns (masked coordinates would
    # make pcolormesh raise)
    ok_x = np.isfinite(lon)
    ok_y = np.isfinite(lat)

    im = ax.pcolormesh(
        lon[ok_x],
        lat[ok_y],
        data[np.ix_(ok_y, ok_x)],
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
    pooled2.append(_to_float(sec_top.values).ravel())

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
    depth_vals = _to_float(sec_top["depth"].values)
    data2d = _to_float(sec_top.transpose("depth", xdim).values)

    # The section longitude coordinate can carry NaN where the 2S-2N mask
    # blanked a whole column; pcolormesh cannot take non-finite x, so drop
    # those columns together with their data.
    lon_vals = _to_float(sec_top[lon_name].values)
    ok_x = np.isfinite(lon_vals)
    ok_z = np.isfinite(depth_vals)

    if ok_x.any() and ok_z.any():
        ax.pcolormesh(
            lon_vals[ok_x],
            depth_vals[ok_z],
            data2d[np.ix_(ok_z, ok_x)],
            shading="auto",
            cmap=cmap2,
            vmin=vmin2,
            vmax=vmax2,
        )

    # nitracline overlain as a line, only where longitude and depth are finite
    nit_vals = _to_float(nit.values)
    if lon_name in nit.coords:
        line_lon = _to_float(nit[lon_name].values)
    else:
        line_lon = _to_float(sec_top[lon_name].values)
    ok_line = np.isfinite(line_lon) & np.isfinite(nit_vals)
    if ok_line.any():
        ax.plot(line_lon[ok_line], nit_vals[ok_line], color="red", linewidth=1.3)

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
