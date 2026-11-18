# ---------------------------------------------------------------------------
# Figure 1: global nitracline-depth maps on one shared sequential scale
# ---------------------------------------------------------------------------
plot_fields = {
    "CESM2": nc_model_reg["CESM2"],
    "GFDL-ESM4": nc_model_reg["GFDL-ESM4"],
    "MPI-ESM1-2-LR": nc_model_reg["MPI-ESM1-2-LR"],
    "WOA23": woa_nc_reg,
}

pooled_nc = []
for da in plot_fields.values():
    vals = np.asarray(da.values).ravel()
    pooled_nc.append(vals[np.isfinite(vals)])
pooled_nc = np.concatenate(pooled_nc)

vmin_nc = 0.0
vmax_nc = float(np.nanpercentile(pooled_nc, 98))

fig1, axes1 = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
axes1 = axes1.ravel()

for ax, (name, da) in zip(axes1, plot_fields.items()):
    lat_name = da.attrs.get("lat_coord", "lat")
    lon_name = da.attrs.get("lon_coord", "lon")
    if lat_name not in da.coords:
        lat_name = "lat"
    if lon_name not in da.coords:
        lon_name = "lon"

    lon = da[lon_name].values
    lat = da[lat_name].values
    mesh1 = ax.pcolormesh(
        lon, lat, da.values,
        shading="auto",
        cmap="viridis_r",
        vmin=vmin_nc,
        vmax=vmax_nc,
    )
    ax.set_title(name)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_aspect("equal", adjustable="box")

cbar1 = fig1.colorbar(
    mesh1, ax=axes1, orientation="horizontal",
    fraction=0.05, pad=0.05, extend="max"
)
cbar1.set_label("Nitracline depth (m)")
fig1.suptitle(
    "Nitracline depth, 1971-2000 mean (shared scale, 98th percentile)",
    fontsize=14,
)

# ---------------------------------------------------------------------------
# Figure 2: equatorial Pacific longitude-depth nitrate sections (2S-2N)
# ---------------------------------------------------------------------------
section_order = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

prepped_sections = {}
for name in section_order:
    conc = section_fields[name]
    lon = section_lon[name]

    lon_vals = np.asarray(lon.values)
    if lon_vals.ndim > 1:
        lon_vals = lon_vals.mean(axis=0)

    lon_cont = lon_vals.copy()
    lon_cont[lon_cont < 0] += 360.0

    keep = (lon_cont >= 120) & (lon_cont <= 280)
    lon_cont = lon_cont[keep]

    order = np.argsort(lon_cont)
    lon_cont = lon_cont[order]

    conc_vals = np.asarray(conc.values)
    dims = conc.dims
    lon_dim = None
    for d in dims:
        if d != "depth":
            lon_dim = d
    if lon_dim is None:
        lon_dim = dims[-1]
    if dims[-1] != lon_dim:
        conc_vals = np.moveaxis(conc_vals, dims.index(lon_dim), -1)

    conc_vals = conc_vals[..., keep][..., order]

    depth_vals = conc["depth"].values if "depth" in conc.coords else np.arange(conc.shape[0])

    nc_vals = np.asarray(section_nc[name].values).ravel()
    nc_vals = nc_vals[keep][order]

    prepped_sections[name] = (lon_cont, depth_vals, conc_vals, nc_vals)

pooled_no3 = []
for name in section_order:
    lon_cont, depth_vals, conc_vals, nc_vals = prepped_sections[name]
    dmask = depth_vals <= 1000
    vals = conc_vals[dmask, :].ravel()
    pooled_no3.append(vals[np.isfinite(vals)])
pooled_no3 = np.concatenate(pooled_no3)

vmin_no3 = float(np.nanpercentile(pooled_no3, 2))
vmax_no3 = float(np.nanpercentile(pooled_no3, 98))

fig2, axes2 = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)
axes2 = axes2.ravel()

for ax, name in zip(axes2, section_order):
    lon_cont, depth_vals, conc_vals, nc_vals = prepped_sections[name]
    dmask = depth_vals <= 1000

    mesh2 = ax.pcolormesh(
        lon_cont,
        depth_vals[dmask],
        conc_vals[dmask, :],
        shading="auto",
        cmap="viridis",
        vmin=vmin_no3,
        vmax=vmax_no3,
    )

    good = np.isfinite(nc_vals) & (nc_vals >= 0) & (nc_vals <= 1000)
    ax.plot(lon_cont[good], nc_vals[good], color="k", linewidth=2.2, label="Nitracline")

    ax.set_xlim(120, 280)
    ax.set_ylim(1000, 0)
    ax.set_title(name)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Depth (m)")
    ax.legend(loc="lower right")

cbar2 = fig2.colorbar(
    mesh2, ax=axes2, orientation="horizontal",
    fraction=0.05, pad=0.05, extend="both"
)
cbar2.set_label("Nitrate (µmol kg$^{-1}$)")
fig2.suptitle(
    "Equatorial Pacific (2°S-2°N) nitrate and nitracline, shared scale",
    fontsize=14,
)
