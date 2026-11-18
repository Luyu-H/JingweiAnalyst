# -----------------------------------------------------------------------------
# Figure 1: historical OMZ thickness at the 60 umol/kg threshold
# -----------------------------------------------------------------------------
gfdl_h = model_results["GFDL-ESM4_historical"]
mpi_h = model_results["MPI-ESM1-2-LR_historical"]

lat_gfdl = gfdl_h["lat2d"]
lon_gfdl = np.mod(gfdl_h["lon2d"], 360.0)
thick_gfdl = gfdl_h["thick_map_60"]

lat_mpi = mpi_h["lat2d"]
lon_mpi = np.mod(mpi_h["lon2d"], 360.0)
thick_mpi = mpi_h["thick_map_60"]

# WOA needs its longitude reordered to 0-360 for pcolormesh
woa_sort = np.argsort(np.mod(woa_lon, 360.0))
woa_lon_plot = np.mod(woa_lon, 360.0)[woa_sort]
woa_lat_plot = woa_lat
woa_thick_plot = metrics_woa["thick_map_60"][:, woa_sort]
lon2d_woa_plot, lat2d_woa_plot = np.meshgrid(woa_lon_plot, woa_lat_plot)

fig, axes = plt.subplots(1, 3, figsize=(18, 4.5))
maps1 = [
    (lat_gfdl, lon_gfdl, thick_gfdl, "GFDL-ESM4 historical\n1971-2000"),
    (lat_mpi, lon_mpi, thick_mpi, "MPI-ESM1-2-LR historical\n1971-2000"),
    (lat2d_woa_plot, lon2d_woa_plot, woa_thick_plot, "WOA23\n1971-2000"),
]
vmax1 = max(np.nanmax(m[2]) for m in maps1)
for ax, (lat, lon, data, title) in zip(axes, maps1):
    pcm = ax.pcolormesh(
        lon, lat, data,
        shading="auto", cmap="viridis", vmin=0, vmax=vmax1
    )
    ax.set_title(title)
    ax.set_xlim(0, 360)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
fig.colorbar(pcm, ax=axes, shrink=0.8, label="OMZ thickness < 60 umol kg$^{-1}$ (m)")

# -----------------------------------------------------------------------------
# Figure 2: ssp245 minus historical OMZ thickness, plus CESM2 ssp245 only
# -----------------------------------------------------------------------------
gfdl_s = model_results["GFDL-ESM4_ssp245"]
mpi_s = model_results["MPI-ESM1-2-LR_ssp245"]
cesm_s = model_results["CESM2_ssp245"]

change_gfdl = gfdl_s["thick_map_60"] - gfdl_h["thick_map_60"]
change_mpi = mpi_s["thick_map_60"] - mpi_h["thick_map_60"]
vmax2 = max(
    np.nanmax(np.abs(change_gfdl)),
    np.nanmax(np.abs(change_mpi)),
)

fig, axes = plt.subplots(1, 3, figsize=(18, 4.5))

ax = axes[0]
pcm0 = ax.pcolormesh(
    lon_gfdl, lat_gfdl, change_gfdl,
    shading="auto", cmap="RdBu_r", vmin=-vmax2, vmax=vmax2
)
ax.set_title("GFDL-ESM4\nssp245 2071-2100 minus historical")
ax.set_xlim(0, 360)
ax.set_ylim(-90, 90)

ax = axes[1]
pcm1 = ax.pcolormesh(
    np.mod(mpi_s["lon2d"], 360.0), mpi_s["lat2d"], change_mpi,
    shading="auto", cmap="RdBu_r", vmin=-vmax2, vmax=vmax2
)
ax.set_title("MPI-ESM1-2-LR\nssp245 2071-2100 minus historical")
ax.set_xlim(0, 360)
ax.set_ylim(-90, 90)

ax = axes[2]
cesm_thick = cesm_s["thick_map_60"]
lon_cesm = np.mod(cesm_s["lon2d"], 360.0)
lat_cesm = cesm_s["lat2d"]
pcm2 = ax.pcolormesh(
    lon_cesm, lat_cesm, cesm_thick,
    shading="auto", cmap="viridis", vmin=0, vmax=np.nanmax(cesm_thick)
)
ax.set_title("CESM2 ssp245 2071-2100\n(no historical o2; change not formable)")
ax.set_xlim(0, 360)
ax.set_ylim(-90, 90)

for ax in axes:
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

fig.colorbar(pcm0, ax=axes[:2], shrink=0.8, label="Change in OMZ thickness (m)")
fig.colorbar(pcm2, ax=axes[2], shrink=0.8, label="OMZ thickness (m)")
