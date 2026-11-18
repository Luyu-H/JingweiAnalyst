# Figure 1: historical OMZ thickness at 60 umol kg-1 on shared sequential scale
fig1, axes1 = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)

hist_stack = np.concatenate([
    np.ravel(gfdl_hist_omz60_thick_plot.values[np.isfinite(gfdl_hist_omz60_thick_plot.values)]),
    np.ravel(mpi_hist_omz60_thick_plot.values[np.isfinite(mpi_hist_omz60_thick_plot.values)]),
    np.ravel(woa_omz60_thickness_plot.values[np.isfinite(woa_omz60_thickness_plot.values)]),
])
vmax_hist = float(np.nanpercentile(hist_stack, 98))

panels1 = [
    (gfdl_hist_omz60_thick_plot, "GFDL-ESM4 historical 1971-2000"),
    (mpi_hist_omz60_thick_plot, "MPI-ESM1-2-LR historical 1971-2000"),
    (woa_omz60_thickness_plot, "WOA23 1971-2000"),
]

mappable1 = None
for ax, (da, title) in zip(axes1, panels1):
    mappable1 = ax.pcolormesh(
        da["lon"], da["lat"], da,
        cmap="viridis", vmin=0, vmax=vmax_hist, shading="auto"
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(float(da.lon.min()), float(da.lon.max()))
    ax.set_ylim(float(da.lat.min()), float(da.lat.max()))

cbar1 = fig1.colorbar(mappable1, ax=axes1, orientation="vertical", extend="max")
cbar1.set_label("OMZ thickness (m)")
fig1.suptitle("Figure 1. OMZ thickness at O$_2$ < 60 µmol kg$^{-1}$; shared scale clipped at pooled 98th percentile")

# Figure 2: scenario change for GFDL/MPI and CESM2 ssp245 thickness
fig2, axes2 = plt.subplots(1, 3, figsize=(17, 4.8), constrained_layout=True)

change_stack = np.concatenate([
    np.ravel(gfdl_change_omz60_plot.values[np.isfinite(gfdl_change_omz60_plot.values)]),
    np.ravel(mpi_change_omz60_plot.values[np.isfinite(mpi_change_omz60_plot.values)]),
])
vmax_change = float(np.nanpercentile(np.abs(change_stack), 98))

cesm_stack = np.ravel(cesm2_ssp245_omz60_thick_plot.values[np.isfinite(cesm2_ssp245_omz60_thick_plot.values)])
vmax_cesm = float(np.nanpercentile(cesm_stack, 98))

m0 = axes2[0].pcolormesh(
    gfdl_change_omz60_plot["lon"], gfdl_change_omz60_plot["lat"], gfdl_change_omz60_plot,
    cmap="RdBu_r", vmin=-vmax_change, vmax=vmax_change, shading="auto"
)
axes2[0].set_title("GFDL-ESM4: ssp245 minus historical")
axes2[0].set_xlabel("Longitude")
axes2[0].set_ylabel("Latitude")

m1 = axes2[1].pcolormesh(
    mpi_change_omz60_plot["lon"], mpi_change_omz60_plot["lat"], mpi_change_omz60_plot,
    cmap="RdBu_r", vmin=-vmax_change, vmax=vmax_change, shading="auto"
)
axes2[1].set_title("MPI-ESM1-2-LR: ssp245 minus historical")
axes2[1].set_xlabel("Longitude")
axes2[1].set_ylabel("Latitude")

m2 = axes2[2].pcolormesh(
    cesm2_ssp245_omz60_thick_plot["lon"], cesm2_ssp245_omz60_thick_plot["lat"], cesm2_ssp245_omz60_thick_plot,
    cmap="viridis", vmin=0, vmax=vmax_cesm, shading="auto"
)
axes2[2].set_title("CESM2 ssp245 thickness\n(no historical o2, so no change field)")
axes2[2].set_xlabel("Longitude")
axes2[2].set_ylabel("Latitude")

cbar2a = fig2.colorbar(m1, ax=axes2[:2], orientation="vertical", extend="both")
cbar2a.set_label("Thickness change (m)")
cbar2b = fig2.colorbar(m2, ax=[axes2[2]], orientation="vertical", extend="max")
cbar2b.set_label("OMZ thickness (m)")
fig2.suptitle("Figure 2. OMZ thickness change at O$_2$ < 60 µmol kg$^{-1}$; change scale clipped at pooled 98th percentile")
