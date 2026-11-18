# -----------------------------------------------------------------------------
# Figure 1: historical OMZ thickness at the 60 umol/kg threshold
# -----------------------------------------------------------------------------
fig1, axes1 = plt.subplots(
    1, 3, figsize=(20, 5), constrained_layout=True
)
panels1 = [
    (thick_hist_plot["GFDL-ESM4"], "GFDL-ESM4 historical"),
    (thick_hist_plot["MPI-ESM1-2-LR"], "MPI-ESM1-2-LR historical"),
    (thick_woa_plot, "WOA23 1971-2000"),
]

vals1 = np.concatenate(
    [np.asarray(da.values).ravel() for da, _ in panels1]
)
vmax1 = float(np.nanpercentile(vals1, 98))
if not np.isfinite(vmax1) or vmax1 <= 0:
    vmax1 = 1.0

for ax, (da, title) in zip(axes1, panels1):
    mesh1 = da.plot(
        ax=ax,
        cmap="viridis",
        vmin=0.0,
        vmax=vmax1,
        add_colorbar=False,
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

fig1.colorbar(
    mesh1,
    ax=axes1,
    orientation="horizontal",
    shrink=0.6,
    pad=0.06,
    label="OMZ thickness (m)",
)
fig1.suptitle(
    "Historical OMZ thickness below 60 µmol kg$^{-1}$ "
    "(colour clipped at 98th percentile)"
)


# -----------------------------------------------------------------------------
# Figure 2: ssp245 minus historical OMZ thickness, plus CESM2 ssp245
# -----------------------------------------------------------------------------
fig2, axes2 = plt.subplots(
    1, 3, figsize=(20, 5), constrained_layout=True
)
change_panels = [
    (thick_change_plot["GFDL-ESM4"], "GFDL-ESM4 change"),
    (thick_change_plot["MPI-ESM1-2-LR"], "MPI-ESM1-2-LR change"),
]

vals2 = np.concatenate(
    [np.asarray(da.values).ravel() for da, _ in change_panels]
)
vmax2 = float(np.nanpercentile(np.abs(vals2), 98))
if not np.isfinite(vmax2) or vmax2 <= 0:
    vmax2 = 1.0

for ax, (da, title) in zip(axes2[:2], change_panels):
    mesh2 = da.plot(
        ax=ax,
        cmap="RdBu_r",
        vmin=-vmax2,
        vmax=vmax2,
        add_colorbar=False,
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

da_ces = thick_ssp_plot["CESM2"]
vals3 = np.asarray(da_ces.values).ravel()
vmax3 = float(np.nanpercentile(vals3, 98))
if not np.isfinite(vmax3) or vmax3 <= 0:
    vmax3 = 1.0

mesh3 = da_ces.plot(
    ax=axes2[2],
    cmap="viridis",
    vmin=0.0,
    vmax=vmax3,
    add_colorbar=False,
)
axes2[2].set_title("CESM2 ssp245 (no historical change possible)")
axes2[2].set_xlabel("Longitude")
axes2[2].set_ylabel("Latitude")
axes2[2].set_xlim(-180, 180)
axes2[2].set_ylim(-90, 90)

fig2.colorbar(
    mesh2,
    ax=axes2[:2],
    orientation="horizontal",
    shrink=0.6,
    pad=0.06,
    label="OMZ thickness change (m)",
)
fig2.colorbar(
    mesh3,
    ax=axes2[2],
    orientation="horizontal",
    shrink=0.6,
    pad=0.06,
    label="OMZ thickness (m)",
)
fig2.suptitle(
    "OMZ thickness change at 60 µmol kg$^{-1}$ "
    "(ssp245 2071-2100 minus historical 1971-2000; "
    "colour clipped at 98th percentile)"
)
