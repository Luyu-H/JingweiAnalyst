# ============================================================================
# Figure 1: Historical 60-threshold OMZ thickness, shared sequential scale
# ============================================================================
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

fields1 = [gfdl_hist_reg, mpi_hist_reg, woa_thick60]
titles1 = [
    "GFDL-ESM4 historical",
    "MPI-ESM1-2-LR historical",
    "WOA23 climatology",
]

vals1 = np.concatenate([np.asarray(f.values).ravel() for f in fields1])
vmax1 = float(np.nanpercentile(vals1, 98))

for ax, f, title in zip(axes1, fields1, titles1):
    pcm1 = ax.pcolormesh(
        f.lon, f.lat, f.values, shading="auto",
        cmap="viridis", vmin=0, vmax=vmax1
    )
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("auto")

fig1.colorbar(
    pcm1, ax=axes1, shrink=0.8, extend="max",
    label="OMZ thickness (m)"
)
fig1.suptitle(
    "Historical OMZ thickness at O2 < 60 µmol kg-1 "
    "(shared scale, clipped at 98th percentile)",
    fontsize=12
)


# ============================================================================
# Figure 2: ssp245-minus-historical OMZ thickness change for GFDL/MPI,
#           and CESM2 ssp245 OMZ thickness (no historical change possible)
# ============================================================================
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4), constrained_layout=True)

change_fields = [change_gfdl, change_mpi]
change_titles = [
    "GFDL-ESM4: ssp245 - historical",
    "MPI-ESM1-2-LR: ssp245 - historical",
]

vals2 = np.concatenate([np.asarray(f.values).ravel() for f in change_fields])
vmax2 = float(np.nanpercentile(np.abs(vals2), 98))

for ax, f, title in zip(axes2[:2], change_fields, change_titles):
    pcm2 = ax.pcolormesh(
        f.lon, f.lat, f.values, shading="auto",
        cmap="RdBu_r", vmin=-vmax2, vmax=vmax2
    )
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_aspect("auto")

fig2.colorbar(
    pcm2, ax=axes2[:2], shrink=0.8, extend="both",
    label="Δ OMZ thickness (m)"
)

f_cesm = cesm_ssp_reg
vmax3 = float(np.nanpercentile(np.asarray(f_cesm.values).ravel(), 98))
pcm3 = axes2[2].pcolormesh(
    f_cesm.lon, f_cesm.lat, f_cesm.values, shading="auto",
    cmap="viridis", vmin=0, vmax=vmax3
)
axes2[2].set_title(
    "CESM2 ssp245 OMZ thickness (m)\n"
    "no change formed: historical r10 unavailable",
    fontsize=9
)
axes2[2].set_xlim(-180, 180)
axes2[2].set_ylim(-90, 90)
axes2[2].set_xlabel("Longitude")
axes2[2].set_ylabel("Latitude")
axes2[2].set_aspect("auto")

fig2.colorbar(
    pcm3, ax=axes2[2], shrink=0.8, extend="max",
    label="OMZ thickness (m)"
)
fig2.suptitle(
    "OMZ thickness at O2 < 60 µmol kg-1: scenario changes and CESM2 ssp245",
    fontsize=12
)
