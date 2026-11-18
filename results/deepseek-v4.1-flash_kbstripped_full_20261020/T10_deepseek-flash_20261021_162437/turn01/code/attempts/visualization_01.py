# ---------------------------------------------------------------------------
# Figure 1: historical OMZ thickness at < 60 umol/kg
# ---------------------------------------------------------------------------
fig1_panels = [
    ("GFDL-ESM4 historical", fig1_fields["GFDL-ESM4 historical"]),
    ("MPI-ESM1-2-LR historical", fig1_fields["MPI-ESM1-2-LR historical"]),
    ("WOA23 1971-2000", fig1_fields["WOA23 1971-2000"]),
]

fig1_data = [np.asarray(da.compute()) for _, da in fig1_panels]
pooled_fig1 = np.concatenate([d[np.isfinite(d)].ravel() for d in fig1_data])
vmax_fig1 = float(np.nanpercentile(pooled_fig1, 98))
vmax_fig1 = max(vmax_fig1, 1e-6)

fig1, axes1 = plt.subplots(1, 3, figsize=(16, 4.2), constrained_layout=True)
for ax, (title, da) in zip(axes1, fig1_panels):
    da_c = da.compute()
    lon = da_c["lon"].values
    lat = da_c["lat"].values
    data = da_c.values
    pcm = ax.pcolormesh(
        lon, lat, data, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax_fig1
    )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cbar1 = fig1.colorbar(pcm, ax=axes1, shrink=0.86, extend="max")
cbar1.set_label("OMZ thickness (m)")
fig1.suptitle(
    "OMZ thickness for O2 < 60 µmol kg$^{-1}$; shared sequential scale, "
    "clipped at pooled 98th percentile",
    fontsize=12,
)

# ---------------------------------------------------------------------------
# Figure 2: ssp245 minus historical thickness, plus CESM2 ssp245 state
# ---------------------------------------------------------------------------
gfdl_change = fig2_change_fields["GFDL-ESM4"].compute()
mpi_change = fig2_change_fields["MPI-ESM1-2-LR"].compute()
cesm2_state = cesm2_ssp_display.compute()

fig2_change_data = [
    np.asarray(gfdl_change),
    np.asarray(mpi_change),
]
pooled_change = np.concatenate(
    [d[np.isfinite(d)].ravel() for d in fig2_change_data]
)
vmax_change = float(np.nanpercentile(np.abs(pooled_change), 98))
vmax_change = max(vmax_change, 1e-6)

cesm2_data = np.asarray(cesm2_state)
pooled_cesm2 = cesm2_data[np.isfinite(cesm2_data)].ravel()
vmax_cesm2 = float(np.nanpercentile(pooled_cesm2, 98))
vmax_cesm2 = max(vmax_cesm2, 1e-6)

fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.2), constrained_layout=True)

# Panel 1 and 2: change fields, shared diverging scale
for ax, (title, da) in zip(
    axes2[:2],
    [
        ("GFDL-ESM4: ssp245 - historical", gfdl_change),
        ("MPI-ESM1-2-LR: ssp245 - historical", mpi_change),
    ],
):
    lon = da["lon"].values
    lat = da["lat"].values
    pcm2 = ax.pcolormesh(
        lon, lat, da.values, shading="auto", cmap="RdBu_r",
        vmin=-vmax_change, vmax=vmax_change
    )
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)

cbar2 = fig2.colorbar(
    pcm2, ax=axes2[:2], shrink=0.86, extend="both"
)
cbar2.set_label("Change in OMZ thickness (m)")

# Panel 3: CESM2 ssp245 state, own sequential scale
ax = axes2[2]
lon = cesm2_state["lon"].values
lat = cesm2_state["lat"].values
pcm3 = ax.pcolormesh(
    lon, lat, cesm2_state.values, shading="auto",
    cmap="viridis", vmin=0.0, vmax=vmax_cesm2
)
ax.set_title("CESM2 ssp245 state (no historical baseline)", fontsize=10)
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)

cbar3 = fig2.colorbar(pcm3, ax=axes2[2], shrink=0.86, extend="max")
cbar3.set_label("OMZ thickness (m)")

fig2.suptitle(
    "OMZ thickness change at O2 < 60 µmol kg$^{-1}$; diverging scale centred on zero, "
    "clipped at pooled 98th percentile",
    fontsize=12,
)
