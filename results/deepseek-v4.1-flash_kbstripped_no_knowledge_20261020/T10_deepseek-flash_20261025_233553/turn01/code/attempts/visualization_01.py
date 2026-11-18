# ---------------------------------------------------------------------
# Figure 1 - historical OMZ thickness at 60 umol/kg (shared sequential)
# ---------------------------------------------------------------------
panels1 = [
    ("GFDL-ESM4  (historical 1971-2000)", th_reg["GFDL-ESM4 hist"]),
    ("MPI-ESM1-2-LR  (historical 1971-2000)", th_reg["MPI-ESM1-2-LR hist"]),
    ("WOA23  (observed climatology)", th_reg["WOA23"]),
]
pooled = np.concatenate([np.asarray(p[1].values).ravel() for p in panels1])
vmax1 = float(np.nanpercentile(pooled, 98))

fig1, axes1 = plt.subplots(1, 3, figsize=(19, 4.6), constrained_layout=True)
for ax, (title, da) in zip(axes1, panels1):
    pm = ax.pcolormesh(da.lon.values, da.lat.values, da.values,
                       cmap="viridis", vmin=0.0, vmax=vmax1, shading="auto")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 90)
fig1.colorbar(pm, ax=axes1, label="OMZ thickness (m)", extend="max")
fig1.suptitle("Figure 1 - Historical OMZ thickness (O$_2$ < 60 $\\mu$mol kg$^{-1}$); "
              "shared scale, upper 2% clipped", fontsize=12)

# ---------------------------------------------------------------------
# Figure 2 - ssp245 minus historical OMZ thickness (+ CESM2 ssp245)
# ---------------------------------------------------------------------
panels2 = [
    ("GFDL-ESM4  (ssp245 2071-2100 minus historical)", diff_gfdl),
    ("MPI-ESM1-2-LR  (ssp245 2071-2100 minus historical)", diff_mpi),
]
pooled2 = np.concatenate([np.asarray(p[1].values).ravel() for p in panels2])
vlim2 = float(np.nanpercentile(np.abs(pooled2), 98))

fig2, axes2 = plt.subplots(1, 3, figsize=(19, 4.6), constrained_layout=True)
for ax, (title, da) in zip(axes2[:2], panels2):
    pm2 = ax.pcolormesh(da.lon.values, da.lat.values, da.values,
                        cmap="RdBu_r", vmin=-vlim2, vmax=vlim2, shading="auto")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 90)
fig2.colorbar(pm2, ax=axes2[:2], label="$\\Delta$OMZ thickness (m)", extend="both")

cesm = th_reg["CESM2 ssp245"]
vmax3 = float(np.nanpercentile(np.asarray(cesm.values).ravel(), 98))
pm3 = axes2[2].pcolormesh(cesm.lon.values, cesm.lat.values, cesm.values,
                          cmap="viridis", vmin=0.0, vmax=vmax3, shading="auto")
axes2[2].set_title("CESM2  (ssp245 2071-2100 only - no historical O$_2$,\n"
                   "so no change can be formed)", fontsize=10)
axes2[2].set_xlabel("Longitude")
axes2[2].set_ylabel("Latitude")
axes2[2].set_xlim(-180, 180)
axes2[2].set_ylim(-80, 90)
fig2.colorbar(pm3, ax=axes2[2], label="OMZ thickness (m)", extend="max")
fig2.suptitle("Figure 2 - OMZ thickness change under ssp245; diverging scale clipped at "
              "the 98th percentile of |change|", fontsize=12)
