# ---------------------------------------------------------------------
# Figure 1: historical OMZ thickness at 60 umol/kg, one shared scale
# ---------------------------------------------------------------------
fig1, axes1 = plt.subplots(1, 3, figsize=(16, 4.4), constrained_layout=True)

pooled1 = np.concatenate([np.asarray(p[1].values, dtype=float).ravel() for p in fig1_panels])
pooled1 = pooled1[np.isfinite(pooled1)]
vmax1 = float(np.nanpercentile(pooled1, 98)) if pooled1.size else 1.0

for ax, (name, fld) in zip(axes1, fig1_panels):
    vals = np.asarray(fld.values, dtype=float)
    pcm1 = ax.pcolormesh(fld["lon"].values, fld["lat"].values, vals,
                         cmap="viridis", vmin=0.0, vmax=vmax1, shading="auto")
    ax.set_title(name, fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 90)
    ax.set_facecolor("0.85")

axes1[0].set_ylabel("Latitude")
fig1.colorbar(pcm1, ax=list(axes1), label="OMZ thickness (m)", shrink=0.85,
              extend="max")
fig1.suptitle("Figure 1  Historical OMZ thickness, O$_2$ < 60 $\\mu$mol kg$^{-1}$ "
              "(shared sequential scale, 98th-pct clip)", fontsize=11)

# ---------------------------------------------------------------------
# Figure 2: ssp245 minus historical thickness (2 models) + CESM2 ssp245
# ---------------------------------------------------------------------
fig2, axes2 = plt.subplots(1, 3, figsize=(16, 4.4), constrained_layout=True)

pooled2 = np.concatenate([np.asarray(p[1].values, dtype=float).ravel() for p in fig2_panels])
pooled2 = pooled2[np.isfinite(pooled2)]
vmax2 = float(np.nanpercentile(np.abs(pooled2), 98)) if pooled2.size else 1.0
vmax2 = max(vmax2, 1e-6)

for ax, (name, fld) in zip(axes2[:2], fig2_panels):
    vals = np.asarray(fld.values, dtype=float)
    pcm2 = ax.pcolormesh(fld["lon"].values, fld["lat"].values, vals,
                         cmap="RdBu_r", vmin=-vmax2, vmax=vmax2, shading="auto")
    ax.set_title(name, fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-80, 90)
    ax.set_facecolor("0.85")

axes2[0].set_ylabel("Latitude")
fig2.colorbar(pcm2, ax=list(axes2[:2]), label="OMZ thickness change (m)",
              shrink=0.85, extend="both")

cesm_vals = np.asarray(cesm2_thick.values, dtype=float)
finite_c = cesm_vals[np.isfinite(cesm_vals)]
vmax3 = float(np.nanpercentile(finite_c, 98)) if finite_c.size else 1.0
pcm3 = axes2[2].pcolormesh(cesm2_thick["lon"].values, cesm2_thick["lat"].values,
                           cesm_vals, cmap="viridis", vmin=0.0, vmax=vmax3,
                           shading="auto")
axes2[2].set_title("CESM2 ssp245 (no historical o$_2$:\nchange cannot be formed)", fontsize=9)
axes2[2].set_xlabel("Longitude")
axes2[2].set_xlim(-180, 180)
axes2[2].set_ylim(-80, 90)
axes2[2].set_facecolor("0.85")
fig2.colorbar(pcm3, ax=axes2[2], label="OMZ thickness (m)", shrink=0.85,
              extend="max")

fig2.suptitle("Figure 2  OMZ thickness change, ssp245 (2071-2100) minus historical "
              "(1971-2000), O$_2$ < 60 $\\mu$mol kg$^{-1}$", fontsize=11)
