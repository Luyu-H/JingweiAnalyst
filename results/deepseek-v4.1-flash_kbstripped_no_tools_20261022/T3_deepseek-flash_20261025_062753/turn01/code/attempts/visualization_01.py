# ============================================================
# Figure 1: 2x2 maps of model-minus-WOA23 SST bias and WOA23 SST
# ============================================================
bias_cesm = np.where(shared, regridded["CESM2"] - woa_sst, np.nan)
bias_gfdl = np.where(shared, regridded["GFDL-ESM4"] - woa_sst, np.nan)
bias_mpi  = np.where(shared, regridded["MPI-ESM1-2-LR"] - woa_sst, np.nan)
woa_plot  = np.where(shared, woa_sst, np.nan)

pooled = np.concatenate([
    bias_cesm[shared],
    bias_gfdl[shared],
    bias_mpi[shared],
])
vmax = float(np.nanpercentile(np.abs(pooled), 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 0.5

fig, axes = plt.subplots(
    2, 2, figsize=(14, 9), constrained_layout=True
)

bias_panels = [
    (axes[0, 0], bias_cesm, "CESM2 - WOA23"),
    (axes[0, 1], bias_gfdl, "GFDL-ESM4 - WOA23"),
    (axes[1, 0], bias_mpi,  "MPI-ESM1-2-LR - WOA23"),
]

for ax, data, title in bias_panels:
    im = ax.pcolormesh(
        woa_lon, woa_lat, data,
        cmap="RdBu_r", vmin=-vmax, vmax=vmax,
        shading="auto"
    )
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    cb = fig.colorbar(im, ax=ax, orientation="vertical", extend="both")
    cb.set_label("SST bias (°C)")

im2 = axes[1, 1].pcolormesh(
    woa_lon, woa_lat, woa_plot,
    cmap="viridis", shading="auto"
)
axes[1, 1].set_title("WOA23 SST (1971-2000, 0 m)")
axes[1, 1].set_xlabel("Longitude")
axes[1, 1].set_ylabel("Latitude")
cb2 = fig.colorbar(im2, ax=axes[1, 1], orientation="vertical")
cb2.set_label("SST (°C)")

fig.suptitle(
    "SST climatology (1971-2000): model biases and WOA23 reference "
    "on shared ocean cells"
)

# ============================================================
# Figure 2: native-grid global mean SST bar chart
# ============================================================
fig2, ax2 = plt.subplots(figsize=(8, 5), constrained_layout=True)

names = models + ["WOA23"]
vals = [native_means[m] for m in models] + [woa_gm]

bars = ax2.bar(
    names, vals,
    color=["C0", "C1", "C2", "C3"]
)
ax2.bar_label(bars, fmt="%.2f", padding=3)
ax2.set_ylabel("SST (°C)")
ax2.set_title("Native-grid global mean SST (1971-2000 climatology)")
ax2.set_ylim(0.0, max(vals) * 1.18)
