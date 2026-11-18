import matplotlib as mpl

# ===========================================================================
# Figure 1 - regression patterns over the PDO domain
# ===========================================================================
def _to360(da):
    """Continuous 0-360 longitude axis so the dateline-crossing domain is one strip."""
    return da.assign_coords(lon=(da["lon"] % 360)).sortby("lon")

panels = [_to360(regression_pattern[k]) for k in order]

pooled = np.concatenate([np.abs(np.asarray(p.values, dtype=float)).ravel()
                         for p in panels])
vmax = float(np.nanpercentile(pooled, 98))
levels = np.linspace(-vmax, vmax, 21)

cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad("0.85")

fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), constrained_layout=True)
axes = axes.ravel()

for ax, k, da in zip(axes, order, panels):
    lon = np.asarray(da["lon"].values)
    lat = np.asarray(da["lat"].values)
    mesh = ax.pcolormesh(lon, lat, np.asarray(da.values, dtype=float),
                         cmap=cmap, vmin=-vmax, vmax=vmax,
                         shading="auto")
    ax.contour(lon, lat, np.asarray(da.values, dtype=float),
               levels=[0.0], colors="k", linewidths=0.6)
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)
    ax.set_xticks([110, 140, 170, 200, 230, 260])
    ax.set_xticklabels(["110E", "140E", "170E", "160W", "130W", "100W"])
    ax.set_yticks([20, 30, 40, 50, 60, 70])
    ax.set_ylabel("Latitude (deg N)")
    ax.set_xlabel("Longitude")

    if k == "ERSSTv5":
        ax.set_title(f"{k}  (var. expl. = {explained_variance[k]:.1f} %)\n"
                     "reference pattern", fontsize=10)
    else:
        ax.set_title(f"{k}  (var. expl. = {explained_variance[k]:.1f} %,\n"
                     f"pattern r vs ERSSTv5 = {pattern_corr[k]:.3f})", fontsize=10)

cbar = fig.colorbar(mesh, ax=axes.tolist(), orientation="vertical",
                    extend="both", shrink=0.9, pad=0.02)
cbar.set_label("SST regression (\u00b0C per SD)", fontsize=9)

fig.suptitle("Figure 1 - PDO regression pattern on the standardized index, "
             "20\u00b0N-70\u00b0N, 110\u00b0E-100\u00b0W, 1950-2014\n"
             f"shared diverging scale clipped at the 98th percentile of |pattern| "
             f"({vmax:.2f} \u00b0C/SD); black contour = zero", fontsize=11)

# ===========================================================================
# Figure 2 - standardized PDO indices + 9-year running mean + transitions
# ===========================================================================
fig2, ax2 = plt.subplots(figsize=(13.5, 6.0), constrained_layout=True)

colors = {"ERSSTv5": "#111111",
          "CESM2": "#d62728",
          "GFDL-ESM4": "#1f77b4",
          "MPI-ESM1-2-LR": "#2ca02c"}

for k in order:
    ax2.plot(np.asarray(pdo_index[k]["time"].values),
             np.asarray(pdo_index[k].values, dtype=float),
             color=colors[k], lw=0.6, alpha=0.75,
             label=f"{k}" + (" (ERSST computed)" if k == "ERSSTv5" else " (model)"))

ax2.plot(np.asarray(pdo_published["time"].values),
         np.asarray(pdo_published.values, dtype=float),
         color="#9467bd", lw=0.9, ls="--", alpha=0.9,
         label="Published PDO (NOAA NCEI)")

ax2.plot(np.asarray(smoothed_ersst["time"].values), sv,
         color="#ff7f0e", lw=2.6, label="ERSSTv5 108-month running mean")

ymin, ymax = -4.0, 4.0
for y in transitions:
    xt = pd.Timestamp(f"{y}-07-01")
    ax2.axvline(xt, color="grey", ls=":", lw=1.0)
    ax2.text(xt, ymax * 0.94, str(y), rotation=90, va="top", ha="right",
             fontsize=8, color="dimgrey")

ax2.set_ylim(ymin, ymax)
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2014-12-31"))
ax2.axhline(0.0, color="grey", lw=0.8)
ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized PDO index")
ax2.set_title("Figure 2 - Standardized PDO index, 1950-2014: models vs ERSSTv5 "
              "and the published index\n"
              "dotted vertical lines = sign changes of the 108-month centred "
              "ERSSTv5 running mean (crossings < 3 yr after the previous one merged)",
              fontsize=11)
ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3,
           fontsize=8.5, frameon=False)
