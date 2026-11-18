fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
vals = []
panels = {}
for n in names:
    p = reg_pattern[n].assign_coords(lon=reg_pattern[n]["lon"] % 360).sortby("lon").sortby("lat")
    p = p.sel(lat=slice(20, 70), lon=slice(110, 260))
    panels[n] = p
    vals.append(np.abs(p.values[np.isfinite(p.values)]))
vmax = float(np.nanpercentile(np.concatenate(vals), 98))

im = None
for ax, n in zip(axes.ravel(), names):
    p = panels[n]
    im = ax.pcolormesh(p["lon"], p["lat"], p.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                       shading="auto")
    txt = "EV = %.1f%%" % expl_var[n]
    if n in pattern_corr:
        txt += "\nr(pattern) = %.3f" % pattern_corr[n]
    ax.text(0.02, 0.04, txt, transform=ax.transAxes, fontsize=10,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
    ax.set_title(n)
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
cb = fig.colorbar(im, ax=axes, orientation="vertical", extend="both", shrink=0.85)
cb.set_label("°C per standard deviation")
fig.suptitle("PDO regression patterns, 1950-2014 (1° bilinear grid; colour scale clipped at 98th percentile)")

fig2, ax2 = plt.subplots(figsize=(14, 5), constrained_layout=True)
colors = {"ERSSTv5": "k", "CESM2": "tab:blue", "GFDL-ESM4": "tab:green", "MPI-ESM1-2-LR": "tab:orange"}
for n in names:
    ax2.plot(pdo_index[n]["time"].values, pdo_index[n].values, color=colors[n],
             lw=1.4 if n == "ERSSTv5" else 0.8, alpha=0.9 if n == "ERSSTv5" else 0.6, label=n)
ax2.plot(pub_std["time"].values, pub_std.values, color="tab:red", ls="--", lw=1.2,
         label="Published PDO (NCEI), standardized")
ax2.plot(smooth["time"].values, smooth.values, color="purple", lw=3, label="ERSSTv5, 9-yr running mean")
for t, y in zip(transition_times, transition_years):
    ax2.axvline(t, color="gray", ls=":", lw=1.5)
    ax2.text(t, ax2.get_ylim()[1] * 0.95, str(y), rotation=90, va="top", ha="right", fontsize=9)
ax2.axhline(0, color="k", lw=0.5)
ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index (unitless)")
ax2.set_title("PDO indices, 1950-2014 (dotted lines: ERSSTv5 phase transitions)")
ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=9)
