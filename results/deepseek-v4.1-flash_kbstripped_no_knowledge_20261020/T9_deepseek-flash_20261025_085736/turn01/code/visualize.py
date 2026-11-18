# Figure 1: PDO regression patterns in a 2x2 panel.
names = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

def _continuous_lon(da):
    if "lon" in da.coords:
        return da.assign_coords(lon=(da.lon % 360)).sortby("lon")
    return da

fig1, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
axes = axes.ravel()

pooled_abs = []
for name in names:
    r = _continuous_lon(results[name]["reg"])
    v = r.values
    pooled_abs.append(np.abs(v[np.isfinite(v)]).ravel())

vmax = float(np.nanpercentile(np.concatenate(pooled_abs), 98))
if not np.isfinite(vmax) or vmax == 0:
    vmax = 1.0

im = None
for ax, name in zip(axes, names):
    r = _continuous_lon(results[name]["reg"])
    im = ax.pcolormesh(
        r.lon, r.lat, r.values,
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="auto",
    )
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    title = f"{name}\nEV = {results[name]['var_expl']:.1f}%"
    if name != "ERSSTv5":
        title += f", r = {pattern_corr[name]:.3f}"
    ax.set_title(title, fontsize=10)

cb = fig1.colorbar(im, ax=axes, shrink=0.8, extend="both")
cb.set_label("SST reg. (°C per SD)")
fig1.suptitle(
    "PDO regression patterns, 1950–2014 (colour clipped at pooled 98th percentile of |value|)",
    fontsize=12,
)

# Figure 2: standardized indices, smoothed ERSSTv5 curve, and transition years.
fig2, ax2 = plt.subplots(figsize=(13, 5), constrained_layout=True)

colors = {
    "ERSSTv5": "tab:blue",
    "CESM2": "tab:orange",
    "GFDL-ESM4": "tab:green",
    "MPI-ESM1-2-LR": "tab:red",
}

for name in names:
    idx = results[name]["index"]
    ax2.plot(
        idx.time, idx.values,
        label=name, color=colors[name], lw=0.7, alpha=0.65,
    )

ax2.plot(
    pdo_pub_std.time, pdo_pub_std.values,
    label="CPC PDO (standardized)",
    color="black", lw=1.0, ls="--", alpha=0.9,
)
ax2.plot(
    smoothed_c.time, smoothed_c.values,
    label="ERSSTv5 9-yr running mean",
    color="black", lw=2.5,
)

ax2.axhline(0, color="k", lw=0.5)
for yr in transition_years:
    t = pd.Timestamp(f"{yr}-07-01")
    ax2.axvline(t, color="gray", ls=":", lw=1, alpha=0.8)
    ymin, ymax = ax2.get_ylim()
    ax2.text(
        t, ymax * 0.95, str(yr),
        rotation=90, va="top", ha="right", fontsize=8, color="gray",
    )

ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index")
ax2.set_title("PDO indices, 1950–2014")
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2015-01-01"))
ax2.legend(
    loc="upper left",
    bbox_to_anchor=(1.01, 1.0),
    frameon=False,
    fontsize=8,
)
