# =============================================================================
# Figure 1 - regression patterns over the PDO domain (2 x 2)
# =============================================================================
order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

def pdo_box(da):
    p = da.sortby("lat").sel(lat=slice(20, 70))
    p = p.assign_coords(lon=p["lon"] % 360).sortby("lon").sel(lon=slice(110, 260))
    return p

boxed = {k: pdo_box(reg_s[k]) for k in order}
pool = np.concatenate([np.asarray(boxed[k].values, dtype=float).ravel() for k in order])
vmax = float(np.nanpercentile(np.abs(pool[np.isfinite(pool)]), 98))
levels = np.linspace(-vmax, vmax, 21)

fig1, axes = plt.subplots(2, 2, figsize=(14.5, 8.5), constrained_layout=True)
mappable = None
for ax, k in zip(axes.ravel(), order):
    b = boxed[k]
    mappable = ax.contourf(b["lon"], b["lat"], b.values, levels=levels,
                           cmap="RdBu_r", extend="both")
    if k == "ERSSTv5":
        ttl = "ERSSTv5 (obs):  var. expl. %.1f %%" % ev_s[k]
    else:
        ttl = "%s:  var. expl. %.1f %%,  r vs obs = %.3f" % (k, ev_s[k], pcorr[k])
    ax.set_title(ttl, fontsize=10)
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.grid(alpha=0.25, lw=0.4)

fig1.colorbar(mappable, ax=axes, label="SST regression (°C per SD)", shrink=0.85)
fig1.suptitle("Figure 1 — PDO regression pattern, 1950–2014 "
              "(colours clipped at the 98th percentile of |regression|, shared scale)",
              fontsize=12)

# =============================================================================
# Figure 2 - standardized indices, smoothed ERSSTv5 curve and transition years
# =============================================================================
fig2, ax2 = plt.subplots(figsize=(13.5, 5.8), constrained_layout=True)
colors = {"ERSSTv5": "black", "CESM2": "tab:red",
          "GFDL-ESM4": "tab:blue", "MPI-ESM1-2-LR": "tab:green"}

for k in order:
    ax2.plot(idx_plot[k]["time"].values, idx_plot[k].values,
             lw=0.7, alpha=0.50, color=colors[k], label=k)
ax2.plot(pdo_pub_c["time"].values, pdo_pub_c.values,
         lw=1.1, ls="--", color="tab:purple", label="Published PDO (NOAA NCEI)")
ax2.plot(sm_da["time"].values, sm_da.values,
         lw=2.6, color="black", label="ERSSTv5 9-yr centred mean")

ylim = ax2.get_ylim()
for y in transitions:
    xt = pd.Timestamp("%d-07-01" % int(y))
    ax2.axvline(xt, color="0.35", ls=":", lw=1.3)
    ax2.text(xt, ylim[1] * 0.86, str(int(y)), fontsize=8, ha="center",
             color="0.25", rotation=90)

ax2.axhline(0, color="0.6", lw=0.8)
ax2.set_xlabel("Year")
ax2.set_ylabel("Standardized index")
ax2.set_title("Figure 2 — Standardized PDO indices 1950–2014 "
              "(dotted verticals: zero crossings of the 9-yr ERSSTv5 curve, "
              "< 3-yr crossings merged)", fontsize=11)
ax2.legend(loc="upper left", ncol=3, fontsize=8, framealpha=0.85)
ax2.grid(alpha=0.25, lw=0.4)
