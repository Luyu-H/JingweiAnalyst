# Figure 1: global nitracline-depth maps on a shared sequential scale
fig1, axes1 = plt.subplots(2, 2, figsize=(14, 7.5), constrained_layout=True)
maps = {"WOA23": nit["WOA23"]}
for m in models:
    maps[m] = nit_reg[m]
pm = None
for ax, s in zip(axes1.ravel(), sources):
    da = maps[s]
    pm = ax.pcolormesh(da.lon.values, da.lat.values, da.values, cmap="YlGnBu",
                       vmin=0, vmax=nit_vmax, shading="auto")
    ax.set_title(f"{s} (global mean {global_mean_rounded[s]} m native)" if s != "WOA23"
                 else f"WOA23 (global mean {global_mean_rounded[s]} m)")
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
fig1.colorbar(pm, ax=axes1, shrink=0.85, extend="max", label="Nitracline depth (m)")
fig1.suptitle("Nitracline depth, 1971-2000 mean state (maps on WOA23 1° grid; scale clipped at 98th percentile)")

# Figure 2: equatorial Pacific nitrate section, 2S-2N
fig2, axes2 = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True, sharex=True, sharey=True)
pm2 = None
for ax, s in zip(axes2.ravel(), sources):
    da = sec[s]
    pm2 = ax.pcolormesh(da.lon.values, da.depth.values, da.values, cmap="viridis",
                        vmin=0, vmax=sec_vmax, shading="auto")
    ax.plot(line[s].lon.values, line[s].values, color="red", lw=2, label="Nitracline depth")
    ax.set_ylim(1000, 0)
    ax.set_xlim(120, 280)
    ax.set_title(s)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Depth (m)")
    ax.legend(loc="lower left", fontsize=8)
fig2.colorbar(pm2, ax=axes2, shrink=0.85, extend="max", label="Nitrate (mmol m$^{-3}$)")
fig2.suptitle("Equatorial Pacific (2°S-2°N) nitrate, 1971-2000; WOA23 per mass converted with constant density, models per volume; scale clipped at 99th pct")
