# Figure 1: reconstruction, CPC, running means, thresholds, event shading, and difference
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True, constrained_layout=True)

ax1.plot(recon_9120.time, recon_9120.values, color="tab:blue", lw=0.7, alpha=0.45,
         label="RECON 1991-2020 monthly")
ax1.plot(cpc_9120.time, cpc_9120.values, color="tab:red", lw=0.7, alpha=0.45,
         label="CPC 1991-2020 monthly")
ax1.plot(recon_9120_sm.time, recon_9120_sm.values, color="tab:blue", lw=2.0,
         label="RECON 3-month running mean")
ax1.plot(cpc_9120_sm.time, cpc_9120_sm.values, color="tab:red", lw=2.0,
         label="CPC 3-month running mean")
ax1.axhline(0.5, color="darkorange", ls="--", lw=1.0, label="±0.5 °C thresholds")
ax1.axhline(-0.5, color="darkorange", ls="--", lw=1.0)
ax1.axhline(0, color="k", lw=0.5, alpha=0.5)

# Warm shading for RECON El Nino events
for _, ev in recon_el_nino_events.iterrows():
    t0, t1 = pd.Timestamp(ev["start"]), pd.Timestamp(ev["end"])
    seg = recon_9120_sm.sel(time=slice(t0, t1))
    ax1.fill_between(seg.time.values, 0.5, seg.values, color="red", alpha=0.25, linewidth=0)

# Cool shading for RECON La Nina events
for _, ev in recon_la_nina_events.iterrows():
    t0, t1 = pd.Timestamp(ev["start"]), pd.Timestamp(ev["end"])
    seg = recon_9120_sm.sel(time=slice(t0, t1))
    ax1.fill_between(seg.time.values, -0.5, seg.values, color="blue", alpha=0.25, linewidth=0)

ax1.set_ylabel("Nino 3.4 anomaly (°C)")
ax1.set_title("Nino 3.4 reconstruction vs CPC published series (1950-2023)")
ax1.legend(loc="best", ncol=2, fontsize=7, framealpha=0.9)
ax1.grid(alpha=0.3)

# Lower row: reconstruction minus CPC
ax2.plot(diff_recon_cpc.time, diff_recon_cpc.values, color="tab:purple", lw=0.9,
         label="RECON minus CPC")
ax2.axhline(0, color="k", lw=0.5)
ax2.set_ylabel("Difference (°C)")
ax2.set_xlabel("Year")
ax2.grid(alpha=0.3)
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2023-12-31"))
ax2.legend(loc="best", fontsize=7)

# Figure 2: square scatter with 1:1 line and statistics
fig2, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)
x = cpc_9120.values
y = recon_9120.values
ax.scatter(x, y, s=8, alpha=0.35, color="tab:blue", edgecolors="none", label="Monthly")

all_vals = np.concatenate([x, y])
vmin = np.floor(np.nanmin(all_vals) * 2.0) / 2.0 - 0.25
vmax = np.ceil(np.nanmax(all_vals) * 2.0) / 2.0 + 0.25
lim = [vmin, vmax]

ax.plot(lim, lim, "k--", lw=1.0, label="1:1 line")
ax.set_xlim(lim)
ax.set_ylim(lim)
ticks = np.arange(np.ceil(lim[0]), np.floor(lim[1]) + 0.5, 1.0)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("CPC 1991-2020 Nino 3.4 anomaly (°C)")
ax.set_ylabel("RECON 1991-2020 Nino 3.4 anomaly (°C)")
ax.set_title(f"Monthly Nino 3.4: r = {corr_value:.6f}, RMSE = {rmse_value:.4f} °C")
ax.grid(alpha=0.3)
ax.legend(loc="upper left", fontsize=8, frameon=False)
