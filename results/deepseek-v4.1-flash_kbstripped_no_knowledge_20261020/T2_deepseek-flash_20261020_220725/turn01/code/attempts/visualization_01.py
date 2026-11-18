# Align series for plotting
recon_plot, cpc_plot = xr.align(recon_9120, cpc_9120, join="inner")
run_plot = run3_recon.sel(time=recon_plot["time"])
diff_plot = diff_series.sel(time=recon_plot["time"])

time_plot = pd.to_datetime(recon_plot["time"].values)
recon_v = recon_plot.values
cpc_v = cpc_plot.values
run_v = run_plot.values
diff_v = diff_plot.values

# Figure 1: monthly anomalies, running mean, thresholds, and event shading
fig1, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True, constrained_layout=True)
ax1, ax2 = axes

ax1.plot(time_plot, recon_v, color="tab:blue", lw=1.0, label="RECON_9120")
ax1.plot(time_plot, cpc_v, color="tab:orange", lw=1.0, alpha=0.85, label="CPC 1991-2020")
ax1.plot(time_plot, run_v, color="black", lw=1.6, label="RECON 3-month mean")
ax1.axhline(0.5, color="red", ls="--", lw=0.9, label="+0.5 °C")
ax1.axhline(-0.5, color="blue", ls="--", lw=0.9, label="-0.5 °C")
ax1.axhline(0, color="gray", lw=0.5)

time_pd = pd.to_datetime(time_plot)
for _, ev in elnino_recon.iterrows():
    start = pd.Timestamp(ev["start"])
    end = pd.Timestamp(ev["end"])
    mask = (time_pd >= start) & (time_pd <= end) & (run_v >= 0.5)
    if np.any(mask):
        ax1.fill_between(
            time_pd, 0.5, run_v, where=mask, color="red", alpha=0.25, interpolate=True
        )

for _, ev in lanina_recon.iterrows():
    start = pd.Timestamp(ev["start"])
    end = pd.Timestamp(ev["end"])
    mask = (time_pd >= start) & (time_pd <= end) & (run_v <= -0.5)
    if np.any(mask):
        ax1.fill_between(
            time_pd, -0.5, run_v, where=mask, color="blue", alpha=0.25, interpolate=True
        )

ax1.set_ylabel("SST anomaly (°C)")
ax1.set_title("Nino 3.4 reconstruction vs CPC; shading marks detected events")
ax1.legend(loc="best", fontsize=8, ncol=2)
ax1.grid(True, alpha=0.25)

ax2.plot(time_plot, diff_v, color="purple", lw=0.9)
ax2.axhline(0, color="gray", lw=0.5)
ax2.set_ylabel("RECON − CPC (°C)")
ax2.set_xlabel("Year")
ax2.grid(True, alpha=0.25)
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2023-12-31"))

# Figure 2: square scatter comparison with 1:1 line
fig2, ax = plt.subplots(figsize=(6.5, 6.5), constrained_layout=True)

ax.scatter(cpc_v, recon_v, s=9, alpha=0.5, edgecolors="none", color="tab:blue")

pooled = np.concatenate([cpc_v[np.isfinite(cpc_v)], recon_v[np.isfinite(recon_v)]])
lim_lo = np.floor(np.nanmin(pooled) * 10.0) / 10.0
lim_hi = np.ceil(np.nanmax(pooled) * 10.0) / 10.0
lims = [lim_lo, lim_hi]
ax.plot(lims, lims, "k--", lw=1.0, label="1:1")
ax.set_xlim(lims)
ax.set_ylim(lims)
ticks = np.linspace(lim_lo, lim_hi, 7)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("CPC 1991-2020 Nino 3.4 anomaly (°C)")
ax.set_ylabel("RECON_9120 Nino 3.4 anomaly (°C)")
ax.set_title("RECON_9120 vs CPC 1991-2020")
ax.grid(True, alpha=0.3)
ax.text(
    0.04,
    0.96,
    f"r = {pearson_r_rounded:.6f}\nRMSE = {rmse_rounded:.4f} °C\nmax |Δ| = {max_abs_error_rounded:.4f} °C",
    transform=ax.transAxes,
    va="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
)
ax.legend(loc="best", fontsize=9)
