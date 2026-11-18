# ---------------------------------------------------------------------
# Figure 1: time series and difference
# ---------------------------------------------------------------------
fig, axes = plt.subplots(
    2, 1, figsize=(14, 8), sharex=True, constrained_layout=True
)

ax = axes[0]
ax.plot(
    recon_9120_series.index, recon_9120_series.values,
    color="tab:blue", lw=0.8, alpha=0.35, label="RECON_9120 monthly"
)
ax.plot(
    cpc91_series.index, cpc91_series.values,
    color="tab:orange", lw=0.8, alpha=0.35, label="CPC 1991-2020 monthly"
)
ax.plot(
    recon_9120_run.index, recon_9120_run.values,
    color="tab:blue", lw=1.8, label="RECON_9120 3-month mean"
)
ax.plot(
    cpc91_run.index, cpc91_run.values,
    color="tab:orange", lw=1.8, label="CPC 3-month mean"
)

ax.axhline(0.5, color="k", ls="--", lw=1)
ax.axhline(-0.5, color="k", ls="--", lw=1)

# Shade warm events between +0.5 C and the RECON running mean.
for ev in recon_elnino_events:
    mask = (
        (recon_9120_run.index >= ev["start"])
        & (recon_9120_run.index <= ev["end"])
    )
    x = recon_9120_run.index[mask]
    y = recon_9120_run.values[mask]
    ax.fill_between(
        x, 0.5, y, where=(y >= 0.5),
        color="red", alpha=0.20, interpolate=True
    )

# Shade cool events between -0.5 C and the RECON running mean.
for ev in recon_lanina_events:
    mask = (
        (recon_9120_run.index >= ev["start"])
        & (recon_9120_run.index <= ev["end"])
    )
    x = recon_9120_run.index[mask]
    y = recon_9120_run.values[mask]
    ax.fill_between(
        x, -0.5, y, where=(y <= -0.5),
        color="blue", alpha=0.20, interpolate=True
    )

ax.set_ylabel("SST anomaly (°C)")
ax.set_title("Nino 3.4: RECON_9120 vs CPC 1991-2020, 1950-2023")
ax.legend(loc="best", ncol=2, fontsize=8)
ax.grid(alpha=0.25)

axes[1].plot(
    (recon_9120_series - cpc91_series).index,
    (recon_9120_series - cpc91_series).values,
    color="k", lw=1.0
)
axes[1].axhline(0.0, color="gray", lw=0.8)
axes[1].set_ylabel("RECON - CPC (°C)")
axes[1].set_xlabel("Year")
axes[1].set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2024-01-01"))
axes[1].grid(alpha=0.25)

# ---------------------------------------------------------------------
# Figure 2: square scatter plot
# ---------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(7, 7), constrained_layout=True)

x = pair["CPC_9120"].values
y = pair["RECON_9120"].values

ax2.scatter(
    x, y, s=10, alpha=0.55, color="tab:blue", edgecolors="none"
)

vmin = float(min(np.nanmin(x), np.nanmin(y)))
vmax = float(max(np.nanmax(x), np.nanmax(y)))
pad = 0.08 * (vmax - vmin)
lims = (vmin - pad, vmax + pad)

ax2.plot(lims, lims, "k--", lw=1.2, label="1:1")
ax2.set_xlim(lims)
ax2.set_ylim(lims)
ticks = np.linspace(lims[0], lims[1], 7)
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)
ax2.set_aspect("equal", adjustable="box")
ax2.set_xlabel("CPC 1991-2020 Nino 3.4 anomaly (°C)")
ax2.set_ylabel("RECON_9120 Nino 3.4 anomaly (°C)")
ax2.set_title("RECON_9120 vs CPC 1991-2020")
ax2.grid(alpha=0.25)
ax2.text(
    0.03, 0.97,
    f"r = {corr:.6f}\nRMSE = {rmse:.4f} °C\nn = {len(pair)}",
    transform=ax2.transAxes,
    va="top", ha="left",
    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9),
)
