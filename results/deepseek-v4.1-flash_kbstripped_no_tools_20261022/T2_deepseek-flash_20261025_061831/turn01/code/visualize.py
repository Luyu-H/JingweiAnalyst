# ---------------------------------------------------------------------------
# Figure 1: monthly anomalies, running mean, detected events, and residual
# ---------------------------------------------------------------------------
fig1, (ax_top, ax_bot) = plt.subplots(
    2, 1, figsize=(15, 8), sharex=True, constrained_layout=True
)

ax_top.plot(
    recon_9120_common.index, recon_9120_common.values,
    color="tab:blue", lw=0.9, label="RECON_9120 monthly"
)
ax_top.plot(
    cpc_9120_common.index, cpc_9120_common.values,
    color="tab:orange", lw=0.9, alpha=0.85, label="CPC ersst5_91-20 monthly"
)
ax_top.plot(
    rm_recon.index, rm_recon.values,
    color="black", lw=1.5, label="RECON_9120 3-month running mean"
)

ax_top.axhline(0.5, color="red", ls="--", lw=1.0, label="±0.5 °C thresholds")
ax_top.axhline(-0.5, color="blue", ls="--", lw=1.0)

x_rm = rm_recon.index.to_numpy()
y_rm = rm_recon.to_numpy(dtype=float)
warm = warm_mask_recon.to_numpy(dtype=bool) & np.isfinite(y_rm)
cold = cold_mask_recon.to_numpy(dtype=bool) & np.isfinite(y_rm)

ax_top.fill_between(
    x_rm, 0.5, y_rm, where=warm,
    color="red", alpha=0.25, label="Detected El Niño"
)
ax_top.fill_between(
    x_rm, -0.5, y_rm, where=cold,
    color="blue", alpha=0.25, label="Detected La Niña"
)

ax_top.set_ylabel("SST anomaly (°C)")
ax_top.grid(True, alpha=0.25)
ax_top.legend(
    loc="upper center", bbox_to_anchor=(0.5, 1.18),
    ncol=4, fontsize=8, frameon=False
)

ax_bot.plot(
    diff_9120_cpc.index, diff_9120_cpc.values,
    color="0.2", lw=0.8
)
ax_bot.axhline(0.0, color="0.6", lw=0.8)
ax_bot.set_ylabel("RECON_9120 - CPC (°C)")
ax_bot.set_xlabel("Year")
ax_bot.grid(True, alpha=0.25)
ax_bot.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2023-12-31"))

fig1.suptitle(
    "Nino 3.4 anomaly: ERSSTv5 reconstruction versus CPC (1950-2023)",
    fontsize=13
)

# ---------------------------------------------------------------------------
# Figure 2: square scatter of RECON_9120 versus CPC ersst5_91-20
# ---------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(7.5, 7.5), constrained_layout=True)

x_scatter = valid_pair["cpc"].to_numpy()
y_scatter = valid_pair["recon"].to_numpy()

ax2.scatter(
    x_scatter, y_scatter,
    s=9, alpha=0.55, color="tab:blue", edgecolors="none"
)

all_values = np.concatenate([x_scatter, y_scatter])
lo = float(np.nanmin(all_values))
hi = float(np.nanmax(all_values))
pad = 0.06 * (hi - lo)
lims = (lo - pad, hi + pad)

ax2.plot(lims, lims, "k--", lw=1.0, label="1:1")
ax2.set_xlim(lims)
ax2.set_ylim(lims)
ax2.set_aspect("equal", adjustable="box")

tick_start = np.ceil(lims[0])
tick_end = np.floor(lims[1])
ticks = np.arange(tick_start, tick_end + 0.5, 1.0)
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)

ax2.set_xlabel("CPC ersst5_91-20 Nino 3.4 anomaly (°C)")
ax2.set_ylabel("RECON_9120 Nino 3.4 anomaly (°C)")
ax2.grid(True, alpha=0.25)

ax2.text(
    0.04, 0.96,
    f"r = {correlation_r:.6f}\nRMSE = {rmse_c:.4f} °C",
    transform=ax2.transAxes,
    ha="left", va="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="0.7")
)

ax2.legend(loc="lower right", fontsize=9)
ax2.set_title("Monthly Nino 3.4 anomalies, 1950-2023")
