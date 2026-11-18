# -----------------------------------------------------------------------------
# Figure 1: reconstructed and CPC monthly anomalies plus 3-month running means
# -----------------------------------------------------------------------------
fig1, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(14, 8), sharex=True,
    gridspec_kw={"height_ratios": [3, 1]}
)

ax1.plot(
    recon_9120_s.index, recon_9120_s.values,
    color="tab:red", lw=0.7, alpha=0.45, label="RECON_9120 monthly"
)
ax1.plot(
    cpc9120_win.index, cpc9120_win.values,
    color="tab:blue", lw=0.7, alpha=0.45, label="CPC 1991-2020 monthly"
)
ax1.plot(
    recon_roll.index, recon_roll.values,
    color="tab:red", lw=2.0, label="RECON 3-month mean"
)
ax1.plot(
    cpc_roll.index, cpc_roll.values,
    color="tab:blue", lw=2.0, label="CPC 3-month mean"
)

ax1.axhline(0.5, color="0.35", ls="--", lw=1.0, label="±0.5°C threshold")
ax1.axhline(-0.5, color="0.35", ls="--", lw=1.0)

for ev in events_recon:
    seg = recon_roll.loc[ev["start"]:ev["end"]]
    ax1.fill_between(
        seg.index, 0.5, seg.values,
        where=(seg.values >= 0.5),
        color="tab:red", alpha=0.20, interpolate=True
    )

for ev in events_recon_ln:
    seg = recon_roll.loc[ev["start"]:ev["end"]]
    ax1.fill_between(
        seg.index, -0.5, seg.values,
        where=(seg.values <= -0.5),
        color="tab:blue", alpha=0.20, interpolate=True
    )

ax1.set_ylabel("Nino 3.4 anomaly (°C)")
ax1.set_title("Reconstructed versus CPC Nino 3.4 index")
ax1.legend(loc="upper left", ncol=2, fontsize=8)
ax1.grid(True, alpha=0.25)

ax2.plot(diff.index, diff.values, color="black", lw=0.8)
ax2.axhline(0, color="0.5", lw=0.8)
ax2.set_ylabel("RECON - CPC (°C)")
ax2.set_xlabel("Year")
ax2.grid(True, alpha=0.25)

ax1.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2023-12-31"))
fig1.tight_layout()

# -----------------------------------------------------------------------------
# Figure 2: square scatter of RECON_9120 versus CPC 1991-2020
# -----------------------------------------------------------------------------
fig2, ax = plt.subplots(figsize=(7, 7))

ax.scatter(
    comp["cpc"].values, comp["recon"].values,
    s=10, alpha=0.6, edgecolors="none"
)

lo = min(comp["cpc"].min(), comp["recon"].min())
hi = max(comp["cpc"].max(), comp["recon"].max())
span = hi - lo
lo = lo - 0.05 * span
hi = hi + 0.05 * span

ax.plot([lo, hi], [lo, hi], "k--", lw=1)
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ax.set_aspect("equal", adjustable="box")

ticks = np.linspace(lo, hi, 9)
ax.set_xticks(ticks)
ax.set_yticks(ticks)

ax.set_xlabel("CPC Nino 3.4 anomaly 1991-2020 (°C)")
ax.set_ylabel("RECON_9120 (°C)")
ax.grid(True, alpha=0.3)

ax.text(
    0.05, 0.95,
    f"r = {corr_recon_cpc_9120:.6f}\nRMSE = {rmse_recon_cpc_9120:.4f} °C",
    transform=ax.transAxes, va="top",
    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9)
)

fig2.tight_layout()
