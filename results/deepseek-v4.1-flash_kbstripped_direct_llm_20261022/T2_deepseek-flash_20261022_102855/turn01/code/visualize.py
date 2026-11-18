# ----------------------------------------------------------------------
# Figure 1: reconstruction, CPC monthly, running mean and event shading
# ----------------------------------------------------------------------
fig1, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(14, 8), sharex=True,
    gridspec_kw={"height_ratios": [3, 1]}
)

t_comp = comp.index.to_timestamp()

ax1.plot(
    t_comp, comp["RECON_9120"].values,
    color="tab:blue", lw=0.7, alpha=0.85,
    label="RECON_9120 monthly"
)
ax1.plot(
    t_comp, comp["CPC_9120"].values,
    color="tab:orange", lw=0.7, alpha=0.85,
    label="CPC 1991-2020 monthly"
)

rm = recon_9120_run3.loc[comp.index]
ax1.plot(
    t_comp, rm.values,
    color="black", lw=1.6,
    label="RECON 3-month running mean"
)

ax1.axhline(0.5, color="red", ls="--", lw=0.9, label="+0.5 °C")
ax1.axhline(-0.5, color="blue", ls="--", lw=0.9, label="-0.5 °C")

for ev in elnino_events_recon:
    seg = recon_9120_run3.loc[ev["start"]:ev["end"]].dropna()
    ax1.fill_between(
        seg.index.to_timestamp(),
        0.5,
        seg.values,
        where=seg.values >= 0.5,
        color="red",
        alpha=0.25,
        interpolate=True,
        linewidth=0,
    )

for ev in lanina_events_recon:
    seg = recon_9120_run3.loc[ev["start"]:ev["end"]].dropna()
    ax1.fill_between(
        seg.index.to_timestamp(),
        -0.5,
        seg.values,
        where=seg.values <= -0.5,
        color="blue",
        alpha=0.25,
        interpolate=True,
        linewidth=0,
    )

ax1.set_ylabel("SST anomaly (°C)")
ax1.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2023-12-31"))
ax1.grid(True, alpha=0.3)
ax1.legend(loc="upper left", ncol=2, fontsize=8)
ax1.set_title("Niño 3.4 reconstruction vs CPC 1991-2020")

ax2.plot(
    t_comp,
    (comp["RECON_9120"] - comp["CPC_9120"]).values,
    color="purple", lw=0.8,
)
ax2.axhline(0.0, color="black", lw=0.6)
ax2.set_ylabel("RECON − CPC (°C)")
ax2.set_xlabel("Time")
ax2.grid(True, alpha=0.3)

# ----------------------------------------------------------------------
# Figure 2: square scatter of reconstruction vs CPC monthly anomalies
# ----------------------------------------------------------------------
fig2, ax = plt.subplots(figsize=(6.5, 6.5))

x = comp["CPC_9120"].values
y = comp["RECON_9120"].values
ax.scatter(x, y, s=8, alpha=0.45, color="tab:green", edgecolors="none")

lo = float(min(np.nanmin(x), np.nanmin(y)))
hi = float(max(np.nanmax(x), np.nanmax(y)))
pad = 0.10
lim = (lo - pad, hi + pad)

ax.plot(lim, lim, "k--", lw=1.0, label="1:1 line")
ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_aspect("equal", adjustable="box")

ticks = np.arange(
    np.floor(lim[0] * 2.0) / 2.0,
    np.ceil(lim[1] * 2.0) / 2.0 + 0.25,
    0.5,
)
ax.set_xticks(ticks)
ax.set_yticks(ticks)

ax.set_xlabel("CPC Niño 3.4 anomaly 1991-2020 (°C)")
ax.set_ylabel("RECON_9120 (°C)")
ax.grid(True, alpha=0.3)
ax.text(
    0.04, 0.96,
    f"Pearson r = {corr_9120_cpc_6:.6f}\nRMSE = {rmse_9120_cpc_4:.4f} °C",
    transform=ax.transAxes,
    va="top",
    fontsize=9,
)
ax.legend(loc="lower right", fontsize=8)
ax.set_title("Monthly Niño 3.4: reconstruction vs CPC")
