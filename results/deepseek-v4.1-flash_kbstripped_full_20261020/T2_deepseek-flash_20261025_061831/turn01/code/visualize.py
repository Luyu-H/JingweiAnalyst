# Figure 1: two-row time series
fig, axes = plt.subplots(
    2, 1, sharex=True, figsize=(12, 7), constrained_layout=True
)

ax = axes[0]
ax.plot(
    recon_9120.time.values,
    recon_9120.values,
    color="tab:blue",
    lw=0.6,
    alpha=0.45,
    label="RECON 1991-2020 monthly",
)
ax.plot(
    cpc_9120.time.values,
    cpc_9120.values,
    color="tab:gray",
    lw=0.6,
    alpha=0.45,
    label="CPC 1991-2020 monthly",
)
ax.plot(
    recon_run3.time.values,
    recon_run3.values,
    color="black",
    lw=1.3,
    label="RECON 3-month running mean",
)
ax.axhline(0.5, color="tab:red", ls="--", lw=1.0, label="+0.5 °C threshold")
ax.axhline(-0.5, color="tab:blue", ls="--", lw=1.0, label="-0.5 °C threshold")
ax.axhline(0.0, color="k", lw=0.5)

# Event shading masks
x_run = pd.to_datetime(recon_run3.time.values)
y_run = recon_run3.values

el_mask = np.zeros(x_run.shape, dtype=bool)
for _, row in el_events.iterrows():
    el_mask |= (x_run >= pd.Timestamp(row["start"])) & (
        x_run <= pd.Timestamp(row["end"])
    )

la_mask = np.zeros(x_run.shape, dtype=bool)
for _, row in la_events.iterrows():
    la_mask |= (x_run >= pd.Timestamp(row["start"])) & (
        x_run <= pd.Timestamp(row["end"])
    )

ax.fill_between(
    x_run,
    0.5,
    y_run,
    where=el_mask & (y_run >= 0.5),
    color="tab:red",
    alpha=0.30,
    step=None,
    label="El Nino shading",
)
ax.fill_between(
    x_run,
    -0.5,
    y_run,
    where=la_mask & (y_run <= -0.5),
    color="tab:blue",
    alpha=0.30,
    step=None,
    label="La Nina shading",
)

ax.set_ylabel("SST anomaly (°C)")
ax.set_title("Nino 3.4 reconstruction vs CPC, 1950-2023")
ax.grid(True, ls=":", alpha=0.4)
ax.legend(
    loc="upper left",
    bbox_to_anchor=(1.01, 1.0),
    fontsize=8,
    frameon=False,
)

ax2 = axes[1]
diff_plot = recon_9120 - cpc_9120
ax2.plot(
    diff_plot.time.values,
    diff_plot.values,
    color="tab:purple",
    lw=0.8,
    label="RECON - CPC",
)
ax2.axhline(0.0, color="k", lw=0.5)
ax2.set_ylabel("RECON - CPC (°C)")
ax2.set_xlabel("Year")
ax2.grid(True, ls=":", alpha=0.4)
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2024-01-01"))
ax2.legend(loc="upper left", fontsize=8)

# Figure 2: square monthly scatter
fig2, ax3 = plt.subplots(figsize=(6.5, 6.5), constrained_layout=True)

xsc = np.asarray(cpc_9120.values).ravel()
ysc = np.asarray(recon_9120.values).ravel()
good = np.isfinite(xsc) & np.isfinite(ysc)
xsc = xsc[good]
ysc = ysc[good]

ax3.scatter(xsc, ysc, s=10, alpha=0.5, color="tab:blue", edgecolors="none")

lim_min = min(np.nanmin(xsc), np.nanmin(ysc))
lim_max = max(np.nanmax(xsc), np.nanmax(ysc))
pad = 0.05 * (lim_max - lim_min)
lims = (lim_min - pad, lim_max + pad)

ax3.plot(lims, lims, "k--", lw=0.9, label="1:1 line")
ax3.set_xlim(lims)
ax3.set_ylim(lims)
ax3.set_aspect("equal", adjustable="box")

tick_min = np.floor(lims[0] * 2.0) / 2.0
tick_max = np.ceil(lims[1] * 2.0) / 2.0
ticks = np.arange(tick_min, tick_max + 0.001, 0.5)
ax3.set_xticks(ticks)
ax3.set_yticks(ticks)

ax3.set_xlabel("CPC Nino 3.4 anomaly (°C)")
ax3.set_ylabel("RECON Nino 3.4 anomaly (°C)")
ax3.set_title("Monthly Nino 3.4 comparison, 1950-2023")
ax3.grid(True, ls=":", alpha=0.4)
ax3.legend(loc="lower right", fontsize=8)

ax3.text(
    0.03,
    0.97,
    f"r = {r:.6f}\nRMSE = {rmse:.4f} °C",
    transform=ax3.transAxes,
    va="top",
    ha="left",
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.7", alpha=0.9),
)
