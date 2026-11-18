# ============================================================
# Figure 1: Reconstruction, CPC reference, running mean,
#           threshold shading, and reconstruction-minus-CPC.
# ============================================================

fig, axes = plt.subplots(
    2, 1, figsize=(13, 8), sharex=True, constrained_layout=True
)

ax = axes[0]

ax.plot(
    recon_9120_s.index,
    recon_9120_s.values,
    color="tab:blue",
    lw=0.8,
    alpha=0.65,
    label="RECON 1991-2020 monthly",
)
ax.plot(
    cpc_9120_s.index,
    cpc_9120_s.values,
    color="tab:orange",
    lw=0.8,
    alpha=0.65,
    label="CPC 1991-2020 monthly",
)
ax.plot(
    run3_recon.index,
    run3_recon.values,
    color="black",
    lw=1.6,
    label="RECON 3-month running mean",
)

ax.axhline(
    0.5, color="red", ls="--", lw=1.0, label="+0.5 degC threshold"
)
ax.axhline(
    -0.5, color="blue", ls="--", lw=1.0, label="-0.5 degC threshold"
)

for k, ev in enumerate(recon_elnino_events):
    seg = run3_recon.loc[ev["start"]:ev["end"]]
    ax.fill_between(
        seg.index,
        0.5,
        seg.values,
        where=seg.values >= 0.5,
        color="red",
        alpha=0.25,
        interpolate=True,
        label="El Nino shading" if k == 0 else None,
    )

for k, ev in enumerate(recon_lanina_events):
    seg = run3_recon.loc[ev["start"]:ev["end"]]
    ax.fill_between(
        seg.index,
        -0.5,
        seg.values,
        where=seg.values <= -0.5,
        color="blue",
        alpha=0.25,
        interpolate=True,
        label="La Nina shading" if k == 0 else None,
    )

ax.set_ylabel("SST anomaly (degC)")
ax.set_title("Nino 3.4 reconstruction and CPC reference, 1950-2023")
ax.grid(alpha=0.25)
ax.legend(
    loc="upper left",
    fontsize=7,
    ncol=2,
    framealpha=0.9,
)

all_upper = np.concatenate(
    [
        recon_9120_s.values,
        cpc_9120_s.values,
        run3_recon.values,
    ]
)
upper_lim = float(np.nanmax(np.abs(all_upper))) * 1.1
ax.set_ylim(-upper_lim, upper_lim)

ax2 = axes[1]
diff = comp["RECON_9120"] - comp["CPC_9120"]
ax2.plot(diff.index, diff.values, color="tab:purple", lw=0.8)
ax2.axhline(0.0, color="grey", lw=0.8)
ax2.set_ylabel("RECON minus CPC (degC)")
ax2.set_xlabel("Year")
ax2.grid(alpha=0.25)

diff_lim = float(np.nanmax(np.abs(diff.values))) * 1.1
ax2.set_ylim(-diff_lim, diff_lim)
ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2023-12-31"))


# ============================================================
# Figure 2: Square scatter with 1:1 line and metrics
# ============================================================

x = comp["CPC_9120"].values
y = comp["RECON_9120"].values

vmax = max(
    float(np.nanmin(x)),
    float(np.nanmax(x)),
    float(np.nanmin(y)),
    float(np.nanmax(y)),
    abs(float(np.nanmin(x))),
    abs(float(np.nanmax(x))),
    abs(float(np.nanmin(y))),
    abs(float(np.nanmax(y))),
)
lim = vmax * 1.08
ticks = np.linspace(-lim, lim, 7)

fig2, ax3 = plt.subplots(figsize=(6.5, 6.5), constrained_layout=True)

ax3.scatter(
    x,
    y,
    s=12,
    alpha=0.55,
    edgecolors="none",
    label="monthly values",
)
ax3.plot(
    [-lim, lim],
    [-lim, lim],
    color="black",
    lw=1.0,
    label="1:1 line",
)

ax3.set_xlim(-lim, lim)
ax3.set_ylim(-lim, lim)
ax3.set_xticks(ticks)
ax3.set_yticks(ticks)
ax3.set_aspect("equal", adjustable="box")

ax3.set_xlabel("CPC 1991-2020 Nino 3.4 anomaly (degC)")
ax3.set_ylabel("RECON 1991-2020 Nino 3.4 anomaly (degC)")
ax3.set_title("Monthly Nino 3.4 anomalies, 1950-2023")
ax3.grid(alpha=0.25)

ax3.text(
    0.03,
    0.97,
    f"r = {corr:.6f}\nRMSE = {rmse:.4f} degC",
    transform=ax3.transAxes,
    va="top",
    ha="left",
    bbox=dict(boxstyle="round", facecolor="white", edgecolor="0.7", alpha=0.9),
)

ax3.legend(loc="lower right", fontsize=8)
