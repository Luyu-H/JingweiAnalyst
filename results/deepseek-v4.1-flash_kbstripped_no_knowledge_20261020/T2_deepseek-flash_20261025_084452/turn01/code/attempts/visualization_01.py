# Figure 1: time series, running mean, event shading, and reconstruction-minus-CPC difference.
fig, axes = plt.subplots(
    2, 1, figsize=(14, 8), sharex=True, constrained_layout=True
)

ax = axes[0]
ax.plot(
    recon_9120.time.values,
    recon_9120.values,
    color="tab:blue",
    lw=1.0,
    label="RECON_9120",
)
ax.plot(
    cpc_9120.time.values,
    cpc_9120.values,
    color="tab:orange",
    lw=1.0,
    alpha=0.85,
    label="CPC 1991-2020",
)
ax.plot(
    recon_run3.time.values,
    recon_run3.values,
    color="black",
    lw=1.5,
    label="RECON 3-month running mean",
)

ax.axhline(0.5, color="red", ls="--", lw=1.0, label="±0.5 °C thresholds")
ax.axhline(-0.5, color="blue", ls="--", lw=1.0)

# Warm shading for detected El Nino events
for _, ev in el_nino_recon.iterrows():
    seg = recon_run3.sel(time=slice(ev["start"], ev["end"]))
    ax.fill_between(
        seg.time.values,
        0.5,
        seg.values,
        where=seg.values >= 0.5,
        color="red",
        alpha=0.25,
        interpolate=True,
    )

# Cool shading for detected La Nina events
for _, ev in la_nina_recon.iterrows():
    seg = recon_run3.sel(time=slice(ev["start"], ev["end"]))
    ax.fill_between(
        seg.time.values,
        -0.5,
        seg.values,
        where=seg.values <= -0.5,
        color="blue",
        alpha=0.25,
        interpolate=True,
    )

ax.set_ylabel("Nino 3.4 anomaly (°C)")
ax.set_title("Nino 3.4 reconstruction vs CPC, 1950-2023")
ax.legend(loc="best", fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)

ax2 = axes[1]
ax2.plot(
    diff_field.time.values,
    diff_field.values,
    color="tab:purple",
    lw=1.0,
)
ax2.axhline(0.0, color="gray", lw=0.8)
ax2.set_ylabel("RECON - CPC (°C)")
ax2.set_xlabel("Year")
ax2.grid(True, alpha=0.3)

# Figure 2: square scatter of monthly RECON_9120 vs CPC 1991-2020.
fig2, axs = plt.subplots(figsize=(7, 7), constrained_layout=True)

x = cpc_9120_al.values
y = recon_9120_al.values
mask = np.isfinite(x) & np.isfinite(y)

axs.scatter(
    x[mask],
    y[mask],
    s=12,
    alpha=0.5,
    color="tab:blue",
    edgecolors="none",
)

lo = float(np.nanmin([np.nanmin(x[mask]), np.nanmin(y[mask])]))
hi = float(np.nanmax([np.nanmax(x[mask]), np.nanmax(y[mask])]))
pad = 0.3
lim_min = np.floor((lo - pad) * 2.0) / 2.0
lim_max = np.ceil((hi + pad) * 2.0) / 2.0
ticks = np.arange(lim_min, lim_max + 0.25, 0.5)

axs.plot([lim_min, lim_max], [lim_min, lim_max], "k--", lw=1.0, label="1:1 line")
axs.set_xlim(lim_min, lim_max)
axs.set_ylim(lim_min, lim_max)
axs.set_xticks(ticks)
axs.set_yticks(ticks)
axs.set_aspect("equal", adjustable="box")
axs.set_xlabel("CPC 1991-2020 Nino 3.4 anomaly (°C)")
axs.set_ylabel("RECON_9120 Nino 3.4 anomaly (°C)")
axs.set_title("Monthly RECON_9120 vs CPC 1991-2020")
axs.grid(True, alpha=0.3)

axs.text(
    0.05,
    0.95,
    f"r = {pearson_r:.6f}\nRMSE = {rmse:.4f} °C",
    transform=axs.transAxes,
    va="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
)
axs.legend(loc="lower right", fontsize=8)
