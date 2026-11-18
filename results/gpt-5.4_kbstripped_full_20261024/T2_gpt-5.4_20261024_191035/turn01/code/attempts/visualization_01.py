# Figure 1: time series comparison and residuals
fig, axes = plt.subplots(
    2, 1, figsize=(14, 8), sharex=True,
    gridspec_kw={"height_ratios": [3, 1]},
    constrained_layout=True
)

time = pd.to_datetime(RECON_9120["time"].values)
recon_vals = RECON_9120.values
cpc_vals = CPC_9120.values
run3_vals = RECON_9120_3mo.values
diff_vals = RECON_minus_CPC.values

ax = axes[0]
ax.plot(time, recon_vals, color="tab:blue", lw=1.2, label="RECON_9120 (ERSSTv5, 1991–2020 baseline)")
ax.plot(time, cpc_vals, color="k", lw=1.0, alpha=0.8, label="CPC Niño 3.4 (ersst5_91-20)")
ax.plot(time, run3_vals, color="tab:red", lw=1.6, label="RECON_9120 3-month centered mean")
ax.axhline(0.5, color="firebrick", lw=1.0, ls="--")
ax.axhline(-0.5, color="royalblue", lw=1.0, ls="--")

# Shade warm events
for _, ev in recon_warm_events.iterrows():
    start = pd.to_datetime(ev["start"])
    end = pd.to_datetime(ev["end"])
    mask = (time >= start) & (time <= end)
    ax.fill_between(
        time[mask],
        0.5,
        run3_vals[mask],
        where=np.isfinite(run3_vals[mask]) & (run3_vals[mask] >= 0.5),
        color="tomato",
        alpha=0.25,
        interpolate=True,
    )

# Shade cool events
for _, ev in recon_cool_events.iterrows():
    start = pd.to_datetime(ev["start"])
    end = pd.to_datetime(ev["end"])
    mask = (time >= start) & (time <= end)
    ax.fill_between(
        time[mask],
        -0.5,
        run3_vals[mask],
        where=np.isfinite(run3_vals[mask]) & (run3_vals[mask] <= -0.5),
        color="cornflowerblue",
        alpha=0.25,
        interpolate=True,
    )

ax.set_ylabel("degC")
ax.set_title("Figure 1. Niño 3.4 reconstruction from ERSSTv5 vs CPC, 1950–2023")
ax.legend(loc="best", fontsize=9)

ax2 = axes[1]
ax2.plot(time, diff_vals, color="purple", lw=1.0)
ax2.axhline(0, color="0.3", lw=0.8)
ax2.set_ylabel("degC")
ax2.set_xlabel("Time")
ax2.set_title("RECON_9120 - CPC ersst5_91-20")

# Figure 2: scatter plot with 1:1 line
fig2, ax = plt.subplots(figsize=(7, 7), constrained_layout=True)

x = CPC_9120.values
y = RECON_9120.values
mask = np.isfinite(x) & np.isfinite(y)
xv = x[mask]
yv = y[mask]

vmin = float(np.nanmin(np.concatenate([xv, yv])))
vmax = float(np.nanmax(np.concatenate([xv, yv])))

ax.scatter(xv, yv, s=12, alpha=0.45, color="tab:blue", edgecolors="none")
ax.plot([vmin, vmax], [vmin, vmax], color="k", lw=1.0, ls="--")
ax.set_xlim(vmin, vmax)
ax.set_ylim(vmin, vmax)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("CPC Niño 3.4 (degC)")
ax.set_ylabel("RECON_9120 (degC)")
ax.set_title("Figure 2. RECON_9120 vs CPC Niño 3.4 monthly anomalies")
ax.text(
    0.03, 0.97,
    f"r = {corr_value:.6f}\nRMSE = {rmse_value:.4f} degC",
    transform=ax.transAxes,
    ha="left", va="top",
    bbox=dict(facecolor="white", alpha=0.8, edgecolor="0.7")
)

plt.show()
