colors = {"ERSSTv5": "black", "CESM2": "tab:blue", "GFDL-ESM4": "tab:red", "MPI-ESM1-2-LR": "tab:green"}
x = np.arange(1, 13)

# Figure 1: seasonal cycle of amplitude
fig1, ax1 = plt.subplots(figsize=(9, 5.5))
for k in src_order:
    s = monthly_std[k]
    lw = 3.5 if k == "ERSSTv5" else 1.8
    ax1.plot(x, s, color=colors[k], lw=lw, marker="o", ms=4, label=k,
             zorder=3 if k == "ERSSTv5" else 2)
    ipk = int(np.argmax(s))
    ax1.plot(x[ipk], s[ipk], marker="*", ms=18, color=colors[k], mec="gold" if k == "ERSSTv5" else "k",
             mew=1.2, zorder=5)
    ax1.annotate(f"{months_names[ipk]} {s[ipk]:.2f}", (x[ipk], s[ipk]), textcoords="offset points",
                 xytext=(0, 10), ha="center", fontsize=8, color=colors[k])
ax1.set_xticks(x)
ax1.set_xticklabels(months_names)
ax1.set_xlabel("Calendar month")
ax1.set_ylabel("Std of Niño 3.4 anomaly (°C)")
ax1.set_title("Seasonal cycle of ENSO amplitude, 1950-2014 (stars mark peak month)")
ax1.grid(alpha=0.3)
ax1.legend(loc="best")
fig1.tight_layout()

# Figure 2: stacked time series
allv = np.concatenate([nino[k].values for k in src_order])
ymax = float(np.ceil(max(abs(allv.min()), abs(allv.max())) * 2) / 2) + 0.25
fig2, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True, sharey=True)
for ax, k in zip(axes, src_order):
    t = nino[k]["time"].values
    ax.plot(t, nino[k].values, color=colors[k], lw=1.0)
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.axhline(-1.0, color="gray", ls="--", lw=1)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_ylim(-ymax, ymax)
    ax.set_ylabel("Niño 3.4 (°C)")
    ax.set_title(k, loc="left", fontsize=10, fontweight="bold")
    ax.text(0.99, 0.95, f"std = {std_all[k]:.2f} °C\n|anom| > 1.0 °C: {exceed[k]['percent']:.1f}%",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="gray", alpha=0.85))
    ax.grid(alpha=0.25)
axes[-1].set_xlabel("Year")
fig2.suptitle("Monthly Niño 3.4 SST anomaly (5S-5N, 170W-120W), 1950-2014; dashed lines at ±1.0 °C")
fig2.tight_layout()
