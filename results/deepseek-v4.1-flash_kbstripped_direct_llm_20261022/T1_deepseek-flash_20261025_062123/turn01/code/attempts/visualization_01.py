# Figure 1: seasonal cycle of Nino 3.4 amplitude
months = np.arange(1, 13)
month_labels = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

colors = {
    "ERSSTv5": "black",
    "CESM2": "tab:blue",
    "GFDL-ESM4": "tab:orange",
    "MPI-ESM1-2-LR": "tab:green",
}

fig1, ax1 = plt.subplots(figsize=(10, 6))
for src in SOURCES:
    lw = 3.0 if src == "ERSSTv5" else 1.8
    ls = "-" if src == "ERSSTv5" else "--"
    ax1.plot(
        months,
        seasonal_std[src],
        color=colors[src],
        lw=lw,
        ls=ls,
        marker="o",
        ms=4,
        label=src,
    )
    ax1.scatter(
        [peak_month[src]],
        [peak_std[src]],
        s=90,
        color=colors[src],
        edgecolor="black",
        zorder=5,
    )

ax1.set_xticks(months)
ax1.set_xticklabels(month_labels)
ax1.set_xlabel("Calendar month")
ax1.set_ylabel("Nino 3.4 anomaly std dev (°C)")
ax1.set_title("Seasonal cycle of Nino 3.4 amplitude, 1950-2014")
ax1.grid(True, alpha=0.3)
ax1.legend(frameon=False)
fig1.tight_layout()


# Figure 2: Nino 3.4 anomaly time series, stacked panels
all_anom = np.concatenate([dfs[src]["anom"].values for src in SOURCES])
ymin = float(np.nanmin(all_anom))
ymax = float(np.nanmax(all_anom))
pad = 0.08 * (ymax - ymin) if ymax > ymin else 0.5
ylim = (ymin - pad, ymax + pad)

fig2, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

for ax, src in zip(axes, SOURCES):
    df = dfs[src]
    t = pd.to_datetime(df[["year", "month"]].assign(day=1))

    ax.plot(t, df["anom"].values, color=colors[src], lw=1.0)
    ax.axhline(1.0, color="red", ls="--", lw=0.8)
    ax.axhline(-1.0, color="red", ls="--", lw=0.8)
    ax.axhline(0.0, color="black", lw=0.5, alpha=0.5)

    ax.set_ylim(*ylim)
    ax.set_ylabel("Anomaly (°C)")
    ax.grid(True, alpha=0.2)

    txt = (
        f"{src}\n"
        f"SD = {std_all[src]:.2f} °C\n"
        f"|anom| > 1.0: {exceed_count[src]}/{n_months[src]} "
        f"({exceed_pct[src]:.1f}%)"
    )
    ax.text(
        0.01,
        0.96,
        txt,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

axes[-1].set_xlabel("Time")
axes[0].set_title("Nino 3.4 SST anomalies, 1950-2014")
fig2.tight_layout()
