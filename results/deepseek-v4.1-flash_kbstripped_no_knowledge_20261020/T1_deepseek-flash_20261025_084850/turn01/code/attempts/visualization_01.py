# ---------------------------------------------------------------------------
# Figure 1 - seasonal cycle of amplitude (calendar-month std, degC)
# ---------------------------------------------------------------------------
colors = {
    "ERSSTv5": "black",
    "CESM2": "#d62728",
    "GFDL-ESM4": "#1f77b4",
    "MPI-ESM1-2-LR": "#2ca02c",
}

fig1, ax1 = plt.subplots(figsize=(10, 6))

for s in sources:
    if s == "ERSSTv5":
        ax1.plot(months_axis, std_month[s], color=colors[s], lw=3.2,
                 marker="o", ms=7, zorder=5, label=f"{s} (obs)")
    else:
        ax1.plot(months_axis, std_month[s], color=colors[s], lw=1.9,
                 ls="--", marker="s", ms=5, zorder=3, label=s)

for s in sources:
    ax1.plot(peak_month[s], peak_value[s], marker="*", ms=17,
             color=colors[s], mec="black", mew=0.6, zorder=6)
    ax1.annotate(f"{month_names[peak_month[s]-1]}\n{peak_value[s]:.2f}",
                 xy=(peak_month[s], peak_value[s]),
                 xytext=(0, 14), textcoords="offset points",
                 ha="center", va="bottom", fontsize=8, color=colors[s])

ax1.set_xticks(months_axis)
ax1.set_xticklabels(month_names)
ax1.set_xlabel("Calendar month")
ax1.set_ylabel("Nino 3.4 anomaly std (degC)")
ax1.set_title("Seasonal cycle of ENSO amplitude, 1950-2014\n"
              "Nino 3.4 (5S-5N, 170W-120W), population std per calendar month")
ax1.grid(alpha=0.3)
ax1.legend(loc="upper left", bbox_to_anchor=(0.0, -0.12), ncol=4, frameon=False)

# ---------------------------------------------------------------------------
# Figure 2 - stacked anomaly time series, shared y axis
# ---------------------------------------------------------------------------
pooled_max = max(float(np.nanmax(np.abs(series_values[s]))) for s in sources)
ymax = np.ceil(pooled_max * 10.0) / 10.0 + 0.1

fig2, axes2 = plt.subplots(4, 1, figsize=(12, 11), sharex=True, sharey=True,
                           constrained_layout=True)

for ax, s in zip(axes2, sources):
    ax.axhline(0.0, color="0.55", lw=0.8, zorder=1)
    ax.axhline(1.0, color="#b2182b", ls=":", lw=1.3, zorder=2)
    ax.axhline(-1.0, color="#2166ac", ls=":", lw=1.3, zorder=2)
    ax.plot(series_time[s], series_values[s], color=colors[s], lw=1.0, zorder=3)
    ax.fill_between(series_time[s], 0.0, series_values[s],
                    where=series_values[s] > 1.0, color="#b2182b", alpha=0.35,
                    interpolate=True, zorder=2)
    ax.fill_between(series_time[s], 0.0, series_values[s],
                    where=series_values[s] < -1.0, color="#2166ac", alpha=0.35,
                    interpolate=True, zorder=2)
    ax.set_ylim(-ymax, ymax)
    ax.set_ylabel("Anomaly (degC)", fontsize=9)
    ax.grid(alpha=0.25)
    ax.text(0.005, 0.94,
            f"{s}\nstd = {std_full[s]:.2f} degC   |anom|>1 degC: "
            f"{exc_count[s]} months ({exc_percent[s]:.1f} %)",
            transform=ax.transAxes, ha="left", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.85))

axes2[-1].set_xlabel("Time")
axes2[0].set_title(
    "Nino 3.4 SST anomaly, 1950-2014 (shared y axis, +/-1 degC dotted); "
    "shading marks |anomaly| > 1 degC", fontsize=11)
