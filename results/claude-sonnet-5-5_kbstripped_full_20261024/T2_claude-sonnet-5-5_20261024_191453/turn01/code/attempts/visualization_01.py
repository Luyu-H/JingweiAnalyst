fig1, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(14, 7),
                                gridspec_kw={"height_ratios": [3, 1.3]},
                                constrained_layout=True)
ax1.plot(t_all, c_vals, color="tab:blue", lw=0.8, alpha=0.7, label="CPC 1991-2020 monthly")
ax1.plot(t_all, r_vals, color="tab:orange", lw=0.8, alpha=0.8, label="RECON_9120 monthly")
ax1.plot(rm_t, rm_vals, color="k", lw=1.4, label="RECON_9120 3-month running mean")
ax1.axhline(0.5, color="firebrick", ls="--", lw=0.9, label="+0.5 °C")
ax1.axhline(-0.5, color="navy", ls="--", lw=0.9, label="-0.5 °C")
ax1.axhline(0, color="grey", lw=0.5)

first_w, first_c = True, True
for r in ev_recon.itertuples():
    m = (rm_t >= np.datetime64(r.start)) & (rm_t <= np.datetime64(r.end))
    ax1.fill_between(rm_t[m], 0.5, rm_vals[m], where=rm_vals[m] >= 0.5, color="red", alpha=0.45,
                     label="El Niño episodes" if first_w else None)
    first_w = False
for r in ln_recon.itertuples():
    m = (rm_t >= np.datetime64(r.start)) & (rm_t <= np.datetime64(r.end))
    ax1.fill_between(rm_t[m], -0.5, rm_vals[m], where=rm_vals[m] <= -0.5, color="dodgerblue", alpha=0.45,
                     label="La Niña episodes" if first_c else None)
    first_c = False
ax1.set_ylabel("Niño 3.4 SST anomaly (°C)")
ax1.set_title("Niño 3.4 index: RECON_9120 (ERSSTv5) vs CPC ersst5_91-20, 1950-2023")
ax1.legend(loc="upper left", ncol=3, fontsize=8, framealpha=0.9)

ax2.plot(t_all, diff_vals, color="purple", lw=0.8)
ax2.axhline(0, color="grey", lw=0.5)
ax2.set_ylabel("RECON_9120 − CPC (°C)")
ax2.set_xlabel("Year")
ax2.set_xlim(t_all[0], t_all[-1])

fig2, ax = plt.subplots(figsize=(6.5, 6.5), constrained_layout=True)
lo = float(np.floor(min(r_vals.min(), c_vals.min()) * 2) / 2)
hi = float(np.ceil(max(r_vals.max(), c_vals.max()) * 2) / 2)
ax.scatter(c_vals, r_vals, s=6, alpha=0.5, color="tab:blue")
ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="1:1 line")
ax.set_xlim(lo, hi)
ax.set_ylim(lo, hi)
ticks = np.arange(lo, hi + 0.01, 0.5)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("CPC 1991-2020 Niño 3.4 anomaly (°C)")
ax.set_ylabel("RECON_9120 Niño 3.4 anomaly (°C)")
ax.set_title("Monthly Niño 3.4: reconstruction vs CPC, 1950-2023")
ax.text(0.04, 0.96, f"r = {pearson_r:.6f}\nRMSE = {rmse:.4f} °C\nn = {n_samples}",
        transform=ax.transAxes, va="top", ha="left",
        bbox=dict(boxstyle="round", fc="white", ec="grey"))
ax.legend(loc="lower right")
