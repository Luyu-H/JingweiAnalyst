# ================= Figure 1: time series + event shading =================
ts = R.index.to_timestamp()

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(16, 8.5), sharex=True,
                               gridspec_kw={"height_ratios": [3, 1]})

ax0.plot(ts, R.values, color="0.55", lw=0.7, label="RECON_9120 monthly")
ax0.plot(ts, C.values, color="tab:orange", lw=0.7, alpha=0.85, label="CPC 1991-2020 monthly")

rm = RECON3.dropna()
rm_ts = rm.index.to_timestamp()
ax0.plot(rm_ts, rm.values, color="tab:blue", lw=1.9, label="RECON_9120 3-month running mean")

ax0.axhline(0.0, color="k", lw=0.6)
ax0.axhline(+0.5, color="tab:red", ls="--", lw=1.1)
ax0.axhline(-0.5, color="tab:blue", ls="--", lw=1.1)

warm = pd.Series(False, index=rm.index)
for e in EN_RECON:
    warm.loc[e["start"]:e["end"]] = True
cool = pd.Series(False, index=rm.index)
for e in LN_RECON:
    cool.loc[e["start"]:e["end"]] = True

ax0.fill_between(rm_ts, 0.5, rm.values, where=warm.values, interpolate=True,
                 color="tab:red", alpha=0.30, label="El Nino (5+ seasons >= +0.5)")
ax0.fill_between(rm_ts, -0.5, rm.values, where=cool.values, interpolate=True,
                 color="tab:blue", alpha=0.30, label="La Nina (5+ seasons <= -0.5)")

ax0.set_ylabel("SST anomaly (degC)")
ax0.set_title("Nino 3.4: RECON_9120 vs CPC 1991-2020 and the 3-month running mean, 1950-2023")
ax0.legend(loc="upper left", ncol=3, fontsize=8, frameon=False)

ax1.plot(ts, DIFF.reindex(R.index).values, color="k", lw=0.8)
ax1.axhline(0.0, color="tab:red", ls=":", lw=0.9)
ax1.set_ylabel("RECON - CPC (degC)")
ax1.set_xlabel("Year")

# ================= Figure 2: 1:1 scatter =================
fig2, ax2 = plt.subplots(figsize=(7.5, 7.5))
ax2.scatter(C[good].values, R[good].values, s=12, alpha=0.6,
            color="tab:blue", edgecolors="none")

lo = float(min(R[good].min(), C[good].min())) - 0.25
hi = float(max(R[good].max(), C[good].max())) + 0.25
ax2.plot([lo, hi], [lo, hi], "k--", lw=1.2, label="1:1")
ax2.set_xlim(lo, hi)
ax2.set_ylim(lo, hi)
ax2.set_aspect("equal", adjustable="box")
ticks = np.arange(np.floor(lo * 2.0) / 2.0, hi + 0.25, 0.5)
ax2.set_xticks(ticks)
ax2.set_yticks(ticks)
ax2.grid(alpha=0.3)
ax2.set_xlabel("CPC 1991-2020 Nino3.4 anomaly (degC)")
ax2.set_ylabel("RECON_9120 Nino3.4 anomaly (degC)")
ax2.set_title("Monthly Nino 3.4 anomaly, 1950-2023")
ax2.text(0.03, 0.97,
         "r = %.6f\nRMSE = %.4f degC\nmax|diff| = %.4f degC" %
         (corr_recon_cpc, rmse_recon_cpc, maxdev_recon_cpc),
         transform=ax2.transAxes, va="top", fontsize=10)
ax2.legend(loc="lower right")
