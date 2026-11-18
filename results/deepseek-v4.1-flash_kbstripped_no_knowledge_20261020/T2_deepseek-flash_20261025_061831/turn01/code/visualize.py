# ---------------------------------------------------------------------------
# Figure 1: two rows sharing the 1950-2023 time axis
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(
    2, 1, sharex=True, figsize=(15, 8.5),
    gridspec_kw={"height_ratios": [3, 1]}, constrained_layout=True
)

t_recon = recon_m["time"].values
y_recon = np.asarray(recon_m.values, dtype=float)
t_cpc = cpc_m["time"].values
y_cpc = np.asarray(cpc_m.values, dtype=float)
t_run = recon_run["time"].values
y_run = np.asarray(recon_run.values, dtype=float)

ax1.plot(t_recon, y_recon, color="C0", lw=0.6, alpha=0.55,
         label="RECON Nino 3.4, 1991-2020 climatology (monthly)")
ax1.plot(t_cpc, y_cpc, color="C1", lw=0.6, alpha=0.55,
         label="CPC Nino 3.4 anomaly, 1991-2020 (monthly)")
ax1.plot(t_run, y_run, color="k", lw=1.6,
         label="RECON 3-month centred running mean")

ax1.axhline(0.5, color="firebrick", ls="--", lw=1.2)
ax1.axhline(-0.5, color="navy", ls="--", lw=1.2)
ax1.axhline(0.0, color="grey", lw=0.6)

ax1.fill_between(t_run, 0.5, y_run, where=(y_run >= 0.5), interpolate=True,
                 color="firebrick", alpha=0.30,
                 label="El Nino shading (>= +0.5 degC)")
ax1.fill_between(t_run, -0.5, y_run, where=(y_run <= -0.5), interpolate=True,
                 color="royalblue", alpha=0.30,
                 label="La Nina shading (<= -0.5 degC)")

ax1.set_ylabel("Nino 3.4 SST anomaly (degC)")
ax1.set_title("Figure 1 - Nino 3.4 rebuilt from ERSSTv5 vs the CPC published series, 1950-2023")
ax1.legend(loc="lower left", ncol=2, fontsize=8, framealpha=0.9)
ax1.grid(alpha=0.25)
ax1.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2023-12-31"))

ax2.plot(t_recon, y_recon - y_cpc, color="dimgrey", lw=0.8)
ax2.axhline(0.0, color="k", lw=0.8)
ax2.set_ylabel("RECON - CPC (degC)")
ax2.set_xlabel("Year")
ax2.grid(alpha=0.25)

# ---------------------------------------------------------------------------
# Figure 2: square scatter, one point per month
# ---------------------------------------------------------------------------
fig2, axs = plt.subplots(figsize=(7.5, 7.5), constrained_layout=True)

axs.scatter(y_cpc, y_recon, s=10, alpha=0.45, color="steelblue", edgecolors="none")

_lo = float(np.nanmin([np.nanmin(y_cpc), np.nanmin(y_recon)]))
_hi = float(np.nanmax([np.nanmax(y_cpc), np.nanmax(y_recon)]))
_pad = 0.15 * (_hi - _lo)
lims = [_lo - _pad, _hi + _pad]

axs.plot(lims, lims, "k--", lw=1.0, label="1:1 line")
axs.set_xlim(lims)
axs.set_ylim(lims)
axs.set_aspect("equal", adjustable="box")
axs.set_xticks(np.arange(np.ceil(lims[0] * 2) / 2, lims[1], 0.5))
axs.set_yticks(np.arange(np.ceil(lims[0] * 2) / 2, lims[1], 0.5))
axs.set_xlabel("CPC Nino 3.4 anomaly, 1991-2020 (degC)")
axs.set_ylabel("RECON Nino 3.4 anomaly, 1991-2020 (degC)")
axs.set_title("Figure 2 - Reconstruction vs CPC index, monthly, 1950-2023")
axs.grid(alpha=0.25)
axs.legend(loc="upper left", fontsize=9)
axs.text(0.97, 0.05,
         f"r = {corr_r:.6f}\nRMSE = {rmse:.4f} degC",
         transform=axs.transAxes, ha="right", va="bottom", fontsize=10,
         bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))
