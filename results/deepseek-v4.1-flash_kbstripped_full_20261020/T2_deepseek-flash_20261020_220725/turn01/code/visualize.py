import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Figure 1 - monthly anomalies, 3-month running means, thresholds and shading
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True,
                         constrained_layout=True,
                         gridspec_kw={"height_ratios": [3, 1]})
ax = axes[0]

ax.plot(recon_9120["time"], recon_9120.values, color="tab:blue",
        lw=0.7, alpha=0.40, label="RECON 1991-2020, monthly")
ax.plot(cpc_9120["time"], cpc_9120.values, color="tab:orange",
        lw=0.7, alpha=0.40, label="CPC 1991-2020, monthly")

ax.plot(rm_recon_9120["time"], rm_recon_9120.values, color="tab:blue",
        lw=2.0, label="RECON, 3-month running mean")
ax.plot(rm_oni["time"], rm_oni.values, color="tab:orange",
        lw=2.0, label="CPC ONI")

# El Nino / La Nina shading from the detected episodes (Step 8)
for _, ev in recon_elnino.iterrows():
    seg = rm_recon_9120.sel(time=slice(ev["start"], ev["end"]))
    ax.fill_between(seg["time"].values, 0.5, seg.values,
                    where=(seg.values >= 0.5), color="red", alpha=0.30,
                    linewidth=0, interpolate=True)
for _, ev in recon_lanina.iterrows():
    seg = rm_recon_9120.sel(time=slice(ev["start"], ev["end"]))
    ax.fill_between(seg["time"].values, -0.5, seg.values,
                    where=(seg.values <= -0.5), color="blue", alpha=0.30,
                    linewidth=0, interpolate=True)

ax.axhline(0.0, color="k", lw=0.6)
ax.axhline(0.5, color="red", ls="--", lw=1.1,
           label="+0.5 °C El Niño threshold")
ax.axhline(-0.5, color="blue", ls="--", lw=1.1,
           label="−0.5 °C La Niña threshold")

ax.set_ylabel("Niño 3.4 SST anomaly (°C)")
ax.set_title(
    "Niño 3.4 rebuilt from ERSSTv5 (1991-2020 base) vs CPC published index, 1950-2023\n"
    "shading: detected El Niño (warm) / La Niña (cool) episodes, ≥5 consecutive seasons "
    "beyond ±0.5 °C"
)
ax.legend(loc="lower left", fontsize=8, ncol=2, framealpha=0.9)
ax.grid(alpha=0.25)

# Lower row: reconstruction minus published series on its own scale
ax2 = axes[1]
ax2.plot(diff_series["time"], diff_series.values, color="k", lw=0.8)
ax2.axhline(0.0, color="gray", lw=0.7)
ax2.set_ylabel("RECON − CPC (°C)", fontsize=9)
ax2.set_xlabel("Year")
ax2.grid(alpha=0.25)

ax2.set_xlim(pd.Timestamp("1950-01-01"), pd.Timestamp("2023-12-31"))

# ---------------------------------------------------------------------------
# Figure 2 - square scatter RECON_9120 (y) vs CPC 1991-2020 (x)
# ---------------------------------------------------------------------------
paired = (
    xr.Dataset({"recon": recon_9120, "cpc": cpc_9120})
    .to_dataframe()
    .dropna()
)

fig2, axs = plt.subplots(figsize=(7.2, 7.2), constrained_layout=True)
axs.scatter(paired["cpc"], paired["recon"], s=9, alpha=0.45,
            color="tab:purple", edgecolors="none")

lo = float(min(paired["cpc"].min(), paired["recon"].min()))
hi = float(max(paired["cpc"].max(), paired["recon"].max()))
pad = 0.10 * (hi - lo)
lo, hi = lo - pad, hi + pad
ticks = np.arange(np.floor(lo * 2) / 2, np.ceil(hi * 2) / 2 + 0.5, 0.5)

axs.plot([lo, hi], [lo, hi], "k--", lw=1.2, label="1:1 reference")
axs.set_xlim(lo, hi)
axs.set_ylim(lo, hi)
axs.set_xticks(ticks)
axs.set_yticks(ticks)
axs.set_aspect("equal", adjustable="box")
axs.set_xlabel("CPC Niño 3.4 anomaly, 1991-2020 base (°C)")
axs.set_ylabel("RECON Niño 3.4 anomaly, 1991-2020 base (°C)")
axs.set_title(
    f"Monthly Niño 3.4, 1950-2023 ({len(paired)} months)\n"
    f"r = {correlation:.6f}    RMSE = {rmse:.4f} °C"
)
axs.text(
    0.03, 0.97,
    f"r = {correlation:.6f}\nRMSE = {rmse:.4f} °C\n"
    f"bias = {bias:+.4f} °C\nmax |Δ| = {max_abs_error:.4f} °C "
    f"({max_abs_error_year}-{max_abs_error_month:02d})",
    transform=axs.transAxes, va="top", ha="left", fontsize=9,
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
)
axs.legend(loc="lower right", fontsize=9)
axs.grid(alpha=0.25)
