import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------
# Figure 1 - historical OMZ thickness at the 60 umol/kg threshold
# ---------------------------------------------------------------------
panels1 = [("GFDL-ESM4 (historical 1971-2000)", MAP1["GFDL-ESM4"]),
           ("MPI-ESM1-2-LR (historical 1971-2000)", MAP1["MPI-ESM1-2-LR"]),
           ("WOA23 (1971-2000)", MAP1["WOA23"])]

pooled = np.concatenate([p[1][2].ravel() for p in panels1])
pooled = pooled[np.isfinite(pooled)]
vmax1 = float(np.nanpercentile(pooled, 98))

fig, axes = plt.subplots(1, 3, figsize=(18, 5.2), constrained_layout=True)
for ax, (title, (tlon, tlat, fld)) in zip(axes, panels1):
    pm = ax.pcolormesh(tlon, tlat, fld, cmap="viridis", vmin=0.0, vmax=vmax1,
                       shading="auto")
    ax.set_xlim(0, 360)
    ax.set_ylim(-80, 90)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
cb = fig.colorbar(pm, ax=axes, shrink=0.85, extend="max")
cb.set_label("OMZ thickness (m)")
fig.suptitle("Figure 1 - Historical OMZ thickness for O$_2$ < 60 umol kg$^{-1}$ "
             "(shared sequential scale, clipped at the pooled 98th percentile = "
             f"{vmax1:.0f} m; each panel on its own native grid)", fontsize=11)
plt.show()

# ---------------------------------------------------------------------
# Figure 2 - scenario change in OMZ thickness + CESM2 ssp245 thickness
# ---------------------------------------------------------------------
chg_pool = np.concatenate([MAP2["GFDL-ESM4"][2].ravel(), MAP2["MPI-ESM1-2-LR"][2].ravel()])
chg_pool = chg_pool[np.isfinite(chg_pool)]
vmax2 = float(np.nanpercentile(np.abs(chg_pool), 98))

cesm_fld = MAP2["CESM2"][2]
vmax3 = float(np.nanpercentile(cesm_fld[np.isfinite(cesm_fld)], 98))

fig = plt.figure(figsize=(18, 5.4), constrained_layout=True)
gs = fig.add_gridspec(1, 3)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

for ax, key, ttl in ((ax1, "GFDL-ESM4", "GFDL-ESM4: ssp245 - historical"),
                     (ax2, "MPI-ESM1-2-LR", "MPI-ESM1-2-LR: ssp245 - historical")):
    tlon, tlat, fld = MAP2[key]
    pm = ax.pcolormesh(tlon, tlat, fld, cmap="RdBu_r", vmin=-vmax2, vmax=vmax2,
                       shading="auto")
    ax.set_xlim(0, 360); ax.set_ylim(-80, 90)
    ax.set_title(ttl + " (m)", fontsize=10)
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")

tlon, tlat, fld = MAP2["CESM2"]
pm3 = ax3.pcolormesh(tlon, tlat, fld, cmap="viridis", vmin=0.0, vmax=vmax3,
                     shading="auto")
ax3.set_xlim(0, 360); ax3.set_ylim(-80, 90)
ax3.set_title("CESM2 ssp245 2071-2100 OMZ thickness (m)", fontsize=10)
ax3.set_xlabel("Longitude"); ax3.set_ylabel("Latitude")
ax3.text(0.02, 0.02, "No change formable:\nssp245 is r10i1p1f1 and no\nr10i1p1f1 historical o2 exists",
         transform=ax3.transAxes, fontsize=8, va="bottom",
         bbox=dict(boxstyle="round", fc="white", alpha=0.85))

cb1 = fig.colorbar(pm, ax=[ax1, ax2], shrink=0.85, extend="both")
cb1.set_label("dOMZ thickness (m)")
cb2 = fig.colorbar(pm3, ax=ax3, shrink=0.85, extend="max")
cb2.set_label("OMZ thickness (m)")

fig.suptitle("Figure 2 - Change in OMZ thickness at O$_2$ < 60 umol kg$^{-1}$ "
             f"(shared diverging scale, clipped at the 98th percentile of |change| = {vmax2:.0f} m); "
             f"right panel sequential scale clipped at {vmax3:.0f} m", fontsize=11)
plt.show()
