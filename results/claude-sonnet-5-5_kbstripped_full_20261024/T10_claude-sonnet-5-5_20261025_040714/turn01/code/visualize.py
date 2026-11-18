import numpy as np
import matplotlib.pyplot as plt

def draw(ax, da, vmin, vmax, cmap, title):
    lon = da["lon"].values
    lat = da["lat"].values
    pm = ax.pcolormesh(lon, lat, da.values, vmin=vmin, vmax=vmax, cmap=cmap, shading="auto")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("Longitude (°E)")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    return pm

# ---- Figure 1: historical OMZ thickness at 60 umol/kg
pooled = np.concatenate([g_hist_m.values.ravel(), m_hist_m.values.ravel(), woa_thk_m.values.ravel()])
vmax1 = float(np.nanpercentile(pooled, 98))
fig1, axes1 = plt.subplots(1, 3, figsize=(18, 4.8), constrained_layout=True, sharey=True)
titles1 = ["GFDL-ESM4 historical 1971-2000", "MPI-ESM1-2-LR historical 1971-2000", "WOA23 observed"]
for ax, da, t in zip(axes1, [g_hist_m, m_hist_m, woa_thk_m], titles1):
    pm = draw(ax, da, 0, vmax1, "viridis", t)
axes1[0].set_ylabel("Latitude (°N)")
cb = fig1.colorbar(pm, ax=axes1, orientation="vertical", extend="max", shrink=0.9)
cb.set_label("OMZ thickness (m)")
fig1.suptitle(f"Fig. 1: OMZ thickness, O2 < 60 µmol/kg (shared 1° mask, {fig1_ncells} cells; colour clipped at 98th pct)")

# ---- Figure 2: change in thickness (ssp245 - hist) and CESM2 ssp245 thickness
pooled_d = np.abs(np.concatenate([g_dthk_m.values.ravel(), m_dthk_m.values.ravel()]))
vmax2 = float(np.nanpercentile(pooled_d, 98))
vmax3 = float(np.nanpercentile(c_thk_m.values, 98))
fig2, axes2 = plt.subplots(1, 3, figsize=(18, 4.8), constrained_layout=True, sharey=True)
pm_a = draw(axes2[0], g_dthk_m, -vmax2, vmax2, "RdBu_r", "GFDL-ESM4: ssp245 (2071-2100) - hist (1971-2000)")
draw(axes2[1], m_dthk_m, -vmax2, vmax2, "RdBu_r", "MPI-ESM1-2-LR: ssp245 (2071-2100) - hist (1971-2000)")
pm_c = draw(axes2[2], c_thk_m, 0, vmax3, "viridis", "CESM2 ssp245 thickness (no change can be formed for CESM2)")
axes2[0].set_ylabel("Latitude (°N)")
cb1 = fig2.colorbar(pm_a, ax=axes2[:2], orientation="vertical", extend="both", shrink=0.9)
cb1.set_label("Thickness change (m)")
cb2 = fig2.colorbar(pm_c, ax=axes2[2], orientation="vertical", extend="max", shrink=0.9)
cb2.set_label("OMZ thickness (m)")
fig2.suptitle(f"Fig. 2: OMZ thickness (<60 µmol/kg) change, shared diverging scale ±98th pct (1° mask, {fig2_ncells} cells)")
plt.show()
