import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# ---------------------------------------------------------------------------
# Figure 1: global nitracline-depth maps, shared sequential colour scale
# ---------------------------------------------------------------------------
map_order = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR", "WOA23"]
map_data = {
    "CESM2": nt_regrid["CESM2"],
    "GFDL-ESM4": nt_regrid["GFDL-ESM4"],
    "MPI-ESM1-2-LR": nt_regrid["MPI-ESM1-2-LR"],
    "WOA23": nt_woa,
}

fig1, axes1 = plt.subplots(
    2, 2, figsize=(12, 8), sharex=True, sharey=True, constrained_layout=True
)
vmin1, vmax1 = 0.0, 1000.0
mesh1 = None

for ax, name in zip(axes1.ravel(), map_order):
    data = np.ma.masked_invalid(map_data[name])
    mesh1 = ax.pcolormesh(
        lon_woa,
        lat_woa,
        data,
        shading="auto",
        cmap="viridis_r",
        vmin=vmin1,
        vmax=vmax1,
    )
    ax.set_title(name)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

cbar1 = fig1.colorbar(
    mesh1, ax=axes1.ravel().tolist(), orientation="vertical", shrink=0.85
)
cbar1.set_label("Nitracline depth (m)")

# ---------------------------------------------------------------------------
# Figure 2: equatorial Pacific longitude-depth nitrate sections
# ---------------------------------------------------------------------------
section_order = ["WOA23", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

all_section_vals = np.concatenate(
    [sections[name]["conc_plot"].ravel() for name in section_order]
)
vmax2 = np.nanpercentile(all_section_vals, 99)
vmin2 = 0.0

fig2, axes2 = plt.subplots(
    4, 1, figsize=(12, 11), sharex=True, sharey=True, constrained_layout=True
)
mesh2 = None

for ax, name in zip(axes2, section_order):
    s = sections[name]
    mesh2 = ax.pcolormesh(
        s["lon"],
        s["depth_plot"],
        s["conc_plot"],
        shading="auto",
        cmap="viridis",
        vmin=vmin2,
        vmax=vmax2,
    )
    ax.plot(s["lon"], s["nt_line"], color="k", linewidth=1.6)
    ax.set_ylim(1000, 0)
    ax.set_title(name)
    ax.set_ylabel("Depth (m)")

axes2[-1].set_xlabel("Longitude")

xticks = [120, 150, 180, 210, 240, 270]
axes2[-1].set_xticks(xticks)


def lon_label(x, pos=None):
    x = x % 360
    if abs(x - 180) < 1e-6:
        return "180"
    if x > 180:
        return f"{360 - x:.0f}W"
    return f"{x:.0f}E"


axes2[-1].xaxis.set_major_formatter(FuncFormatter(lon_label))

cbar2 = fig2.colorbar(
    mesh2, ax=axes2.tolist(), orientation="vertical", shrink=0.85
)
cbar2.set_label(r"Nitrate ($\mu$mol kg$^{-1}$)")
