# Figure 1: 2x2 regression-pattern panels on shared diverging color scale
fields_for_scale = [ersst_reg_np.values] + [model_reg_np[m].values for m in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]]
pooled = np.concatenate([f[np.isfinite(f)] for f in fields_for_scale])
vmax = float(np.nanpercentile(np.abs(pooled), 98))

fig1, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
panel_order = ["ERSSTv5", "CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
panel_data = {
    "ERSSTv5": ersst_reg_np,
    "CESM2": model_reg_np["CESM2"],
    "GFDL-ESM4": model_reg_np["GFDL-ESM4"],
    "MPI-ESM1-2-LR": model_reg_np["MPI-ESM1-2-LR"],
}
panel_expl = {
    "ERSSTv5": ersst_expl_var_value,
    "CESM2": model_expl_var_value["CESM2"],
    "GFDL-ESM4": model_expl_var_value["GFDL-ESM4"],
    "MPI-ESM1-2-LR": model_expl_var_value["MPI-ESM1-2-LR"],
}
panel_pcorr = {
    "CESM2": pattern_corr["CESM2"],
    "GFDL-ESM4": pattern_corr["GFDL-ESM4"],
    "MPI-ESM1-2-LR": pattern_corr["MPI-ESM1-2-LR"],
}

mappable = None
for ax, name in zip(axes.ravel(), panel_order):
    da = panel_data[name]
    p = da.assign_coords(lon=da.lon % 360).sortby("lon")
    mappable = ax.pcolormesh(
        p["lon"].values,
        p["lat"].values,
        p.values,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        shading="auto",
    )
    ax.set_xlim(110, 260)
    ax.set_ylim(20, 70)
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    title = f"{name}\nExplained variance: {panel_expl[name]:.1f}%"
    if name != "ERSSTv5":
        title += f"; r={panel_pcorr[name]:.3f}"
    ax.set_title(title)

cbar = fig1.colorbar(mappable, ax=axes, shrink=0.92, extend="both")
cbar.set_label("PDO regression (°C per std)")
fig1.suptitle("Figure 1. PDO regression patterns, clipped at pooled 98th percentile", y=1.02)

# Figure 2: PDO indices plus published PDO, with 9-year smoothed ERSST and transition years
fig2, ax = plt.subplots(figsize=(13, 5), constrained_layout=True)

ax.plot(pd.to_datetime(ersst_pdo_idx.time.values), ersst_pdo_idx.values, label="ERSSTv5-derived PDO", lw=1.2, color="black")
ax.plot(pd.to_datetime(model_pdo_idx["CESM2"].time.values), model_pdo_idx["CESM2"].values, label="CESM2", lw=1.0)
ax.plot(pd.to_datetime(model_pdo_idx["GFDL-ESM4"].time.values), model_pdo_idx["GFDL-ESM4"].values, label="GFDL-ESM4", lw=1.0)
ax.plot(pd.to_datetime(model_pdo_idx["MPI-ESM1-2-LR"].time.values), model_pdo_idx["MPI-ESM1-2-LR"].values, label="MPI-ESM1-2-LR", lw=1.0)
ax.plot(pd.to_datetime(published_pdo.time.values), published_pdo.values, label="Published PDO", lw=1.0, ls="--", color="tab:purple")
ax.plot(pd.to_datetime(ersst_pdo_9yr.time.values), ersst_pdo_9yr.values, label="ERSSTv5 9-year smooth", lw=2.5, color="tab:red")

for yr in transition_summary["transition_years"]:
    ax.axvline(pd.Timestamp(f"{yr}-01-01"), color="0.5", ls=":", lw=1)

ax.axhline(0, color="0.4", lw=0.8)
ax.set_xlabel("Year")
ax.set_ylabel("Standardized index")
ax.set_title("Figure 2. Monthly PDO indices and smoothed ERSSTv5 PDO")
ax.legend(loc="best")
