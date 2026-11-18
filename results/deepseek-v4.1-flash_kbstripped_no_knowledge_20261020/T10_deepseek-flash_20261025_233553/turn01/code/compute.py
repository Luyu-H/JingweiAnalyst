# =====================================================================
# 1. WOA23 observed climatology (oxygen, saturation, AOU) + grid info
# =====================================================================
o2_woa  = load_climatology_data(dataset="woa23", variable="oxygen")
po2_woa = load_climatology_data(dataset="woa23", variable="percent_oxygen_saturation")
aou_woa = load_climatology_data(dataset="woa23", variable="apparent_oxygen_utilization")

lev_woa, bnds_woa = load_depth_levels(dataset="woa23", variable="oxygen")
dz_woa = xr.DataArray(bnds_woa[:, 1] - bnds_woa[:, 0], dims=["depth"], coords={"depth": lev_woa})
area_woa = load_areas(dataset="woa23")
vol_woa  = area_woa * dz_woa                       # m3 per cell (ocean masked by load_areas)
landsea_woa = load_landsea_mask()

# =====================================================================
# 2. CMIP6 model fields: day-weighted 1971-2000 / 2071-2100 means
# =====================================================================
def model_oxygen(dataset, model, t0, t1):
    """Return (o2 mean in mol m-3, o2 mean in umol kg-1) for one model/experiment."""
    o2 = load_regular_data(dataset=dataset, variable="o2",     model=model, time_range=[t0, t1])
    th = load_regular_data(dataset=dataset, variable="thetao", model=model, time_range=[t0, t1])
    sa = load_regular_data(dataset=dataset, variable="so",     model=model, time_range=[t0, t1])

    o2m = temporal_aggregate(o2, freq="mean", op="mean").assign_attrs(**o2.attrs)
    thm = temporal_aggregate(th, freq="mean", op="mean")
    sam = temporal_aggregate(sa, freq="mean", op="mean")

    sigma = compute_potential_density_anomaly(temp_data=thm, sal_data=sam, ref_pressure=None)

    o2kg = transform_unit(data=o2m,
                          from_unit=o2m.attrs.get("units", "mol m-3"),
                          to_unit="umol/kg",
                          density=sigma)
    o2kg = o2kg.assign_attrs(**o2.attrs)
    o2kg.attrs["units"] = "umol kg-1"
    return o2m, o2kg

gfdl_hist_mol, gfdl_hist = model_oxygen("cmip6_historical", "GFDL-ESM4",      "1971-01", "2000-12")
mpi_hist_mol,  mpi_hist  = model_oxygen("cmip6_historical", "MPI-ESM1-2-LR",  "1971-01", "2000-12")
cesm_ssp_mol,  cesm_ssp  = model_oxygen("cmip6_ssp245",    "CESM2",           "2071-01", "2100-12")
gfdl_ssp_mol,  gfdl_ssp  = model_oxygen("cmip6_ssp245",    "GFDL-ESM4",       "2071-01", "2100-12")
mpi_ssp_mol,   mpi_ssp   = model_oxygen("cmip6_ssp245",    "MPI-ESM1-2-LR",   "2071-01", "2100-12")

# model cell volumes / areas / vertical grids (requested; attrs already let the tools rebuild them)
model_grid = {}
for ds, mdl in [("cmip6_historical", "GFDL-ESM4"),
                ("cmip6_historical", "MPI-ESM1-2-LR"),
                ("cmip6_ssp245", "CESM2"),
                ("cmip6_ssp245", "GFDL-ESM4"),
                ("cmip6_ssp245", "MPI-ESM1-2-LR")]:
    model_grid[(ds, mdl)] = dict(
        volcello=load_volumes_for_cmip6(dataset=ds, model=mdl),
        areacello=load_areas(dataset=ds, model=mdl),
        levels=load_depth_levels(dataset=ds, model=mdl, variable="o2"),
    )

# =====================================================================
# 3. OMZ volumes (strict <20 and <60 umol/kg), maxima and thickness maps
# =====================================================================
def omz_volumes(field, thr):
    r = compute_exceedance_fraction(field, threshold=thr, comparison="less",
                                    dim=None, weighting="volume")
    return float(r["exceeding_amount"]), float(r["total_amount"])

sources = [
    ("WOA23",              o2_woa),
    ("GFDL-ESM4 hist",     gfdl_hist),
    ("MPI-ESM1-2-LR hist", mpi_hist),
    ("CESM2 ssp245",       cesm_ssp),
    ("GFDL-ESM4 ssp245",   gfdl_ssp),
    ("MPI-ESM1-2-LR ssp245", mpi_ssp),
]

raw_volumes     = {}   # m3
volumes_1e15    = {}   # 1e15 m3
thickness_fields = {}  # (depth-reduced) OMZ thickness at 60 umol/kg
thickness_max   = {}

for name, fld in sources:
    v20, tot20 = omz_volumes(fld, 20.0)
    v60, tot60 = omz_volumes(fld, 60.0)
    raw_volumes[name] = {"lt20": v20, "lt60": v60, "total": tot60}
    volumes_1e15[name] = {"lt20": v20 / 1e15, "lt60": v60 / 1e15, "total": tot60 / 1e15}

    th = identify_oxygen_minimum_zone(fld, threshold=60.0)[0].compute()
    th = th.assign_attrs(**fld.attrs)
    thickness_fields[name] = th
    thickness_max[name] = locate_extreme(th, op="max")

# =====================================================================
# 4. Historical -> ssp245 percentage change (unrounded volumes)
# =====================================================================
def pct_change(hist_name, ssp_name, thr):
    v0 = raw_volumes[hist_name][thr]
    v1 = raw_volumes[ssp_name][thr]
    return 100.0 * (v1 - v0) / v0

pct = {}
for mdl in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    pct[mdl] = {
        "lt20": pct_change(f"{mdl} hist", f"{mdl} ssp245", "lt20"),
        "lt60": pct_change(f"{mdl} hist", f"{mdl} ssp245", "lt60"),
    }

# =====================================================================
# 5. Total dissolved oxygen inventory (moles)
# =====================================================================
def oxygen_inventory(o2_mol_mean):
    content = convert_concentration_to_content(o2_mol_mean, density=1025.0)
    return float(compute_basic_stat(content, op="sum", dim=None).compute())

inventory = {
    "GFDL-ESM4 hist":       oxygen_inventory(gfdl_hist_mol),
    "MPI-ESM1-2-LR hist":   oxygen_inventory(mpi_hist_mol),
    "CESM2 ssp245":         oxygen_inventory(cesm_ssp_mol),
    "GFDL-ESM4 ssp245":     oxygen_inventory(gfdl_ssp_mol),
    "MPI-ESM1-2-LR ssp245": oxygen_inventory(mpi_ssp_mol),
}

# =====================================================================
# 6. Basin-wise OMZ volume at 60 umol/kg (Pacific, Indian, Atlantic)
# =====================================================================
basin_names = ["Pacific", "Indian", "Atlantic"]

def basin_omz(field, thr=60.0):
    out = {}
    for b in basin_names:
        fb = select_basin(field, basin=b, include_marginal_seas=True, depth_dependent=False)
        v, t = omz_volumes(fb, thr)
        out[b] = {"omz": v / 1e15, "total": t / 1e15}
    return out

basin_volumes = {
    "WOA23":              basin_omz(o2_woa),
    "GFDL-ESM4 hist":     basin_omz(gfdl_hist),
    "MPI-ESM1-2-LR hist": basin_omz(mpi_hist),
}

# =====================================================================
# 7. WOA23 AOU and O2 saturation inside the <60 umol/kg OMZ
# =====================================================================
omz_mask_woa = (o2_woa < 60.0)

aou_omz_mean = float(compute_basic_stat(aou_woa.where(omz_mask_woa),
                                        op="volume_weighted_mean", dim=None).compute())
po2_omz_mean = float(compute_basic_stat(po2_woa.where(omz_mask_woa),
                                        op="volume_weighted_mean", dim=None).compute())

# =====================================================================
# 8. Regular-grid versions for the figures
# =====================================================================
def to_regular(da, name):
    if name == "WOA23":
        return da
    return regrid_horizontal(da, target="woa23", method="bilinear").compute()

th_reg = {name: to_regular(th, name) for name, th in thickness_fields.items()}

diff_gfdl = th_reg["GFDL-ESM4 ssp245"] - th_reg["GFDL-ESM4 hist"]
diff_mpi  = th_reg["MPI-ESM1-2-LR ssp245"] - th_reg["MPI-ESM1-2-LR hist"]

# =====================================================================
# 9. Report
# =====================================================================
for name in volumes_1e15:
    e = thickness_max[name]
    print(f"{name:22s}  <20: {volumes_1e15[name]['lt20']:.3g}e15 m3   "
          f"<60: {volumes_1e15[name]['lt60']:.3g}e15 m3   "
          f"max thickness {e['value']:.1f} m at "
          f"{e['coordinates'].get('latitude')}, {e['coordinates'].get('longitude')}")
for mdl in pct:
    print(f"{mdl}: dV20 = {pct[mdl]['lt20']:+.2f}%   dV60 = {pct[mdl]['lt60']:+.2f}%")
for k, v in inventory.items():
    print(f"Inventory {k:22s}: {v:.3g} mol")
for k, d in basin_volumes.items():
    print(k, {b: f"{d[b]['omz']:.3g}e15 m3" for b in basin_names})
print(f"WOA23 OMZ AOU = {aou_omz_mean:.1f} umol/kg, O2 saturation = {po2_omz_mean:.1f} %")
