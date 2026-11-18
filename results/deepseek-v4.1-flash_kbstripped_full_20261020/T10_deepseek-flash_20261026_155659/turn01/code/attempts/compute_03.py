# =====================================================================
# 1.  Observations: WOA23 oxygen, AOU, O2 saturation, masks, geometry
# =====================================================================
woa_o2    = load_climatology_data(dataset="woa23", variable="oxygen", time_slice="annual")
woa_aou   = load_climatology_data(dataset="woa23", variable="apparent_oxygen_utilization", time_slice="annual")
woa_o2sat = load_climatology_data(dataset="woa23", variable="percent_oxygen_saturation", time_slice="annual")

landsea    = load_landsea_mask()
basin_mask = load_basin_mask()

lev_woa, bnd_woa = load_depth_levels(dataset="woa23", variable="oxygen", time_slice="annual")
area_woa         = load_areas(dataset="woa23")

# NOTE (product difference): the WOA23 file in this archive publishes oxygen /
# AOU / percent-saturation as a single all-period composite normal, although the
# question labels them 1971-2000.  We use the available normal and record it.
WOA_O2_PERIOD = "all-period composite (question labels it 1971-2000)"

def mask_keep_attrs(da, mask, template=None):
    """Apply a boolean mask but keep the provenance attrs the weighting /
    thickness tools read (masking with .where alone can strip them)."""
    out = da.where(mask)
    attrs = dict(template.attrs) if template is not None else {}
    attrs.update(da.attrs)
    out.attrs = attrs
    return out

# =====================================================================
# 2.  CMIP6 period means (day-weighted), in-situ density, unit change
# =====================================================================
HIST = ["GFDL-ESM4", "MPI-ESM1-2-LR"]          # CESM2 has no historical o2
SSP  = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

def load_period_mean(dataset, model):
    """Average-first: day-weighted long-term mean over the 30-yr window."""
    tr = ["1971-01", "2000-12"] if dataset == "cmip6_historical" else ["2071-01", "2100-12"]
    o2 = load_regular_data(dataset=dataset, variable="o2",     model=model, time_range=tr)
    th = load_regular_data(dataset=dataset, variable="thetao", model=model, time_range=tr)
    sa = load_regular_data(dataset=dataset, variable="so",     model=model, time_range=tr)

    o2m = temporal_aggregate(o2, freq="mean", op="mean")   # day weighted
    thm = temporal_aggregate(th, freq="mean", op="mean")
    sam = temporal_aggregate(sa, freq="mean", op="mean")

    # in-situ density anomaly (sigma), evaluated at each cell's own pressure
    sigma = compute_potential_density_anomaly(temp_data=thm, sal_data=sam,
                                              ref_pressure=None, temperature_type="in_situ")
    # mol m-3  ->  micromoles per kilogram, cell by cell, using the density
    o2kg = transform_unit(data=o2m, from_unit="mol m-3", to_unit="umol/kg", density=sigma)

    vol = load_volumes_for_cmip6(dataset=dataset, model=model, variable="thetao")
    if "time" in vol.dims:                 # a model publishing monthly cell volume
        vol = vol.mean("time")

    return {"o2_native": o2m, "o2_kg": o2kg, "sigma": sigma, "vol": vol}

cmip = {}
for m in HIST:
    cmip[("historical", m)] = load_period_mean("cmip6_historical", m)
for m in SSP:
    cmip[("ssp245", m)]     = load_period_mean("cmip6_ssp245", m)

DENSITY_ASSUMPTION = {k: v["o2_kg"].attrs.get("density_assumption", "in-situ density field")
                      for k, v in cmip.items()}

# per-model vertical geometry (metres, positive down) - used for reporting only
depth_levels = {}
for m in HIST:
    depth_levels[("historical", m)] = load_depth_levels(dataset="cmip6_historical", model=m, variable="thetao")
for m in SSP:
    depth_levels[("ssp245", m)] = load_depth_levels(dataset="cmip6_ssp245", model=m, variable="thetao")

# =====================================================================
# 3.  The six source-and-experiment combinations
# =====================================================================
labels = ["WOA23"] \
       + [f"{m} historical" for m in HIST] \
       + [f"{m} ssp245" for m in SSP]

o2_native = {"WOA23": woa_o2}
for m in HIST:
    o2_native[f"{m} historical"] = cmip[("historical", m)]["o2_kg"]
for m in SSP:
    o2_native[f"{m} ssp245"]     = cmip[("ssp245", m)]["o2_kg"]

# =====================================================================
# 4.  OMZ volume (<20 and <60 umol/kg) and column thickness (60) on
#     each source's OWN grid / layer geometry
# =====================================================================
def omz_volume(field, threshold):
    """Volume in m3 of water strictly below the threshold."""
    r = compute_exceedance_fraction(data=field, threshold=float(threshold),
                                    comparison="less", dim=None, weighting="volume")
    return float(r["exceeding_amount"]), float(r["fraction"]), float(r["total_amount"])

omz = {}
for lab in labels:
    f = o2_native[lab]
    v20, f20, tot20 = omz_volume(f, 20.0)
    v60, f60, tot60 = omz_volume(f, 60.0)
    thick, mask60 = identify_oxygen_minimum_zone(o2_data=f, threshold=60.0)
    thick = thick.compute()
    omz[lab] = {"v20": v20 / 1e15, "v60": v60 / 1e15,
                "f20": f20, "f60": f60, "valid_vol": tot60,
                "thick_native": thick, "mask60": mask60,
                "max60": locate_extreme(thick, op="max")}

# =====================================================================
# 5.  Global dissolved-oxygen inventory (native mol m-3, volcello)
# =====================================================================
inventory = {}
for m in HIST:
    c = convert_concentration_to_content(concentration=cmip[("historical", m)]["o2_native"],
                                         density=1025.0)
    inventory[f"{m} historical"] = float(c.sum())
for m in SSP:
    c = convert_concentration_to_content(concentration=cmip[("ssp245", m)]["o2_native"],
                                         density=1025.0)
    inventory[f"{m} ssp245"] = float(c.sum())

# =====================================================================
# 6.  Common 1-degree grid + one comparison mask
# =====================================================================
o2_common = {"WOA23": woa_o2}
for m in HIST:
    o2_common[f"{m} historical"] = regrid_horizontal(cmip[("historical", m)]["o2_kg"],
                                                     target="woa23", method="conservative")
for m in SSP:
    o2_common[f"{m} ssp245"] = regrid_horizontal(cmip[("ssp245", m)]["o2_kg"],
                                                 target="woa23", method="conservative")

mask_fields = [o2_common["WOA23"],
               o2_common["GFDL-ESM4 historical"], o2_common["MPI-ESM1-2-LR historical"],
               o2_common["GFDL-ESM4 ssp245"],     o2_common["MPI-ESM1-2-LR ssp245"]]
cmp_mask = build_comparison_mask(fields=mask_fields, ocean_only=True, reduce_other_dims="all")
N_MASK_CELLS = int(np.asarray(cmp_mask.values).sum())
N_MASK_TOTAL = int(np.asarray(cmp_mask.notnull().values).sum())

o2_common = {k: mask_keep_attrs(v, cmp_mask) for k, v in o2_common.items()}

# --- column OMZ thickness on the common grid -------------------------------
# identify_oxygen_minimum_zone builds the layer thickness from the field's OWN
# vertical axis and its provenance attrs.  A horizontally regridded model field
# no longer belongs to a single (variable, time-slice) vertical axis, so the
# zone is identified on each source's native grid (step 4) and only the
# resulting per-column thickness MAP is moved onto the common 1-degree grid.
thick_common = {}
for lab in labels:
    if lab == "WOA23":
        t = omz["WOA23"]["thick_native"]
    else:
        t = regrid_horizontal(omz[lab]["thick_native"], target="woa23", method="bilinear")
    thick_common[lab] = mask_keep_attrs(t, cmp_mask, template=o2_common[lab])

# =====================================================================
# 7.  Basin volumes at the 60 threshold (WOA23 basin labels)
#     computed on each source's OWN grid / cell volumes
# =====================================================================
BASINS = ["Pacific", "Indian", "Atlantic"]
basin_vol = {}
for lab in ["WOA23", "GFDL-ESM4 historical", "MPI-ESM1-2-LR historical"]:
    for b in BASINS:
        sub = select_basin(data=o2_native[lab], basin=b,
                           include_marginal_seas=True, depth_dependent=False)
        r = compute_exceedance_fraction(data=sub, threshold=60.0, comparison="less",
                                        dim=None, weighting="volume")
        basin_vol[(lab, b)] = float(r["exceeding_amount"]) / 1e15

# =====================================================================
# 8.  WOA23 within-OMZ (60 umol/kg) AOU and O2 saturation
# =====================================================================
woa_omz_mask = omz["WOA23"]["mask60"]

aou_in_omz   = mask_keep_attrs(woa_aou,   woa_omz_mask)
o2sat_in_omz = mask_keep_attrs(woa_o2sat, woa_omz_mask)

aou_omz_mean = float(compute_basic_stat(data=aou_in_omz, op="volume_weighted_mean",
                                        dim=["depth", "lat", "lon"]))
o2sat_omz_mean = float(compute_basic_stat(data=o2sat_in_omz, op="volume_weighted_mean",
                                          dim=["depth", "lat", "lon"]))

# =====================================================================
# 9.  Scenario change, model-vs-WOA skill at 60 umol/kg
# =====================================================================
scenario_change = {}
for m in HIST:
    for thr, key in [("20", "v20"), ("60", "v60")]:
        h = omz[f"{m} historical"][key]
        s = omz[f"{m} ssp245"][key]
        scenario_change[(m, thr)] = 100.0 * (s - h) / h       # + = expansion

volume_bias = {}
for m in HIST:
    for key in ["v20", "v60"]:
        d = omz[f"{m} historical"][key] - omz["WOA23"][key]
        volume_bias[(m, key, "abs")] = d
        volume_bias[(m, key, "rel")] = 100.0 * d / omz["WOA23"][key]

diff_metrics = {}
for m in HIST:
    diff_metrics[m] = compute_difference_metrics(
        data=thick_common[f"{m} historical"],
        reference=thick_common["WOA23"],
        dim=["lat", "lon"], weighting="area")

# =====================================================================
# 10.  Consolidated result objects (materialised numbers)
# =====================================================================
report_volume = pd.DataFrame(
    [{"source": lab,
      "V_lt20 (1e15 m3)": float(f"{omz[lab]['v20']:.3g}"),
      "V_lt60 (1e15 m3)": float(f"{omz[lab]['v60']:.3g}"),
      "max_thick60 (m)": float(f"{omz[lab]['max60']['value']:.1f}"),
      "max_lat": float(f"{omz[lab]['max60']['coordinates'].get('latitude', np.nan):.2f}"),
      "max_lon": float(f"{omz[lab]['max60']['coordinates'].get('longitude', np.nan):.2f}")}
     for lab in labels])

report_inventory = {k: float(f"{v:.3g}") for k, v in inventory.items()}
report_scenario  = {f"{m} {t} umol/kg (%)": float(f"{scenario_change[(m, t)]:.2f}")
                    for (m, t) in scenario_change}
report_basin     = {f"{lab} {b}": float(f"{basin_vol[(lab, b)]:.3g}")
                    for (lab, b) in basin_vol}

fig1_panels = [("GFDL-ESM4 historical", thick_common["GFDL-ESM4 historical"]),
               ("MPI-ESM1-2-LR historical", thick_common["MPI-ESM1-2-LR historical"]),
               ("WOA23", thick_common["WOA23"])]

d_thick_gfdl = (thick_common["GFDL-ESM4 ssp245"] - thick_common["GFDL-ESM4 historical"])
d_thick_mpi  = (thick_common["MPI-ESM1-2-LR ssp245"] - thick_common["MPI-ESM1-2-LR historical"])
fig2_panels  = [("GFDL-ESM4 ssp245 - historical", d_thick_gfdl),
                ("MPI-ESM1-2-LR ssp245 - historical", d_thick_mpi)]
cesm2_thick  = thick_common["CESM2 ssp245"]
