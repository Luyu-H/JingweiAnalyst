models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
DS = "cmip6_historical"
T2002 = ["2002-01", "2002-12"]
TBASE = ["1971-01", "2000-12"]
PROF = [0, 100, 500, 1000, 2000]

missing = {}      # (model, variable) -> reason, only for talk / ph, which a model may not publish
provenance = {}


def load_model(m, v, tr, dr, optional=False):
    if not optional:
        return load_regular_data(dataset=DS, variable=v, model=m, time_range=tr, depth_range=dr)
    try:
        return load_regular_data(dataset=DS, variable=v, model=m, time_range=tr, depth_range=dr)
    except Exception as e:  # the task asks for a missing variable to be reported, not substituted
        missing[(m, v)] = repr(e)
        return None


def to_umolkg(conc, thetao, so):
    """mol m-3 -> umol kg-1 with the in-situ density anomaly (ref_pressure=None)."""
    rho = compute_potential_density_anomaly(thetao, so, ref_pressure=None)
    out = transform_unit(conc, conc.attrs.get("units", "mol m-3"), "umol/kg", density=rho)
    return out, rho


def vol_mean(da):
    return compute_basic_stat(da, op="volume_weighted_mean", dim=["depth", "lat", "lon"])


def thickness_weighted_map(da, thk):
    """Each cell's mean over its own levels, weighted by layer thickness, valid cells only."""
    if thk.sizes.get("depth") == da.sizes["depth"]:
        thk = thk.assign_coords(depth=da["depth"])
    num = (da * thk).sum("depth", skipna=True, min_count=1)
    den = thk.where(da.notnull()).sum("depth")
    out = num / den
    return out.assign_attrs(units=da.attrs.get("units", ""))


# ------------------------------------------------------------------ GLODAP
gl_dic = load_climatology_data(dataset="glodap", variable="TCO2", depth_range=[0, 50])
gl_talk = load_climatology_data(dataset="glodap", variable="TAlk", depth_range=[0, 50])
gl_ph = load_climatology_data(dataset="glodap", variable="pHtsinsitutp", depth_range=[0, 50])
gl_levels, gl_bounds = load_depth_levels(dataset="glodap", depth_range=[0, 50], variable="TCO2")
gl_thk_vals = np.asarray(gl_bounds)[:, 1] - np.asarray(gl_bounds)[:, 0]
assert len(gl_thk_vals) == gl_dic.sizes["depth"]
gl_thk = xr.DataArray(gl_thk_vals, dims="depth", coords={"depth": gl_dic["depth"].values},
                      attrs={"units": "m"})
gl_area = load_areas()   # shared 1-degree grid
provenance["GLODAP"] = {"area_source": gl_area.attrs.get("area_source"),
                        "levels_0_50m": [float(x) for x in np.asarray(gl_levels)],
                        "layer_bounds_m": np.asarray(gl_bounds).tolist()}

gl_glob = {
    "DIC": float(vol_mean(gl_dic).compute()),
    "TAlk": float(vol_mean(gl_talk).compute()),
    "pH": float(vol_mean(gl_ph).compute()),
}
gl_dic_map = thickness_weighted_map(gl_dic, gl_thk).compute()

# ------------------------------------------------------------------ models, 0-50 m, 2002
sh = {}       # shallow 2002 annual means (converted)
glob = {}     # global 0-50 m volume weighted means
dens_assump = {}
thk_m = {}
for m in models:
    dic = load_model(m, "dissic", T2002, [0, 50])
    thetao = load_model(m, "thetao", T2002, [0, 50])
    so = load_model(m, "so", T2002, [0, 50])
    talk = load_model(m, "talk", T2002, [0, 50], optional=True)
    ph = load_model(m, "ph", T2002, [0, 50], optional=True)

    dic_a = temporal_aggregate(dic, freq="mean")
    thetao_a = temporal_aggregate(thetao, freq="mean")
    so_a = temporal_aggregate(so, freq="mean")
    dic_u, rho = to_umolkg(dic_a, thetao_a, so_a)
    dens_assump[m] = dic_u.attrs.get("density_assumption")
    entry = {"dic": dic_u}
    glob[m] = {"DIC": float(vol_mean(dic_u).compute())}
    if talk is not None:
        talk_u, _ = to_umolkg(temporal_aggregate(talk, freq="mean"), thetao_a, so_a)
        entry["talk"] = talk_u
        glob[m]["TAlk"] = float(vol_mean(talk_u).compute())
    else:
        glob[m]["TAlk"] = np.nan
    if ph is not None:
        ph_a = temporal_aggregate(ph, freq="mean")
        entry["ph"] = ph_a
        glob[m]["pH"] = float(vol_mean(ph_a).compute())
    else:
        glob[m]["pH"] = np.nan
    sh[m] = entry

    # static grid fields: thickness = volcello / areacello
    vol = load_volumes_for_cmip6(dataset=DS, model=m, depth_range=[0, 50])
    if "time" in vol.dims:
        vol = vol.sel(time=slice("2002-01", "2002-12")).mean("time")
    area = load_areas(dataset=DS, model=m)
    thk_m[m] = (vol / area).assign_attrs(units="m")
    provenance[m] = {"area_source": area.attrs.get("area_source"),
                     "volume_source": vol.attrs.get("volume_source"),
                     "density_assumption": dens_assump[m]}

diff_glob = {m: {k: glob[m][k] - gl_glob[k] for k in ["DIC", "TAlk", "pH"]} for m in models}
glob_rounded = {m: {"DIC": round(glob[m]["DIC"]), "TAlk": (round(glob[m]["TAlk"]) if np.isfinite(glob[m]["TAlk"]) else np.nan),
                    "pH": round(glob[m]["pH"], 3)} for m in models}
gl_glob_rounded = {"DIC": round(gl_glob["DIC"]), "TAlk": round(gl_glob["TAlk"]), "pH": round(gl_glob["pH"], 3)}
diff_glob_rounded = {m: {"DIC": round(diff_glob[m]["DIC"]), "TAlk": (round(diff_glob[m]["TAlk"]) if np.isfinite(diff_glob[m]["TAlk"]) else np.nan),
                         "pH": round(diff_glob[m]["pH"], 3)} for m in models}

# ------------------------------------------------------------------ 0-50 m DIC map on GLODAP grid
reg_maps = {}
for m in models:
    mp = thickness_weighted_map(sh[m]["dic"], thk_m[m]).compute()
    reg_maps[m] = regrid_horizontal(mp, target="glodap", method="conservative").compute()

map_mask = build_comparison_mask([reg_maps[m] for m in models] + [gl_dic_map], ocean_only=False)
map_mask_counts = {"common_cells": int(map_mask.sum()), "attrs": dict(map_mask.attrs)}

gl_map_c = gl_dic_map.where(map_mask)
map_metrics = {}
bias_maps = {}
for m in models:
    mm = reg_maps[m].where(map_mask)
    met = compute_difference_metrics(mm, gl_map_c)
    corr = analyze_correlation(mm, gl_map_c, dim=["lat", "lon"], weighting="area")
    r = float(corr["correlation"].compute())
    map_metrics[m] = {"bias": float(met["bias"]), "rmse": float(met["rmse"]), "corr": r,
                      "mean_model": float(met["mean_data"]), "mean_glodap": float(met["mean_reference"]),
                      "n_valid": met["n_valid"]}
    bias_maps[m] = met["difference"].compute()
map_metrics_rounded = {m: {"bias": round(v["bias"]), "rmse": round(v["rmse"]), "corr": round(v["corr"], 3)}
                       for m, v in map_metrics.items()}

# ------------------------------------------------------------------ DIC bias profile
gl_dic_deep = load_climatology_data(dataset="glodap", variable="TCO2", depth_range=[0, 2000])
gl_prof = regrid_vertical(gl_dic_deep, target=PROF).compute()

prof_bias = {}
prof_counts = {}
for m in models:
    dic = load_model(m, "dissic", T2002, [0, 2500])
    thetao = load_model(m, "thetao", T2002, [0, 2500])
    so = load_model(m, "so", T2002, [0, 2500])
    thetao_a = temporal_aggregate(thetao, freq="mean")
    so_a = temporal_aggregate(so, freq="mean")
    dic_u, _ = to_umolkg(temporal_aggregate(dic, freq="mean"), thetao_a, so_a)
    dic_p = regrid_vertical(dic_u, target=PROF)
    dic_pg = regrid_horizontal(dic_p, target="glodap", method="conservative").compute()
    b, c = [], []
    for i, d in enumerate(PROF):
        md = dic_pg.isel(depth=i)
        gd = gl_prof.isel(depth=i)
        mk = build_comparison_mask([md, gd], ocean_only=False)
        met = compute_difference_metrics(md.where(mk), gd.where(mk))
        b.append(float(met["bias"]))
        c.append(int(mk.sum()))
    prof_bias[m] = b
    prof_counts[m] = c
prof_bias_rounded = {m: [round(x) for x in v] for m, v in prof_bias.items()}

# ------------------------------------------------------------------ period change (0-50 m DIC)
period_change = {}
for m in models:
    dic = load_model(m, "dissic", TBASE, [0, 50])
    thetao = load_model(m, "thetao", TBASE, [0, 50])
    so = load_model(m, "so", TBASE, [0, 50])
    dic_b, _ = to_umolkg(temporal_aggregate(dic, freq="mean"),
                         temporal_aggregate(thetao, freq="mean"),
                         temporal_aggregate(so, freq="mean"))
    base_mean = float(vol_mean(dic_b).compute())
    change = glob[m]["DIC"] - base_mean
    ratio = change / diff_glob[m]["DIC"]
    period_change[m] = {"mean_1971_2000": base_mean, "mean_2002": glob[m]["DIC"],
                        "change": change, "change_rounded": round(change),
                        "ratio_to_bias": round(ratio, 2),
                        "baseline": "model's own 1971-2000 climatology (does not remove model bias)"}

# ------------------------------------------------------------------ synthesis
def rank_vals(vals, higher_better=False):
    v = np.array(vals, dtype=float)
    key = -v if higher_better else v
    r = np.full(len(v), np.nan)
    ok = np.isfinite(key)
    order = np.argsort(key[ok])
    rr = np.empty(ok.sum())
    rr[order] = np.arange(1, ok.sum() + 1)
    r[ok] = rr
    return r

criteria = {
    "|dDIC 0-50m|": rank_vals([abs(diff_glob[m]["DIC"]) for m in models]),
    "|dTAlk 0-50m|": rank_vals([abs(diff_glob[m]["TAlk"]) if np.isfinite(diff_glob[m]["TAlk"]) else np.nan for m in models]),
    "|dpH 0-50m|": rank_vals([abs(diff_glob[m]["pH"]) if np.isfinite(diff_glob[m]["pH"]) else np.nan for m in models]),
    "map RMSE": rank_vals([map_metrics[m]["rmse"] for m in models]),
    "map corr": rank_vals([map_metrics[m]["corr"] for m in models], higher_better=True),
}
mean_rank = {m: float(np.nanmean([criteria[k][i] for k in criteria])) for i, m in enumerate(models)}
best_model = min(mean_rank, key=mean_rank.get)
sign_summary = {k: {m: ("+" if diff_glob[m][k] > 0 else "-") if np.isfinite(diff_glob[m][k]) else "missing" for m in models}
                for k in ["DIC", "TAlk", "pH"]}
signs_consistent = {k: len({s for s in v.values() if s != "missing"}) == 1 for k, v in sign_summary.items()}

results = {
    "glodap_global_0_50m": gl_glob_rounded,
    "model_global_0_50m": glob_rounded,
    "model_minus_glodap_0_50m": diff_glob_rounded,
    "map_mask_counts": map_mask_counts,
    "map_metrics": map_metrics_rounded,
    "map_metrics_same_sample": map_metrics,
    "profile_depths_m": PROF,
    "profile_bias": prof_bias_rounded,
    "profile_cell_counts": prof_counts,
    "period_change": period_change,
    "mean_rank": mean_rank,
    "best_model": best_model,
    "sign_summary": sign_summary,
    "signs_consistent": signs_consistent,
    "missing_variables": missing,
    "provenance": provenance,
}
for k, v in results.items():
    print(k, ":", v)
