# ---------------------------------------------------------------
# Evaluation of CMIP6 ocean carbon fields against GLODAPv2.2016b
# ---------------------------------------------------------------
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
depths_profile = [0.0, 100.0, 500.0, 1000.0, 2000.0]

def _scalar(x):
    """Materialise a scalar DataArray / numpy value as a Python float."""
    if hasattr(x, "compute"):
        x = x.compute()
    return float(x)

# ----------------------------- GLODAP ---------------------------
glodap_tco2 = load_climatology_data(dataset="glodap", variable="TCO2")
glodap_talk = load_climatology_data(dataset="glodap", variable="TAlk")
glodap_ph = load_climatology_data(dataset="glodap", variable="pHtsinsitutp")

glodap_levels, glodap_bounds = load_depth_levels(
    dataset="glodap", depth_range=[0, 50]
)
glodap_areas = load_areas(dataset="glodap")

# Thickness-weighted 0-50 m mean for GLODAP, then area-weighted globally.
def glodap_0_50_field(da):
    da0 = da.sel(depth=slice(0, 50)).sel(depth=glodap_levels)
    thick = xr.DataArray(
        glodap_bounds[:, 1] - glodap_bounds[:, 0],
        dims="depth",
        coords={"depth": glodap_levels},
    )
    thick = thick.sel(depth=da0.depth)
    num = (da0 * thick).sum("depth", skipna=True)
    den = thick.where(da0.notnull()).sum("depth", skipna=True)
    out = num / den
    out = out.assign_attrs(**da.attrs)
    return out

glodap_dic_field = glodap_0_50_field(glodap_tco2)
glodap_talk_field = glodap_0_50_field(glodap_talk)
glodap_ph_field = glodap_0_50_field(glodap_ph)

# Unweighted level mean for maps (each cell's own 0-50 m levels).
glodap_dic_map = compute_basic_stat(
    glodap_tco2.sel(depth=slice(0, 50)), op="mean", dim="depth"
)
glodap_talk_map = compute_basic_stat(
    glodap_talk.sel(depth=slice(0, 50)), op="mean", dim="depth"
)
glodap_ph_map = compute_basic_stat(
    glodap_ph.sel(depth=slice(0, 50)), op="mean", dim="depth"
)

glodap_dic_mean = _scalar(
    compute_basic_stat(glodap_dic_field, op="area_weighted_mean", dim=["lat", "lon"])
)
glodap_talk_mean = _scalar(
    compute_basic_stat(glodap_talk_field, op="area_weighted_mean", dim=["lat", "lon"])
)
glodap_ph_mean = _scalar(
    compute_basic_stat(glodap_ph_field, op="area_weighted_mean", dim=["lat", "lon"])
)

global_means = {
    "GLODAP": {
        "DIC": glodap_dic_mean,
        "TAlk": glodap_talk_mean,
        "pH": glodap_ph_mean,
        "DIC_bias": 0.0,
        "TAlk_bias": 0.0,
        "pH_bias": 0.0,
    }
}

# Store small / lazy products for plotting and reporting.
dic_bias_maps = {}
map_metrics = {}
profile_bias = {}
period_changes = {}
model_areas = {}
model_volumes = {}

# ----------------------------- models ---------------------------
for model in models:
    # The loader supplies the model's own areacello and volcello provenance.
    model_areas[model] = load_areas(dataset="cmip6_historical", model=model)
    model_volumes[model] = load_volumes_for_cmip6(
        dataset="cmip6_historical", model=model, variable="thetao"
    )

    # ------------------- 2002 monthly fields -------------------
    dissic = load_regular_data(
        dataset="cmip6_historical",
        variable="dissic",
        model=model,
        time_range=["2002-01", "2002-12"],
    )
    talk = load_regular_data(
        dataset="cmip6_historical",
        variable="talk",
        model=model,
        time_range=["2002-01", "2002-12"],
    )
    ph = load_regular_data(
        dataset="cmip6_historical",
        variable="ph",
        model=model,
        time_range=["2002-01", "2002-12"],
    )
    thetao = load_regular_data(
        dataset="cmip6_historical",
        variable="thetao",
        model=model,
        time_range=["2002-01", "2002-12"],
    )
    so = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=model,
        time_range=["2002-01", "2002-12"],
    )

    sigma = compute_potential_density_anomaly(
        temp_data=thetao, sal_data=so, ref_pressure=0
    )

    dissic_umol = transform_unit(
        dissic,
        from_unit=dissic.attrs.get("units", "mol m-3"),
        to_unit="umol kg-1",
        density=sigma,
    )
    talk_umol = transform_unit(
        talk,
        from_unit=talk.attrs.get("units", "mol m-3"),
        to_unit="umol kg-1",
        density=sigma,
    )

    # Day-weighted 2002 annual means.
    dissic_ann = temporal_aggregate(dissic_umol, freq="mean", op="mean")
    talk_ann = temporal_aggregate(talk_umol, freq="mean", op="mean")
    ph_ann = temporal_aggregate(ph, freq="mean", op="mean")

    dissic_0_50 = dissic_ann.sel(depth=slice(0, 50))
    talk_0_50 = talk_ann.sel(depth=slice(0, 50))
    ph_0_50 = ph_ann.sel(depth=slice(0, 50))

    # Volume-weighted global means over the model's own 0-50 m levels.
    dic_mean_da = compute_basic_stat(
        dissic_0_50,
        op="volume_weighted_mean",
        dim=["depth", "lat", "lon"],
    )
    talk_mean_da = compute_basic_stat(
        talk_0_50,
        op="volume_weighted_mean",
        dim=["depth", "lat", "lon"],
    )
    ph_mean_da = compute_basic_stat(
        ph_0_50,
        op="volume_weighted_mean",
        dim=["depth", "lat", "lon"],
    )

    dic_mean_val = _scalar(dic_mean_da)
    talk_mean_val = _scalar(talk_mean_da)
    ph_mean_val = _scalar(ph_mean_da)

    global_means[model] = {
        "DIC": dic_mean_val,
        "TAlk": talk_mean_val,
        "pH": ph_mean_val,
        "DIC_bias": dic_mean_val - glodap_dic_mean,
        "TAlk_bias": talk_mean_val - glodap_talk_mean,
        "pH_bias": ph_mean_val - glodap_ph_mean,
    }

    # ---------------- 0-50 m DIC map comparison ----------------
    dic_map_model = compute_basic_stat(dissic_0_50, op="mean", dim="depth")
    dic_map_regrid = regrid_horizontal(
        dic_map_model, target="glodap", method="conservative"
    )

    mask_map = build_comparison_mask(
        [dic_map_regrid, glodap_dic_map], ocean_only=True
    )
    model_map_masked = dic_map_regrid.where(mask_map)
    glodap_map_masked = glodap_dic_map.where(mask_map)

    metrics_map = compute_difference_metrics(
        data=model_map_masked,
        reference=glodap_map_masked,
        dim=["lat", "lon"],
        weighting="area",
    )
    corr_map = analyze_correlation(
        var1=model_map_masked,
        var2=glodap_map_masked,
        dim=["lat", "lon"],
        method="pearson",
        weighting="area",
    )

    map_metrics[model] = {
        "bias": _scalar(metrics_map["bias"]),
        "rmse": _scalar(metrics_map["rmse"]),
        "correlation": _scalar(corr_map["correlation"]),
    }
    dic_bias_maps[model] = metrics_map["difference"].compute()

    # ---------------- vertical DIC bias profile ----------------
    dic_vert = regrid_vertical(dissic_ann, target="glodap", method="linear")
    dic_vert = dic_vert.sel(depth=depths_profile)
    dic_vert = regrid_horizontal(dic_vert, target="glodap", method="conservative")

    prof = []
    for dep in depths_profile:
        model_d = dic_vert.sel(depth=dep)
        glodap_d = glodap_tco2.sel(depth=dep)
        mask_d = build_comparison_mask([model_d, glodap_d], ocean_only=True)
        metrics_d = compute_difference_metrics(
            data=model_d.where(mask_d),
            reference=glodap_d.where(mask_d),
            dim=["lat", "lon"],
            weighting="area",
        )
        prof.append((dep, _scalar(metrics_d["bias"])))
    profile_bias[model] = prof

    # --------------- 1971-2000 DIC reference mean ---------------
    dissic_hist = load_regular_data(
        dataset="cmip6_historical",
        variable="dissic",
        model=model,
        time_range=["1971-01", "2000-12"],
    )
    thetao_hist = load_regular_data(
        dataset="cmip6_historical",
        variable="thetao",
        model=model,
        time_range=["1971-01", "2000-12"],
    )
    so_hist = load_regular_data(
        dataset="cmip6_historical",
        variable="so",
        model=model,
        time_range=["1971-01", "2000-12"],
    )

    sigma_hist = compute_potential_density_anomaly(
        temp_data=thetao_hist, sal_data=so_hist, ref_pressure=0
    )
    dissic_hist_umol = transform_unit(
        dissic_hist,
        from_unit=dissic_hist.attrs.get("units", "mol m-3"),
        to_unit="umol kg-1",
        density=sigma_hist,
    )
    dissic_hist_ann = temporal_aggregate(
        dissic_hist_umol, freq="mean", op="mean"
    )
    dissic_hist_0_50 = dissic_hist_ann.sel(depth=slice(0, 50))
    dic_hist_mean_da = compute_basic_stat(
        dissic_hist_0_50,
        op="volume_weighted_mean",
        dim=["depth", "lat", "lon"],
    )

    change_val = _scalar(dic_mean_da - dic_hist_mean_da)
    bias_val = global_means[model]["DIC_bias"]
    ratio_val = change_val / bias_val if bias_val != 0 else np.nan
    period_changes[model] = {
        "change_2002_minus_1971_2000": change_val,
        "ratio_to_glodap_bias": ratio_val,
    }

# ----------------------------- summary --------------------------
print("Global 0-50 m means and differences from GLODAP")
print("source        DIC(umol/kg)  TAlk(umol/kg)  pH       DIC_bias  TAlk_bias  pH_bias")
for src in ["GLODAP"] + models:
    g = global_means[src]
    print(
        f"{src:12s} {g['DIC']:12.1f} {g['TAlk']:13.1f} {g['pH']:8.3f} "
        f"{g['DIC_bias']:9.1f} {g['TAlk_bias']:10.1f} {g['pH_bias']:8.3f}"
    )

print("\n0-50 m DIC map comparison on the GLODAP 1-degree grid")
print("model          bias(umol/kg)  RMSE(umol/kg)  correlation")
for model in models:
    m = map_metrics[model]
    print(
        f"{model:12s} {m['bias']:12.1f} {m['rmse']:13.1f} {m['correlation']:12.3f}"
    )

print("\nVertical global-mean DIC bias (model minus GLODAP), umol kg-1")
print("model          " + "  ".join(f"{d:>6.0f} m" for d in depths_profile))
for model in models:
    vals = [v for _, v in profile_bias[model]]
    print(f"{model:12s} " + "  ".join(f"{v:8.1f}" for v in vals))

print("\nDIC period change and ratio to 0-50 m GLODAP bias")
for model in models:
    p = period_changes[model]
    print(
        f"{model:12s} change={p['change_2002_minus_1971_2000']:8.1f} umol kg-1, "
        f"ratio={p['ratio_to_glodap_bias']:6.2f}"
    )
