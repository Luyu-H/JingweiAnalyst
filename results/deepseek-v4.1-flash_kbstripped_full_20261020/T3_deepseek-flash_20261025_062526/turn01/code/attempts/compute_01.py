models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]

def to_degc(da):
    u = str(da.attrs.get("units", "")).strip()
    if u == "":
        return da
    ul = u.lower().replace(" ", "")
    if ul in {"degc", "degreecelsius", "degreescelsius", "celsius", "°c", "c"}:
        return da
    if ul in {"k", "kelvin", "degreesk"}:
        return transform_unit(da, from_unit="K", to_unit="degC")
    return transform_unit(da, from_unit=u, to_unit="degC")

model_clim = {}
model_area = {}
native_gm = {}
native_area = {}

for model in models:
    tos = load_regular_data(
        dataset="cmip6_historical", variable="tos", model=model,
        time_range=["1971-01", "2000-12"]
    )
    tos = to_degc(tos)

    clim = temporal_aggregate(data=tos, freq="mean")
    area = load_areas(dataset="cmip6_historical", model=model, variable="tos")

    gm = compute_basic_stat(data=clim, op="area_weighted_mean", dim=["lat", "lon"])
    gm_val = float(gm.compute())

    valid = clim.notnull()
    total_area = float(area.where(valid).sum().compute())

    model_clim[model] = clim
    model_area[model] = area
    native_gm[model] = gm_val
    native_area[model] = total_area

woa = load_climatology_data(
    dataset="woa23", variable="temperature",
    time_slice="annual", period="1971-2000",
    depth_range=[0, 0]
)
woa_sst = to_degc(woa)
if "depth" in woa_sst.dims:
    woa_sst = woa_sst.isel(depth=0, drop=True)
if "depth" in woa_sst.dims:
    woa_sst = woa_sst.squeeze("depth", drop=True)

woa_area = load_areas(dataset="woa23", variable="temperature")
woa_gm = compute_basic_stat(data=woa_sst, op="area_weighted_mean", dim=["lat", "lon"])
woa_global_gm = float(woa_gm.compute())
woa_global_area = float(woa_area.where(woa_sst.notnull()).sum().compute())

model_regrid = {}
for model in models:
    model_regrid[model] = regrid_horizontal(
        model_clim[model], target="woa23", method="conservative"
    )

landsea = load_landsea_mask()
ocean_mask = landsea > 1

fields_for_mask = [model_regrid[m] for m in models] + [woa_sst]
shared_mask = build_comparison_mask(
    fields=fields_for_mask, ocean_only=True, extra_mask=ocean_mask
)
n_comp = int(shared_mask.sum().compute())

model_masked = {m: model_regrid[m].where(shared_mask) for m in models}
woa_masked = woa_sst.where(shared_mask)

metrics = {}
for model in models:
    m = compute_difference_metrics(
        data=model_masked[model], reference=woa_masked,
        dim=None, weighting="auto"
    )
    c = analyze_correlation(
        var1=model_masked[model], var2=woa_masked,
        dim=["lat", "lon"], method="pearson", weighting="area"
    )

    metrics[model] = {
        "mean_data": float(m["mean_data"]),
        "mean_reference": float(m["mean_reference"]),
        "bias": float(m["bias"]),
        "rmse": float(m["rmse"]),
        "max_abs_error": float(m["max_abs_error"]),
        "max_abs_error_at": m["max_abs_error_at"],
        "difference": m["difference"],
        "correlation": float(c["correlation"]),
        "pvalue": float(c["pvalue"]),
    }

abs_delta = {m: abs(metrics[m]["bias"]) for m in models}
ranking = sorted(models, key=lambda m: abs(metrics[m]["bias"]))

bias_map = {m: (model_masked[m] - woa_masked).compute() for m in models}
woa_ref_map = woa_masked.compute()

print("Global native-grid SST means (degC):")
for m in models:
    print(f"  {m}: {native_gm[m]:.2f}  (area {native_area[m]:.3e} m2)")
print(f"  WOA23: {woa_global_gm:.2f}  (area {woa_global_area:.3e} m2)")
print(f"Comparison cells: {n_comp}")
print("Shared-sample metrics:")
for m in models:
    r = metrics[m]
    print(
        f"  {m}: mean={r['mean_data']:.2f}, WOA={r['mean_reference']:.2f}, "
        f"bias={r['bias']:.2f}, RMSE={r['rmse']:.3f}, "
        f"r={r['correlation']:.3f}, max|err|={r['max_abs_error']:.2f} "
        f"at {r['max_abs_error_at']}"
    )
print("Ranking by |bias|:", ranking)
