# Computation code only
# Open the raw CMIP6 files, compute the density-criterion MLD from each model's
# own thetao/so, compare it with the model's published mlotst, and store metrics.

import os
import glob
from collections import Counter
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import gsw

root = "/data1/lorraine-20419/benchmark_cmip6/historical"

models = {
    "CESM2": {
        "grid_type": "curvilinear",
        "grid_label": "gn",
        "y_dim": "nlat",
        "x_dim": "nlon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 0.01,
        "calendar": "noleap",
        "ocean_mask_from": "sftof",
    },
    "GFDL-ESM4": {
        "grid_type": "regular",
        "grid_label": "gr",
        "y_dim": "lat",
        "x_dim": "lon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 1.0,
        "calendar": "noleap",
    },
    "MPI-ESM1-2-LR": {
        "grid_type": "curvilinear",
        "grid_label": "gn",
        "y_dim": "j",
        "x_dim": "i",
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "depth_coord": "lev",
        "depth_scale": 1.0,
        "calendar": "proleptic_gregorian",
    },
}

model_order = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
months = ["March", "September"]


def files_for(model, var):
    """Return files for one model/variable, using the variant parsed from filenames."""
    d = os.path.join(root, model, var)
    files = sorted(glob.glob(os.path.join(d, f"{var}_*.nc")))
    if not files:
        raise FileNotFoundError(f"No files found for {model} {var} under {d}")

    variants = []
    for f in files:
        parts = os.path.basename(f).split("_")
        if len(parts) < 5:
            raise ValueError(f"Cannot parse variant from filename: {f}")
        variants.append(parts[4])

    # Do not assume r1i1p1f1.  Pick the variant actually present.
    variant = Counter(variants).most_common(1)[0][0]
    files = [f for f, v in zip(files, variants) if v == variant]
    return files, variant


def open_omon(model, var):
    files, variant = files_for(model, var)
    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        use_cftime=True,
        decode_times=True,
        chunks={"time": 1},
    )
    return ds, variant


def select_1971_2000(ds):
    t = ds.time
    sel = (t.dt.year >= 1971) & (t.dt.year <= 2000)
    idx = np.where(sel.values)[0]
    return ds.isel(time=idx)


def open_ofx(model, var):
    files, variant = files_for(model, var)
    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True, decode_times=True)
    else:
        ds = xr.open_mfdataset(
            files, combine="by_coords", use_cftime=True, decode_times=True
        )
    return ds, variant


def get_lat_lon_2d(ds, meta):
    lat_da, lon_da = xr.broadcast(
        ds[meta["lat_coord"]], ds[meta["lon_coord"]]
    )
    lat2d = lat_da.transpose(meta["y_dim"], meta["x_dim"]).values.astype("float64")
    lon2d = lon_da.transpose(meta["y_dim"], meta["x_dim"]).values.astype("float64")
    return lat2d, lon2d


def compute_operator_mld(thetao, so, depth_m, lat2d, lon2d):
    """
    Density-criterion MLD from one monthly thetao/so profile.

    thetao, so : (nlev, ny, nx), thetao is potential temperature and so is
                 practical salinity.  Density is potential density referenced
                 to the shallowest level in each column.  The crossing of
                 rho = rho(ref) + 0.03 kg m-3 is linearly interpolated.
    """
    nlev, ny, nx = thetao.shape

    # GSW expects longitude in a conventional range.
    lon_gsw = ((lon2d + 180.0) % 360.0) - 180.0

    z = -np.broadcast_to(depth_m[:, None, None], (nlev, ny, nx))
    lat3 = np.broadcast_to(lat2d[None, :, :], (nlev, ny, nx))
    lon3 = np.broadcast_to(lon_gsw[None, :, :], (nlev, ny, nx))

    p = gsw.p_from_z(z, lat3)

    # so is practical salinity (units "0.001", values ~35).
    SA = gsw.SA_from_SP(so, p, lon3, lat3)
    CT = gsw.CT_from_pt(SA, thetao)

    # Potential density referenced to the shallowest level in each column.
    p_ref = p[0, :, :]
    p_ref3 = np.broadcast_to(p_ref[None, :, :], (nlev, ny, nx))
    rho = gsw.rho(SA, CT, p_ref3)

    ref_rho = rho[0, :, :]
    threshold = ref_rho + 0.03

    cond = rho > threshold[None, :, :]
    has_crossing = np.any(cond, axis=0)
    idx = np.argmax(cond, axis=0)

    idx_exp = idx[None, :, :]
    rho_idx = np.take_along_axis(rho, idx_exp, axis=0)[0]

    i_prev = np.maximum(idx - 1, 0)
    rho_prev = np.take_along_axis(rho, i_prev[None, :, :], axis=0)[0]

    z3 = np.broadcast_to(depth_m[:, None, None], (nlev, ny, nx))
    z_idx = np.take_along_axis(z3, idx_exp, axis=0)[0]
    z_prev = np.take_along_axis(z3, i_prev[None, :, :], axis=0)[0]

    den = rho_idx - rho_prev
    safe_den = np.where(np.abs(den) > 1e-12, den, np.nan)

    with np.errstate(invalid="ignore", divide="ignore"):
        cross = z_prev + (threshold - rho_prev) * (z_idx - z_prev) / safe_den

    cross = np.where(idx == 0, z_idx, cross)
    cross = np.where(has_crossing, cross, np.nan)

    good_surface = np.isfinite(thetao[0, :, :]) & np.isfinite(so[0, :, :])
    cross = np.where(good_surface, cross, np.nan)

    return cross.astype("float32")


def build_comparison_mask(op, mlt):
    """True only where both fields are finite in every collected month."""
    return np.all(np.isfinite(op) & np.isfinite(mlt), axis=0)


def compute_difference_metrics(op_clim, mlt_clim, mask, area, total_ocean_area):
    diff = op_clim - mlt_clim
    valid = (
        mask
        & np.isfinite(op_clim)
        & np.isfinite(mlt_clim)
        & np.isfinite(diff)
        & np.isfinite(area)
    )

    w = area[valid].astype("float64")
    x = op_clim[valid].astype("float64")
    y = mlt_clim[valid].astype("float64")
    d = diff[valid].astype("float64")

    W = np.sum(w)
    bias = np.sum(w * d) / W
    rms = np.sqrt(np.sum(w * d * d) / W)

    mx = np.sum(w * x) / W
    my = np.sum(w * y) / W
    cov = np.sum(w * (x - mx) * (y - my)) / W
    vx = np.sum(w * (x - mx) ** 2) / W
    vy = np.sum(w * (y - my) ** 2) / W
    corr = cov / np.sqrt(vx * vy)

    exceed_area = np.sum(w[np.abs(d) > 20.0])
    exceed_frac = 100.0 * exceed_area / W

    shared_area = W
    shared_pct = 100.0 * W / total_ocean_area

    return {
        "bias": float(bias),
        "rms": float(rms),
        "corr": float(corr),
        "exceed_area": float(exceed_area),
        "exceed_frac": float(exceed_frac),
        "shared_area": float(shared_area),
        "shared_pct": float(shared_pct),
        "cell_count": int(np.sum(mask)),
        "n_valid": int(np.sum(valid)),
    }


def summarize_exceedance(diff, mask, area, lat2d, lon2d):
    """Area-weighted exceedance fractions by latitude band and named basin."""
    valid = mask & np.isfinite(diff) & np.isfinite(area)
    ad = np.abs(diff)
    out = {}

    for low in np.arange(-90.0, 90.0, 10.0):
        high = low + 10.0
        sel = valid & (lat2d >= low) & (lat2d < high)
        w = np.sum(area[sel])
        if w > 0:
            frac = 100.0 * np.sum(area[sel & (ad > 20.0)]) / w
            out[f"lat {low:+.0f} to {high:+.0f}"] = float(frac)

    lon = lon2d % 360.0
    basins = {
        "North Atlantic": (lat2d >= 20.0) & (lat2d <= 80.0) & (lon >= 280.0),
        "North Pacific": (lat2d >= 20.0)
        & (lat2d <= 80.0)
        & (lon >= 100.0)
        & (lon < 280.0),
        "Southern Ocean": lat2d < -40.0,
        "Equatorial band": (lat2d >= -20.0) & (lat2d <= 20.0),
    }
    for name, sel0 in basins.items():
        sel = valid & sel0
        w = np.sum(area[sel])
        if w > 0:
            frac = 100.0 * np.sum(area[sel & (ad > 20.0)]) / w
            out[name] = float(frac)

    return out


def conservative_regrid(field, area, lat2d, lon2d, target_lon, target_lat):
    """
    Conservative area-weighted cell-centre binning to a regular lon/lat grid.
    Source cells are binned by their centres into 1-degree target cells and
    averaged with areacello as the weight.  Used for display only.
    """
    ny, nx = field.shape
    nlat = len(target_lat)
    nlon = len(target_lon)

    lon = lon2d % 360.0
    lat = lat2d

    dlon = target_lon[1] - target_lon[0]
    dlat = target_lat[1] - target_lat[0]
    lon0 = target_lon[0] - 0.5 * dlon
    lat0 = target_lat[0] - 0.5 * dlat

    i = np.floor((lon - lon0) / dlon).astype(int)
    j = np.floor((lat - lat0) / dlat).astype(int)

    valid = (
        np.isfinite(field)
        & np.isfinite(area)
        & (i >= 0)
        & (i < nlon)
        & (j >= 0)
        & (j < nlat)
    )

    out_sum = np.zeros((nlat, nlon), dtype="float64")
    out_w = np.zeros((nlat, nlon), dtype="float64")

    np.add.at(
        out_sum,
        (j[valid], i[valid]),
        field[valid].astype("float64") * area[valid].astype("float64"),
    )
    np.add.at(out_w, (j[valid], i[valid]), area[valid].astype("float64"))

    out = np.full((nlat, nlon), np.nan, dtype="float64")
    good = out_w > 0
    out[good] = out_sum[good] / out_w[good]
    return out


# Common 1-degree display grid.
target_lon_edges = np.arange(0.0, 361.0, 1.0)
target_lat_edges = np.arange(-90.0, 91.0, 1.0)
target_lon = 0.5 * (target_lon_edges[:-1] + target_lon_edges[1:])
target_lat = 0.5 * (target_lat_edges[:-1] + target_lat_edges[1:])

results = {}

for model_name in model_order:
    meta = models[model_name]
    print(f"\nProcessing {model_name} ...")

    ds_t, var_t = open_omon(model_name, "thetao")
    ds_t = select_1971_2000(ds_t)

    ds_s, var_s = open_omon(model_name, "so")
    ds_s = select_1971_2000(ds_s)

    ds_m, var_m = open_omon(model_name, "mlotst")
    ds_m = select_1971_2000(ds_m)

    if len(ds_t.time) != 360:
        raise ValueError(f"{model_name}: expected 360 months, got {len(ds_t.time)}")

    # Align all variables to the thetao time axis.
    ds_s = ds_s.sel(time=ds_t.time)
    ds_m = ds_m.sel(time=ds_t.time)

    lev_vals = ds_t[meta["depth_coord"]].values.astype("float64")
    order = np.argsort(lev_vals)
    depth_m = lev_vals[order] * meta["depth_scale"]

    lat2d, lon2d = get_lat_lon_2d(ds_t, meta)

    ds_a, var_a = open_ofx(model_name, "areacello")
    area_da = ds_a["areacello"].squeeze().transpose(meta["y_dim"], meta["x_dim"])
    area = area_da.values.astype("float64")
    area = np.where(area > 1.0e10, np.nan, area)

    if meta.get("ocean_mask_from") == "sftof":
        ds_sf, var_sf = open_ofx(model_name, "sftof")
        sftof_da = ds_sf["sftof"].squeeze().transpose(meta["y_dim"], meta["x_dim"])
        sftof_vals = sftof_da.values.astype("float64")
        sftof_vals = np.where(sftof_vals < 1.0e10, sftof_vals, 0.0)
        area = np.where(sftof_vals > 0.0, area, np.nan)

    total_ocean_area = float(np.nansum(area))

    time_month = ds_t.time.dt.month.values
    idx_march = np.where(time_month == 3)[0]
    idx_sept = np.where(time_month == 9)[0]

    if len(idx_march) != 30 or len(idx_sept) != 30:
        raise ValueError(
            f"{model_name}: expected 30 Marches and 30 Septembers, got "
            f"{len(idx_march)} and {len(idx_sept)}"
        )

    def collect_operator(indices):
        out = []
        for it in indices:
            th = (
                ds_t["thetao"]
                .isel(time=it)
                .transpose(meta["depth_coord"], meta["y_dim"], meta["x_dim"])
                .values.astype("float64")
            )
            so = (
                ds_s["so"]
                .isel(time=it)
                .transpose(meta["depth_coord"], meta["y_dim"], meta["x_dim"])
                .values.astype("float64")
            )

            th = np.where(np.abs(th) > 1.0e10, np.nan, th)
            so = np.where(np.abs(so) > 1.0e10, np.nan, so)

            th = th[order, :, :]
            so = so[order, :, :]

            op = compute_operator_mld(th, so, depth_m, lat2d, lon2d)
            out.append(op)
        return np.stack(out, axis=0)

    def collect_mlotst(indices):
        out = []
        for it in indices:
            ml = (
                ds_m["mlotst"]
                .isel(time=it)
                .transpose(meta["y_dim"], meta["x_dim"])
                .values.astype("float64")
            )
            ml = np.where(np.abs(ml) > 1.0e10, np.nan, ml)
            out.append(ml.astype("float32"))
        return np.stack(out, axis=0)

    op_by_month = {}
    mlt_by_month = {}
    mask_by_month = {}
    op_clim_by_month = {}
    mlt_clim_by_month = {}
    diff_by_month = {}
    metrics_by_month = {}
    regions_by_month = {}
    diff_maps = {}
    scatter = {}

    for month_name, indices in [("March", idx_march), ("September", idx_sept)]:
        print(f"  {model_name} {month_name}: computing operator MLD for 30 months ...")
        op = collect_operator(indices)
        mlt = collect_mlotst(indices)

        mask = build_comparison_mask(op, mlt)

        op_masked = np.where(mask[None, :, :], op, np.nan)
        mlt_masked = np.where(mask[None, :, :], mlt, np.nan)

        op_clim = np.nanmean(op_masked, axis=0)
        mlt_clim = np.nanmean(mlt_masked, axis=0)
        diff = op_clim - mlt_clim

        met = compute_difference_metrics(
            op_clim, mlt_clim, mask, area, total_ocean_area
        )
        reg = summarize_exceedance(diff, mask, area, lat2d, lon2d)

        op_by_month[month_name] = op
        mlt_by_month[month_name] = mlt
        mask_by_month[month_name] = mask
        op_clim_by_month[month_name] = op_clim
        mlt_clim_by_month[month_name] = mlt_clim
        diff_by_month[month_name] = diff
        metrics_by_month[month_name] = met
        regions_by_month[month_name] = reg

        diff_maps[month_name] = conservative_regrid(
            diff.astype("float64"),
            area,
            lat2d,
            lon2d,
            target_lon,
            target_lat,
        )

        valid_scatter = (
            mask
            & np.isfinite(op_clim)
            & np.isfinite(mlt_clim)
            & np.isfinite(area)
            & (op_clim > 0.0)
            & (mlt_clim > 0.0)
        )
        scatter[month_name] = (
            op_clim[valid_scatter].astype("float64"),
            mlt_clim[valid_scatter].astype("float64"),
        )

    results[model_name] = {
        "variant_thetao": var_t,
        "variant_so": var_s,
        "variant_mlotst": var_m,
        "variant_areacello": var_a,
        "depth_m": depth_m,
        "lat2d": lat2d,
        "lon2d": lon2d,
        "area": area,
        "total_ocean_area": total_ocean_area,
        "op": op_by_month,
        "mlotst": mlt_by_month,
        "mask": mask_by_month,
        "op_clim": op_clim_by_month,
        "mlotst_clim": mlt_clim_by_month,
        "diff": diff_by_month,
        "metrics": metrics_by_month,
        "regions": regions_by_month,
        "diff_maps": diff_maps,
        "scatter": scatter,
    }

    print(f"  {model_name}: total ocean area = {total_ocean_area:.3e} m2")

# Ranking from the mean of the March and September RMS differences.
ranking = sorted(
    model_order,
    key=lambda m: 0.5
    * (
        results[m]["metrics"]["March"]["rms"]
        + results[m]["metrics"]["September"]["rms"]
    ),
)

print("\n=== Metrics on each model's native grid and shared mask ===")
for model_name in model_order:
    print(f"\n{model_name}:")
    for month_name in months:
        met = results[model_name]["metrics"][month_name]
        print(
            f"  {month_name:9s}  "
            f"bias = {met['bias']:8.1f} m   "
            f"RMS = {met['rms']:8.1f} m   "
            f"r = {met['corr']:.3f}   "
            f"|diff|>20 m = {met['exceed_frac']:5.1f}%   "
            f"shared area = {met['shared_area']:.3e} m2 "
            f"({met['shared_pct']:.1f}% of ocean)   "
            f"cells = {met['cell_count']}"
        )

print("\n=== Internal-consistency ranking ===")
print("Criterion: mean of March and September RMS(operator MLD - mlotst), native grid; smallest = most consistent.")
for rank, model_name in enumerate(ranking, start=1):
    rms_mean = 0.5 * (
        results[model_name]["metrics"]["March"]["rms"]
        + results[model_name]["metrics"]["September"]["rms"]
    )
    print(f"  {rank}. {model_name:14s} mean RMS = {rms_mean:.1f} m")

print("\n=== Regions where |operator MLD - mlotst| > 20 m is concentrated ===")
for model_name in model_order:
    print(f"\n{model_name}:")
    for month_name in months:
        reg = results[model_name]["regions"][month_name]
        lat_items = [(k, v) for k, v in reg.items() if k.startswith("lat ")]
        lat_items.sort(key=lambda kv: kv[1], reverse=True)
        basin_items = [(k, v) for k, v in reg.items() if not k.startswith("lat ")]
        basin_items.sort(key=lambda kv: kv[1], reverse=True)
        print(f"  {month_name}:")
        print("    top latitude bands by exceedance fraction:")
        for k, v in lat_items[:4]:
            print(f"      {k:18s} {v:5.1f}%")
        print("    named regions:")
        for k, v in basin_items:
            print(f"      {k:18s} {v:5.1f}%")

# Figure 1 display scale from pooled regridded differences.
all_diff_maps = np.concatenate(
    [
        results[m]["diff_maps"][mon].ravel()
        for m in model_order
        for mon in months
    ]
)
vmax = float(np.nanpercentile(np.abs(all_diff_maps), 98))
if not np.isfinite(vmax) or vmax <= 0.0:
    vmax = 1.0

# Figure 2 common logarithmic limits from all native shared-mask cells.
all_scatter_vals = []
for m in model_order:
    for mon in months:
        x, y = results[m]["scatter"][mon]
        all_scatter_vals.append(x)
        all_scatter_vals.append(y)
all_scatter_vals = np.concatenate(
    [v[np.isfinite(v) & (v > 0.0)] for v in all_scatter_vals if len(v) > 0]
)
scatter_vmin = float(np.nanmin(all_scatter_vals))
scatter_vmax = float(np.nanmax(all_scatter_vals))
scatter_vmin = max(scatter_vmin, 0.1)
exp_min = np.floor(np.log10(scatter_vmin))
exp_max = np.ceil(np.log10(scatter_vmax))
scatter_ticks = 10.0 ** np.arange(exp_min, exp_max + 1.0, 1.0)
scatter_lim_min = float(scatter_ticks[0])
scatter_lim_max = float(scatter_ticks[-1])
