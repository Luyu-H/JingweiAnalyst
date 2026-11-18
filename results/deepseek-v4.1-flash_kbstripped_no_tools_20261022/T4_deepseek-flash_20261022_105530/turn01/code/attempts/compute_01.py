import os
import glob
import re
import numpy as np
import xarray as xr
import gsw

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
root = "/data1/lorraine-20419/benchmark_cmip6/historical"

models = {
    "CESM2": {
        "grid_label": "gn",
        "y_dim": "nlat",
        "x_dim": "nlon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 0.01,
    },
    "GFDL-ESM4": {
        "grid_label": "gr",
        "y_dim": "lat",
        "x_dim": "lon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "grid_label": "gn",
        "y_dim": "j",
        "x_dim": "i",
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "depth_coord": "lev",
        "depth_scale": 1.0,
    },
}

month_names = {3: "March", 9: "September"}


# -----------------------------------------------------------------------------
# File discovery and opening
# -----------------------------------------------------------------------------
def get_variant(model):
    grid = models[model]["grid_label"]
    pattern = os.path.join(
        root, model, "thetao",
        f"thetao_Omon_{model}_historical_*_{grid}*.nc"
    )
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No thetao files found for {model}: {pattern}")
    variants = []
    for f in files:
        parts = os.path.basename(f).split("_")
        if len(parts) >= 5:
            variants.append(parts[4])
    if not variants:
        raise RuntimeError(f"Could not parse variant from {files}")
    return sorted(set(variants))[0]


def open_var(model, var, variant, table="Omon"):
    grid = models[model]["grid_label"]
    pattern = os.path.join(
        root, model, var,
        f"{var}_{table}_{model}_historical_{variant}_{grid}*.nc"
    )
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found for {model} {var}: {pattern}")
    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True)
    else:
        ds = xr.open_mfdataset(files, combine="by_coords", use_cftime=True)
    return ds


def subset_time(ds):
    idx = (ds.indexes["time"].year >= 1971) & (ds.indexes["time"].year <= 2000)
    return ds.isel(time=np.where(idx)[0])


# -----------------------------------------------------------------------------
# Mask helper
# -----------------------------------------------------------------------------
def build_comparison_mask(*fields, ocean_only=False):
    """Common finite mask for the supplied fields."""
    mask = xr.ones_like(fields[0], dtype=bool)
    for f in fields:
        mask = mask & np.isfinite(f)
    if ocean_only:
        # Not used here: WOA23 is not part of this comparison.
        pass
    return mask


# -----------------------------------------------------------------------------
# Operator MLD from thetao/so
# -----------------------------------------------------------------------------
def compute_operator_mld(thetao_t, so_t, lat2d, lon2d, z, threshold=0.03):
    """
    Compute mixed-layer depth from a single monthly thetao/so profile.

    Parameters
    ----------
    thetao_t, so_t : ndarray (nlev, ny, nx)
        Potential temperature and practical salinity.
    lat2d, lon2d : ndarray (ny, nx)
    z : ndarray (nlev,)
        Positive-down model level depths in metres.

    Returns
    -------
    mld : ndarray (ny, nx)
        Operator MLD, NaN where no crossing is found.
    no_cross : ndarray (ny, nx)
        True where the surface profile is valid but the threshold is never
        reached.
    """
    # Treat unphysically large values as fill.
    th = np.where(np.isfinite(thetao_t) & (np.abs(thetao_t) < 1e10), thetao_t, np.nan)
    sp = np.where(np.isfinite(so_t) & (np.abs(so_t) < 1e10), so_t, np.nan)

    nlev, ny, nx = th.shape

    # Pressure and TEOS-10 conversions.
    z3 = np.broadcast_to(-z[:, None, None], (nlev, ny, nx))
    lat3 = np.broadcast_to(lat2d[None, :, :], (nlev, ny, nx))
    lon3 = np.broadcast_to(lon2d[None, :, :], (nlev, ny, nx))

    p = gsw.p_from_z(z3, lat3)
    SA = gsw.SA_from_SP(sp, p, lon3, lat3)
    CT = gsw.CT_from_pt(SA, th)
    sig0 = gsw.sigma0(SA, CT)

    # Reference = shallowest model level.
    ref = sig0[0, :, :]
    diff = sig0 - ref[None, :, :]

    finite = np.isfinite(diff)
    cond = (diff >= threshold) & finite

    # Do not cross a missing level: require no NaNs before the candidate level.
    cum_nan = np.cumsum(~finite, axis=0)
    cond = cond & (cum_nan == 0)

    has = cond.any(axis=0)
    k = np.argmax(cond, axis=0)
    has = has & (k > 0)

    idx_prev = np.maximum(k - 1, 0)

    sig_prev = np.take_along_axis(sig0, idx_prev[None, :, :], axis=0)[0]
    sig_k = np.take_along_axis(sig0, k[None, :, :], axis=0)[0]

    z_prev = z[idx_prev]
    z_k = z[k]

    denom = sig_k - sig_prev
    frac = np.where(denom != 0.0, (threshold - (sig_prev - ref)) / denom, 0.0)
    mld = z_prev + frac * (z_k - z_prev)
    mld = np.where(has, mld, np.nan)

    surface_valid = np.isfinite(th[0]) & np.isfinite(sp[0])
    no_cross = surface_valid & ~has

    return mld, no_cross


# -----------------------------------------------------------------------------
# Weighted statistics
# -----------------------------------------------------------------------------
def weighted_pearson(x, y, w):
    sw = w.sum()
    if sw <= 0:
        return np.nan
    xm = (w * x).sum() / sw
    ym = (w * y).sum() / sw
    cov = (w * (x - xm) * (y - ym)).sum() / sw
    vx = (w * (x - xm) ** 2).sum() / sw
    vy = (w * (y - ym) ** 2).sum() / sw
    if vx <= 0 or vy <= 0:
        return np.nan
    return cov / np.sqrt(vx * vy)


def compute_metrics(op_clim, ml_clim, area, total_ocean_area):
    mask = build_comparison_mask(op_clim, ml_clim, ocean_only=False)
    op_m = op_clim.where(mask)
    ml_m = ml_clim.where(mask)
    D = (op_m - ml_m).where(mask)
    w = area.where(mask)

    sw = float(w.sum())
    mean = float((D * w).sum() / sw)
    rmse = float(np.sqrt((D ** 2 * w).sum() / sw))

    x = op_m.values.ravel()
    y = ml_m.values.ravel()
    ww = w.values.ravel()
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(ww) & (ww > 0)
    corr = weighted_pearson(x[valid], y[valid], ww[valid])

    exceed_area = float(w.where(np.abs(D) > 20.0).sum())
    exceed_pct = 100.0 * exceed_area / sw
    coverage_pct = 100.0 * sw / total_ocean_area
    n_cells = int(mask.sum())

    metrics = {
        "mean": mean,
        "rmse": rmse,
        "corr": corr,
        "exceed_pct": exceed_pct,
        "coverage_pct": coverage_pct,
        "n_cells": n_cells,
        "common_area": sw,
    }
    return metrics, mask, D, op_m, ml_m


def regional_stats(D, area, lat2d):
    bands = [
        (-90.0, -60.0),
        (-60.0, -40.0),
        (-40.0, -20.0),
        (-20.0, 0.0),
        (0.0, 20.0),
        (20.0, 40.0),
        (40.0, 60.0),
        (60.0, 90.0),
    ]
    out = []
    for lo, hi in bands:
        m = (lat2d >= lo) & (lat2d < hi)
        w = area.where(m)
        d = D.where(m)
        sw = float(w.sum())
        if sw > 0:
            mean = float((d * w).sum() / sw)
            rmse = float(np.sqrt((d ** 2 * w).sum() / sw))
        else:
            mean = np.nan
            rmse = np.nan
        out.append({"lat_lo": lo, "lat_hi": hi, "mean": mean, "rmse": rmse, "area": sw})
    return out


# -----------------------------------------------------------------------------
# Conservative area-weighted remapping for display only
# -----------------------------------------------------------------------------
def remap_conservative(field, area, lat2d, lon2d, dlon=2.0, dlat=2.0):
    lon = ((lon2d + 180.0) % 360.0) - 180.0
    lat = lat2d

    lon_edges = np.arange(-180.0, 180.0 + dlon, dlon)
    lat_edges = np.arange(-90.0, 90.0 + dlat, dlat)
    nlon = len(lon_edges) - 1
    nlat = len(lat_edges) - 1

    ix = np.clip(np.floor((lon - lon_edges[0]) / dlon).astype(int), 0, nlon - 1)
    iy = np.clip(np.floor((lat - lat_edges[0]) / dlat).astype(int), 0, nlat - 1)

    valid = np.isfinite(field) & np.isfinite(area) & (area > 0)
    flat_idx = iy[valid] * nlon + ix[valid]

    sum_w = np.bincount(
        flat_idx, weights=area[valid], minlength=nlat * nlon
    ).reshape(nlat, nlon)
    sum_fw = np.bincount(
        flat_idx, weights=area[valid] * field[valid], minlength=nlat * nlon
    ).reshape(nlat, nlon)

    out = np.full((nlat, nlon), np.nan)
    np.divide(sum_fw, sum_w, out=out, where=sum_w > 0)

    lon_centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lat_centers = 0.5 * (lat_edges[:-1] + lat_edges[1:])

    R = 6371000.0
    target_area = (
        R ** 2
        * np.deg2rad(dlon)
        * (np.sin(np.deg2rad(lat_edges[1:])) - np.sin(np.deg2rad(lat_edges[:-1])))[:, None]
        * np.ones((1, nlon))
    )

    return out, lon_centers, lat_centers, sum_w, target_area


# -----------------------------------------------------------------------------
# Main loop over models
# -----------------------------------------------------------------------------
results = {}

for model in models:
    print(f"Processing {model} ...")
    cfg = models[model]
    variant = get_variant(model)

    thetao_ds = subset_time(open_var(model, "thetao", variant, table="Omon"))
    so_ds = subset_time(open_var(model, "so", variant, table="Omon"))
    mlotst_ds = subset_time(open_var(model, "mlotst", variant, table="Omon"))

    thetao_ds, so_ds, mlotst_ds = xr.align(
        thetao_ds, so_ds, mlotst_ds, join="inner"
    )

    area_ds = open_var(model, "areacello", variant, table="Ofx")
    area = area_ds["areacello"].squeeze(drop=True)
    if model == "CESM2":
        sftof_ds = open_var(model, "sftof", variant, table="Ofx")
        sftof = sftof_ds["sftof"].squeeze(drop=True)
        sftof = sftof.where(np.isfinite(sftof) & (sftof < 1e10))
        area = area.where(sftof > 0)

    area = area.where(np.isfinite(area) & (area > 0) & (area < 1e18))
    area = area.transpose(cfg["y_dim"], cfg["x_dim"])

    lat_da = thetao_ds[cfg["lat_coord"]]
    lon_da = thetao_ds[cfg["lon_coord"]]
    if lat_da.ndim == 1 and lon_da.ndim == 1:
        lat2d, lon2d = xr.broadcast(lat_da, lon_da)
    else:
        lat2d, lon2d = lat_da, lon_da
    lat2d = lat2d.transpose(cfg["y_dim"], cfg["x_dim"])
    lon2d = lon2d.transpose(cfg["y_dim"], cfg["x_dim"])

    lev_name = cfg["depth_coord"]
    z = thetao_ds[lev_name].values.astype(float) * cfg["depth_scale"]

    template = thetao_ds["thetao"].isel(time=0, lev=0, drop=True)

    ny, nx = lat2d.shape

    month_sums = {}
    for m in (3, 9):
        month_sums[m] = {
            "op": np.zeros((ny, nx), dtype=float),
            "ml": np.zeros((ny, nx), dtype=float),
            "both": np.zeros((ny, nx), dtype=int),
            "nocross": np.zeros((ny, nx), dtype=int),
        }

    times = thetao_ds.time.values
    for i, t in enumerate(times):
        m = t.month
        if m not in (3, 9):
            continue

        thetao_t = thetao_ds["thetao"].isel(time=i).load().values
        so_t = so_ds["so"].isel(time=i).load().values
        ml_t = mlotst_ds["mlotst"].isel(time=i).load().values

        ml_t = np.where(np.isfinite(ml_t) & (np.abs(ml_t) < 1e10), ml_t, np.nan)

        op_mld, no_cross = compute_operator_mld(
            thetao_t, so_t, lat2d.values, lon2d.values, z
        )

        valid = np.isfinite(op_mld) & np.isfinite(ml_t)
        month_sums[m]["both"] += valid.astype(int)
        month_sums[m]["op"] += np.where(valid, op_mld, 0.0)
        month_sums[m]["ml"] += np.where(valid, ml_t, 0.0)
        month_sums[m]["nocross"] += no_cross.astype(int)

    area2 = area.broadcast_like(template)
    area2 = area2.where(np.isfinite(area2) & (area2 > 0))
    total_ocean_area = float(area2.sum())

    # Display grid for Figure 1 (computed on native grid and remapped only here).
    ones = np.ones_like(area2.values)
    _, lon_disp, lat_disp, sum_w_all, target_area = remap_conservative(
        ones, area2.values, lat2d.values, lon2d.values, dlon=2.0, dlat=2.0
    )
    ocean_frac_disp = sum_w_all / target_area

    results[model] = {
        "variant": variant,
        "area": area2,
        "total_ocean_area": total_ocean_area,
        "lat2d": lat2d,
        "lon2d": lon2d,
        "display": {
            "lon": lon_disp,
            "lat": lat_disp,
            "ocean_frac": ocean_frac_disp,
        },
    }

    for m in (3, 9):
        both = month_sums[m]["both"]
        mask30 = both == 30

        op_clim = template.copy(
            data=np.where(mask30, month_sums[m]["op"] / 30.0, np.nan)
        )
        ml_clim = template.copy(
            data=np.where(mask30, month_sums[m]["ml"] / 30.0, np.nan)
        )
        op_clim.name = "operator_mld"
        op_clim.attrs = {"units": "m", "long_name": "operator MLD"}
        ml_clim.name = "mlotst"
        ml_clim.attrs = {"units": "m", "long_name": "published mlotst"}

        metrics, mask, D, op_m, ml_m = compute_metrics(
            op_clim, ml_clim, area2, total_ocean_area
        )

        nocross_any = month_sums[m]["nocross"] > 0
        nocross_area = float(area2.where(nocross_any).sum())
        nocross_pct = 100.0 * nocross_area / total_ocean_area
        nocross_cells = int(nocross_any.sum())

        reg_stats = regional_stats(D, area2, lat2d.values)

        D_disp, _, _, _, _ = remap_conservative(
            D.values, area2.values, lat2d.values, lon2d.values, dlon=2.0, dlat=2.0
        )

        results[model][m] = {
            "op_clim": op_clim,
            "ml_clim": ml_clim,
            "mask": mask,
            "D": D,
            "op_m": op_m,
            "ml_m": ml_m,
            "metrics": metrics,
            "regional": reg_stats,
            "nocross_area": nocross_area,
            "nocross_pct": nocross_pct,
            "nocross_cells": nocross_cells,
            "D_display": D_disp,
        }

# -----------------------------------------------------------------------------
# Ranking: primary = mean of March/September area-weighted RMSE
# -----------------------------------------------------------------------------
ranking = []
for model in models:
    rmse_avg = np.nanmean(
        [results[model][3]["metrics"]["rmse"], results[model][9]["metrics"]["rmse"]]
    )
    abs_mean_avg = np.nanmean(
        [
            abs(results[model][3]["metrics"]["mean"]),
            abs(results[model][9]["metrics"]["mean"]),
        ]
    )
    exceed_avg = np.nanmean(
        [
            results[model][3]["metrics"]["exceed_pct"],
            results[model][9]["metrics"]["exceed_pct"],
        ]
    )
    ranking.append(
        {
            "model": model,
            "rmse_avg": rmse_avg,
            "abs_mean_avg": abs_mean_avg,
            "exceed_avg": exceed_avg,
        }
    )

ranking = sorted(ranking, key=lambda d: d["rmse_avg"])

# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------
print("\nOperator MLD minus published mlotst, 1971-2000 climatologies")
print("=" * 78)
for model in models:
    for m in (3, 9):
        r = results[model][m]["metrics"]
        print(
            f"{model:14s} {month_names[m]:10s} "
            f"mean={r['mean']:8.1f} m  "
            f"RMSE={r['rmse']:7.1f} m  "
            f"r={r['corr']:7.3f}  "
            f"|D|>20 m={r['exceed_pct']:5.1f}%  "
            f"coverage={r['coverage_pct']:5.1f}%  "
            f"cells={r['n_cells']:7d}"
        )
        print(
            f"{'':14s} {'':10s} "
            f"non-crossing exclusion area={results[model][m]['nocross_pct']:5.1f}% "
            f"of total ocean ({results[model][m]['nocross_cells']} native cells)"
        )

print("\nRanking (primary criterion: mean area-weighted RMSE over March and September; lower is more consistent)")
for rank, d in enumerate(ranking, start=1):
    print(
        f"{rank}. {d['model']:14s} "
        f"mean RMSE={d['rmse_avg']:6.1f} m  "
        f"mean |bias|={d['abs_mean_avg']:6.1f} m  "
        f"mean |D|>20 m area={d['exceed_avg']:5.1f}%"
    )

print("\nRegional area-weighted mean and RMSE of D (operator - mlotst), by latitude band")
for model in models:
    for m in (3, 9):
        print(f"\n{model} {month_names[m]}")
        for row in results[model][m]["regional"]:
            print(
                f"  {row['lat_lo']:6.0f} to {row['lat_hi']:6.0f}: "
                f"mean={row['mean']:8.1f} m  RMSE={row['rmse']:8.1f} m"
            )
