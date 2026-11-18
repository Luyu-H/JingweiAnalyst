import os
import glob
import numpy as np
import pandas as pd
import xarray as xr

# ---------------------------------------------------------------------------
# CMIP6 historical thetao / so / mlotst, 1971-2000
# ---------------------------------------------------------------------------
ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
START_YEAR = 1971
END_YEAR = 2000
MONTHS = [("March", 3), ("September", 9)]

MODELS = {
    "CESM2": dict(
        y="nlat", x="nlon", lat="lat", lon="lon", lev="lev",
        depth_scale=0.01, ocean_mask_from="sftof",
    ),
    "GFDL-ESM4": dict(
        y="lat", x="lon", lat="lat", lon="lon", lev="lev",
        depth_scale=1.0, ocean_mask_from=None,
    ),
    "MPI-ESM1-2-LR": dict(
        y="j", x="i", lat="latitude", lon="longitude", lev="lev",
        depth_scale=1.0, ocean_mask_from=None,
    ),
}


def _find_files(model, var):
    return sorted(glob.glob(os.path.join(ROOT, model, var, "*.nc")))


def open_var(model, var):
    files = _find_files(model, var)
    if not files:
        raise FileNotFoundError(f"No files for {model}/{var} under {ROOT}")
    ds = xr.open_mfdataset(
        files, combine="by_coords", use_cftime=True, chunks={"time": 1}
    )
    da = ds[var]

    # Subset 1971-2000 using cftime/calendar-safe year access.
    yrs = da.time.dt.year.values
    msk = (yrs >= START_YEAR) & (yrs <= END_YEAR)
    da = da.isel(time=np.where(msk)[0])

    # Drop duplicate time labels, if any.
    tv = da.time.values
    if len(tv) > 1:
        _, idx = np.unique(tv, return_index=True)
        if len(idx) < len(tv):
            da = da.isel(time=np.sort(idx))
    return da


def open_fixed(model, var):
    files = _find_files(model, var)
    if not files:
        return None
    if len(files) == 1:
        ds = xr.open_dataset(files[0])
    else:
        ds = xr.open_mfdataset(files, combine="by_coords")
    if var in ds:
        return ds[var]
    for v in ds.data_vars:
        return ds[v]
    return None


# ---------------------------------------------------------------------------
# EOS-80 density at p=0 (potential density) and MLD from a density threshold
# ---------------------------------------------------------------------------
def eos80_rho(T, S):
    """UNESCO/EOS-80 in-situ density at zero pressure, T in degC, S in PSU."""
    T = np.asarray(T, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)

    # CMIP6 salinity is often stored with units "0.001".  If a file stores
    # practical salinity as ~0.035 instead of ~35, convert it.
    finite = np.isfinite(S)
    if np.any(finite):
        if np.nanmax(S[finite]) < 2.0:
            S = S * 1000.0

    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    T5 = T4 * T

    rho_w = (
        999.842594
        + 6.793952e-2 * T
        - 9.095290e-3 * T2
        + 1.001685e-4 * T3
        - 1.120083e-6 * T4
        + 6.536332e-9 * T5
    )
    A = 8.24493e-1 - 4.0899e-3 * T + 7.724e-5 * T2 - 8.114e-7 * T3
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T2
    C = 4.8314e-4

    S = np.maximum(S, 0.0)
    return rho_w + A * S + B * S ** 1.5 + C * S * S


def mld_from_profile(T, S, z, delta=0.03):
    """
    MLD as the first depth where rho(z) - rho(surface) >= 0.03 kg m-3,
    linearly interpolated between the two bracketing level centres.
    If no crossing exists in the column, use the deepest valid level.
    """
    T = np.asarray(T, dtype=np.float64)
    S = np.asarray(S, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)

    rho = eos80_rho(T, S)
    valid = np.isfinite(rho)
    rho0 = rho[0]
    target = rho0 + delta

    exceed = valid & (rho >= target)
    has = np.any(exceed, axis=0)
    count = np.sum(valid, axis=0)

    first = np.argmax(exceed, axis=0)
    first = np.where(has, first, 1)  # avoid index -1 when there is no crossing
    k = first

    rho_k = np.take_along_axis(rho, k[None, ...], axis=0)[0]
    rho_km1 = np.take_along_axis(rho, (k - 1)[None, ...], axis=0)[0]
    z_k = z[k]
    z_km1 = z[k - 1]

    denom = rho_k - rho_km1
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = (target - rho_km1) / denom
    mld_cross = z_km1 + frac * (z_k - z_km1)

    last = np.maximum(count - 1, 0)
    z_last = z[last]
    mld = np.where(has, mld_cross, z_last)
    mld = np.where(valid[0] & (count > 0), mld, np.nan)
    return mld


# ---------------------------------------------------------------------------
# Regridding helper for display only (native -> regular lon/lat)
# ---------------------------------------------------------------------------
LON_EDGES = np.arange(-180.0, 181.0, 2.0)
LAT_EDGES = np.arange(-90.0, 91.0, 2.0)
LON_REG = 0.5 * (LON_EDGES[:-1] + LON_EDGES[1:])
LAT_REG = 0.5 * (LAT_EDGES[:-1] + LAT_EDGES[1:])


def regrid_to_regular(lon2d, lat2d, field, area, lon_edges, lat_edges):
    lon = ((np.asarray(lon2d) + 180.0) % 360.0) - 180.0
    lat = np.asarray(lat2d)
    field = np.asarray(field)

    good = np.isfinite(field) & np.isfinite(lon) & np.isfinite(lat)
    lonf = lon[good]
    latf = lat[good]
    valf = field[good]

    if area is not None:
        wf = np.asarray(area)[good].astype(np.float64)
    else:
        wf = np.ones_like(valf, dtype=np.float64)
    wf = np.where(np.isfinite(wf), wf, 0.0)

    nlon = len(lon_edges) - 1
    nlat = len(lat_edges) - 1
    ix = np.digitize(lonf, lon_edges) - 1
    iy = np.digitize(latf, lat_edges) - 1

    keep = (ix >= 0) & (ix < nlon) & (iy >= 0) & (iy < nlat) & (wf > 0)
    ix = ix[keep]
    iy = iy[keep]
    valf = valf[keep]
    wf = wf[keep]

    idx = iy * nlon + ix
    sumw = np.bincount(idx, weights=wf, minlength=nlat * nlon).reshape(nlat, nlon)
    sumv = np.bincount(idx, weights=wf * valf, minlength=nlat * nlon).reshape(nlat, nlon)

    with np.errstate(divide="ignore", invalid="ignore"):
        out = sumv / sumw
    out[sumw == 0] = np.nan
    return out


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------
results = {}

for model, meta in MODELS.items():
    thetao = open_var(model, "thetao")
    so = open_var(model, "so")
    mlotst = open_var(model, "mlotst")

    thetao, so, mlotst = xr.align(thetao, so, mlotst, join="inner")

    # Native horizontal coordinates.
    lon_da = thetao[meta["lon"]]
    lat_da = thetao[meta["lat"]]
    if lon_da.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon_da.values, lat_da.values)
    else:
        lon2d = lon_da.transpose(meta["y"], meta["x"]).values
        lat2d = lat_da.transpose(meta["y"], meta["x"]).values

    # Vertical coordinate, positive downward, metres.
    lev = thetao[meta["lev"]].values.astype(np.float64) * meta["depth_scale"]
    lev = np.abs(lev)
    flip_lev = lev[0] > lev[-1]
    z = lev[::-1] if flip_lev else lev

    # Cell areas and ocean mask.
    area_da = open_fixed(model, "areacello")
    if area_da is None:
        area = np.cos(np.deg2rad(lat2d))
    else:
        area = area_da.transpose(meta["y"], meta["x"]).values.astype(np.float64)

    sftof_da = open_fixed(model, "sftof")
    if sftof_da is not None:
        sftof = sftof_da.transpose(meta["y"], meta["x"]).values
        ocean_mask = (
            np.isfinite(area) & (area > 0)
            & np.isfinite(sftof) & (sftof > 0)
        )
    else:
        ocean_mask = np.isfinite(area) & (area > 0)

    total_ocean_area = float(np.nansum(np.where(ocean_mask, area, np.nan)))

    yrs = thetao.time.dt.year.values
    mons = thetao.time.dt.month.values

    model_res = dict(
        lon2d=lon2d,
        lat2d=lat2d,
        area=area,
        ocean_mask=ocean_mask,
        total_ocean_area=total_ocean_area,
        months={},
    )

    for mname, month in MONTHS:
        idx = np.where(
            (yrs >= START_YEAR) & (yrs <= END_YEAR) & (mons == month)
        )[0]

        if len(idx) != 30:
            print(f"Warning: {model} {mname} has {len(idx)} slices, expected 30")

        op_list = []
        ml_list = []

        for i in idx:
            T = thetao.isel(time=int(i)).transpose(
                meta["lev"], meta["y"], meta["x"]
            ).values
            S = so.isel(time=int(i)).transpose(
                meta["lev"], meta["y"], meta["x"]
            ).values

            if flip_lev:
                T = T[::-1, ...]
                S = S[::-1, ...]

            op_mld = mld_from_profile(T, S, z)
            ml = mlotst.isel(time=int(i)).transpose(
                meta["y"], meta["x"]
            ).values

            op_list.append(op_mld.astype(np.float32))
            ml_list.append(np.asarray(ml, dtype=np.float32))

        op_stack = np.stack(op_list, axis=0)
        ml_stack = np.stack(ml_list, axis=0)

        # Cells where both fields are defined in all 30 instances.
        valid = (
            np.all(np.isfinite(op_stack), axis=0)
            & np.all(np.isfinite(ml_stack), axis=0)
            & ocean_mask
        )

        op_clim = np.full(valid.shape, np.nan, dtype=np.float64)
        ml_clim = np.full(valid.shape, np.nan, dtype=np.float64)
        if np.any(valid):
            op_clim[valid] = np.mean(op_stack[:, valid], axis=0)
            ml_clim[valid] = np.mean(ml_stack[:, valid], axis=0)

        diff_clim = np.full(valid.shape, np.nan, dtype=np.float64)
        diff_clim[valid] = op_clim[valid] - ml_clim[valid]

        n_valid = int(np.sum(valid))
        bias_w = rmse_w = corr_w = np.nan
        bias_unw = rmse_unw = corr_unw = np.nan
        pct_gt20 = pct_area = np.nan
        W = np.nan

        if n_valid > 0:
            x = op_clim[valid].astype(np.float64)
            y = ml_clim[valid].astype(np.float64)
            diff = x - y
            w = area[valid].astype(np.float64)
            W = float(np.sum(w))

            if W > 0:
                bias_w = float(np.sum(w * diff) / W)
                rmse_w = float(np.sqrt(np.sum(w * diff * diff) / W))

                mx = np.sum(w * x) / W
                my = np.sum(w * y) / W
                cov = np.sum(w * (x - mx) * (y - my)) / W
                vx = np.sum(w * (x - mx) ** 2) / W
                vy = np.sum(w * (y - my) ** 2) / W
                corr_w = float(cov / np.sqrt(vx * vy)) if vx > 0 and vy > 0 else np.nan

                pct_gt20 = float(100.0 * np.sum(w[np.abs(diff) > 20.0]) / W)
                pct_area = float(100.0 * W / total_ocean_area)

            bias_unw = float(np.mean(diff))
            rmse_unw = float(np.sqrt(np.mean(diff * diff)))
            corr_unw = (
                float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 else np.nan
            )

        diff_reg = regrid_to_regular(
            lon2d, lat2d, diff_clim, area, LON_EDGES, LAT_EDGES
        )

        # Coarse latitude-band diagnostics of where disagreement lives.
        band_abs = {}
        for lo, hi in [(-90, -60), (-60, -30), (-30, 0),
                       (0, 30), (30, 60), (60, 90)]:
            m = valid & (lat2d >= lo) & (lat2d < hi)
            if np.any(m):
                ww = area[m].astype(np.float64)
                dd = np.abs(diff_clim[m])
                band_abs[f"{lo}..{hi}"] = float(np.sum(ww * dd) / np.sum(ww))
            else:
                band_abs[f"{lo}..{hi}"] = np.nan

        top_bands = sorted(
            [(k, v) for k, v in band_abs.items() if np.isfinite(v)],
            key=lambda kv: kv[1],
            reverse=True,
        )[:3]

        model_res["months"][mname] = dict(
            op_clim=op_clim,
            ml_clim=ml_clim,
            valid=valid,
            diff_clim=diff_clim,
            diff_reg=diff_reg,
            bias_w=bias_w,
            rmse_w=rmse_w,
            corr_w=corr_w,
            bias_unw=bias_unw,
            rmse_unw=rmse_unw,
            corr_unw=corr_unw,
            pct_gt20=pct_gt20,
            pct_area=pct_area,
            n_valid=n_valid,
            W=W,
            band_abs=band_abs,
            top_bands=top_bands,
        )

    results[model] = model_res


# ---------------------------------------------------------------------------
# Summary table, ranking and caveat
# ---------------------------------------------------------------------------
summary_rows = []
for model in MODELS:
    for mname, _ in MONTHS:
        r = results[model]["months"][mname]
        summary_rows.append(
            dict(
                model=model,
                month=mname,
                bias_m=r["bias_w"],
                rmse_m=r["rmse_w"],
                corr=r["corr_w"],
                pct_absdiff_gt_20m=r["pct_gt20"],
                pct_area_both_defined=r["pct_area"],
                n_cells=r["n_valid"],
            )
        )
summary_table = pd.DataFrame(summary_rows)

rank_scores = {}
for model in MODELS:
    vals = [results[model]["months"][mn]["rmse_w"] for mn, _ in MONTHS]
    rank_scores[model] = float(np.nanmean(vals))

model_ranking = sorted(
    MODELS,
    key=lambda m: (not np.isfinite(rank_scores[m]), rank_scores[m]),
)
ranking_criterion = (
    "Mean area-weighted RMSE of operator MLD minus mlotst over March and "
    "September (lower RMSE = more internally consistent)."
)

disagreement_regions = {
    model: {
        mname: results[model]["months"][mname]["top_bands"]
        for mname, _ in MONTHS
    }
    for model in MODELS
}

agreement_caveat = (
    "mlotst is a model-diagnosed field produced by each model's own mixed-layer "
    "algorithm, reference-depth convention, vertical discretization and temporal "
    "averaging. The offline operator here uses EOS-80 potential density at p=0, "
    "a surface reference density, first crossing of the 0.03 kg m-3 threshold, and "
    "linear interpolation between level centres. It also uses monthly mean thetao/so "
    "and therefore cannot reproduce sub-monthly mixing events. Nonzero differences "
    "are therefore expected even for the same simulation; they diagnose algorithm "
    "and discretization differences, not model physics."
)
