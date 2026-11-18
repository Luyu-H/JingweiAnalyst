# Computation code only
# Open the raw CMIP6 files, compute the offline MLD from thetao/so, and compare
# it with each model's published mlotst on the common mask.

import os
import glob
import re
import numpy as np
import pandas as pd
import xarray as xr
import gsw
from scipy.interpolate import griddata

ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
MONTHS = [3, 9]
MONTH_NAMES = {3: "March", 9: "September"}

MODEL_INFO = {
    "CESM2": {
        "grid_label": "gn",
        "y_dim": "nlat",
        "x_dim": "nlon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 0.01,   # lev is in cm for CESM2
        "n_levels": 60,
        "shape": (384, 320),
    },
    "GFDL-ESM4": {
        "grid_label": "gr",
        "y_dim": "lat",
        "x_dim": "lon",
        "lat_coord": "lat",
        "lon_coord": "lon",
        "depth_coord": "lev",
        "depth_scale": 1.0,
        "n_levels": 35,
        "shape": (180, 360),
    },
    "MPI-ESM1-2-LR": {
        "grid_label": "gn",
        "y_dim": "j",
        "x_dim": "i",
        "lat_coord": "latitude",
        "lon_coord": "longitude",
        "depth_coord": "lev",
        "depth_scale": 1.0,
        "n_levels": 40,
        "shape": (220, 256),
    },
}


def variant_from_filename(path):
    m = re.search(r"_(r\d+i\d+p\d+f\d+)_", os.path.basename(path))
    return m.group(1) if m else None


def find_files(model, var, table, variant=None):
    d = os.path.join(ROOT, model, var)
    if not os.path.isdir(d):
        return []
    fs = sorted(glob.glob(os.path.join(d, f"{var}_{table}_{model}_*.nc")))
    if not fs:
        fs = sorted(glob.glob(os.path.join(d, f"{var}_*.nc")))
    if variant is not None:
        fs = [f for f in fs if f"_{variant}_" in os.path.basename(f)]
    return fs


def open_var(model, var, table, variant=None, chunks=None):
    fs = find_files(model, var, table, variant)
    if not fs:
        return None, []
    if len(fs) == 1:
        ds = xr.open_dataset(
            fs[0], use_cftime=True, decode_times=True, chunks=chunks
        )
    else:
        ds = xr.open_mfdataset(
            fs,
            combine="by_coords",
            use_cftime=True,
            decode_times=True,
            chunks=chunks,
        )
    return ds, fs


def compute_mld(thetao, so, z, lon2d, lat2d, crit=3.0e-4):
    """
    Offline mixed-layer depth from potential temperature and practical salinity.

    The criterion is displaced-parcel buoyancy:
        dB = g * (rho_ambient(z) - rho_parcel_from_reference(z)) / rho_ref
        dB > 3e-4 m s-2  (~0.03 kg m-3)
    The crossing is linearly interpolated between the bracketing levels.
    """
    thetao = np.asarray(thetao, dtype=np.float64)
    so = np.asarray(so, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)

    nlev, ny, nx = thetao.shape

    # Local pressure at each level from depth and latitude.
    lat3 = np.broadcast_to(lat2d[None, :, :], thetao.shape)
    z3 = np.broadcast_to(-z[:, None, None], thetao.shape)  # gsw expects z positive up
    p = gsw.p_from_z(z3, lat3)

    lon3 = np.broadcast_to(lon2d[None, :, :], thetao.shape)

    # Practical salinity -> Absolute Salinity; potential temperature -> CT.
    SA = gsw.SA_from_SP(so, p, lon3, lat3)
    CT = gsw.CT_from_pt(SA, thetao)
    rho = gsw.rho(SA, CT, p)

    SA_ref = SA[0, :, :]
    CT_ref = CT[0, :, :]
    rho_ref = rho[0, :, :]

    # Reference parcel density evaluated at each target level's local pressure.
    rho_parcel = gsw.rho(SA_ref[None, :, :], CT_ref[None, :, :], p)
    dB = 9.81 * (rho - rho_parcel) / rho_ref[None, :, :]

    exceed = np.isfinite(dB) & (dB > crit)

    # Only columns with a valid shallowest level can define a reference parcel.
    surface_valid = np.isfinite(thetao[0, :, :]) & np.isfinite(so[0, :, :])
    exceed[:, ~surface_valid] = False

    has = np.any(exceed, axis=0)
    mld = np.full((ny, nx), np.nan, dtype=np.float64)

    if np.any(has):
        kidx = np.argmax(exceed, axis=0)
        idx = np.where(has)
        kk = kidx[idx]
        prev = kk - 1

        dB_prev = dB[prev, idx[0], idx[1]]
        dB_k = dB[kk, idx[0], idx[1]]
        z_prev = z[prev]
        z_k = z[kk]

        mld[idx] = z_prev + (crit - dB_prev) / (dB_k - dB_prev) * (z_k - z_prev)

    return mld


def weighted_stats(op, ml, area, mask, total_area):
    valid = (
        mask
        & np.isfinite(op)
        & np.isfinite(ml)
        & np.isfinite(area)
        & (area > 0.0)
    )
    n = int(np.sum(valid))
    if n == 0:
        return {
            "bias": np.nan,
            "rmsd": np.nan,
            "corr": np.nan,
            "pct_gt20": np.nan,
            "common_area_pct": np.nan,
            "n": 0,
        }

    w = area[valid]
    x = op[valid]
    y = ml[valid]
    d = x - y
    sw = float(np.sum(w))

    bias = float(np.sum(w * d) / sw)
    rmsd = float(np.sqrt(np.sum(w * d * d) / sw))

    mx = float(np.sum(w * x) / sw)
    my = float(np.sum(w * y) / sw)
    vx = float(np.sum(w * (x - mx) ** 2) / sw)
    vy = float(np.sum(w * (y - my) ** 2) / sw)
    cov = float(np.sum(w * (x - mx) * (y - my)) / sw)
    corr = cov / np.sqrt(vx * vy) if (vx > 0.0 and vy > 0.0) else np.nan

    pct_gt20 = float(np.sum(w[np.abs(d) > 20.0]) / sw * 100.0)
    common_area_pct = float(sw / total_area * 100.0) if total_area > 0.0 else np.nan

    return {
        "bias": bias,
        "rmsd": rmsd,
        "corr": corr,
        "pct_gt20": pct_gt20,
        "common_area_pct": common_area_pct,
        "n": n,
    }


def basin_masks(lat2d, lon2d):
    lon = np.mod(lon2d, 360.0)
    southern = lat2d < -50.0
    arctic = lat2d > 60.0
    nonpolar = ~(southern | arctic)

    pacific = nonpolar & (lon >= 120.0) & (lon < 280.0)
    atlantic = nonpolar & ((lon >= 280.0) | (lon < 20.0))
    indian = nonpolar & (lon >= 20.0) & (lon < 120.0)

    return {
        "Southern Ocean": southern,
        "Arctic": arctic,
        "Pacific": pacific,
        "Atlantic": atlantic,
        "Indian": indian,
    }


def interp_native_to_regular(lon2d, lat2d, field, valid, ocean_mask_native,
                             reg_lon, reg_lat):
    """
    Interpolate a native curvilinear/regular field to a common regular
    longitude-latitude grid for display only.
    """
    LON, LAT = np.meshgrid(reg_lon, reg_lat)
    lon = np.mod(lon2d, 360.0)

    good = (
        valid
        & np.isfinite(field)
        & np.isfinite(lon)
        & np.isfinite(lat2d)
    )

    out = np.full(LON.shape, np.nan, dtype=np.float64)

    if np.any(good):
        plon = lon[good]
        plat = lat2d[good]
        vals = field[good]

        # Duplicate across the 0/360 seam to avoid an artificial gap.
        plon3 = np.concatenate([plon, plon - 360.0, plon + 360.0])
        plat3 = np.concatenate([plat, plat, plat])
        vals3 = np.concatenate([vals, vals, vals])

        out = griddata((plon3, plat3), vals3, (LON, LAT), method="linear")

    # Land-sea mask for display, interpolated from all native grid points.
    allgood = np.isfinite(lon) & np.isfinite(lat2d)
    om = np.asarray(ocean_mask_native, dtype=float)

    if np.any(allgood):
        plon = lon[allgood]
        plat = lat2d[allgood]
        vals = om[allgood]

        plon3 = np.concatenate([plon, plon - 360.0, plon + 360.0])
        plat3 = np.concatenate([plat, plat, plat])
        vals3 = np.concatenate([vals, vals, vals])

        mask_reg = griddata((plon3, plat3), vals3, (LON, LAT), method="nearest")
        out = np.where(mask_reg >= 0.5, out, np.nan)
    else:
        mask_reg = np.zeros(LON.shape, dtype=float)

    return LON, LAT, out, mask_reg


# ---------------------------------------------------------------------------
# Discover variant labels from the actual filenames (never assume r1i1p1f1).
# ---------------------------------------------------------------------------
variants = {}
for model in MODELS:
    fs = find_files(model, "thetao", "Omon")
    variants[model] = variant_from_filename(fs[0]) if fs else None

# ---------------------------------------------------------------------------
# Open all required variables.
# ---------------------------------------------------------------------------
missing = []
datasets = {}

for model in MODELS:
    variant = variants[model]

    ds_th, f_th = open_var(
        model, "thetao", "Omon", variant, chunks={"time": 12}
    )
    ds_so, f_so = open_var(
        model, "so", "Omon", variant, chunks={"time": 12}
    )
    ds_ml, f_ml = open_var(
        model, "mlotst", "Omon", variant, chunks={"time": 12}
    )
    ds_ar, f_ar = open_var(
        model, "areacello", "Ofx", variant, chunks=None
    )

    ds_sf = None
    if model == "CESM2":
        ds_sf, f_sf = open_var(
            model, "sftof", "Ofx", variant, chunks=None
        )
        if not f_sf:
            missing.append((model, "sftof"))

    for name, ds, fs in [
        ("thetao", ds_th, f_th),
        ("so", ds_so, f_so),
        ("mlotst", ds_ml, f_ml),
        ("areacello", ds_ar, f_ar),
    ]:
        if not fs:
            missing.append((model, name))

    datasets[model] = {
        "thetao": ds_th,
        "so": ds_so,
        "mlotst": ds_ml,
        "areacello": ds_ar,
        "sftof": ds_sf,
        "variant": variant,
    }

processed_models = []
for model in MODELS:
    d = datasets[model]
    if (
        d["thetao"] is not None
        and d["so"] is not None
        and d["mlotst"] is not None
        and d["areacello"] is not None
        and (model != "CESM2" or d["sftof"] is not None)
    ):
        processed_models.append(model)

missing_report = pd.DataFrame(missing, columns=["model", "variable"])

# Common plotting grid.
reg_lon = np.arange(0.5, 360.0, 1.0)
reg_lat = np.arange(-89.5, 90.0, 1.0)

# ---------------------------------------------------------------------------
# Process each model and each calendar month.
# ---------------------------------------------------------------------------
results = []
disagreement_rows = []
plot_data = {}
scatter_data = {}

for model in processed_models:
    info = MODEL_INFO[model]
    ds_th = datasets[model]["thetao"]
    ds_so = datasets[model]["so"]
    ds_ml = datasets[model]["mlotst"]
    ds_ar = datasets[model]["areacello"]
    ds_sf = datasets[model]["sftof"]

    ny, nx = info["shape"]

    # Depth coordinate in metres.  Each coordinate is converted according to
    # its own declared/model-specific units; lev_bnds is not used here.
    z = np.asarray(ds_th[info["depth_coord"]].values, dtype=np.float64) * info["depth_scale"]

    # Native latitude/longitude coordinates.
    lat_da = ds_th[info["lat_coord"]]
    lon_da = ds_th[info["lon_coord"]]

    if lat_da.ndim == 1 and lon_da.ndim == 1:
        lon2d, lat2d = np.meshgrid(
            np.asarray(lon_da.values), np.asarray(lat_da.values)
        )
    else:
        lat2d = np.asarray(lat_da.values)
        lon2d = np.asarray(lon_da.values)
        if lat2d.ndim == 1:
            lat2d = np.broadcast_to(lat2d[:, None], (ny, nx))
        if lon2d.ndim == 1:
            lon2d = np.broadcast_to(lon2d[None, :], (ny, nx))

    lat2d = np.asarray(lat2d, dtype=np.float64)
    lon2d = np.asarray(lon2d, dtype=np.float64)

    # Cell areas.  CESM2 publishes a global areacello, so apply sftof.
    area_da = ds_ar["areacello"].transpose(info["y_dim"], info["x_dim"])

    if model == "CESM2":
        sftof_da = ds_sf["sftof"].transpose(info["y_dim"], info["x_dim"])
        area_ocean_da = area_da.where((sftof_da > 0.0) & (area_da > 0.0))
    else:
        area_ocean_da = area_da.where(area_da > 0.0)

    area_ocean_da = area_ocean_da.where(area_ocean_da < 1.0e20)
    area_ocean = np.asarray(area_ocean_da.values, dtype=np.float64)

    if area_ocean.shape != (ny, nx):
        area_ocean = np.asarray(
            area_ocean_da.transpose(info["y_dim"], info["x_dim"]).values,
            dtype=np.float64,
        )

    ocean_mask_native = np.isfinite(area_ocean) & (area_ocean > 0.0)
    total_area = float(np.nansum(area_ocean))

    years = ds_th.time.dt.year.values
    months = ds_th.time.dt.month.values

    for month in MONTHS:
        idxs = np.where(
            (years >= 1971) & (years <= 2000) & (months == month)
        )[0]

        if len(idxs) != 30:
            raise ValueError(
                f"{model} {MONTH_NAMES[month]}: expected 30 records, got {len(idxs)}"
            )

        years_so = ds_so.time.dt.year.values
        months_so = ds_so.time.dt.month.values
        idxs_so = np.where(
            (years_so >= 1971) & (years_so <= 2000) & (months_so == month)
        )[0]

        years_ml = ds_ml.time.dt.year.values
        months_ml = ds_ml.time.dt.month.values
        idxs_ml = np.where(
            (years_ml >= 1971) & (years_ml <= 2000) & (months_ml == month)
        )[0]

        if not np.array_equal(idxs, idxs_so) or not np.array_equal(idxs, idxs_ml):
            raise ValueError(
                f"{model} {MONTH_NAMES[month]}: time axes are not aligned"
            )

        sum_op = np.zeros((ny, nx), dtype=np.float64)
        sum_ml = np.zeros((ny, nx), dtype=np.float64)
        common_count = np.zeros((ny, nx), dtype=np.int16)
        excluded_fracs = []

        for it in idxs:
            theta = (
                ds_th["thetao"]
                .isel(time=it)
                .transpose(
                    info["depth_coord"], info["y_dim"], info["x_dim"]
                )
                .values.astype(np.float64)
            )
            so = (
                ds_so["so"]
                .isel(time=it)
                .transpose(
                    info["depth_coord"], info["y_dim"], info["x_dim"]
                )
                .values.astype(np.float64)
            )
            ml = (
                ds_ml["mlotst"]
                .isel(time=it)
                .transpose(info["y_dim"], info["x_dim"])
                .values.astype(np.float64)
            )

            op = compute_mld(theta, so, z, lon2d, lat2d)

            valid = (
                np.isfinite(op)
                & np.isfinite(ml)
                & np.isfinite(area_ocean)
                & (area_ocean > 0.0)
            )

            common_count += valid.astype(np.int16)
            sum_op += np.where(valid, op, 0.0)
            sum_ml += np.where(valid, ml, 0.0)

            # Cells with a valid profile but no threshold crossing are the
            # step-2 excluded cells.
            profile_valid = np.isfinite(theta[0, :, :]) & np.isfinite(so[0, :, :])
            no_cross = profile_valid & ~np.isfinite(op)
            excluded_fracs.append(
                float(np.nansum(area_ocean[no_cross]) / total_area * 100.0)
            )

        common = common_count == len(idxs)

        op_clim = np.full((ny, nx), np.nan, dtype=np.float64)
        ml_clim = np.full((ny, nx), np.nan, dtype=np.float64)

        op_clim[common] = sum_op[common] / float(len(idxs))
        ml_clim[common] = sum_ml[common] / float(len(idxs))

        stats = weighted_stats(
            op_clim, ml_clim, area_ocean, common, total_area
        )

        excl_mean = float(np.mean(excluded_fracs))
        excl_max = float(np.max(excluded_fracs))

        results.append(
            {
                "model": model,
                "month": month,
                "bias": stats["bias"],
                "rmsd": stats["rmsd"],
                "corr": stats["corr"],
                "pct_gt20": stats["pct_gt20"],
                "common_area_pct": stats["common_area_pct"],
                "excluded_area_pct": excl_mean,
                "excluded_area_pct_max": excl_max,
                "n_common_cells": stats["n"],
            }
        )

        plot_data[(model, month)] = {
            "lon2d": lon2d,
            "lat2d": lat2d,
            "diff": op_clim - ml_clim,
            "common": common,
            "area": area_ocean,
            "ocean_mask": ocean_mask_native,
        }

        scatter_data[(model, month)] = {
            "op": op_clim,
            "ml": ml_clim,
            "common": common,
        }

        # Tabulate where |difference| > 20 m concentrates.
        diff = op_clim - ml_clim
        basins = basin_masks(lat2d, lon2d)

        lat_edges = [-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0]

        for ib in range(len(lat_edges) - 1):
            lo = lat_edges[ib]
            hi = lat_edges[ib + 1]

            if hi < 90.0:
                band = (lat2d >= lo) & (lat2d < hi)
            else:
                band = (lat2d >= lo) & (lat2d <= hi)

            for bname, bmask in basins.items():
                sel = common & band & bmask & np.isfinite(diff)
                area_sel = float(np.nansum(area_ocean[sel]))

                if area_sel > 0.0:
                    bad = sel & (np.abs(diff) > 20.0)
                    area_bad = float(np.nansum(area_ocean[bad]))

                    disagreement_rows.append(
                        {
                            "model": model,
                            "month": month,
                            "lat_band": f"{lo:g} to {hi:g}",
                            "basin": bname,
                            "pct_gt20_in_region": area_bad / area_sel * 100.0,
                            "region_area_pct_of_total": area_sel / total_area * 100.0,
                        }
                    )

# ---------------------------------------------------------------------------
# Summaries.
# ---------------------------------------------------------------------------
results_df = pd.DataFrame(results)

if len(results_df) > 0:
    results_df = results_df[
        [
            "model",
            "month",
            "bias",
            "rmsd",
            "corr",
            "pct_gt20",
            "common_area_pct",
            "excluded_area_pct",
            "excluded_area_pct_max",
            "n_common_cells",
        ]
    ]

    rank_base = (
        results_df.groupby("model", as_index=False)
        .agg(
            rmsd_mean=("rmsd", "mean"),
            bias_mean=("bias", "mean"),
            corr_mean=("corr", "mean"),
            pct_gt20_mean=("pct_gt20", "mean"),
            excluded_area_pct_mean=("excluded_area_pct", "mean"),
        )
        .sort_values("rmsd_mean", ascending=True)
        .reset_index(drop=True)
    )
    rank_base["rank"] = np.arange(1, len(rank_base) + 1)
    rank_df = rank_base
else:
    rank_df = pd.DataFrame()

disagreement_df = pd.DataFrame(disagreement_rows)
excluded_df = (
    results_df[["model", "month", "excluded_area_pct", "excluded_area_pct_max"]].copy()
    if len(results_df) > 0
    else pd.DataFrame()
)

print("Missing variables (reported as missing):")
if len(missing_report) > 0:
    print(missing_report.to_string(index=False))
else:
    print("  none")

print("\nInternal-consistency scores: operator MLD minus published mlotst")
print("(common mask requires both fields to be defined in all 30 records of that month)")
if len(results_df) > 0:
    print(
        results_df.to_string(
            index=False,
            formatters={
                "bias": "{:.1f}".format,
                "rmsd": "{:.1f}".format,
                "corr": "{:.3f}".format,
                "pct_gt20": "{:.1f}".format,
                "common_area_pct": "{:.1f}".format,
                "excluded_area_pct": "{:.1f}".format,
                "excluded_area_pct_max": "{:.1f}".format,
            },
        )
    )

print("\nRanking: lowest mean area-weighted RMSD over March and September is most internally consistent")
if len(rank_df) > 0:
    print(
        rank_df.to_string(
            index=False,
            formatters={
                "rmsd_mean": "{:.1f}".format,
                "bias_mean": "{:.1f}".format,
                "corr_mean": "{:.3f}".format,
                "pct_gt20_mean": "{:.1f}".format,
                "excluded_area_pct_mean": "{:.1f}".format,
            },
        )
    )

if len(disagreement_df) > 0:
    print("\nLargest |MLD_op - mlotst| > 20 m fractions by latitude band and basin:")
    print(
        disagreement_df.sort_values(
            "pct_gt20_in_region", ascending=False
        )
        .head(20)
        .to_string(
            index=False,
            formatters={
                "pct_gt20_in_region": "{:.1f}".format,
                "region_area_pct_of_total": "{:.1f}".format,
            },
        )
    )

print(
    "\nResidual causes of operator-vs-mlotst disagreement (not model error):\n"
    "  1. reference-level choice: the operator uses the model's shallowest level, "
    "whereas mlotst is referenced to the ocean surface;\n"
    "  2. displaced-parcel buoyancy at local pressure vs a sigma-0 density criterion;\n"
    "  3. instantaneous diagnostic MLD averaged monthly vs an MLD formed from monthly-mean "
    "profiles, with the offline diagnostic expected to be shallower.\n"
    "  The common-mask statistics are on each model's native grid with its own areacello "
    "ocean total; they are not a cell-by-cell cross-model comparison."
)
