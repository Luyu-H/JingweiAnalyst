import glob
import os
import numpy as np
import pandas as pd
import xarray as xr
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Paths and helper functions
# ---------------------------------------------------------------------------
WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASK_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"
CMIP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"


def read_msk_numeric(path):
    """Read a whitespace/comma-delimited WOA .msk table, returning numeric rows."""
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except ValueError:
                    pass
            if len(vals) >= 2:
                rows.append(vals)
    return rows


def nitracline_depth(no3, depth):
    """
    Nitracline depth = depth of maximum vertical nitrate gradient within the
    top 1000 m.  The gradient at the deepest top-1000-m level still uses the
    next level below 1000 m because the full profile is passed in.
    """
    arr = np.asarray(no3, dtype=float)
    z = np.asarray(depth, dtype=float)
    out = np.full(arr.shape[:-1], np.nan, dtype=float)
    out_flat = out.reshape(-1)
    flat = arr.reshape(-1, arr.shape[-1])
    z_ok = np.isfinite(z)

    for i, prof in enumerate(flat):
        valid = np.isfinite(prof) & z_ok
        if valid.sum() < 2:
            continue

        idx = np.where(valid)[0]
        i0, i1 = idx[0], idx[-1]
        if i1 - i0 + 1 < 2:
            continue

        zz = z[i0:i1 + 1]
        yy = prof[i0:i1 + 1]

        # Remove any internal gaps if they occur (normally only below-bottom NaNs).
        if not np.all(np.isfinite(yy)):
            m = np.isfinite(yy)
            if m.sum() < 2:
                continue
            zz = zz[m]
            yy = yy[m]

        if zz[0] > zz[-1]:
            zz = zz[::-1]
            yy = yy[::-1]

        try:
            g = np.gradient(yy, zz)
        except Exception:
            continue

        cand = (zz <= 1000.0) & np.isfinite(g) & np.isfinite(yy)
        if not np.any(cand):
            continue

        j = np.argmax(np.where(cand, g, -np.inf))
        out_flat[i] = zz[j]

    return out


def build_folded_basin(raw_basin, lat2d, lon2d):
    """
    Rebuild the WOA basin mask on the 1-degree grid and fold any non-major
    basin codes into the nearest major basin.  Output codes are:
      1 Atlantic, 2 Pacific, 3 Indian, 4 Arctic, 5 Southern.
    """
    folded = np.full(raw_basin.shape, np.nan, dtype=float)
    vals = raw_basin[np.isfinite(raw_basin)]
    vals = vals[vals > 0]
    if vals.size == 0:
        return folded

    vals_int = np.rint(vals).astype(int)
    uniq, counts = np.unique(vals_int, return_counts=True)

    if all(c in uniq for c in [1, 2, 3, 4, 5]):
        major_codes = np.array([1, 2, 3, 4, 5])
    else:
        order = np.argsort(counts)[::-1]
        major_codes = uniq[order[:min(5, len(uniq))]]

    name_code = {"Atlantic": 1, "Pacific": 2, "Indian": 3,
                 "Arctic": 4, "Southern": 5}

    code_to_name = {}
    for code in major_codes:
        m = np.rint(raw_basin) == code
        if not np.any(m):
            continue
        lat_c = np.mean(lat2d[m])
        lon_rad = np.deg2rad(lon2d[m])
        lon_c = np.rad2deg(np.angle(np.mean(np.exp(1j * lon_rad))))

        if lat_c < -55:
            name = "Southern"
        elif lat_c > 55:
            name = "Arctic"
        elif -80 <= lon_c <= 20:
            name = "Atlantic"
        elif 20 < lon_c <= 120:
            name = "Indian"
        else:
            name = "Pacific"

        code_to_name[code] = name

    # Fallback to the standard WOA13 numbering if classification is ambiguous.
    if len(set(code_to_name.values())) < len(code_to_name) or len(code_to_name) < 5:
        if all(c in uniq for c in [1, 2, 3, 4, 5]):
            code_to_name = {1: "Atlantic", 2: "Pacific", 3: "Indian",
                            4: "Arctic", 5: "Southern"}
            major_codes = np.array([1, 2, 3, 4, 5])

    for code, name in code_to_name.items():
        folded[np.rint(raw_basin) == code] = name_code[name]

    major_mask = (
        np.isfinite(raw_basin)
        & (np.rint(raw_basin) > 0)
        & np.isin(np.rint(raw_basin), major_codes)
    )
    need = (
        np.isfinite(raw_basin)
        & (np.rint(raw_basin) > 0)
        & ~np.isin(np.rint(raw_basin), major_codes)
    )

    if np.any(need):
        pts = np.column_stack([lat2d[major_mask], lon2d[major_mask]])
        vals_m = np.rint(raw_basin[major_mask]).astype(int)

        pts_wrap = np.vstack([
            pts,
            np.column_stack([pts[:, 0], pts[:, 1] + 360.0]),
            np.column_stack([pts[:, 0], pts[:, 1] - 360.0]),
        ])
        vals_wrap = np.tile(vals_m, 3)

        tree = cKDTree(pts_wrap)
        q = np.column_stack([lat2d[need], lon2d[need]])
        dist, idx = tree.query(q, k=1)
        nearest = vals_wrap[idx]
        folded[need] = [name_code[code_to_name[c]] for c in nearest]

    return folded


# ---------------------------------------------------------------------------
# WOA23 grid and masks
# ---------------------------------------------------------------------------
woa_lat = np.arange(-89.5, 89.6, 1.0)
woa_lon = np.arange(-179.5, 179.6, 1.0)
lat2d_w, lon2d_w = np.meshgrid(woa_lat, woa_lon, indexing="ij")
area_woa = np.cos(np.deg2rad(lat2d_w))

# Land-sea mask: Bottom_Standard_Level == 1 is land; > 1 is ocean.
ls_rows = read_msk_numeric(os.path.join(WOA_MASK_ROOT, "landsea_01.msk"))
bottom = np.full((len(woa_lat), len(woa_lon)), np.nan)
for vals in ls_rows:
    if len(vals) < 3:
        continue
    lat, lon, b = vals[0], vals[1], vals[2]
    iy = int(round(lat - woa_lat[0]))
    ix = int(round(lon - woa_lon[0]))
    if 0 <= iy < len(woa_lat) and 0 <= ix < len(woa_lon):
        bottom[iy, ix] = b
woa_ocean_mask = bottom > 1

# Basin mask: sparse table; first basin value after lon/lat is the surface basin.
basin_rows = read_msk_numeric(os.path.join(WOA_MASK_ROOT, "basinmask_01.msk"))
raw_basin = np.full((len(woa_lat), len(woa_lon)), np.nan)
for vals in basin_rows:
    if len(vals) < 3:
        continue
    lat, lon, b = vals[0], vals[1], vals[2]
    iy = int(round(lat - woa_lat[0]))
    ix = int(round(lon - woa_lon[0]))
    if 0 <= iy < len(woa_lat) and 0 <= ix < len(woa_lon):
        raw_basin[iy, ix] = b

folded_basin = build_folded_basin(raw_basin, lat2d_w, lon2d_w)
basin_names = {"Atlantic": 1, "Pacific": 2, "Indian": 3}

# ---------------------------------------------------------------------------
# WOA23 nitrate annual climatology
# ---------------------------------------------------------------------------
woa_file = os.path.join(WOA_ROOT, "nitrate", "all_time_periods", "woa23_all_n00_01.nc")
ds_woa = xr.open_dataset(woa_file, decode_times=False)


def pick_woa_var(ds, field="an"):
    for v in ds.data_vars:
        if v.endswith("_" + field) or v == field:
            return v
    for v in ds.data_vars:
        sn = ds[v].attrs.get("standard_name", "")
        if "nitrate" in sn:
            return v
    return list(ds.data_vars)[0]


woa_var = pick_woa_var(ds_woa, "an")
woa_da = ds_woa[woa_var]
if "time" in woa_da.dims:
    woa_da = woa_da.isel(time=0)
woa_da = woa_da.squeeze()

dims = woa_da.dims
depth_name = [d for d in dims if "depth" in d.lower()][0]
lat_name = [d for d in dims if "lat" in d.lower()][0]
lon_name = [d for d in dims if "lon" in d.lower()][0]

woa_no3 = woa_da.transpose(lat_name, lon_name, depth_name).values.astype(float)
woa_depth = ds_woa[depth_name].values.astype(float)
woa_lat = ds_woa[lat_name].values.astype(float)
woa_lon = ds_woa[lon_name].values.astype(float)

if woa_depth[0] > woa_depth[-1]:
    woa_depth = woa_depth[::-1]
    woa_no3 = woa_no3[..., ::-1]

lat2d_w, lon2d_w = np.meshgrid(woa_lat, woa_lon, indexing="ij")
area_woa = np.cos(np.deg2rad(lat2d_w))
woa_surf_ocean = np.isfinite(woa_no3[:, :, 0])

woa_nit = nitracline_depth(woa_no3, woa_depth)

valid_woa = np.isfinite(woa_nit)
woa_global_nit_depth = (
    np.nansum(woa_nit[valid_woa] * area_woa[valid_woa])
    / np.nansum(area_woa[valid_woa])
)

woa_basin_nit_depth = {}
for name, code in basin_names.items():
    m = valid_woa & (folded_basin == code)
    woa_basin_nit_depth[name] = (
        np.nansum(woa_nit[m] * area_woa[m]) / np.nansum(area_woa[m])
    )

# ---------------------------------------------------------------------------
# CMIP6 model configuration and loading
# ---------------------------------------------------------------------------
model_cfg = {
    "CESM2": {
        "y": "nlat", "x": "nlon", "lat": "lat", "lon": "lon",
        "depth": "lev", "depth_scale": 0.01,
    },
    "GFDL-ESM4": {
        "y": "lat", "x": "lon", "lat": "lat", "lon": "lon",
        "depth": "lev", "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "y": "j", "x": "i", "lat": "latitude", "lon": "longitude",
        "depth": "lev", "depth_scale": 1.0,
    },
}
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]


def load_model_data(model):
    cfg = model_cfg[model]

    files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, "no3", "*.nc")))
    if not files:
        raise FileNotFoundError(f"No CMIP6 no3 files for {model}")

    ds = xr.open_mfdataset(
        files, combine="by_coords", use_cftime=True, decode_times=True
    )
    da = ds["no3"]

    time = ds["time"]
    sel = (time.dt.year >= 1971) & (time.dt.year <= 2000)
    da = da.isel(time=sel)
    clim = da.mean(dim="time", skipna=True)

    clim = clim.transpose(cfg["y"], cfg["x"], cfg["depth"])
    no3 = clim.values.astype(float)

    depth = ds[cfg["depth"]].values.astype(float) * cfg["depth_scale"]

    lat = ds[cfg["lat"]].values.astype(float)
    lon = ds[cfg["lon"]].values.astype(float)
    if lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon, lat)
    else:
        lat2d = lat
        lon2d = lon
    lon2d = (lon2d + 180.0) % 360.0 - 180.0

    if depth[0] > depth[-1]:
        depth = depth[::-1]
        no3 = no3[..., ::-1]

    # CESM2 areacello is global; apply its sftof mask if available.
    if model == "CESM2":
        sftof_files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, "sftof", "*.nc")))
        if sftof_files:
            ds_s = xr.open_dataset(sftof_files[0], use_cftime=True)
            sftof = ds_s["sftof"].squeeze().transpose(cfg["y"], cfg["x"]).values
            no3 = np.where((sftof > 0)[..., None], no3, np.nan)

    # Native areacello.
    area_files = sorted(glob.glob(os.path.join(CMIP_ROOT, model, "areacello", "*.nc")))
    if not area_files:
        raise FileNotFoundError(f"No areacello files for {model}")

    if len(area_files) == 1:
        ds_area = xr.open_dataset(area_files[0], use_cftime=True)
    else:
        ds_area = xr.open_mfdataset(
            area_files, combine="by_coords", use_cftime=True
        )

    area_da = ds_area["areacello"]
    if "time" in area_da.dims:
        area_da = area_da.isel(time=0)
    area_da = area_da.squeeze()
    area = area_da.transpose(cfg["y"], cfg["x"]).values.astype(float)

    return no3, depth, lat2d, lon2d, area


# ---------------------------------------------------------------------------
# Model nitracline depths, basin means, regridded comparison, and sections
# ---------------------------------------------------------------------------
target_woa_pts = np.column_stack([lat2d_w.ravel(), lon2d_w.ravel()])

# Basin lookup tree on WOA ocean points.
basin_valid = np.isfinite(folded_basin) & np.isin(folded_basin, [1, 2, 3, 4, 5])
basin_pts = np.column_stack([lat2d_w[basin_valid], lon2d_w[basin_valid]])
basin_vals = folded_basin[basin_valid].astype(int)
basin_pts_wrap = np.vstack([
    basin_pts,
    np.column_stack([basin_pts[:, 0], basin_pts[:, 1] + 360.0]),
    np.column_stack([basin_pts[:, 0], basin_pts[:, 1] - 360.0]),
])
basin_vals_wrap = np.tile(basin_vals, 3)
tree_basin = cKDTree(basin_pts_wrap)

# Equatorial Pacific section indices on the WOA 1-degree grid.
lon360_w = (woa_lon + 360.0) % 360.0
sec_lon_idx = np.where((lon360_w >= 120.0) & (lon360_w <= 280.0))[0]
sec_order = np.argsort(lon360_w[sec_lon_idx])
sec_lon = lon360_w[sec_lon_idx][sec_order]
sec_lon_idx = sec_lon_idx[sec_order]

lat_band_idx = np.where((woa_lat >= -2.0) & (woa_lat <= 2.0))[0]
depth_idx = np.where(woa_depth <= 1000.0)[0]
sec_depth = woa_depth[depth_idx]

woa_sub = woa_no3[lat_band_idx, :, :][:, sec_lon_idx, :][:, :, depth_idx]
woa_section_no3 = np.nanmean(woa_sub, axis=0)
woa_section_nit_line = np.nanmean(
    woa_nit[lat_band_idx, :][:, sec_lon_idx], axis=0
)

model_global_nit_depth = {}
model_basin_nit_depth = {}
model_nit_on_woa = {}
model_corr = {}
model_rmse = {}
model_section_no3 = {}
model_section_nit_line = {}

for model in models:
    no3, depth, lat2d, lon2d, area = load_model_data(model)
    nit = nitracline_depth(no3, depth)

    # Global area-weighted mean.
    valid = np.isfinite(nit) & np.isfinite(area) & (area > 0)
    model_global_nit_depth[model] = (
        np.sum(nit[valid] * area[valid]) / np.sum(area[valid])
    )

    # Assign each model cell to a WOA basin by nearest WOA ocean basin cell.
    q = np.column_stack([lat2d.ravel(), lon2d.ravel()])
    dist, bidx = tree_basin.query(q, k=1)
    model_basin = basin_vals_wrap[bidx].reshape(lat2d.shape)

    model_basin_nit_depth[model] = {}
    for name, code in basin_names.items():
        m = valid & (model_basin == code)
        if np.any(m):
            model_basin_nit_depth[model][name] = (
                np.sum(nit[m] * area[m]) / np.sum(area[m])
            )
        else:
            model_basin_nit_depth[model][name] = np.nan

    # Regrid nitracline depth to WOA grid for correlation/RMSE and maps.
    valid_h = np.isfinite(nit) & np.isfinite(lat2d) & np.isfinite(lon2d)
    pts = np.column_stack([lat2d[valid_h], lon2d[valid_h]])
    vals = nit[valid_h]

    pts_wrap = np.vstack([
        pts,
        np.column_stack([pts[:, 0], pts[:, 1] + 360.0]),
        np.column_stack([pts[:, 0], pts[:, 1] - 360.0]),
    ])
    vals_wrap = np.tile(vals, 3)
    tree = cKDTree(pts_wrap)
    dist, idx = tree.query(target_woa_pts, k=1)
    nit_on_woa = vals_wrap[idx].reshape(woa_nit.shape)
    model_nit_on_woa[model] = nit_on_woa

    common = np.isfinite(woa_nit) & np.isfinite(nit_on_woa)
    model_corr[model] = np.corrcoef(
        woa_nit[common], nit_on_woa[common]
    )[0, 1]
    model_rmse[model] = np.sqrt(
        np.mean((woa_nit[common] - nit_on_woa[common]) ** 2)
    )

    # Build model nitrate on the WOA grid for the equatorial Pacific section.
    surf_valid = (
        np.isfinite(no3[:, :, 0])
        & np.isfinite(lat2d)
        & np.isfinite(lon2d)
    ).ravel()
    pts3 = np.column_stack([lat2d.ravel()[surf_valid], lon2d.ravel()[surf_valid]])
    vals3 = no3.reshape(-1, no3.shape[-1])[surf_valid]

    pts3_wrap = np.vstack([
        pts3,
        np.column_stack([pts3[:, 0], pts3[:, 1] + 360.0]),
        np.column_stack([pts3[:, 0], pts3[:, 1] - 360.0]),
    ])
    vals3_wrap = np.tile(vals3, (3, 1))

    tree3 = cKDTree(pts3_wrap)
    dist3, idx3 = tree3.query(target_woa_pts, k=1)
    no3_on_woa = vals3_wrap[idx3].reshape(
        woa_nit.shape[0], woa_nit.shape[1], no3.shape[-1]
    )

    # Convert model mol m-3 to approx umol kg-1 using rho = 1025 kg m-3.
    no3_on_woa_umol = no3_on_woa * 1.0e6 / 1025.0
    no3_on_woa_umol[~woa_surf_ocean, :] = np.nan

    model_sub = no3_on_woa_umol[lat_band_idx, :, :][:, sec_lon_idx, :]
    prof = np.nanmean(model_sub, axis=0)

    model_sec = np.full((len(sec_lon_idx), len(sec_depth)), np.nan)
    z = depth
    for j in range(len(sec_lon_idx)):
        y = prof[j, :]
        validz = np.isfinite(y) & np.isfinite(z)
        if validz.sum() >= 2:
            zv = z[validz]
            yv = y[validz]
            if zv[0] > zv[-1]:
                zv = zv[::-1]
                yv = yv[::-1]
            model_sec[j, :] = np.interp(
                sec_depth, zv, yv, left=yv[0], right=yv[-1]
            )

    model_section_no3[model] = model_sec
    model_section_nit_line[model] = np.nanmean(
        nit_on_woa[lat_band_idx, :][:, sec_lon_idx], axis=0
    )

# ---------------------------------------------------------------------------
# Summary quantities
# ---------------------------------------------------------------------------
global_nit_depth_mean = {"WOA23": woa_global_nit_depth}
global_nit_depth_mean.update(model_global_nit_depth)
global_nit_depth_mean_round = {
    k: int(round(v)) for k, v in global_nit_depth_mean.items()
}

global_nit_depth_diff = {
    m: model_global_nit_depth[m] - woa_global_nit_depth for m in models
}
global_nit_depth_diff_round = {
    m: int(round(v)) for m, v in global_nit_depth_diff.items()
}

basin_nit_depth_mean = {"WOA23": woa_basin_nit_depth.copy()}
for m in models:
    basin_nit_depth_mean[m] = model_basin_nit_depth[m].copy()

basin_nit_depth_mean_round = {
    source: {b: int(round(v)) for b, v in vals.items()}
    for source, vals in basin_nit_depth_mean.items()
}

basin_nit_depth_diff = {}
for m in models:
    basin_nit_depth_diff[m] = {
        b: model_basin_nit_depth[m][b] - woa_basin_nit_depth[b]
        for b in basin_names
    }

basin_nit_depth_diff_round = {
    m: {b: int(round(v)) for b, v in vals.items()}
    for m, vals in basin_nit_depth_diff.items()
}

model_corr_3 = {m: round(model_corr[m], 3) for m in models}
model_rmse_1 = {m: round(model_rmse[m], 1) for m in models}

best_model = min(models, key=lambda m: (model_rmse[m], -model_corr[m]))

woa_basin_order = sorted(
    basin_names.keys(), key=lambda b: woa_basin_nit_depth[b], reverse=True
)
model_basin_order = {
    m: sorted(
        basin_names.keys(),
        key=lambda b: model_basin_nit_depth[m][b],
        reverse=True,
    )
    for m in models
}
basin_order_match = {
    m: model_basin_order[m] == woa_basin_order for m in models
}

section_data = {
    "lon": sec_lon,
    "depth": sec_depth,
    "WOA23": woa_section_no3,
    "CESM2": model_section_no3["CESM2"],
    "GFDL-ESM4": model_section_no3["GFDL-ESM4"],
    "MPI-ESM1-2-LR": model_section_no3["MPI-ESM1-2-LR"],
    "nit_line": {
        "WOA23": woa_section_nit_line,
        "CESM2": model_section_nit_line["CESM2"],
        "GFDL-ESM4": model_section_nit_line["GFDL-ESM4"],
        "MPI-ESM1-2-LR": model_section_nit_line["MPI-ESM1-2-LR"],
    },
}

summary_df = pd.DataFrame(
    {
        "global_mean_nitracline_m": [
            global_nit_depth_mean_round[k] for k in ["WOA23"] + models
        ],
        "Atlantic_m": [
            basin_nit_depth_mean_round[k]["Atlantic"] for k in ["WOA23"] + models
        ],
        "Pacific_m": [
            basin_nit_depth_mean_round[k]["Pacific"] for k in ["WOA23"] + models
        ],
        "Indian_m": [
            basin_nit_depth_mean_round[k]["Indian"] for k in ["WOA23"] + models
        ],
    },
    index=["WOA23"] + models,
)
