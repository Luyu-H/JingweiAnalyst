import os
import glob
import numpy as np
import pandas as pd
import xarray as xr

# ---------------------------------------------------------------------------
# Paths / constants
# ---------------------------------------------------------------------------
WOA_CLIM_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASK_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/masks"
CMIP_HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
RHO_REF = 1025.0

MODEL_META = {
    "CESM2": {
        "y_dim": "nlat", "x_dim": "nlon",
        "lat_coord": "lat", "lon_coord": "lon",
        "depth_coord": "lev", "depth_scale": 0.01,
    },
    "GFDL-ESM4": {
        "y_dim": "lat", "x_dim": "lon",
        "lat_coord": "lat", "lon_coord": "lon",
        "depth_coord": "lev", "depth_scale": 1.0,
    },
    "MPI-ESM1-2-LR": {
        "y_dim": "j", "x_dim": "i",
        "lat_coord": "latitude", "lon_coord": "longitude",
        "depth_coord": "lev", "depth_scale": 1.0,
    },
}
MODEL_NAMES = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
BASIN_IDS = {"Pacific": 2, "Atlantic": 1, "Indian": 3}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def find_var(ds, candidates):
    for c in candidates:
        if c in ds.data_vars:
            return c
    for v in ds.data_vars:
        for c in candidates:
            if c.lower() in v.lower():
                return v
    raise KeyError(f"None of {candidates} found in {list(ds.data_vars)}")


def nitracline_depth(conc, depth, max_depth=1000.0, absolute=False):
    """
    Depth of maximum vertical nitrate gradient, searching only intervals whose
    upper level is <= max_depth.  The gradient at the deepest searched level
    uses the level immediately below it.

    conc : array with vertical axis first, shape (nlev, ...)
    depth: 1-D vertical coordinate, positive down
    """
    conc = np.asarray(conc, dtype=float)
    depth = np.asarray(depth, dtype=float)

    # ensure increasing depth
    if depth.size > 1 and depth[0] > depth[-1]:
        depth = depth[::-1]
        conc = conc[::-1, ...]

    idx = np.where(depth[:-1] <= max_depth)[0]
    out_shape = conc.shape[1:]
    if idx.size == 0:
        return np.full(out_shape, np.nan)

    dz = depth[idx + 1] - depth[idx]
    dz_shape = (len(idx),) + (1,) * (conc.ndim - 1)
    grad = (conc[idx + 1] - conc[idx]) / dz.reshape(dz_shape)

    if absolute:
        grad = np.abs(grad)

    grad_filled = np.where(np.isfinite(grad), grad, -np.inf)
    imax = np.argmax(grad_filled, axis=0)
    maxg = np.take_along_axis(grad_filled, imax[None, ...], axis=0)[0]
    bad = ~np.isfinite(maxg)

    depth_upper = depth[idx]
    out = np.full(out_shape, np.nan)
    out[~bad] = depth_upper[imax[~bad]]
    return out


def latlon_to_xyz(lat, lon):
    lat = np.deg2rad(np.asarray(lat, dtype=float).ravel())
    lon = np.deg2rad(np.mod(np.asarray(lon, dtype=float).ravel(), 360.0))
    return np.column_stack((
        np.cos(lat) * np.cos(lon),
        np.cos(lat) * np.sin(lon),
        np.sin(lat),
    ))


def nearest_indices(src_lat, src_lon, tgt_lat, tgt_lon):
    src_xyz = latlon_to_xyz(src_lat, src_lon)
    tgt_xyz = latlon_to_xyz(tgt_lat, tgt_lon)

    if src_xyz.shape[0] == 0:
        raise ValueError("No source points for nearest-neighbor mapping")
    if tgt_xyz.shape[0] == 0:
        return np.empty(0, dtype=int)

    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(src_xyz)
        _, ind = tree.query(tgt_xyz)
        return ind
    except Exception:
        # Fallback brute-force nearest neighbour, chunked in target space.
        ind = np.empty(tgt_xyz.shape[0], dtype=int)
        chunk = 100
        for start in range(0, tgt_xyz.shape[0], chunk):
            stop = min(start + chunk, tgt_xyz.shape[0])
            d2 = ((tgt_xyz[start:stop, None, :] - src_xyz[None, :, :]) ** 2).sum(axis=-1)
            ind[start:stop] = np.argmin(d2, axis=1)
        return ind


def read_landsea(path, lat_axis, lon_axis):
    """Read WOA landsea_01.msk; ocean where Bottom_Standard_Level > 1."""
    nlat, nlon = len(lat_axis), len(lon_axis)
    ocean = np.zeros((nlat, nlon), dtype=bool)
    if nlat < 2 or nlon < 2:
        return ocean

    dlat = lat_axis[1] - lat_axis[0]
    dlon = lon_axis[1] - lon_axis[0]

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            if len(parts) < 3:
                continue
            try:
                lat = float(parts[0])
                lon = float(parts[1])
                bottom = float(parts[2])
            except ValueError:
                continue
            i = int(round((lat - lat_axis[0]) / dlat))
            j = int(round((lon - lon_axis[0]) / dlon))
            if 0 <= i < nlat and 0 <= j < nlon:
                ocean[i, j] = bottom > 1.0
    return ocean


def read_basin_mask(path, lat_axis, lon_axis, nlev=None):
    """
    Read basinmask_01.msk.  The file may contain either one basin number per
    ocean point, or a sequence of basin numbers on standard levels.  This
    routine extracts the shallowest valid basin code at each point and returns
    a 2-D array.
    """
    nlat, nlon = len(lat_axis), len(lon_axis)
    basin2d = np.full((nlat, nlon), np.nan)
    if nlat < 2 or nlon < 2:
        return basin2d, None

    dlat = lat_axis[1] - lat_axis[0]
    dlon = lon_axis[1] - lon_axis[0]

    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.replace(",", " ").split()
            if len(parts) < 3:
                continue
            try:
                lat = float(parts[0])
                lon = float(parts[1])
            except ValueError:
                continue

            vals = []
            for p in parts[2:]:
                try:
                    v = float(p)
                except ValueError:
                    continue
                if np.isfinite(v):
                    vals.append(v)
            if not vals:
                continue

            # Drop a possible leading level count / bottom index.
            start = 0
            if len(vals) > 1 and vals[0] >= 1 and abs(vals[0] - (len(vals) - 1)) < 1e-6:
                start = 1

            code = None
            for v in vals[start:]:
                if 1 <= v <= 11:
                    code = v
                    break
            if code is None:
                for v in vals[start:]:
                    if v > 0:
                        code = v
                        break
            if code is None:
                continue

            i = int(round((lat - lat_axis[0]) / dlat))
            j = int(round((lon - lon_axis[0]) / dlon))
            if 0 <= i < nlat and 0 <= j < nlon:
                basin2d[i, j] = code

    return basin2d, None


def fold_basin(b):
    """Fold WOA13 marginal-sea basin numbers into their parent ocean."""
    out = np.array(b, dtype=float, copy=True)
    # WOA13 marginal seas: Mediterranean, Black Sea, Baltic, Hudson Bay -> Atlantic
    out[np.isin(out, [6, 7, 8, 11])] = 1.0
    # Red Sea and Persian Gulf -> Indian
    out[np.isin(out, [9, 10])] = 3.0
    return out


def load_woa_nitrate():
    path = os.path.join(
        WOA_CLIM_ROOT, "nitrate", "all_time_periods", "woa23_all_n00_01.nc"
    )
    ds = xr.open_dataset(path, decode_times=False)
    var = find_var(ds, ["n_an", "n", "nitrate", "no3"])
    da = ds[var]

    # annual climatology: remove singleton time axis
    for d in list(da.dims):
        if d.lower() in ("time", "month") and da.sizes[d] == 1:
            da = da.isel({d: 0}, drop=True)

    depth_name = None
    for c in ["depth", "lev", "st_level"]:
        if c in ds.coords or c in ds.dims:
            depth_name = c
            break
    if depth_name is None:
        depth_name = da.dims[0]

    lat_name = "lat" if ("lat" in ds.coords or "lat" in ds.dims) else "latitude"
    lon_name = "lon" if ("lon" in ds.coords or "lon" in ds.dims) else "longitude"

    da = da.transpose(depth_name, lat_name, lon_name)
    depth = ds[depth_name].values.astype(float)
    lat = ds[lat_name].values.astype(float)
    lon = ds[lon_name].values.astype(float)
    N = da.values.astype(float)

    order = np.argsort(depth)
    depth = depth[order]
    N = N[order, ...]

    if lat.size > 1 and lat[0] > lat[-1]:
        lat = lat[::-1]
        N = N[:, ::-1, :]
    if lon.size > 1 and lon[0] > lon[-1]:
        lon = lon[::-1]
        N = N[:, :, ::-1]

    return N, depth, lat, lon


def load_cmip_no3(model, meta):
    files = sorted(glob.glob(os.path.join(CMIP_HIST_ROOT, model, "no3", "*.nc")))
    if not files:
        raise FileNotFoundError(f"No no3 files for {model}")

    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        use_cftime=True,
        chunks={"time": 60},
    )
    var = find_var(ds, ["no3"])
    da_all = ds[var]

    t = ds.indexes["time"]
    if hasattr(t[0], "year"):
        years = np.array([x.year for x in t])
    else:
        years = np.asarray(pd.to_datetime(t).year)
    sel = (years >= 1971) & (years <= 2000)
    if sel.sum() == 0:
        raise ValueError(f"No 1971-2000 data for {model}")

    da = da_all.isel(time=np.where(sel)[0]).mean("time")

    zname = meta["depth_coord"]
    z = ds[zname].values.astype(float) * float(meta["depth_scale"])

    lat = ds[meta["lat_coord"]]
    lon = ds[meta["lon_coord"]]
    if lat.ndim == 1 and lon.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon.values, lat.values)
    else:
        lat2d = lat.values
        lon2d = lon.values

    da = da.transpose(zname, meta["y_dim"], meta["x_dim"])
    N = da.compute().values.astype(float)
    units = da_all.attrs.get("units", "")

    order = np.argsort(z)
    z = z[order]
    N = N[order, ...]

    return N, z, lat2d, lon2d, units


def load_areacello(model, meta):
    files = sorted(glob.glob(os.path.join(CMIP_HIST_ROOT, model, "areacello", "*.nc")))
    if not files:
        raise FileNotFoundError(f"No areacello file for {model}")
    ds = xr.open_dataset(files[0])
    var = find_var(ds, ["areacello"])
    da = ds[var]
    if "time" in da.dims and da.sizes["time"] == 1:
        da = da.isel(time=0, drop=True)

    try:
        da = da.transpose(meta["y_dim"], meta["x_dim"])
    except Exception:
        dims = [d for d in da.dims if d not in ("time", "lev", "depth")]
        if len(dims) >= 2:
            da = da.transpose(dims[-2], dims[-1])
    return da.values.astype(float)


def make_lon_depth_section(N, depth, lat2d, lon2d, lon_edges, lat_min=-2.0, lat_max=2.0):
    """
    Average a 3-D field (depth, y, x) over a latitude band and bin by longitude.
    Longitudes are treated modulo 360 so 120E-80W crosses the dateline.
    """
    N = np.asarray(N, dtype=float)
    nlev = N.shape[0]
    centers = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    sec = np.full((nlev, centers.size), np.nan)

    lat_band = (lat2d >= lat_min) & (lat2d <= lat_max)
    lon360 = np.mod(lon2d, 360.0)

    Nflat = N.reshape(nlev, -1)
    latflat = lat_band.ravel()
    lonflat = lon360.ravel()

    for j in range(centers.size):
        m = latflat & (lonflat >= lon_edges[j]) & (lonflat < lon_edges[j + 1])
        if not np.any(m):
            continue
        vals = Nflat[:, m]
        cnt = np.sum(np.isfinite(vals), axis=1)
        sums = np.nansum(vals, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            sec[:, j] = np.where(cnt > 0, sums / cnt, np.nan)

    return centers, sec


# ---------------------------------------------------------------------------
# WOA23 observed nitrate, nitracline, masks
# ---------------------------------------------------------------------------
N_woa, depth_woa, lat_woa, lon_woa = load_woa_nitrate()
lat2d_woa, lon2d_woa = np.meshgrid(lat_woa, lon_woa, indexing="ij")

nt_woa = nitracline_depth(N_woa, depth_woa, max_depth=1000.0)

landsea_path = os.path.join(WOA_MASK_ROOT, "landsea_01.msk")
basin_path = os.path.join(WOA_MASK_ROOT, "basinmask_01.msk")

ocean_woa = read_landsea(landsea_path, lat_woa, lon_woa)
ocean_from_data = np.isfinite(N_woa).any(axis=0)

if ocean_woa.sum() == 0:
    ocean_woa = ocean_from_data
else:
    ocean_woa = ocean_woa & ocean_from_data
    if ocean_woa.sum() < 0.5 * ocean_from_data.sum():
        ocean_woa = ocean_from_data

basin2d_woa, basin3d_woa = read_basin_mask(
    basin_path, lat_woa, lon_woa, nlev=depth_woa.size
)

if basin3d_woa is not None:
    valid3 = np.isfinite(basin3d_woa) & (basin3d_woa > 0)
    has = valid3.any(axis=0)
    first_idx = np.argmax(valid3, axis=0)
    basin2d_woa = np.full(valid3.shape[1:], np.nan)
    basin2d_woa[has] = np.take_along_axis(
        basin3d_woa, first_idx[None, :, :], axis=0
    )[0][has]

basin2d_woa = fold_basin(basin2d_woa)
basin2d_woa[~ocean_woa] = np.nan

weights_woa = np.cos(np.deg2rad(lat2d_woa))
valid_woa = ocean_woa & np.isfinite(nt_woa)

# WOA global and basin means
global_raw = {}
basin_raw = {"WOA23": {}}

m = valid_woa
global_raw["WOA23"] = (
    np.nansum(nt_woa[m] * weights_woa[m]) / np.nansum(weights_woa[m])
)

for bname, bid in BASIN_IDS.items():
    m = valid_woa & (basin2d_woa == bid)
    denom = np.nansum(weights_woa[m])
    if denom > 0:
        basin_raw["WOA23"][bname] = (
            np.nansum(nt_woa[m] * weights_woa[m]) / denom
        )
    else:
        basin_raw["WOA23"][bname] = np.nan

# ---------------------------------------------------------------------------
# WOA section for Figure 2
# ---------------------------------------------------------------------------
lon_edges = np.arange(120.0, 281.0, 1.0)

lon_centers_sec, sec_woa = make_lon_depth_section(
    N_woa, depth_woa, lat2d_woa, lon2d_woa, lon_edges
)
nt_line_woa = nitracline_depth(sec_woa, depth_woa, max_depth=1000.0)
mask_depth_woa = depth_woa <= 1000.0

sections = {
    "WOA23": {
        "lon": lon_centers_sec,
        "depth_plot": depth_woa[mask_depth_woa],
        "conc_plot": sec_woa[mask_depth_woa, :],
        "nt_line": nt_line_woa,
    }
}

# Precompute WOA ocean-basin lookup for mapping model native grids to basins
woa_basin_ocean = (
    ocean_woa
    & np.isfinite(basin2d_woa)
    & (basin2d_woa > 0)
)
if not np.any(woa_basin_ocean):
    woa_basin_ocean = np.isfinite(basin2d_woa) & (basin2d_woa > 0)
if not np.any(woa_basin_ocean):
    raise RuntimeError("No WOA ocean basin points available for model-basin mapping")

src_lat_b = lat2d_woa[woa_basin_ocean]
src_lon_b = lon2d_woa[woa_basin_ocean]
src_basin_b = basin2d_woa[woa_basin_ocean]

# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
nt_native = {}
nt_regrid = {}
corr_raw = {}
rmse_raw = {}

for model in MODEL_NAMES:
    meta = MODEL_META[model]
    N_model, depth_model, lat2d_model, lon2d_model, units_model = load_cmip_no3(
        model, meta
    )
    area_model = load_areacello(model, meta)

    if area_model.shape != N_model.shape[1:]:
        if area_model.T.shape == N_model.shape[1:]:
            area_model = area_model.T
        else:
            area_model = np.squeeze(area_model)
            if area_model.shape != N_model.shape[1:]:
                raise ValueError(
                    f"areacello shape {area_model.shape} does not match "
                    f"{model} no3 spatial shape {N_model.shape[1:]}"
                )

    nt_model = nitracline_depth(N_model, depth_model, max_depth=1000.0)
    nt_native[model] = nt_model

    valid_model = (
        np.isfinite(nt_model)
        & np.isfinite(area_model)
        & (area_model > 0)
    )

    # Native-grid global mean weighted by the model's own areacello
    denom = np.nansum(area_model[valid_model])
    if denom > 0:
        global_raw[model] = (
            np.nansum(nt_model[valid_model] * area_model[valid_model])
            / denom
        )
    else:
        global_raw[model] = np.nan

    # Map model native points to WOA basins by nearest WOA ocean point
    valid_b = np.isfinite(nt_model)
    basin_model = np.full(nt_model.shape, np.nan)
    if np.any(valid_b):
        inds_b = nearest_indices(
            src_lat_b,
            src_lon_b,
            lat2d_model[valid_b],
            lon2d_model[valid_b],
        )
        basin_model[valid_b] = src_basin_b[inds_b]

    basin_raw[model] = {}
    for bname, bid in BASIN_IDS.items():
        m = valid_model & (basin_model == bid)
        denom = np.nansum(area_model[m])
        if denom > 0:
            basin_raw[model][bname] = (
                np.nansum(nt_model[m] * area_model[m]) / denom
            )
        else:
            basin_raw[model][bname] = np.nan

    # Regrid native nitracline to the WOA 1-degree grid by nearest valid model
    # point (used for maps, spatial correlation and RMSE).
    valid_r = valid_model
    if np.any(valid_r):
        inds_r = nearest_indices(
            lat2d_model[valid_r],
            lon2d_model[valid_r],
            lat2d_woa.ravel(),
            lon2d_woa.ravel(),
        )
        arr = nt_model[valid_r][inds_r].reshape(lat2d_woa.shape)
    else:
        arr = np.full(lat2d_woa.shape, np.nan)
    arr[~ocean_woa] = np.nan
    nt_regrid[model] = arr

    m = ocean_woa & np.isfinite(nt_woa) & np.isfinite(arr)
    x = arr[m]
    y = nt_woa[m]
    if x.size > 1:
        corr_raw[model] = float(np.corrcoef(x, y)[0, 1])
        rmse_raw[model] = float(np.sqrt(np.mean((x - y) ** 2)))
    else:
        corr_raw[model] = np.nan
        rmse_raw[model] = np.nan

    # Section: convert CMIP6 mol m-3 to approximate micromol kg-1 using a
    # constant reference density.  The nitracline depth itself is invariant
    # to this constant multiplicative factor.
    units_norm = str(units_model).lower().replace(" ", "")
    if "mol" in units_norm and ("m-3" in units_norm or "m^-3" in units_norm):
        conv = 1.0e6 / RHO_REF
    else:
        conv = 1.0
    N_model_plot = N_model * conv

    lon_centers_m, sec_m = make_lon_depth_section(
        N_model_plot, depth_model, lat2d_model, lon2d_model, lon_edges
    )
    nt_line_m = nitracline_depth(sec_m, depth_model, max_depth=1000.0)
    mask_depth_m = depth_model <= 1000.0

    sections[model] = {
        "lon": lon_centers_m,
        "depth_plot": depth_model[mask_depth_m],
        "conc_plot": sec_m[mask_depth_m, :],
        "nt_line": nt_line_m,
    }

# ---------------------------------------------------------------------------
# Rounded result variables for the answer
# ---------------------------------------------------------------------------
global_mean_nitracline_raw = global_raw.copy()
global_mean_nitracline = {
    k: int(round(v)) for k, v in global_mean_nitracline_raw.items()
}
global_diff_from_woa = {
    k: int(round(global_mean_nitracline_raw[k] - global_mean_nitracline_raw["WOA23"]))
    for k in MODEL_NAMES
}

basin_mean_nitracline_raw = {
    source: {bname: basin_raw[source][bname] for bname in BASIN_IDS}
    for source in ["WOA23", *MODEL_NAMES]
}
basin_mean_nitracline = {
    source: {
        bname: int(round(basin_mean_nitracline_raw[source][bname]))
        for bname in BASIN_IDS
    }
    for source in ["WOA23", *MODEL_NAMES]
}
basin_diff_from_woa = {
    model: {
        bname: int(
            round(
                basin_mean_nitracline_raw[model][bname]
                - basin_mean_nitracline_raw["WOA23"][bname]
            )
        )
        for bname in BASIN_IDS
    }
    for model in MODEL_NAMES
}

spatial_corr = {m: round(corr_raw[m], 3) for m in MODEL_NAMES}
spatial_rmse = {m: round(rmse_raw[m], 1) for m in MODEL_NAMES}

best_model = min(MODEL_NAMES, key=lambda m: rmse_raw[m])
best_model_by_corr = max(MODEL_NAMES, key=lambda m: corr_raw[m])


def ordered_basins(d):
    return [b for b, _ in sorted(d.items(), key=lambda kv: kv[1], reverse=True)]


basin_order_woa = ordered_basins(basin_mean_nitracline_raw["WOA23"])
basin_order_models = {
    m: ordered_basins(basin_mean_nitracline_raw[m]) for m in MODEL_NAMES
}
basin_order_matches_woa = {
    m: (basin_order_models[m] == basin_order_woa) for m in MODEL_NAMES
}
