import os
import glob
import re
import warnings
warnings.filterwarnings("ignore")

WOA_ROOT = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASKS = "/data1/lorraine-20419/ocean_data/WOA23/masks"
CMIP6_HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
CMIP6_SSP_ROOT = "/data1/lorraine-20419/benchmark_cmip6/ssp245"

R_EARTH = 6371000.0
BASIN_NAMES = {1: "Atlantic", 2: "Pacific", 3: "Indian", 4: "Southern", 5: "Arctic"}

# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------
def parse_variant(fname):
    b = os.path.basename(fname)
    m = re.search(r"_(r\d+i\d+p\d+f\d+)_", b)
    return m.group(1) if m else None

def nearest_index_1d(vals, centers):
    vals = np.asarray(vals, dtype=float)
    centers = np.asarray(centers, dtype=float)
    idx = np.searchsorted(centers, vals)
    idx = np.clip(idx, 1, len(centers) - 1)
    left = centers[idx - 1]
    right = centers[idx]
    choose = np.where(np.abs(vals - left) <= np.abs(vals - right), idx - 1, idx)
    return choose.astype(int)

def compute_dz(depth):
    depth = np.asarray(depth, dtype=float)
    if depth.size == 1:
        return np.array([max(depth[0], 1.0)])
    edges = np.empty(depth.size + 1, dtype=float)
    edges[0] = 0.0
    edges[1:-1] = 0.5 * (depth[:-1] + depth[1:])
    edges[-1] = depth[-1] + 0.5 * (depth[-1] - depth[-2])
    return np.diff(edges)

def compute_area(lat, lon, R=R_EARTH):
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    if lat.size > 1:
        lat_edges = np.empty(lat.size + 1)
        lat_edges[1:-1] = 0.5 * (lat[:-1] + lat[1:])
        lat_edges[0] = lat[0] - 0.5 * (lat[1] - lat[0])
        lat_edges[-1] = lat[-1] + 0.5 * (lat[-1] - lat[-2])
    else:
        lat_edges = np.array([lat[0] - 0.5, lat[0] + 0.5])
    if lon.size > 1:
        lon_edges = np.empty(lon.size + 1)
        lon_edges[1:-1] = 0.5 * (lon[:-1] + lon[1:])
        lon_edges[0] = lon[0] - 0.5 * (lon[1] - lon[0])
        lon_edges[-1] = lon[-1] + 0.5 * (lon[-1] - lon[-2])
    else:
        lon_edges = np.array([lon[0] - 0.5, lon[0] + 0.5])
    dphi = np.deg2rad(np.diff(lon_edges))
    dsin = np.diff(np.sin(np.deg2rad(lat_edges)))
    return R * R * np.outer(dsin, dphi)

def read_msk_table(path):
    with open(path, "r") as f:
        lines = f.readlines()
    start = 0
    for i, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("%"):
            continue
        tok = s.replace(",", " ").split()[0]
        try:
            float(tok)
            start = i
            break
        except ValueError:
            pass
    # Preferred parser
    try:
        df = pd.read_csv(
            path, delim_whitespace=True, comment="#", header=None, skiprows=start,
            na_values=["", "NaN", "nan", "****"]
        )
        arr = df.to_numpy(dtype=float)
        if arr.ndim == 2 and arr.shape[1] >= 3:
            return arr
    except Exception:
        pass
    # Fallback parser, comma or whitespace
    rows = []
    for line in lines[start:]:
        s = line.strip()
        if not s or s.startswith("#") or s.startswith("%"):
            continue
        s = s.replace(",", " ")
        toks = s.split()
        if not toks:
            continue
        try:
            row = [float(t) for t in toks]
        except ValueError:
            continue
        rows.append(row)
    if not rows:
        return np.empty((0, 0))
    maxlen = max(len(r) for r in rows)
    arr = np.full((len(rows), maxlen), np.nan)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = r
    return arr

def build_woa_land_mask(land_path, woa_lat, woa_lon, nz):
    try:
        arr = read_msk_table(land_path)
        if arr.shape[1] < 3:
            return None
        latvals = arr[:, 0]
        lonvals = arr[:, 1]
        bottom = arr[:, 2]
        iy = nearest_index_1d(latvals, woa_lat)
        lon180 = ((lonvals + 180.0) % 360.0) - 180.0
        ix = nearest_index_1d(lon180, woa_lon)
        bottom_grid = np.full((len(woa_lat), len(woa_lon)), np.nan)
        for a, b, c in zip(iy, ix, bottom):
            if np.isfinite(c):
                bottom_grid[a, b] = c
        ocean = np.zeros((nz, len(woa_lat), len(woa_lon)), dtype=bool)
        for k in range(nz):
            ocean[k] = (bottom_grid > 1) & ((bottom_grid - 1) >= k)
        return ocean
    except Exception:
        return None

def build_woa_basin(basin_path, woa_lat, woa_lon, nz):
    arr = read_msk_table(basin_path)
    if arr.shape[1] < 3:
        return None, None
    latvals = arr[:, 0]
    lonvals = arr[:, 1]
    levs = arr[:, 2:]
    iy = nearest_index_1d(latvals, woa_lat)
    lon180 = ((lonvals + 180.0) % 360.0) - 180.0
    ix = nearest_index_1d(lon180, woa_lon)
    nzb = min(nz, levs.shape[1])
    basin3d = np.full((nz, len(woa_lat), len(woa_lon)), np.nan)
    for k in range(nzb):
        vals = levs[:, k]
        ok = np.isfinite(vals) & (vals > 0)
        basin3d[k, iy[ok], ix[ok]] = vals[ok]
    basin2d = np.full((len(woa_lat), len(woa_lon)), np.nan)
    for k in range(nzb - 1, -1, -1):
        vals = basin3d[k]
        ok = np.isfinite(vals) & (vals > 0)
        basin2d[ok] = vals[ok]
    return basin3d, basin2d

def map_basin_to_grid(basin2d, lat2d, lon2d, woa_lat, woa_lon):
    iy = nearest_index_1d(lat2d.ravel(), woa_lat).reshape(lat2d.shape)
    lon180 = ((lon2d + 180.0) % 360.0) - 180.0
    ix = nearest_index_1d(lon180.ravel(), woa_lon).reshape(lon2d.shape)
    return basin2d[iy, ix]

def open_woa_var(prefix, subdir, tag, nn=0):
    path = os.path.join(WOA_ROOT, subdir, f"woa23_{tag}_{prefix}{nn:02d}_01.nc")
    ds = xr.open_dataset(path, decode_times=False)
    candidates = [
        f"{prefix}_an", f"{prefix}an", f"{prefix.upper()}_an", f"{prefix.lower()}_an", "an"
    ]
    vname = None
    for c in candidates:
        if c in ds.data_vars:
            vname = c
            break
    if vname is None:
        for n in ds.data_vars:
            if "an" in n.lower() and not n.endswith("_bnds"):
                vname = n
                break
    da = ds[vname]
    if "time" in da.dims:
        da = da.isel(time=0)
    da = da.squeeze(drop=True)
    return da, ds

def density_kg_m3(theta_c, salt, depth_m, lon2d, lat2d):
    # theta_c: (nz, ny, nx), salt: same, depth_m: (nz,)
    try:
        import gsw
        p3 = np.broadcast_to(depth_m[:, None, None], theta_c.shape)
        lon3 = np.broadcast_to(lon2d[None, :, :], theta_c.shape)
        lat3 = np.broadcast_to(lat2d[None, :, :], theta_c.shape)
        SA = gsw.SA_from_SP(salt, p3, lon3, lat3)
        rho = gsw.rho(SA, theta_c, p3)
        return rho
    except Exception:
        T = theta_c
        S = salt
        rho_w = (
            999.842594 + 6.793952e-2 * T - 9.095290e-3 * T**2
            + 1.001685e-4 * T**3 - 1.120083e-6 * T**4 + 6.536332e-9 * T**5
        )
        A = (
            0.824493 - 4.0899e-3 * T + 7.6438e-5 * T**2
            - 8.2467e-7 * T**3 + 5.3875e-9 * T**4
        )
        B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
        C = 4.8314e-4
        rho = rho_w + A * S + B * np.power(S, 1.5) + C * S**2
        p3 = np.broadcast_to(depth_m[:, None, None], T.shape)
        rho = rho + 0.0045 * p3
        return rho

# ---------------------------------------------------------------------------
# WOA23 oxygen, AOU, saturation
# ---------------------------------------------------------------------------
woa_o2_da, woa_o2_ds = open_woa_var(
    "o", "oxygen/1971-2000_climate_normal", "decav71A0", nn=0
)
woa_aou_da, _ = open_woa_var(
    "A", "apparent_oxygen_utilization/1971-2000_climate_normal", "decav71A0", nn=0
)
woa_sat_da, _ = open_woa_var(
    "O", "percent_oxygen_saturation/1971-2000_climate_normal", "decav71A0", nn=0
)

woa_depth = np.asarray(woa_o2_ds["depth"].values, dtype=float)
woa_lat = np.asarray(woa_o2_ds["lat"].values, dtype=float)
woa_lon = np.asarray(woa_o2_ds["lon"].values, dtype=float)
woa_nz = len(woa_depth)

woa_land_mask = build_woa_land_mask(
    os.path.join(WOA_MASKS, "landsea_01.msk"), woa_lat, woa_lon, woa_nz
)
woa_basin3d, woa_basin2d = build_woa_basin(
    os.path.join(WOA_MASKS, "basinmask_01.msk"), woa_lat, woa_lon, woa_nz
)

def analyze_woa(o2_da, depth_vals, lat_vals, lon_vals, land_mask, basin3d):
    o2 = o2_da.transpose("depth", "lat", "lon").values.astype(float)
    nz, ny, nx = o2.shape
    depth_vals = np.asarray(depth_vals, dtype=float)[:nz]
    dz = compute_dz(depth_vals)
    area = compute_area(lat_vals, lon_vals)
    vol = area[None, :, :] * dz[:, None, None]
    ocean = np.isfinite(o2)
    if land_mask is not None:
        ocean &= land_mask[:nz, :, :]

    vol20 = np.nansum(np.where(ocean & (o2 < 20.0), vol, 0.0))
    vol60 = np.nansum(np.where(ocean & (o2 < 60.0), vol, 0.0))

    thick60 = np.nansum(
        np.where(ocean & (o2 < 60.0), dz[:, None, None], 0.0), axis=0
    )
    idx = np.nanargmax(thick60)
    iy, ix = np.unravel_index(idx, thick60.shape)
    lon2d, lat2d = np.meshgrid(lon_vals, lat_vals)
    lat_max = float(lat2d[iy, ix])
    lon_max = float(((lon2d[iy, ix] + 180.0) % 360.0) - 180.0)

    basin_vols = {}
    if basin3d is not None:
        b3 = basin3d[:nz, :, :]
        for code in [1, 2, 3]:
            basin_vols[code] = (
                np.nansum(
                    np.where(ocean & (o2 < 60.0) & (b3 == code), vol, 0.0)
                ) / 1e15
            )

    return dict(
        vol_below20_1e15=vol20 / 1e15,
        vol_below60_1e15=vol60 / 1e15,
        max_thickness60_m=float(np.nanmax(thick60)),
        max_lat=lat_max,
        max_lon=lon_max,
        basin_vol60_1e15=basin_vols,
        thickness60_m=thick60,
        lat2d=lat2d,
        lon2d=lon2d,
        o2_umolkg=o2,
        vol_m3=vol,
        ocean_mask=ocean,
    )

woa_res = analyze_woa(
    woa_o2_da, woa_depth, woa_lat, woa_lon, woa_land_mask, woa_basin3d
)
woa_res["basin_vol60_named"] = {
    BASIN_NAMES[k]: v for k, v in woa_res["basin_vol60_1e15"].items()
}

aou_arr = woa_aou_da.transpose("depth", "lat", "lon").values.astype(float)
sat_arr = woa_sat_da.transpose("depth", "lat", "lon").values.astype(float)
mask60_woa = woa_res["ocean_mask"] & (woa_res["o2_umolkg"] < 60.0)
w_omz = woa_res["vol_m3"] * mask60_woa
woa_res["omz_aou_umolkg"] = float(
    np.nansum(aou_arr * w_omz) / np.nansum(w_omz)
)
woa_res["omz_o2sat_percent"] = float(
    np.nansum(sat_arr * w_omz) / np.nansum(w_omz)
)

omz_results = {"WOA23": woa_res}
hist_maps = {
    "WOA23": {
        "thickness60": woa_res["thickness60_m"],
        "lat": woa_res["lat2d"],
        "lon": woa_res["lon2d"],
    }
}
ssp_maps = {}
change_maps = {}

# ---------------------------------------------------------------------------
# CMIP6 loaders
# ---------------------------------------------------------------------------
MODEL_CFG = {
    "CESM2": dict(
        y_dim="nlat", x_dim="nlon", lat_coord="lat", lon_coord="lon",
        depth_coord="lev", depth_scale=0.01
    ),
    "GFDL-ESM4": dict(
        y_dim="lat", x_dim="lon", lat_coord="lat", lon_coord="lon",
        depth_coord="lev", depth_scale=1.0
    ),
    "MPI-ESM1-2-LR": dict(
        y_dim="j", x_dim="i", lat_coord="latitude", lon_coord="longitude",
        depth_coord="lev", depth_scale=1.0
    ),
}

def select_files(root, model, var, exp, preferred_variant=None, fallback_root=None):
    files = sorted(glob.glob(os.path.join(root, model, var, "*.nc")))
    used_root = root
    if not files and fallback_root is not None:
        files = sorted(glob.glob(os.path.join(fallback_root, model, var, "*.nc")))
        used_root = fallback_root
    if not files:
        return None, None, None
    groups = {}
    for f in files:
        v = parse_variant(f) or "unknown"
        groups.setdefault(v, []).append(f)
    if preferred_variant is not None and preferred_variant in groups:
        v = preferred_variant
    else:
        v = max(groups, key=lambda k: len(groups[k]))
    return sorted(groups[v]), v, used_root

def open_cmip_da(
    root, model, var, exp, y0=None, y1=None,
    preferred_variant=None, fallback_root=None, mean_time=False
):
    files, variant, _ = select_files(
        root, model, var, exp, preferred_variant, fallback_root
    )
    if not files:
        return None, None
    if var in ["areacello", "deptho", "sftof"] and len(files) > 1:
        files = files[:1]
    if len(files) == 1:
        ds = xr.open_dataset(files[0], use_cftime=True, chunks="auto")
    else:
        ds = xr.open_mfdataset(
            files, combine="by_coords", use_cftime=True, chunks="auto"
        )
    vname = var if var in ds else None
    if vname is None:
        for n in ds.data_vars:
            if n.lower() == var.lower():
                vname = n
                break
    if vname is None:
        candidates = [
            n for n in ds.data_vars
            if not n.endswith("_bnds")
            and n not in ["time_bnds", "lat_bnds", "lon_bnds", "lev_bnds"]
        ]
        if not candidates:
            return None, variant
        vname = candidates[0]
    da = ds[vname]
    if "time" in da.dims and y0 is not None:
        yrs = da.time.dt.year
        mask = (yrs >= y0) & (yrs <= y1)
        if bool(mask.any().compute()):
            da = da.where(mask, drop=True)
    if mean_time and "time" in da.dims:
        da = da.mean("time")
    return da, variant

def to_3d(da, cfg):
    da = da.squeeze(drop=True)
    dims = [cfg["depth_coord"], cfg["y_dim"], cfg["x_dim"]]
    da = da.transpose(*dims)
    return np.asarray(da.values, dtype=float)

def to_2d(da, cfg):
    da = da.squeeze(drop=True)
    da = da.transpose(cfg["y_dim"], cfg["x_dim"])
    return np.asarray(da.values, dtype=float)

def get_latlon_2d(da, cfg):
    lat = da[cfg["lat_coord"]]
    lon = da[cfg["lon_coord"]]
    if lat.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon.values, lat.values)
    else:
        lat2d = np.asarray(lat.values)
        lon2d = np.asarray(lon.values)
    return lat2d, lon2d

def load_model_fields(root, hist_root, model, exp, y0, y1):
    cfg = MODEL_CFG[model]

    o2_da, o2_variant = open_cmip_da(
        root, model, "o2", exp, y0, y1, mean_time=False
    )
    if o2_da is None:
        return None
    theta_da, _ = open_cmip_da(
        root, model, "thetao", exp, y0, y1,
        preferred_variant=o2_variant, mean_time=False
    )
    so_da, _ = open_cmip_da(
        root, model, "so", exp, y0, y1,
        preferred_variant=o2_variant, mean_time=False
    )
    if theta_da is None or so_da is None:
        return None

    o2_da, theta_da, so_da = xr.align(o2_da, theta_da, so_da, join="inner")
    if o2_da.sizes.get("time", 0) == 0:
        return None

    o2_mean = o2_da.mean("time").compute()
    theta_mean = theta_da.mean("time").compute()
    so_mean = so_da.mean("time").compute()

    vol_da, _ = open_cmip_da(
        root, model, "volcello", exp, y0, y1,
        preferred_variant=o2_variant, fallback_root=hist_root, mean_time=True
    )
    if vol_da is None:
        return None

    area_da, _ = open_cmip_da(
        root, model, "areacello", exp, None, None,
        fallback_root=hist_root, mean_time=True
    )
    if area_da is None:
        return None

    o2_arr = to_3d(o2_mean, cfg)
    theta_arr = to_3d(theta_mean, cfg)
    so_arr = to_3d(so_mean, cfg)
    vol_arr = to_3d(vol_da, cfg)
    area_arr = to_2d(area_da, cfg)
    lat2d, lon2d = get_latlon_2d(o2_mean, cfg)
    depth = np.asarray(o2_mean[cfg["depth_coord"]].values, dtype=float)
    depth = depth * cfg["depth_scale"]

    o2_units = str(o2_mean.attrs.get("units", "")).lower()
    o2_is_mass = ("kg" in o2_units)

    theta_units = str(theta_mean.attrs.get("units", ""))
    if (("K" in theta_units and "deg" not in theta_units.lower())
            or (np.isfinite(np.nanmean(theta_arr)) and np.nanmean(theta_arr) > 100)):
        theta_arr = theta_arr - 273.15

    return dict(
        model=model, exp=exp,
        o2_arr=o2_arr, o2_is_mass=o2_is_mass,
        theta_arr=theta_arr, so_arr=so_arr,
        vol_arr=vol_arr, area_arr=area_arr,
        depth=depth, lat2d=lat2d, lon2d=lon2d,
        o2_units=o2_units
    )

def analyze_model(fields, basin2d_woa=None, woa_lat=None, woa_lon=None):
    o2 = fields["o2_arr"]
    theta = fields["theta_arr"]
    salt = fields["so_arr"]
    vol = fields["vol_arr"]
    area = fields["area_arr"]
    depth = fields["depth"]
    lat2d = fields["lat2d"]
    lon2d = fields["lon2d"]

    rho = density_kg_m3(theta, salt, depth, lon2d, lat2d)

    if fields["o2_is_mass"]:
        o2_umolkg = o2.copy()
        o2_mol = o2_umolkg * rho / 1e6
    else:
        o2_umolkg = o2 * 1e6 / rho
        med = np.nanmedian(o2_umolkg)
        if np.isfinite(med):
            if med < 1.0:
                o2_umolkg = o2_umolkg * 1000.0
            elif med > 10000.0:
                # Defensive fallback for files whose numerical values are already per-mass.
                o2_umolkg = o2.copy()
        o2_mol = o2_umolkg * rho / 1e6

    ocean = np.isfinite(o2_umolkg) & np.isfinite(vol) & (vol > 0)

    vol20 = np.nansum(np.where(ocean & (o2_umolkg < 20.0), vol, 0.0))
    vol60 = np.nansum(np.where(ocean & (o2_umolkg < 60.0), vol, 0.0))

    area3 = np.broadcast_to(area[None, :, :], vol.shape)
    thick_lev = np.divide(
        vol, area3, out=np.zeros_like(vol, dtype=float), where=(area3 > 0)
    )
    thick60 = np.nansum(
        np.where(ocean & (o2_umolkg < 60.0), thick_lev, 0.0), axis=0
    )
    idx = np.nanargmax(thick60)
    iy, ix = np.unravel_index(idx, thick60.shape)
    lat_max = float(lat2d[iy, ix])
    lon_max = float(((lon2d[iy, ix] + 180.0) % 360.0) - 180.0)

    inventory_mol = np.nansum(np.where(ocean, o2_mol * vol, 0.0))

    basin_vols = {}
    if basin2d_woa is not None:
        basin = map_basin_to_grid(basin2d_woa, lat2d, lon2d, woa_lat, woa_lon)
        for code in [1, 2, 3]:
            basin_vols[code] = (
                np.nansum(
                    np.where(
                        ocean & (o2_umolkg < 60.0) & (basin[None, :, :] == code),
                        vol, 0.0
                    )
                ) / 1e15
            )

    return dict(
        vol_below20_1e15=vol20 / 1e15,
        vol_below60_1e15=vol60 / 1e15,
        max_thickness60_m=float(np.nanmax(thick60)),
        max_lat=lat_max,
        max_lon=lon_max,
        inventory_mol=float(inventory_mol),
        basin_vol60_1e15=basin_vols,
        thickness60_m=thick60,
        lat2d=lat2d,
        lon2d=lon2d,
        o2_umolkg=o2_umolkg,
        vol_m3=vol,
        ocean_mask=ocean,
    )

# ---------------------------------------------------------------------------
# Historical CMIP6: GFDL-ESM4 and MPI-ESM1-2-LR, 1971-2000
# ---------------------------------------------------------------------------
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    fields = load_model_fields(
        CMIP6_HIST_ROOT, None, model, "historical", 1971, 2000
    )
    if fields is None:
        continue
    res = analyze_model(fields, woa_basin2d, woa_lat, woa_lon)
    res["basin_vol60_named"] = {
        BASIN_NAMES[k]: v for k, v in res["basin_vol60_1e15"].items()
    }
    key = f"{model}_historical"
    omz_results[key] = res
    hist_maps[model] = {
        "thickness60": res["thickness60_m"],
        "lat": res["lat2d"],
        "lon": res["lon2d"],
    }

# ---------------------------------------------------------------------------
# ssp245 CMIP6: CESM2, GFDL-ESM4, MPI-ESM1-2-LR, 2071-2100
# ---------------------------------------------------------------------------
for model in ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]:
    fields = load_model_fields(
        CMIP6_SSP_ROOT, CMIP6_HIST_ROOT, model, "ssp245", 2071, 2100
    )
    if fields is None:
        continue
    res = analyze_model(fields, woa_basin2d, woa_lat, woa_lon)
    res["basin_vol60_named"] = {
        BASIN_NAMES[k]: v for k, v in res["basin_vol60_1e15"].items()
    }
    key = f"{model}_ssp245"
    omz_results[key] = res
    ssp_maps[model] = {
        "thickness60": res["thickness60_m"],
        "lat": res["lat2d"],
        "lon": res["lon2d"],
    }

# ---------------------------------------------------------------------------
# Changes and summary tables
# ---------------------------------------------------------------------------
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    if model in hist_maps and model in ssp_maps:
        h = hist_maps[model]
        s = ssp_maps[model]
        change_maps[model] = {
            "thickness60_change": s["thickness60"] - h["thickness60"],
            "lat": s["lat"],
            "lon": s["lon"],
        }

omz_change_pct = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    hkey = f"{model}_historical"
    skey = f"{model}_ssp245"
    if hkey in omz_results and skey in omz_results:
        omz_change_pct[model] = {}
        for thr in [20, 60]:
            hval = omz_results[hkey][f"vol_below{thr}_1e15"]
            sval = omz_results[skey][f"vol_below{thr}_1e15"]
            omz_change_pct[model][str(thr)] = (
                (sval - hval) / hval * 100.0 if hval != 0 else np.nan
            )

woa_vol60 = omz_results["WOA23"]["vol_below60_1e15"]
hist_model_bias_60 = {}
for model in ["GFDL-ESM4", "MPI-ESM1-2-LR"]:
    key = f"{model}_historical"
    if key in omz_results:
        hist_model_bias_60[model] = (
            omz_results[key]["vol_below60_1e15"] - woa_vol60
        )
better_model_60 = (
    min(hist_model_bias_60, key=lambda m: abs(hist_model_bias_60[m]))
    if hist_model_bias_60 else None
)

summary_rows = []
for key, res in omz_results.items():
    summary_rows.append(dict(
        source=key,
        vol_below20_1e15=res["vol_below20_1e15"],
        vol_below60_1e15=res["vol_below60_1e15"],
        max_thickness60_m=res["max_thickness60_m"],
        max_lat=res["max_lat"],
        max_lon=res["max_lon"],
        inventory_mol=res.get("inventory_mol", np.nan),
    ))
omz_summary = pd.DataFrame(summary_rows)

basin_rows = []
for key in ["WOA23", "GFDL-ESM4_historical", "MPI-ESM1-2-LR_historical"]:
    if key in omz_results:
        for bname in ["Pacific", "Indian", "Atlantic"]:
            basin_rows.append(dict(
                source=key,
                basin=bname,
                vol60_1e15=omz_results[key]["basin_vol60_named"].get(bname, np.nan),
            ))
basin_vol_table = pd.DataFrame(basin_rows)

woa_omz_aou_umolkg = woa_res["omz_aou_umolkg"]
woa_omz_o2sat_percent = woa_res["omz_o2sat_percent"]
