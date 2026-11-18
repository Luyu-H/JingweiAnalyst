import os, glob, re, warnings
warnings.filterwarnings('ignore')

try:
    import numpy as np
except Exception:
    pass
try:
    import xarray as xr
except Exception:
    pass
try:
    import pandas as pd
except Exception:
    pass

# -----------------------------------------------------------------------------
# Constants and configuration
# -----------------------------------------------------------------------------
RHO = 1025.0
R_EARTH = 6371000.0
CMIP_ROOT = '/data1/lorraine-20419/benchmark_cmip6/historical'
WOA_ROOT = '/data1/lorraine-20419/ocean_data/WOA23/climatology'
WOA_MASK = '/data1/lorraine-20419/ocean_data/WOA23/masks/basinmask_01.msk'

SOURCES = ['WOA23', 'CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']
MODELS = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']

MODEL_CFG = {
    'CESM2': dict(lat='lat', lon='lon', depth='lev', depth_scale=0.01),
    'GFDL-ESM4': dict(lat='lat', lon='lon', depth='lev', depth_scale=1.0),
    'MPI-ESM1-2-LR': dict(lat='latitude', lon='longitude', depth='lev', depth_scale=1.0),
}

def parse_variant(fn):
    m = re.search(r'(r\d+i\d+p\d+f\d+)', os.path.basename(fn))
    return m.group(1) if m else None

def files_in_dir(d):
    return sorted(glob.glob(os.path.join(d, '*.nc')))

def choose_variant_files(d, desired=None):
    fs = files_in_dir(d)
    if not fs:
        raise FileNotFoundError(d)
    vs = {parse_variant(f) for f in fs}
    vs.discard(None)
    if desired is not None and desired in vs:
        return [f for f in fs if parse_variant(f) == desired]
    if 'r1i1p1f1' in vs:
        return [f for f in fs if parse_variant(f) == 'r1i1p1f1']
    if vs:
        counts = {v: sum(parse_variant(f) == v for f in fs) for v in vs}
        v = max(counts, key=counts.get)
        return [f for f in fs if parse_variant(f) == v]
    return fs

def open_nc_maybe_mf(files):
    if len(files) == 1:
        try:
            return xr.open_dataset(files[0], use_cftime=True, chunks={'time': 12})
        except Exception:
            return xr.open_dataset(files[0], use_cftime=True)
    try:
        return xr.open_mfdataset(files, combine='by_coords', use_cftime=True,
                                  chunks={'time': 12})
    except Exception:
        return xr.open_mfdataset(files, combine='by_coords', use_cftime=True)

def mean_period(da, y0=1971, y1=2000):
    """Robust 1971-2000 mean that works with cftime indexes."""
    if 'time' not in da.dims:
        return da
    yr = None
    try:
        tidx = da.indexes.get('time', None)
        if tidx is not None and hasattr(tidx, 'year'):
            yr = np.asarray(tidx.year)
    except Exception:
        yr = None
    if yr is None:
        try:
            yr = np.asarray(da.time.dt.year.values)
        except Exception:
            yr = None
    if yr is not None:
        sel = (yr >= y0) & (yr <= y1)
        if np.any(sel):
            return da.isel(time=np.where(sel)[0]).mean('time')
    return da.mean('time')

def find_dim(da, names):
    for n in names:
        if n in da.dims:
            return n
    for d in da.dims:
        ld = d.lower()
        for n in names:
            if n in ld:
                return d
    return None

def get_latlon_2d(ds, cfg, fallback=None):
    try:
        lat = ds[cfg['lat']]
        lon = ds[cfg['lon']]
    except Exception:
        ds = fallback
        lat = ds[cfg['lat']]
        lon = ds[cfg['lon']]
    if lat.ndim == 1 and lon.ndim == 1:
        lat2d, lon2d = xr.broadcast(lat, lon)
    else:
        lat2d, lon2d = lat, lon
    return np.asarray(lat2d.values, dtype=float), np.asarray(lon2d.values, dtype=float)

def unit_umol_per_kg(units, rho=RHO):
    """Return multiplication factor to convert a concentration to umol kg-1."""
    if units is None:
        return None
    u = str(units).strip().lower()
    u = (u.replace('μ', 'u').replace('µ', 'u')
           .replace('−', '-').replace('⁻', '-')
           .replace('³', '3').replace('²', '2'))
    u = u.replace(' ', '').replace('**', '^')
    per_vol = ('m-3' in u) or ('m^-3' in u) or ('/m3' in u) or ('m3-1' in u)
    per_kg = ('kg-1' in u) or ('kg^-1' in u) or ('/kg' in u)
    if u.startswith('nmol'):
        amt = 1.0e-3
    elif u.startswith('umol') or u.startswith('umole'):
        amt = 1.0
    elif u.startswith('mmol'):
        amt = 1.0e3
    elif u.startswith('mol'):
        amt = 1.0e6
    else:
        return None
    if per_vol:
        return amt / rho
    if per_kg:
        return amt
    return None

# -----------------------------------------------------------------------------
# WOA readers
# -----------------------------------------------------------------------------
def load_woa_field(varkey, prefix):
    if varkey == 'nitrate':
        path = os.path.join(WOA_ROOT, 'nitrate', 'all_time_periods', 'woa23_all_n00_01.nc')
    elif varkey == 'phosphate':
        path = os.path.join(WOA_ROOT, 'phosphate', 'all_time_periods', 'woa23_all_p00_01.nc')
    else:
        raise ValueError(varkey)
    ds = xr.open_dataset(path, decode_times=False)

    cands = []
    for v in ds.data_vars:
        vl = v.lower()
        if ('_an' in vl or vl.endswith('an') or vl == 'an') and vl.startswith(prefix.lower()):
            cands.append(v)
    if not cands:
        cands = [v for v in ds.data_vars if '_an' in v.lower() or v.lower().endswith('an')]
    if not cands:
        cands = [v for v in ds.data_vars if ds[v].ndim >= 3]
    var = cands[0]
    da = ds[var]
    if 'time' in da.dims:
        da = da.isel(time=0)
    da = da.squeeze(drop=True)
    return ds, da

def woa_layer_fields(n, p):
    depth_name = find_dim(n, ['depth', 'lev', 'z'])
    lat_name = find_dim(n, ['lat', 'latitude'])
    lon_name = find_dim(n, ['lon', 'longitude'])
    if depth_name is None or lat_name is None or lon_name is None:
        raise ValueError('Could not identify WOA dimensions')

    n = n.transpose(depth_name, lat_name, lon_name)
    p = p.transpose(depth_name, lat_name, lon_name)

    depth = np.asarray(n[depth_name].values, dtype=float)
    lat = np.asarray(n[lat_name].values, dtype=float)
    lon = np.asarray(n[lon_name].values, dtype=float)

    if depth[0] > depth[-1]:
        n = n.isel({depth_name: slice(None, None, -1)})
        p = p.isel({depth_name: slice(None, None, -1)})
        depth = depth[::-1]

    if len(depth) > 1:
        edges = np.empty(len(depth) + 1)
        edges[1:-1] = 0.5 * (depth[:-1] + depth[1:])
        edges[0] = max(0.0, depth[0] - 0.5 * (depth[1] - depth[0]))
        edges[-1] = depth[-1] + 0.5 * (depth[-1] - depth[-2])
    else:
        edges = np.array([0.0, depth[0] + 1.0])

    thick_layer = np.clip(
        np.minimum(edges[1:], 1000.0) - np.maximum(edges[:-1], 200.0),
        0.0, None
    )

    dlat = np.deg2rad(np.median(np.abs(np.diff(lat)))) if len(lat) > 1 else np.deg2rad(1.0)
    dlon = np.deg2rad(np.median(np.abs(np.diff(lon)))) if len(lon) > 1 else np.deg2rad(1.0)
    area = R_EARTH**2 * np.cos(np.deg2rad(lat)) * dlat * dlon
    area2d = np.broadcast_to(area[:, None], (len(lat), len(lon)))

    w = thick_layer[:, None, None] * area2d[None, :, :]
    nv = np.asarray(n.values, dtype=float)
    pv = np.asarray(p.values, dtype=float)

    valid = np.isfinite(nv) & np.isfinite(pv) & (w > 0)
    wj = np.where(valid, w, 0.0)
    vol_layer = wj.sum(axis=0)

    no3_layer = np.divide(
        (np.where(valid, nv, 0.0) * wj).sum(axis=0),
        vol_layer,
        out=np.full_like(vol_layer, np.nan),
        where=vol_layer > 0
    )
    po4_layer = np.divide(
        (np.where(valid, pv, 0.0) * wj).sum(axis=0),
        vol_layer,
        out=np.full_like(vol_layer, np.nan),
        where=vol_layer > 0
    )

    lat2d, lon2d = np.meshgrid(lat, lon, indexing='ij')
    return no3_layer, po4_layer, vol_layer, lat2d, lon2d, lat, lon

# -----------------------------------------------------------------------------
# WOA basin mask
# -----------------------------------------------------------------------------
_NUM_RE = re.compile(r'[-+]?(?:(?:\d+\.\d*)|(?:\.\d+)|(?:\d+))(?:[EeDd][-+]?\d+)?')

def _parse_numeric_line(s):
    """Extract all numeric tokens from a line, tolerating labels/missing markers."""
    tokens = [m.group(0) for m in _NUM_RE.finditer(s)]
    vals = []
    for t in tokens:
        tt = t.replace('D', 'E').replace('d', 'e')
        try:
            v = float(tt)
        except ValueError:
            continue
        vals.append((v, t))
    # If a short row contains one very long integer token, it is almost
    # certainly a run of concatenated single-digit basin codes.
    if 0 < len(vals) <= 8:
        expanded = []
        for v, t in vals:
            if re.fullmatch(r'\d{15,}', t):
                expanded.extend(float(ch) for ch in t)
            else:
                expanded.append(v)
        return expanded
    return [v for v, _ in vals]

def load_basin_mask(path, nlev, lat_centers, lon_centers):
    rows = []
    with open(path, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s[0] in '#!%':
                continue
            vals = _parse_numeric_line(s)
            if len(vals) >= 3:
                rows.append(vals)

    if not rows:
        raise ValueError(f'No numeric rows found in {path}')

    lats = np.asarray(lat_centers, dtype=float)
    lons = np.asarray(lon_centers, dtype=float)
    nlat, nlon = len(lats), len(lons)
    basin_h = np.zeros((nlat, nlon), dtype=int)

    def idx(v, arr):
        if len(arr) <= 1:
            return 0
        step = arr[1] - arr[0]
        if step == 0:
            return 0
        return int(np.clip(int(round((v - arr[0]) / step)), 0, len(arr) - 1))

    # --- Identify the latitude/longitude columns robustly -------------------
    minlen = min(len(r) for r in rows)
    ncol = min(8, max(3, minlen))
    use = [r for r in rows if len(r) >= ncol]
    if not use:
        use = rows
        ncol = max(3, min(len(r) for r in rows))

    stats = {}
    for c in range(ncol):
        col = np.array([r[c] for r in use], dtype=float)
        fin = col[np.isfinite(col)]
        if fin.size == 0:
            stats[c] = None
        else:
            stats[c] = (float(fin.min()), float(fin.max()),
                        int(len(np.unique(np.round(fin, 4)))))

    lat_cands = [c for c in range(ncol)
                 if stats.get(c) is not None
                 and -90.0 <= stats[c][0] and stats[c][1] <= 90.0
                 and stats[c][2] > 20
                 and (stats[c][1] - stats[c][0]) > 20.0]
    lon_cands = [c for c in range(ncol)
                 if stats.get(c) is not None
                 and -360.0 <= stats[c][0] and stats[c][1] <= 360.0
                 and stats[c][2] > 20
                 and (stats[c][1] - stats[c][0]) > 20.0]

    latcol = lat_cands[0] if lat_cands else 0
    lon_cands = [c for c in lon_cands if c != latcol]
    if lon_cands:
        loncol = max(lon_cands, key=lambda c: stats[c][1] - stats[c][0])
    else:
        loncol = 1 if latcol == 0 else 0

    # --- Assign the surface/first valid basin code at each point ------------
    start = max(latcol, loncol) + 1
    for r in rows:
        if len(r) <= max(latcol, loncol):
            continue
        la, lo = r[latcol], r[loncol]
        if not (np.isfinite(la) and np.isfinite(lo)):
            continue

        code_seq = []
        for v in r[start:]:
            if not np.isfinite(v):
                continue
            if abs(v - round(v)) < 1e-6 and v >= 0.0:
                iv = int(round(v))
                if iv < 50:
                    code_seq.append(iv)
                else:
                    # concatenated run of single-digit codes
                    code_seq.extend(int(ch) for ch in str(iv))
        code = next((c for c in code_seq if 0 < c < 50), 0)
        if code > 0:
            i = idx(la, lats)
            j = idx(lo, lons)
            basin_h[i, j] = code
    return basin_h

def fill_basin_nearest(b):
    b = b.astype(int).copy()
    for _ in range(60):
        if not np.any(b == 0):
            break
        neigh = []
        for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1),
                       (-1, -1), (-1, 1), (1, -1), (1, 1)]:
            neigh.append(np.roll(b, shift=(di, dj), axis=(0, 1)))
        fill = np.zeros_like(b)
        for s in neigh:
            fill = np.where((fill == 0) & (s > 0), s, fill)
        if not np.any(fill > 0):
            break
        b = np.where(b == 0, fill, b)
    return b

def identify_basins(basin_h, lat2d, lon2d):
    """Classify mask codes into Arctic/Southern/Pacific/Atlantic/Indian.

    The five largest codes are treated as the principal basins. Any remaining
    (marginal-sea) codes are folded into the nearest non-polar principal basin.
    """
    codes = [int(c) for c in np.unique(basin_h) if int(c) > 0]
    info = {}
    area = {}
    lat2d = np.asarray(lat2d, dtype=float)
    lon2d = np.asarray(lon2d, dtype=float)

    for c in codes:
        m = basin_h == c
        la = lat2d[m]
        lo = lon2d[m]
        if la.size == 0:
            continue
        w = np.cos(np.deg2rad(la))
        if np.sum(w) == 0:
            w = np.ones_like(la)
        mlat = float(np.average(la, weights=w))
        ang = np.deg2rad(lo)
        mlon = float(np.rad2deg(np.arctan2(
            np.average(np.sin(ang), weights=w),
            np.average(np.cos(ang), weights=w)
        )))
        info[c] = (mlat, mlon)
        area[c] = float(np.sum(w))

    names = {}
    usable = [c for c in codes if c in info]
    if not usable:
        return names, info

    main = sorted(usable, key=lambda c: area.get(c, 0.0), reverse=True)[:5]
    main_set = set(main)

    for c in main:
        mlat, mlon = info[c]
        if mlat > 60.0:
            names[c] = 'Arctic'
        elif mlat < -60.0:
            names[c] = 'Southern'
        else:
            lon360 = mlon % 360.0
            if 30.0 <= lon360 < 110.0:
                names[c] = 'Indian'
            elif 110.0 <= lon360 < 280.0:
                names[c] = 'Pacific'
            else:
                names[c] = 'Atlantic'

    parents = [c for c in main if names.get(c) in ('Atlantic', 'Pacific', 'Indian')]

    def gc_dist(a, b):
        la1, lo1 = np.deg2rad(a[0]), np.deg2rad(a[1])
        la2, lo2 = np.deg2rad(b[0]), np.deg2rad(b[1])
        dlon = lo2 - lo1
        d = (np.sin((la2 - la1) / 2.0) ** 2
             + np.cos(la1) * np.cos(la2) * np.sin(dlon / 2.0) ** 2)
        return float(np.arcsin(np.sqrt(np.clip(d, 0.0, 1.0))))

    for c in usable:
        if c in names:
            continue
        if parents:
            p = min(parents, key=lambda x: gc_dist(info[c], info[x]))
            names[c] = names[p]
        else:
            mlat, mlon = info[c]
            lon360 = mlon % 360.0
            if 30.0 <= lon360 < 110.0:
                names[c] = 'Indian'
            elif 110.0 <= lon360 < 280.0:
                names[c] = 'Pacific'
            else:
                names[c] = 'Atlantic'

    return names, info

def map_basin(lat2d, lon2d, basin_h, lat_woa, lon_woa):
    lon_w = ((lon2d + 180.0) % 360.0) - 180.0
    step_lat = lat_woa[1] - lat_woa[0] if len(lat_woa) > 1 else 1.0
    step_lon = lon_woa[1] - lon_woa[0] if len(lon_woa) > 1 else 1.0
    i = np.rint((lat2d - lat_woa[0]) / step_lat).astype(int)
    j = np.rint((lon_w - lon_woa[0]) / step_lon).astype(int)
    i = np.clip(i, 0, len(lat_woa) - 1)
    j = np.clip(j, 0, len(lon_woa) - 1)
    return basin_h[i, j]

# -----------------------------------------------------------------------------
# CMIP6 readers
# -----------------------------------------------------------------------------
def model_layer_fields(model, cfg):
    dir_n = os.path.join(CMIP_ROOT, model, 'no3')
    dir_p = os.path.join(CMIP_ROOT, model, 'po4')
    dir_v = os.path.join(CMIP_ROOT, model, 'volcello')

    fs_n_all = files_in_dir(dir_n)
    fs_p_all = files_in_dir(dir_p)
    vs_n = {parse_variant(f) for f in fs_n_all}
    vs_p = {parse_variant(f) for f in fs_p_all}
    vs_n.discard(None)
    vs_p.discard(None)
    common = vs_n & vs_p

    if common:
        desired = 'r1i1p1f1' if 'r1i1p1f1' in common else sorted(common)[0]
        fs_n = [f for f in fs_n_all if parse_variant(f) == desired]
        fs_p = [f for f in fs_p_all if parse_variant(f) == desired]
    else:
        desired = None
        fs_n, fs_p = fs_n_all, fs_p_all

    fs_v = choose_variant_files(dir_v, desired)

    ds_n = open_nc_maybe_mf(fs_n)
    ds_p = open_nc_maybe_mf(fs_p)
    ds_v = open_nc_maybe_mf(fs_v)

    no3_tm = mean_period(ds_n['no3'])
    po4_tm = mean_period(ds_p['po4'])
    vol_tm = mean_period(ds_v['volcello'])

    no3_tm, po4_tm, vol_tm = xr.align(no3_tm, po4_tm, vol_tm, join='inner')

    lev_name = cfg['depth']
    if lev_name not in no3_tm.dims:
        lev_name = find_dim(no3_tm, ['lev', 'depth'])

    depth_da = no3_tm[lev_name]
    depth_m = depth_da * cfg['depth_scale']
    depth_mask = (depth_m >= 200.0) & (depth_m <= 1000.0)

    valid3 = (
        np.isfinite(no3_tm) &
        np.isfinite(po4_tm) &
        np.isfinite(vol_tm) &
        (vol_tm > 0.0) &
        depth_mask
    )

    w3 = xr.where(valid3, vol_tm, 0.0)
    vol_layer = w3.sum(dim=lev_name)
    num_n = xr.where(valid3, no3_tm * vol_tm, 0.0).sum(dim=lev_name)
    num_p = xr.where(valid3, po4_tm * vol_tm, 0.0).sum(dim=lev_name)

    no3_layer = (num_n / vol_layer).where(vol_layer > 0.0)
    po4_layer = (num_p / vol_layer).where(vol_layer > 0.0)

    out = xr.Dataset({
        'no3': no3_layer,
        'po4': po4_layer,
        'vol': vol_layer,
    }).compute()

    # Convert model concentrations to umol kg-1 using the declared units.
    u_n = ds_n['no3'].attrs.get('units')
    u_p = ds_p['po4'].attrs.get('units')
    fac_n = unit_umol_per_kg(u_n)
    fac_p = unit_umol_per_kg(u_p)
    if fac_n is None:
        fac_n = 1.0e6 / RHO
    if fac_p is None:
        fac_p = 1.0e6 / RHO

    no3 = np.asarray(out['no3'].values, dtype=float) * fac_n
    po4 = np.asarray(out['po4'].values, dtype=float) * fac_p
    vol = np.asarray(out['vol'].values, dtype=float)

    lat2d, lon2d = get_latlon_2d(out, cfg, fallback=ds_n)
    return no3, po4, vol, lat2d, lon2d

# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------
def valid_mask(no3, po4, vol):
    return (
        np.isfinite(no3) & np.isfinite(po4) &
        (no3 >= 0.0) & (po4 > 0.0) &
        np.isfinite(vol) & (vol > 0.0)
    )

def compute_source_metrics(name, no3, po4, vol, basin_codes, basin_names):
    no3 = np.asarray(no3, dtype=float)
    po4 = np.asarray(po4, dtype=float)
    vol = np.asarray(vol, dtype=float)
    basin_codes = np.asarray(basin_codes)

    valid = valid_mask(no3, po4, vol)
    ratio = np.full(no3.shape, np.nan, dtype=float)
    ratio[valid] = no3[valid] / po4[valid]

    w = vol[valid]
    r = ratio[valid]
    x = po4[valid]
    y = no3[valid]
    total_w = float(np.sum(w))

    # Headline means: volume-weighted across water parcels.
    np_mean_cell = float(np.sum(w * r) / total_w) if total_w > 0 else np.nan
    np_ratio_means = float(np.sum(w * y) / np.sum(w * x)) if np.sum(w * x) > 0 else np.nan

    # Equal-cell versions retained for transparency.
    np_mean_cell_unw = float(np.nanmean(r)) if len(r) > 0 else np.nan
    np_ratio_means_unw = (
        float(np.nanmean(y) / np.nanmean(x))
        if len(x) > 0 and np.nanmean(x) > 0 else np.nan
    )

    if len(x) > 1:
        slope, intercept = np.polyfit(x, y, 1)
        yhat = slope * x + intercept
        ss_res = float(np.sum((y - yhat) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    else:
        slope = intercept = r2 = np.nan

    pct_14_18 = (
        100.0 * float(np.sum(w[(r > 14.0) & (r < 18.0)]) / total_w)
        if total_w > 0 else np.nan
    )
    pct_lt12 = (
        100.0 * float(np.sum(w[r < 12.0]) / total_w)
        if total_w > 0 else np.nan
    )

    # basenames can contain several mask codes (main basin + folded marginal seas)
    codes_by_name = {}
    for code, nm in basin_names.items():
        codes_by_name.setdefault(nm, []).append(int(code))

    out = {
        'source': name,
        'np_mean_cell': np_mean_cell,
        'np_ratio_means': np_ratio_means,
        'np_mean_cell_unweighted': np_mean_cell_unw,
        'np_ratio_means_unweighted': np_ratio_means_unw,
        'slope': float(slope),
        'intercept': float(intercept),
        'r2': float(r2),
        'pct_14_18': pct_14_18,
        'pct_lt12_global': pct_lt12,
    }

    for bname in ['Pacific', 'Atlantic', 'Indian']:
        codes_b = codes_by_name.get(bname, [])
        key_np = 'np_' + bname.lower()
        key_pct = 'pct_lt12_' + bname.lower()
        if not codes_b:
            out[key_np] = np.nan
            out[key_pct] = np.nan
            continue
        m = valid & np.isin(basin_codes, np.asarray(codes_b, dtype=int))
        if not np.any(m):
            out[key_np] = np.nan
            out[key_pct] = np.nan
            continue
        wb = vol[m]
        rb = ratio[m]
        tw = float(np.sum(wb))
        out[key_np] = float(np.sum(wb * rb) / tw) if tw > 0 else np.nan
        out[key_pct] = 100.0 * float(np.sum(wb[rb < 12.0]) / tw) if tw > 0 else np.nan

    if np.isfinite(np_mean_cell):
        if np_mean_cell < 15.95:
            out['redfield_position'] = 'below'
        elif np_mean_cell > 16.05:
            out['redfield_position'] = 'above'
        else:
            out['redfield_position'] = 'at'
    else:
        out['redfield_position'] = np.nan

    return out

# -----------------------------------------------------------------------------
# Read WOA23 nitrate and phosphate
# -----------------------------------------------------------------------------
dsn, n_woa = load_woa_field('nitrate', 'n')
dsp, p_woa = load_woa_field('phosphate', 'p')
n_woa, p_woa = xr.align(n_woa, p_woa, join='inner')

no3_woa, po4_woa, vol_woa, lat2d_woa, lon2d_woa, lat_woa, lon_woa = woa_layer_fields(n_woa, p_woa)

depth_name_woa = find_dim(n_woa, ['depth', 'lev', 'z'])
nlev_woa = int(n_woa.sizes[depth_name_woa])

basin_raw = load_basin_mask(WOA_MASK, nlev_woa, lat_woa, lon_woa)
basin_names, basin_info = identify_basins(basin_raw, lat2d_woa, lon2d_woa)
basin_h = fill_basin_nearest(basin_raw)
basin_code_by_name = {v: k for k, v in basin_names.items()}

# -----------------------------------------------------------------------------
# Build fields, scatter data and metrics
# -----------------------------------------------------------------------------
fields = {}
scatter_data = {}
metric_rows = []

# WOA23
m_woa = valid_mask(no3_woa, po4_woa, vol_woa)
ratio_woa = np.full(no3_woa.shape, np.nan, dtype=float)
ratio_woa[m_woa] = no3_woa[m_woa] / po4_woa[m_woa]

fields['WOA23'] = {
    'no3': no3_woa,
    'po4': po4_woa,
    'ratio': ratio_woa,
    'vol': vol_woa,
    'lat': lat2d_woa,
    'lon': lon2d_woa,
    'basin': basin_raw,
}
scatter_data['WOA23'] = (po4_woa[m_woa], no3_woa[m_woa])
metric_rows.append(compute_source_metrics(
    'WOA23', no3_woa, po4_woa, vol_woa, basin_raw, basin_names
))

# CMIP6 models
for model in MODELS:
    no3_m, po4_m, vol_m, lat_m, lon_m = model_layer_fields(model, MODEL_CFG[model])
    basin_m = map_basin(lat_m, lon_m, basin_h, lat_woa, lon_woa)

    mm = valid_mask(no3_m, po4_m, vol_m)
    ratio_m = np.full(no3_m.shape, np.nan, dtype=float)
    ratio_m[mm] = no3_m[mm] / po4_m[mm]

    fields[model] = {
        'no3': no3_m,
        'po4': po4_m,
        'ratio': ratio_m,
        'vol': vol_m,
        'lat': lat_m,
        'lon': lon_m,
        'basin': basin_m,
    }
    scatter_data[model] = (po4_m[mm], no3_m[mm])
    metric_rows.append(compute_source_metrics(
        model, no3_m, po4_m, vol_m, basin_m, basin_names
    ))

metrics_df = pd.DataFrame(metric_rows)

# Add explicit model-vs-WOA differences.
woa_row = metrics_df.loc[metrics_df['source'] == 'WOA23'].iloc[0]
for col in ['np_mean_cell', 'np_ratio_means',
            'pct_lt12_global', 'pct_lt12_pacific', 'pct_lt12_indian',
            'np_pacific', 'np_atlantic', 'np_indian']:
    if col in metrics_df.columns:
        metrics_df['delta_' + col + '_vs_WOA'] = metrics_df[col] - woa_row[col]

# Useful summary objects left for the answer.
basin_metrics = metrics_df[[
    'source', 'np_pacific', 'np_atlantic', 'np_indian',
    'pct_lt12_global', 'pct_lt12_pacific', 'pct_lt12_indian'
]].copy()

model_low_np_assessment = metrics_df[[
    'source',
    'pct_lt12_global', 'delta_pct_lt12_global_vs_WOA',
    'pct_lt12_pacific', 'delta_pct_lt12_pacific_vs_WOA',
    'pct_lt12_indian', 'delta_pct_lt12_indian_vs_WOA',
    'np_pacific', 'delta_np_pacific_vs_WOA',
    'np_indian', 'delta_np_indian_vs_WOA'
]].copy()
