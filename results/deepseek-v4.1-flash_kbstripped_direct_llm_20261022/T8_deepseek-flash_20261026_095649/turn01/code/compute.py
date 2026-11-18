import os, glob
from collections import Counter
import numpy as np
import pandas as pd
import xarray as xr

CMIP_ROOT   = "/data1/lorraine-20419/benchmark_cmip6/historical"
GLODAP_ROOT = "/data1/lorraine-20419/ocean_data/GLODAPv2.2016b_MappedClimatologies"

MODELS = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
INFO = {
    "CESM2":         dict(y="nlat", x="nlon", lat="lat",       lon="lon",       lev="lev", scale=0.01),
    "GFDL-ESM4":     dict(y="lat",  x="lon",  lat="lat",       lon="lon",       lev="lev", scale=1.0),
    "MPI-ESM1-2-LR": dict(y="j",    x="i",    lat="latitude",  lon="longitude", lev="lev", scale=1.0),
}

# ----------------------------------------------------------------------------- helpers
_DC = {}

def _files(model, var):
    fs = sorted(glob.glob(os.path.join(CMIP_ROOT, model, var, "*.nc")))
    if not fs:
        return []
    def ens(f):
        p = os.path.basename(f)[:-3].split("_")
        return p[4] if len(p) > 4 else ""
    keep = Counter(ens(f) for f in fs).most_common(1)[0][0]
    return [f for f in fs if ens(f) == keep]

def open_var(model, var):
    """Open every time-split file of one model/variable, keeping one ensemble member."""
    k = (model, var)
    if k not in _DC:
        fs = _files(model, var)
        if not fs:
            _DC[k] = None
        else:
            try:
                _DC[k] = xr.open_mfdataset(fs, combine="by_coords", use_cftime=True, chunks={})
            except Exception:
                _DC[k] = xr.open_mfdataset(fs, combine="by_coords",
                                           decode_times=True, use_cftime=True)
    return _DC[k]

def _bounds(z):
    z = np.asarray(z, float)
    b = np.empty(z.size + 1)
    b[0] = 0.0
    b[1:-1] = 0.5 * (z[:-1] + z[1:])
    b[-1] = z[-1] + 0.5 * (z[-1] - z[-2])
    return b

def band_thickness(z, z0=0.0, z1=50.0):
    """Thickness of every level that falls inside the z0-z1 band (0 outside)."""
    b = _bounds(z)
    return np.clip(np.minimum(b[1:], z1) - np.maximum(b[:-1], z0), 0.0, None)

def rho_sw(T, S, z):
    """Simple in-situ density (kg m-3) from theta, practical salinity and depth (m)."""
    sig = 28.106 - 0.0735 * T - 0.00469 * T * T + (0.802 - 0.002 * T) * (S - 35.0)
    r0 = 1000.0 + sig
    return r0 * (1.0 + np.abs(z) * 1.0e4 / 2.2e9)      # ~compressibility term

def model_lev(model):
    ds = open_var(model, "thetao")
    return np.asarray(ds[INFO[model]["lev"]].values, float) * INFO[model]["scale"]

def model_area(model):
    ds = open_var(model, "areacello")
    return np.asarray(ds["areacello"].values, float)

def model_coords2d(model):
    ds = open_var(model, "thetao")
    la, lo = ds[INFO[model]["lat"]], ds[INFO[model]["lon"]]
    if la.ndim == 1:
        LAT, LON = np.meshgrid(np.asarray(la.values, float),
                               np.asarray(lo.values, float), indexing="ij")
    else:
        LAT = np.asarray(la.values, float)
        LON = np.asarray(lo.values, float)
    return LAT, LON

def model_slab(model, period, varlist=("dissic", "talk", "ph", "thetao", "so")):
    """Time-mean of the requested variables over the 0-50 m levels."""
    z = model_lev(model)
    dz = np.diff(_bounds(z))
    dzb = band_thickness(z)
    idx = np.where(dzb > 0)[0]
    out = {}
    for v in varlist:
        ds = open_var(model, v)
        if ds is None:
            out[v] = None
            continue
        da = ds[v]
        if period == "2002":
            da = da.sel(time=slice("2002-01-01", "2002-12-31"))
        else:
            da = da.sel(time=slice("1971-01-01", "2000-12-31"))
        da = da.isel(**{INFO[model]["lev"]: idx}).mean("time")
        out[v] = np.asarray(da.compute().values, float)
    return out, z[idx], dzb[idx], np.clip(dzb[idx] / dz[idx], 0.0, 1.0)

def model_band_stats(model, period="2002"):
    varlist = ("dissic", "talk", "ph", "thetao", "so")
    slab, z, dzb, frac = model_slab(model, period, varlist)
    A = model_area(model)
    Z = z[:, None, None]
    rho = rho_sw(slab["thetao"], slab["so"], Z)
    dic = slab["dissic"] * 1e6 / rho            # mol m-3 -> umol kg-1
    alk = slab["talk"]    * 1e6 / rho
    ph  = slab["ph"]
    W = dzb[:, None, None] * A[None, :, :]       # volume inside the 0-50 m band

    res = {}
    for nm, arr in (("DIC", dic), ("TAlk", alk), ("pH", ph)):
        m = np.isfinite(arr) & np.isfinite(W) & (W > 0)
        res[nm] = float(np.sum(np.where(m, arr * W, 0.0)) / np.sum(np.where(m, W, 0.0)))

    # per-cell 0-50 m mean of DIC
    Wv = np.broadcast_to(dzb[:, None, None], dic.shape).astype(float)
    ok = np.isfinite(dic)
    num = np.sum(np.where(ok, dic * Wv, 0.0), axis=0)
    den = np.sum(np.where(ok, Wv, 0.0), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        res["map_dic"] = np.where(den > 0, num / den, np.nan)
    res["area"] = A
    return res

def model_at(model, v, zt, period="2002"):
    """Linear interpolation of a 3-D CMIP6 field to depth zt (m), annual mean."""
    ds = open_var(model, v)
    z = np.asarray(ds[INFO[model]["lev"]].values, float) * INFO[model]["scale"]
    lev = INFO[model]["lev"]
    da = ds[v]
    if period == "2002":
        da = da.sel(time=slice("2002-01-01", "2002-12-31"))
    else:
        da = da.sel(time=slice("1971-01-01", "2000-12-31"))
    da = da.mean("time")
    i = int(np.searchsorted(z, zt))
    if i <= 0:
        return np.asarray(da.isel(**{lev: 0}).compute().values, float), float(z[0])
    if i >= z.size:
        i = z.size - 1
    z0, z1 = z[i - 1], z[i]
    w = 0.0 if z1 == z0 else float(np.clip((zt - z0) / (z1 - z0), 0.0, 1.0))
    a0 = np.asarray(da.isel(**{lev: i - 1}).compute().values, float)
    a1 = np.asarray(da.isel(**{lev: i}).compute().values, float)
    return (1.0 - w) * a0 + w * a1, float(z0 + w * (z1 - z0))

def to_glodap(vals, lat2d, lon2d, area, glat, glon):
    """Area-weighted nearest-cell binning of a native-grid field onto the 1-deg GLODAP grid."""
    v  = np.asarray(vals, float).ravel()
    la = np.asarray(lat2d, float).ravel()
    lo = np.mod(np.asarray(lon2d, float).ravel(), 360.0)
    ar = np.asarray(area, float).ravel()
    ok = np.isfinite(v) & np.isfinite(ar) & (ar > 0) & np.isfinite(la) & np.isfinite(lo)
    dlat, dlon = glat[1] - glat[0], glon[1] - glon[0]
    ilat = np.rint((la[ok] - glat[0]) / dlat).astype(int)
    ilon = np.rint(np.mod(lo[ok] - glon[0], 360.0) / dlon).astype(int)
    good = (ilat >= 0) & (ilat < glat.size) & (ilon >= 0) & (ilon < glon.size)
    flat = (ilat * glon.size + ilon)[good]
    w = ar[ok][good]; x = v[ok][good]
    n = glat.size * glon.size
    num = np.bincount(flat, weights=w * x, minlength=n)
    den = np.bincount(flat, weights=w,       minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return out.reshape(glat.size, glon.size)

def wstats(x, y, w):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    x, y, w = x[m], y[m], w[m]
    W = w.sum()
    mx = (w * x).sum() / W
    my = (w * y).sum() / W
    dx, dy = x - mx, y - my
    cov = (w * dx * dy).sum() / W
    sx = np.sqrt((w * dx * dx).sum() / W)
    sy = np.sqrt((w * dy * dy).sum() / W)
    return ((w * (x - y)).sum() / W,
            np.sqrt((w * (x - y) ** 2).sum() / W),
            cov / (sx * sy))

# ----------------------------------------------------------------------------- GLODAP
_g0 = xr.open_dataset(os.path.join(GLODAP_ROOT, "GLODAPv2.2016b.TCO2.nc"))
GLAT = np.asarray(_g0["lat"].values, float) if "lat" in _g0 else np.asarray(_g0["latitude"].values, float)
GLON = np.asarray(_g0["lon"].values, float) if "lon" in _g0 else np.asarray(_g0["longitude"].values, float)
GZ   = np.asarray(_g0["Depth"].values, float)
GDZ  = band_thickness(GZ)
GIDX = np.where(GDZ > 0)[0]
GLATW = np.cos(np.deg2rad(GLAT))[:, None] * np.ones((1, GLON.size))

def glodap_sel(name, zt=None):
    ds = xr.open_dataset(os.path.join(GLODAP_ROOT, "GLODAPv2.2016b.%s.nc" % name))
    if name in ds.data_vars:
        v = ds[name]
    else:
        v = ds[[k for k in ds.data_vars if ds[k].ndim == 3][0]]
    a = np.asarray(v.values, float)
    if zt is None:
        return a
    i = int(np.argmin(np.abs(GZ - zt)))
    return a[i]

# GLODAP 0-50 m (volume-weighted global means and per-cell maps)
GLODAP_MEAN, GLODAP_MAP = {}, {}
for nm, var in (("DIC", "TCO2"), ("TAlk", "TAlk"), ("pH", "pHtsinsitutp")):
    a = glodap_sel(var)
    sub = a[GIDX]
    Wg = GDZ[GIDX][:, None, None] * GLATW[None, :, :]
    m = np.isfinite(sub) & (Wg > 0)
    GLODAP_MEAN[nm] = float(np.sum(np.where(m, sub * Wg, 0.0)) / np.sum(np.where(m, Wg, 0.0)))
    Wv = np.broadcast_to(GDZ[GIDX][:, None, None], sub.shape).astype(float)
    ok = np.isfinite(sub)
    num = np.sum(np.where(ok, sub * Wv, 0.0), axis=0)
    den = np.sum(np.where(ok, Wv, 0.0), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        GLODAP_MAP[nm] = np.where(den > 0, num / den, np.nan)

# ----------------------------------------------------------------------------- CMIP6 2002
RES_2002, MAPS, AREA, LAT2D, LON2D = {}, {}, {}, {}, {}
for model in MODELS:
    r = model_band_stats(model, "2002")
    RES_2002[model] = r
    MAPS[model] = r["map_dic"]
    AREA[model] = r["area"]
    LAT2D[model], LON2D[model] = model_coords2d(model)

BIAS_2002 = {m: {k: (RES_2002[m][k] - GLODAP_MEAN[k]) for k in ("DIC", "TAlk", "pH")} for m in MODELS}

# ----------------------------------------------------------------------------- map comparison (0-50 m DIC)
MAP_STATS, BIAS_MAPS = {}, {}
for model in MODELS:
    mg = to_glodap(MAPS[model], LAT2D[model], LON2D[model], AREA[model], GLAT, GLON)
    BIAS_MAPS[model] = mg - GLODAP_MAP["DIC"]
    b, rmse, corr = wstats(mg, GLODAP_MAP["DIC"], GLATW)
    MAP_STATS[model] = dict(bias=b, rmse=rmse, corr=corr)

# ----------------------------------------------------------------------------- vertical bias profile
DEPTHS = [0.0, 100.0, 500.0, 1000.0, 2000.0]
PROFILE = {m: {} for m in MODELS}
for model in MODELS:
    for zt in DEPTHS:
        d, zd = model_at(model, "dissic", zt)
        t, _  = model_at(model, "thetao", zt)
        s, _  = model_at(model, "so", zt)
        rho = rho_sw(t, s, max(zd, 0.0))
        dic = d * 1e6 / rho
        mg = to_glodap(dic, LAT2D[model], LON2D[model], AREA[model], GLAT, GLON)
        g  = glodap_sel("TCO2", zt)
        both = np.isfinite(mg) & np.isfinite(g)
        w = np.where(both, GLATW, 0.0)
        PROFILE[model][zt] = float(np.sum(np.where(both, (mg - g) * w, 0.0)) / np.sum(w))

# ----------------------------------------------------------------------------- 1971-2000 and period change
DIC_7100 = {}
for model in MODELS:
    slab, z, dzb, frac = model_slab(model, "7100", ("dissic", "thetao", "so"))
    A = model_area(model)
    Z = z[:, None, None]
    rho = rho_sw(slab["thetao"], slab["so"], Z)
    dic = slab["dissic"] * 1e6 / rho
    W = dzb[:, None, None] * A[None, :, :]
    m = np.isfinite(dic) & np.isfinite(W) & (W > 0)
    DIC_7100[model] = float(np.sum(np.where(m, dic * W, 0.0)) / np.sum(np.where(m, W, 0.0)))

PERIOD_CHANGE = {m: RES_2002[m]["DIC"] - DIC_7100[m] for m in MODELS}
PERIOD_RATIO  = {m: PERIOD_CHANGE[m] / BIAS_2002[m]["DIC"] for m in MODELS}

# ----------------------------------------------------------------------------- summary table
SUMMARY = pd.DataFrame(
    [dict(model=m,
          DIC=RES_2002[m]["DIC"], DIC_bias=BIAS_2002[m]["DIC"],
          TAlk=RES_2002[m]["TAlk"], TAlk_bias=BIAS_2002[m]["TAlk"],
          pH=RES_2002[m]["pH"], pH_bias=BIAS_2002[m]["pH"],
          map_bias=MAP_STATS[m]["bias"], map_rmse=MAP_STATS[m]["rmse"],
          map_corr=MAP_STATS[m]["corr"],
          DIC_1971_2000=DIC_7100[m], DIC_change=PERIOD_CHANGE[m],
          change_over_bias=PERIOD_RATIO[m])
     for m in MODELS],
    index=MODELS)
SUMMARY.loc["GLODAP"] = dict(DIC=GLODAP_MEAN["DIC"], TAlk=GLODAP_MEAN["TAlk"], pH=GLODAP_MEAN["pH"])

print("0-50 m global means (umol/kg, pH dimensionless)")
print(SUMMARY[["DIC", "DIC_bias", "TAlk", "TAlk_bias", "pH", "pH_bias"]].round(3).to_string())
print()
print("0-50 m DIC map comparison on the GLODAP 1-deg grid")
for m in MODELS:
    print("  %-14s bias %8.2f  rmse %8.2f  corr %.4f"
          % (m, MAP_STATS[m]["bias"], MAP_STATS[m]["rmse"], MAP_STATS[m]["corr"]))
print()
print("Global mean DIC bias (model - GLODAP) vs depth")
for m in MODELS:
    print("  %-14s " % m + "  ".join("%6.0fm:%8.2f" % (z, PROFILE[m][z]) for z in DEPTHS))
print()
print("DIC 0-50 m drift 2002 minus 1971-2000, and its ratio to the GLODAP offset")
for m in MODELS:
    print("  %-14s change %7.2f  bias %8.2f  ratio %6.2f"
          % (m, PERIOD_CHANGE[m], BIAS_2002[m]["DIC"], PERIOD_RATIO[m]))
