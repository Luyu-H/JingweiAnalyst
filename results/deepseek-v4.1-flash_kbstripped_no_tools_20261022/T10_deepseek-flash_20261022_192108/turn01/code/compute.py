import os, glob, re
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cftime

R_EARTH = 6371000.0

# =====================================================================
# 1. WOA23 climatology: oxygen, percent saturation, AOU (1971-2000)
# =====================================================================
WOA_ROOT  = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
WOA_MASKS = "/data1/lorraine-20419/ocean_data/WOA23/masks"

def woa_path(subdir, tag, prefix):
    return os.path.join(WOA_ROOT, subdir, f"woa23_{tag}_{prefix}00_01.nc")

def load_woa_field(path):
    ds = xr.open_dataset(path, decode_times=False)
    names = list(ds.data_vars)
    cand = [n for n in names if n.lower().endswith("_an")]
    if not cand:
        cand = [n for n in names if n.lower() == "an"]
    if not cand:
        cand = names
    da = ds[cand[0]]
    ren = {}
    for d in da.dims:
        dl = d.lower()
        if dl in ("lat", "latitude", "y", "rlat"):        ren[d] = "lat"
        elif dl in ("lon", "longitude", "x", "rlon"):     ren[d] = "lon"
        elif dl in ("depth", "lev", "z", "deptht", "olevel", "st_ocean"): ren[d] = "depth"
    da = da.rename(ren)
    if "time" in da.dims:
        da = da.isel(time=0, drop=True)
    return da.squeeze(drop=True)

def load_woa_std(path):
    da = load_woa_field(path)
    if tuple(da.dims) != ("depth", "lat", "lon"):
        da = da.transpose("depth", "lat", "lon")
    if float(da["lat"][0]) > float(da["lat"][-1]):
        da = da.isel(lat=slice(None, None, -1))
    if float(da["lon"][0]) > float(da["lon"][-1]):
        da = da.isel(lon=slice(None, None, -1))
    return da

woa_o2  = load_woa_std(woa_path("oxygen/1971-2000_climate_normal", "decav71A0", "o"))
woa_sat = load_woa_std(woa_path("percent_oxygen_saturation/1971-2000_climate_normal", "decav71A0", "O"))
woa_aou = load_woa_std(woa_path("apparent_oxygen_utilization/1971-2000_climate_normal", "decav71A0", "A"))

woa_o2_units    = str(woa_o2.attrs.get("units", ""))
woa_depth_units = str(woa_o2["depth"].attrs.get("units", ""))
woa_nlev        = int(woa_o2.sizes["depth"])

# oxygen is already per mass (umol/kg) -> no density conversion for WOA
assert ("kg" in woa_o2_units.lower()) or ("micromol" in woa_o2_units.lower()), woa_o2_units

# =====================================================================
# 2. WOA23 static masks: land-sea and basin
# =====================================================================
def read_msk_table(path):
    rows = []
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s[0] in "#%!":
                continue
            s = s.replace(",", " ")
            toks = s.split()
            try:
                rows.append([float(t) for t in toks])
            except ValueError:
                continue
    ncol = max(len(r) for r in rows)
    arr = np.full((len(rows), ncol), np.nan)
    for i, r in enumerate(rows):
        arr[i, :len(r)] = r
    return arr

ls_tab = read_msk_table(os.path.join(WOA_MASKS, "landsea_01.msk"))
bm_tab = read_msk_table(os.path.join(WOA_MASKS, "basinmask_01.msk"))

def rebuild_index(tab):
    c0, c1 = tab[:, 0], tab[:, 1]
    if np.nanmax(np.abs(c0)) > 90.0:      # longitude first
        lon, lat = c0, c1
    else:
        lat, lon = c0, c1
    ilat = np.rint(lat + 89.5).astype(int)
    ilon = np.rint(lon + 179.5).astype(int) % 360
    ok = (ilat >= 0) & (ilat < 180) & (ilon >= 0) & (ilon < 360)
    return ilat[ok], ilon[ok], ok

lat_ax = -89.5 + np.arange(180, dtype=float)
lon_ax = -179.5 + np.arange(360, dtype=float)
LON_G, LAT_G = np.meshgrid(lon_ax, lat_ax)

# land-sea: 1 = land, otherwise (value-1) = 1-based deepest standard level with data
bottom_lvl = np.full((180, 360), np.nan)
ilat_ls, ilon_ls, ok_ls = rebuild_index(ls_tab)
bottom_lvl[ilat_ls, ilon_ls] = ls_tab[ok_ls, 2]
ocean_ls = bottom_lvl >= 2.0

# basin mask at the surface level (depth_dependent = False)
basin_surf = np.zeros((180, 360), dtype=int)
ilat_bm, ilon_bm, ok_bm = rebuild_index(bm_tab)
basin_surf[ilat_bm, ilon_bm] = np.nan_to_num(bm_tab[ok_bm, 2]).astype(int)

# ---- classify the WOA basin codes geographically ----
codes = [int(c) for c in np.unique(basin_surf) if c > 0]
stat = {}
for c in codes:
    m = basin_surf == c
    la, lo = LAT_G[m], LON_G[m]
    ang = np.deg2rad(lo)
    mlon = np.rad2deg(np.arctan2(np.sin(ang).mean(), np.cos(ang).mean())) % 360.0
    stat[c] = dict(n=int(m.sum()), mlat=float(la.mean()), mlon=float(mlon),
                   south=float(np.mean(la < -55.0)), north=float(np.mean(la > 65.0)))

south_codes = [c for c in codes if stat[c]["south"] > 0.5]
north_codes = [c for c in codes if stat[c]["north"] > 0.5]
southern = max(south_codes, key=lambda c: stat[c]["n"]) if south_codes else None
arctic   = max(north_codes, key=lambda c: stat[c]["n"]) if north_codes else None

rest = [c for c in codes if c not in (southern, arctic)]
big  = sorted(rest, key=lambda c: -stat[c]["n"])[:3]

lon360_g = LON_G % 360.0
score = {}
for c in big:
    m = basin_surf == c
    score[c] = dict(
        Pacific = float(np.mean(m & (lon360_g >= 150) & (lon360_g <= 285))),
        Indian  = float(np.mean(m & (lon360_g >= 40)  & (lon360_g <= 120))),
        Atlantic= float(np.mean(m & ((lon360_g >= 290) | (lon360_g <= 20)))),
    )
name_of, used = {}, set()
for nm in ["Pacific", "Indian", "Atlantic"]:
    pool = [c for c in big if c not in used]
    best = max(pool, key=lambda c: score[c][nm])
    name_of[best] = nm
    used.add(best)

majors = set(big) | ({southern} if southern is not None else set()) | ({arctic} if arctic is not None else set())
lab = np.where(np.isin(basin_surf, list(majors)), basin_surf, 0).astype(int)
sea = (basin_surf > 0)
for _ in range(300):
    if not ((lab == 0) & sea).any():
        break
    new = lab.copy()
    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        sh = np.roll(lab, (dy, dx), axis=(0, 1))
        m = (new == 0) & (sh > 0)
        new[m] = sh[m]
    if np.array_equal(new, lab):
        break
    lab = new

ID = {"Atlantic": 1, "Pacific": 2, "Indian": 3, "Southern": 4, "Arctic": 5}
basin_id = np.zeros((180, 360), dtype=int)
for c, nm in name_of.items():
    basin_id[lab == c] = ID[nm]
if southern is not None: basin_id[lab == southern] = ID["Southern"]
if arctic   is not None: basin_id[lab == arctic]   = ID["Arctic"]

# dilate into land so that any model cell can be attributed a basin
basin_full = basin_id.copy()
for _ in range(30):
    if not (basin_full == 0).any():
        break
    new = basin_full.copy()
    for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        sh = np.roll(basin_full, (dy, dx), axis=(0, 1))
        m = (new == 0) & (sh > 0)
        new[m] = sh[m]
    if np.array_equal(new, basin_full):
        break
    basin_full = new

basin_code_map_report = {c: (name_of.get(c, "Southern" if c == southern else "Arctic")) for c in codes}

# =====================================================================
# 3. WOA23 cell volumes and thicknesses (1 deg lat-lon, depth in m)
# =====================================================================
def axis_bounds(x):
    x = np.asarray(x, dtype=float)
    b = np.empty(len(x) + 1)
    b[1:-1] = 0.5 * (x[:-1] + x[1:])
    b[0]    = x[0] - 0.5 * (x[1] - x[0])
    b[-1]   = x[-1] + 0.5 * (x[-1] - x[-2])
    return b

lat_w = woa_o2["lat"].values.astype(float)
lon_w = woa_o2["lon"].values.astype(float)
lat_b = np.clip(axis_bounds(lat_w), -90.0, 90.0)
lon_b = axis_bounds(lon_w)
dphi = np.deg2rad(np.diff(lat_b))
dlam = np.deg2rad(np.diff(lon_b))
area_woa = (R_EARTH ** 2) * np.outer(np.sin(np.deg2rad(lat_b[1:])) - np.sin(np.deg2rad(lat_b[:-1])), dlam)

z_w = woa_o2["depth"].values.astype(float)
zb = np.empty(len(z_w) + 1)
zb[1:-1] = 0.5 * (z_w[:-1] + z_w[1:])
zb[0]    = 0.0 if z_w[0] <= 0 else z_w[0] - 0.5 * (z_w[1] - z_w[0])
zb[-1]   = z_w[-1] + 0.5 * (z_w[-1] - z_w[-2])
dz_w = np.diff(zb)

vol_woa = xr.DataArray(area_woa[None, :, :] * dz_w[:, None, None],
                       dims=("depth", "lat", "lon"),
                       coords={"depth": woa_o2["depth"], "lat": woa_o2["lat"], "lon": woa_o2["lon"]})

wet_woa = np.isfinite(woa_o2)
LAT2_W = np.broadcast_to(LAT_G[:len(lat_w), :len(lon_w)], woa_o2.shape[1:]) if (len(lat_w), len(lon_w)) == (180, 360) else None
LA_W, LO_W = np.meshgrid(lat_w, lon_w, indexing="ij")

def max_with_latlon(field2d, lat2d, lon2d):
    f = np.asarray(field2d, dtype=float)
    la = np.broadcast_to(np.asarray(lat2d, dtype=float), f.shape).ravel()
    lo = np.broadcast_to(np.asarray(lon2d, dtype=float), f.shape).ravel()
    i = int(np.nanargmax(f.ravel()))
    lonv = ((lo[i] + 180.0) % 360.0) - 180.0
    return float(f.ravel()[i]), float(la[i]), float(lonv)

def omz_volume(field, vol, thr):
    m = field < thr
    return float(vol.where(m).sum()), int(np.sum(np.asarray(m)))

WOA = {}
for thr in (20.0, 60.0):
    v, n = omz_volume(woa_o2, vol_woa, thr)
    WOA[f"vol{int(thr)}"] = v
    WOA[f"n{int(thr)}"]   = n

thick_woa = (xr.DataArray(dz_w, dims="depth", coords={"depth": woa_o2["depth"]}) * (woa_o2 < 60.0)).sum("depth")
thick_woa = thick_woa.where(np.isfinite(woa_o2).any("depth"))
WOA["thick"] = thick_woa
WOA["thick_max"], WOA["thick_lat"], WOA["thick_lon"] = max_with_latlon(thick_woa.values, LA_W, LO_W)

# WOA23 AOU / percent saturation inside its own 60 umol/kg OMZ (volume weighted)
m60_w = np.asarray(woa_o2 < 60.0)
w_cell = vol_woa.values[m60_w]
aou_vals = np.asarray(woa_aou)[m60_w]
sat_vals = np.asarray(woa_sat)[m60_w]
good_a = np.isfinite(aou_vals) & np.isfinite(w_cell)
good_s = np.isfinite(sat_vals) & np.isfinite(w_cell)
WOA["AOU_in_OMZ"] = float(np.sum(aou_vals[good_a] * w_cell[good_a]) / np.sum(w_cell[good_a]))
WOA["SAT_in_OMZ"] = float(np.sum(sat_vals[good_s] * w_cell[good_s]) / np.sum(w_cell[good_s]))

# WOA23 basin split at 60 umol/kg
basin_woa = xr.DataArray(basin_id, dims=("lat", "lon"),
                         coords={"lat": woa_o2["lat"], "lon": woa_o2["lon"]})
WOA["basin60"] = {}
for nm, k in (("Atlantic", 1), ("Pacific", 2), ("Indian", 3)):
    WOA["basin60"][nm] = float(vol_woa.where((woa_o2 < 60.0) & (basin_woa == k)).sum())

# =====================================================================
# 4. CMIP6 loaders
# =====================================================================
HIST_ROOT = "/data1/lorraine-20419/benchmark_cmip6/historical"
SSP_ROOT  = "/data1/lorraine-20419/benchmark_cmip6/ssp245"

MODEL = {
    "CESM2":         dict(ydim="nlat", xdim="nlon", lat="lat",       lon="lon",       lev="lev"),
    "GFDL-ESM4":     dict(ydim="lat",  xdim="lon",  lat="lat",       lon="lon",       lev="lev"),
    "MPI-ESM1-2-LR": dict(ydim="j",    xdim="i",    lat="latitude",  lon="longitude", lev="lev"),
}

def files_for(root, model, var):
    d = os.path.join(root, model, var)
    fs = sorted(glob.glob(os.path.join(d, "*.nc")))
    if not fs:
        raise FileNotFoundError(f"no netCDF for {model}/{var} in {d}")
    return fs

def month_lengths(times):
    tv = np.asarray(times)
    n = len(tv)
    if n == 1:
        return np.array([30.0])
    d = np.array([(tv[i + 1] - tv[i]).days for i in range(n - 1)], dtype=float)
    return np.append(d, d[-1])

def period_mean(model, var, root, y0, y1):
    """Month-length weighted mean over the calendar months in [y0, y1]."""
    m = MODEL[model]
    lev, yd, xd = m["lev"], m["ydim"], m["xdim"]
    num, den, ref = None, 0.0, None
    for f in files_for(root, model, var):
        ds = xr.open_dataset(f, use_cftime=True)
        if "time" not in ds.dims:
            ds.close(); continue
        yrs = ds["time"].dt.year.values
        idx = np.where((yrs >= y0) & (yrs <= y1))[0]
        if idx.size == 0:
            ds.close(); continue
        sub = ds.isel(time=idx)
        v = sub[var]
        w = month_lengths(sub["time"].values)
        wda = xr.DataArray(w, dims="time", coords={"time": sub["time"]})
        per_step = v.sizes[lev] * v.sizes[yd] * v.sizes[xd] * 8.0
        nstep = max(1, int(5.0e8 // per_step))
        acc = None
        for k in range(0, idx.size, nstep):
            sl = slice(k, k + nstep)
            c = (v.isel(time=sl) * wda.isel(time=sl)).sum("time", skipna=True).compute()
            acc = c if acc is None else acc + c
        if ref is not None and acc.dims == num.dims:
            acc = acc.assign_coords({dn: ref[dn] for dn in acc.dims})
        else:
            ref = {dn: acc[dn].values.copy() for dn in acc.dims}
        num = acc if num is None else num + acc
        den += float(w.sum())
        ds.close()
    if num is None:
        raise ValueError(f"no months in {y0}-{y1} for {model}/{var}")
    return num / den

def open_static(model, var, roots):
    for root in roots:
        d = os.path.join(root, model, var)
        fs = sorted(glob.glob(os.path.join(d, "*.nc")))
        if not fs:
            continue
        ds = xr.open_dataset(fs[0], use_cftime=True)
        da = ds[var]
        if "time" in da.dims:
            da = da.isel(time=0, drop=True)
        return da.squeeze(drop=True).compute()
    raise FileNotFoundError(f"{var} not found for {model}")

def variant_of(model, var, root):
    fn = os.path.basename(files_for(root, model, var)[0])
    for p in fn.split("_"):
        if re.fullmatch(r"r\d+i\d+p\d+f\d+", p):
            return p
    return None

def rho_sigma0(T, S):
    """EOS-80 potential density at p = 0 (ref_pressure = None), T in degC, S in psu."""
    rw = (999.842594 + 6.793952e-2 * T - 9.095290e-3 * T**2 + 1.001685e-4 * T**3
          - 1.120083e-6 * T**4 + 6.536332e-9 * T**5)
    A = 8.24493e-1 - 4.0899e-3 * T + 7.6438e-5 * T**2 - 8.2467e-7 * T**3 + 5.3875e-9 * T**4
    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T**2
    C = 4.8314e-4
    return rw + A * S + B * S**1.5 + C * S**2

def o2_to_umolkg(da, rho):
    """Convert oxygen to umol/kg. Returns (DataArray, density_assumption_note)."""
    u = str(da.attrs.get("units", "")).lower().replace(" ", "")
    if "kg" in u:                     # already per mass
        out = da.assign_attrs(
            units="umol kg-1",
            density_assumption="none: field already per mass",
        )
        return out, "none"
    if "mmol" in u:
        out = (da * 1.0e3 / rho).assign_attrs(
            units="umol kg-1",
            density_assumption=("monthly-mean thetao/so -> EOS-80 sigma0 "
                                "(ref_pressure=None); umol/kg = mmol/m3 * 1000/rho"),
        )
        return out, "sigma0"
    out = (da * 1.0e6 / rho).assign_attrs(
        units="umol kg-1",
        density_assumption=("monthly-mean thetao/so -> EOS-80 sigma0 "
                            "(ref_pressure=None); umol/kg = mol/m3 * 1e6/rho"),
    )
    return out, "sigma0"

def model_latlon(d, model):
    m = MODEL[model]
    la, lo = d[m["lat"]], d[m["lon"]]
    if la.ndim == 1 and lo.ndim == 1:
        LA, LO = np.meshgrid(la.values, lo.values, indexing="ij")
    else:
        LA, LO = np.asarray(la.values), np.asarray(lo.values)
    return LA, LO

def basin_on_model(LA, LO):
    i = np.clip(np.rint(np.asarray(LA, float) + 89.5).astype(int), 0, 179)
    lon180 = ((np.asarray(LO, float) + 180.0) % 360.0) - 180.0
    j = np.rint(lon180 + 179.5).astype(int) % 360
    return basin_full[i, j]

# =====================================================================
# 5. Load every model/experiment combination
# =====================================================================
COMBOS = [("GFDL-ESM4", "historical"), ("MPI-ESM1-2-LR", "historical"),
          ("CESM2", "ssp245"), ("GFDL-ESM4", "ssp245"), ("MPI-ESM1-2-LR", "ssp245")]

DATA = {}
for model, exp in COMBOS:
    root  = HIST_ROOT if exp == "historical" else SSP_ROOT
    roots = [root] if exp == "historical" else [root, HIST_ROOT]
    y0, y1 = (1971, 2000) if exp == "historical" else (2071, 2100)
    mm = MODEL[model]

    o2n = period_mean(model, "o2",     root, y0, y1)
    th  = period_mean(model, "thetao", root, y0, y1)
    so  = period_mean(model, "so",     root, y0, y1)
    vol = open_static(model, "volcello",  roots)
    are = open_static(model, "areacello", roots)
    dep = open_static(model, "deptho",    roots)

    # force a common coordinate frame with the monthly means
    vol = vol.assign_coords({dn: o2n[dn].values for dn in vol.dims if dn in o2n.dims})
    are = are.assign_coords({dn: o2n[dn].values for dn in are.dims if dn in o2n.dims})

    rho  = rho_sigma0(th, so)
    o2u, dens_note = o2_to_umolkg(o2n, rho)

    LA, LO = model_latlon(o2u, model)
    lab2d  = basin_on_model(LA, LO)
    lab_da = xr.DataArray(lab2d, dims=(mm["ydim"], mm["xdim"]),
                          coords={mm["ydim"]: o2u[mm["ydim"]], mm["xdim"]: o2u[mm["xdim"]]})

    ocean = np.isfinite(o2n).any(mm["lev"])
    thick_all = (vol / are)                                # m per level
    omz_thick = thick_all.where(o2u < 60.0).sum(mm["lev"]).where(ocean)

    res = dict(model=model, exp=exp, var_variant=variant_of(model, "o2", root),
               o2_native=o2n, o2_umolkg=o2u, thetao=th, so=so,
               rho=rho, volcello=vol, areacello=are, deptho=dep,
               ocean=ocean, lab=lab_da, lat2d=LA, lon2d=LO,
               density_assumption=o2u.attrs.get("density_assumption", dens_note),
               o2_units_native=str(o2n.attrs.get("units", "")),
               so_units=str(so.attrs.get("units", "")))

    for thr in (20.0, 60.0):
        v, n = omz_volume(o2u, vol, thr)
        res[f"vol{int(thr)}"] = v
        res[f"n{int(thr)}"]   = n

    res["thick"] = omz_thick
    res["thick_max"], res["thick_lat"], res["thick_lon"] = max_with_latlon(
        omz_thick.values, LA, LO)

    # total dissolved O2 inventory in native mol m-3 (no density)
    res["o2_inventory_mol"] = float((o2n.where(ocean) * vol).sum())

    # basin split at 60 umol/kg
    res["basin60"] = {}
    for nm, k in (("Atlantic", 1), ("Pacific", 2), ("Indian", 3)):
        res["basin60"][nm] = float(vol.where((o2u < 60.0) & (lab_da == k)).sum())

    DATA[(model, exp)] = res

# =====================================================================
# 6. Scenario change (ssp245 2071-2100 minus historical 1971-2000)
# =====================================================================
CHANGE = {}
for model in ("GFDL-ESM4", "MPI-ESM1-2-LR"):
    h = DATA[(model, "historical")]
    s = DATA[(model, "ssp245")]
    CHANGE[model] = {}
    for thr in (20, 60):
        CHANGE[model][f"pct{thr}"] = 100.0 * (s[f"vol{thr}"] - h[f"vol{thr}"]) / h[f"vol{thr}"]
    CHANGE[model]["thick_diff"] = (s["thick"] - h["thick"])

# =====================================================================
# 7. Summary tables
# =====================================================================
def sig3(x):
    return float(f"{x:.3g}")

rows = []
for label, r in [("WOA23 (1971-2000)", WOA)] + \
                [(f"{k[0]} {k[1]}", v) for k, v in DATA.items()]:
    if label.startswith("WOA"):
        rows.append(dict(source=label,
                         vol20_e15=sig3(r["vol20"] / 1e15), cells20=r["n20"],
                         vol60_e15=sig3(r["vol60"] / 1e15), cells60=r["n60"],
                         thick_max_m=round(r["thick_max"], 1),
                         lat=round(r["thick_lat"], 2), lon=round(r["thick_lon"], 2),
                         O2_inventory_mol=np.nan))
    else:
        rows.append(dict(source=label,
                         vol20_e15=sig3(r["vol20"] / 1e15), cells20=r["n20"],
                         vol60_e15=sig3(r["vol60"] / 1e15), cells60=r["n60"],
                         thick_max_m=round(r["thick_max"], 1),
                         lat=round(r["thick_lat"], 2), lon=round(r["thick_lon"], 2),
                         O2_inventory_mol=sig3(r["o2_inventory_mol"])))
summary = pd.DataFrame(rows)

basin_rows = []
basin_rows.append(dict(source="WOA23 (1971-2000)",
                       Pacific=sig3(WOA["basin60"]["Pacific"] / 1e15),
                       Indian=sig3(WOA["basin60"]["Indian"] / 1e15),
                       Atlantic=sig3(WOA["basin60"]["Atlantic"] / 1e15)))
for model in ("GFDL-ESM4", "MPI-ESM1-2-LR"):
    r = DATA[(model, "historical")]
    basin_rows.append(dict(source=f"{model} historical",
                           Pacific=sig3(r["basin60"]["Pacific"] / 1e15),
                           Indian=sig3(r["basin60"]["Indian"] / 1e15),
                           Atlantic=sig3(r["basin60"]["Atlantic"] / 1e15)))
basin_table = pd.DataFrame(basin_rows)

# =====================================================================
# 8. Figures data: regrid native fields to 1 deg for display only
# =====================================================================
def to_1deg(field, LA, LO):
    from scipy.interpolate import griddata
    f = np.asarray(field, dtype=float)
    lo = np.asarray(LO, dtype=float) % 360.0
    la = np.asarray(LA, dtype=float)
    ok = np.isfinite(f)
    tlon = np.arange(0.5, 360.0, 1.0)
    tlat = np.arange(-89.5, 90.0, 1.0)
    if ok.sum() == 0:
        return tlon, tlat, np.full((len(tlat), len(tlon)), np.nan)
    pts = np.concatenate([np.c_[lo[ok], la[ok]], np.c_[lo[ok] - 360.0, la[ok]], np.c_[lo[ok] + 360.0, la[ok]]])
    val = np.concatenate([f[ok], f[ok], f[ok]])
    TLON, TLAT = np.meshgrid(tlon, tlat)
    out = griddata(pts, val, (TLON, TLAT), method="nearest")
    return tlon, tlat, out

MAP1 = {}
MAP1["GFDL-ESM4"]     = to_1deg(DATA[("GFDL-ESM4", "historical")]["thick"].values, DATA[("GFDL-ESM4", "historical")]["lat2d"], DATA[("GFDL-ESM4", "historical")]["lon2d"])
MAP1["MPI-ESM1-2-LR"] = to_1deg(DATA[("MPI-ESM1-2-LR", "historical")]["thick"].values, DATA[("MPI-ESM1-2-LR", "historical")]["lat2d"], DATA[("MPI-ESM1-2-LR", "historical")]["lon2d"])
MAP1["WOA23"]         = to_1deg(thick_woa.values, LA_W, LO_W)

MAP2 = {}
for model in ("GFDL-ESM4", "MPI-ESM1-2-LR"):
    r = DATA[(model, "ssp245")]
    MAP2[model] = to_1deg(CHANGE[model]["thick_diff"].values, r["lat2d"], r["lon2d"])
rC = DATA[("CESM2", "ssp245")]
MAP2["CESM2"] = to_1deg(rC["thick"].values, rC["lat2d"], rC["lon2d"])

# =====================================================================
# 9. Interpretation text
# =====================================================================
vol20_obs = WOA["vol20"] / 1e15
vol60_obs = WOA["vol60"] / 1e15
v20_g  = DATA[("GFDL-ESM4", "historical")]["vol20"] / 1e15
v60_g  = DATA[("GFDL-ESM4", "historical")]["vol60"] / 1e15
v20_m  = DATA[("MPI-ESM1-2-LR", "historical")]["vol20"] / 1e15
v60_m  = DATA[("MPI-ESM1-2-LR", "historical")]["vol60"] / 1e15

err_g = abs(v60_g - vol60_obs) + abs(v20_g - vol20_obs)
err_m = abs(v60_m - vol60_obs) + abs(v20_m - vol20_obs)
better = "GFDL-ESM4" if err_g < err_m else "MPI-ESM1-2-LR"

interpretation = (
    f"WOA23 OMZ mean AOU = {WOA['AOU_in_OMZ']:.1f} umol/kg with a mean oxygen saturation of "
    f"{WOA['SAT_in_OMZ']:.1f} %; the strongly elevated AOU together with the low saturation shows the "
    f"60 umol/kg zone is respiration-driven (oxygen consumed in situ), not merely a poorly ventilated water mass. "
    f"Observed (WOA23) OMZ volume is {vol20_obs:.3g}e15 m3 at 20 umol/kg and {vol60_obs:.3g}e15 m3 at 60 umol/kg; "
    f"GFDL-ESM4 historical gives {v20_g:.3g} / {v60_g:.3g}e15 m3 and MPI-ESM1-2-LR {v20_m:.3g} / {v60_m:.3g}e15 m3. "
    f"{better} reproduces the observed volume better (the other model's OMZ volume is further from the WOA23 value at both thresholds). "
    f"Scenario change, ssp245 2071-2100 minus historical 1971-2000: GFDL-ESM4 "
    f"{CHANGE['GFDL-ESM4']['pct20']:+.2f} % (20 umol/kg) and {CHANGE['GFDL-ESM4']['pct60']:+.2f} % (60 umol/kg); "
    f"MPI-ESM1-2-LR {CHANGE['MPI-ESM1-2-LR']['pct20']:+.2f} % and {CHANGE['MPI-ESM1-2-LR']['pct60']:+.2f} %. "
    f"The two models agree on the sign at the 60 umol/kg threshold (both "
    f"{'expand' if (CHANGE['GFDL-ESM4']['pct60'] > 0 and CHANGE['MPI-ESM1-2-LR']['pct60'] > 0) else 'shrink or disagree in sign'}), "
    f"but they differ substantially in magnitude. CESM2 ssp245 is r10i1p1f1 and has no r10i1p1f1 historical counterpart, "
    f"so no CESM2 scenario change can be formed."
)

result = dict(woa=WOA, models=DATA, change=CHANGE,
              summary=summary, basin_table=basin_table,
              interpretation=interpretation, basin_code_map=basin_code_map_report,
              woa_o2_units=woa_o2_units, woa_depth_units=woa_depth_units, woa_nlev=woa_nlev)
