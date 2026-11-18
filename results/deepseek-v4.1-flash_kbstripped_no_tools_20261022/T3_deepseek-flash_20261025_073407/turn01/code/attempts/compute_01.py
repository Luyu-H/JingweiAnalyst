import os
import glob

# ======================================================================
# 1) WOA23 -- temperature, 1971-2000 normal, annual, objectively
#    analysed field ("an"), 0 m, 1 degree grid  (decode_times=False)
# ======================================================================
woa_root = "/data1/lorraine-20419/ocean_data/WOA23/climatology"
woa_file = os.path.join(woa_root, "temperature", "1971-2000_climate_normal",
                        "woa23_decav71A0_t00_01.nc")
ds_woa = xr.open_dataset(woa_file, decode_times=False)

an_name = sorted([v for v in ds_woa.data_vars if v.endswith("_an")])[0]   # 't_an'
woa_an = ds_woa[an_name]
if "time" in woa_an.dims:                       # annual slice -> index 0
    woa_an = woa_an.isel(time=0)
woa_depth0 = float(np.ravel(woa_an["depth"].values)[0])   # must be 0 m
woa_an = woa_an.isel(depth=0).squeeze(drop=True)
woa_an = woa_an.where(np.abs(woa_an) < 1e10)              # guard against undecoded fills
woa_sst = woa_an.rename("woa_sst").assign_attrs(
    units="degrees_celsius",
    long_name="WOA23 1971-2000 annual objectively analysed temperature at 0 m")

woa_lat = np.asarray(woa_sst["lat"].values, dtype=float)
woa_lon = np.asarray(woa_sst["lon"].values, dtype=float)
if woa_lat[0] > woa_lat[-1]:
    woa_sst = woa_sst.isel(lat=slice(None, None, -1)); woa_lat = woa_lat[::-1]
if woa_lon[0] > woa_lon[-1]:
    woa_sst = woa_sst.isel(lon=slice(None, None, -1)); woa_lon = woa_lon[::-1]
NLAT, NLON = woa_sst.sizes["lat"], woa_sst.sizes["lon"]
assert (NLAT, NLON) == (180, 360)

# ---- WOA23 land-sea mask (text table, 2 header rows, 64800 ocean+land pts) --
mask_file = "/data1/lorraine-20419/ocean_data/WOA23/masks/landsea_01.msk"
_rows = []
with open(mask_file, "r") as fh:
    for line in fh:
        p = line.replace(",", " ").split()
        if len(p) < 3:
            continue
        try:
            _rows.append((float(p[0]), float(p[1]), float(p[2])))
        except ValueError:
            continue
tab = np.array(_rows, dtype=float)
if tab.shape[0] > 180 * 360:            # drop any numeric-looking header line
    tab = tab[-180 * 360:]
assert tab.shape[0] == 180 * 360

lat_v = tab[:, 0]
lon_v = ((tab[:, 1] + 180.0) % 360.0) - 180.0
m_v = tab[:, 2]
li = np.clip(np.round((lat_v + 89.5) / 1.0).astype(int), 0, 179)
oi = np.clip(np.round((lon_v + 179.5) / 1.0).astype(int), 0, 359)
_msk = np.full((180, 360), np.nan)
_msk[li, oi] = m_v
landsea = xr.DataArray(_msk, dims=("lat", "lon"),
                       coords={"lat": woa_lat, "lon": woa_lon}, name="landsea")

woa_ocean = (landsea > 1)          # WOA23: ocean when mask value > 1

# ---- requested whole-1-degree-grid, UNWEIGHTED mean of non-missing cells ----
woa_global_mean_unweighted = float(woa_sst.mean(skipna=True))
woa_valid_count = int(woa_sst.notnull().sum())
woa_ocean_count = int(woa_ocean.sum())
print(f"WOA23 0 m depth level = {woa_depth0} m")
print(f"WOA23 whole-1deg-grid unweighted mean SST = {woa_global_mean_unweighted:.2f} degC "
      f"({woa_valid_count} valid cells; landsea ocean cells = {woa_ocean_count})")

# ======================================================================
# 2) CMIP6 historical `tos` 1971-2000 climatology on each native grid
# ======================================================================
cmip_root = "/data1/lorraine-20419/benchmark_cmip6/historical"
models = ["CESM2", "GFDL-ESM4", "MPI-ESM1-2-LR"]
coords_of = {"CESM2": ("lat", "lon"),
             "GFDL-ESM4": ("lat", "lon"),
             "MPI-ESM1-2-LR": ("latitude", "longitude")}


def cmip_open(model, var):
    """Open all historical files of one model/variable; variant is parsed
    from the file name, never assumed."""
    fs = sorted(glob.glob(os.path.join(cmip_root, model, var, "*.nc")))
    if not fs:
        raise FileNotFoundError(f"no files under {cmip_root}/{model}/{var}")
    latn, lonn = coords_of[model]
    parts = []
    for f in fs:
        d = xr.open_dataset(f, use_cftime=True, decode_times=True)
        keep = [var] + [c for c in (latn, lonn)
                        if c in d.variables and c not in d.dims]
        keep = list(dict.fromkeys([k for k in keep if k in d.variables]))
        parts.append(d[keep])
    if len(parts) == 1:
        ds = parts[0]
    else:
        ds = xr.concat(parts, dim="time", data_vars="minimal",
                       coords="minimal", compat="override")
    if "time" in ds.dims:
        ds = ds.sortby("time")
        ti = ds.indexes["time"]
        if ti.has_duplicates:
            ds = ds.isel(time=~ti.duplicated())
    return ds, fs


def drop_singletons(da, keep):
    drop = [d for d in da.dims if d not in keep and da.sizes[d] == 1]
    if drop:
        da = da.isel({d: 0 for d in drop})
    return da


def conservative_regrid(vals2d, lon2d, lat2d, area2d, lon_t, lat_t):
    """First-order (donor-cell) conservative remapping: every native cell
    donates its full area*value to the 1 degree target cell containing its
    centre, so the area-weighted integral is preserved."""
    nlat_t, nlon_t = len(lat_t), len(lon_t)
    v = np.asarray(vals2d, float)
    a = np.asarray(area2d, float)
    la = np.asarray(lat2d, float)
    lo = np.asarray(lon2d, float)
    good = np.isfinite(v) & np.isfinite(a) & (a > 0) & np.isfinite(la) & np.isfinite(lo)

    li_t = np.zeros(la.shape, dtype=int)
    oi_t = np.zeros(lo.shape, dtype=int)
    li_t[good] = np.clip(np.floor(la[good] + 90.0).astype(int), 0, nlat_t - 1)
    oi_t[good] = np.mod(np.floor(np.mod(lo[good], 360.0) + 180.0).astype(int), nlon_t)

    w = np.where(good, a, 0.0)
    f = np.where(good, v, 0.0)
    flat = (li_t * nlon_t + oi_t).ravel()
    num = np.bincount(flat, weights=(f * w).ravel(), minlength=nlat_t * nlon_t)
    den = np.bincount(flat, weights=w.ravel(), minlength=nlat_t * nlon_t)
    out = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    return out.reshape(nlat_t, nlon_t)


native_mean = {}
native_area = {}
model_area_total = {}
n_months = {}
variants = {}
tos_clim_native = {}
regridded = {}

for m in models:
    ds, fs = cmip_open(m, "tos")
    variants[m] = os.path.basename(fs[0]).split("_")[4]     # e.g. r1i1p1f1

    tos = drop_singletons(ds["tos"], keep=("time",))
    tos = tos.where(np.abs(tos) < 1e19)                     # undecoded fill guard

    years = ds["time"].dt.year.values
    idx = np.where((years >= 1971) & (years <= 2000))[0]
    tos_p = tos.isel(time=idx)
    n_months[m] = int(tos_p.sizes["time"])

    tos_clim = tos_p.mean(dim="time", skipna=True)          # cell-wise climatology
    tos_clim.load()
    tos_clim = tos_clim.assign_attrs(units="degrees_celsius")
    tos_clim_native[m] = tos_clim

    dsa, _ = cmip_open(m, "areacello")
    area = drop_singletons(dsa["areacello"], keep=("lat", "lon", "nlat", "nlon", "j", "i"))
    assert area.shape == tos_clim.shape, (m, area.shape, tos_clim.shape)

    tv = np.asarray(tos_clim.values, float)
    av = np.asarray(area.values, float)
    valid = np.isfinite(tv) & np.isfinite(av) & (av > 0)
    w = np.where(valid, av, 0.0)
    native_mean[m] = float(np.sum(np.where(valid, tv, 0.0) * w) / np.sum(w))
    native_area[m] = float(np.sum(w))

    # ---- conservative regridding onto the WOA23 1 degree grid ----
    latn, lonn = coords_of[m]
    latc = np.asarray(tos_clim[latn].values, float)
    lonc = np.asarray(tos_clim[lonn].values, float)
    if latc.ndim == 1:
        lon2d, lat2d = np.meshgrid(lonc, latc)
    else:
        lat2d, lon2d = latc, lonc

    reg = conservative_regrid(tv, lon2d, lat2d, av, woa_lon, woa_lat)
    regridded[m] = xr.DataArray(reg, dims=("lat", "lon"),
                                coords={"lat": woa_lat, "lon": woa_lon},
                                name="tos").assign_attrs(units="degrees_celsius")

    print(f"{m:14s} variant={variants[m]}  n_months={n_months[m]}  "
          f"native area = {native_area[m]/1e12:.4f} x10^12 m2  "
          f"native-grid mean SST = {native_mean[m]:.2f} degC")

# ======================================================================
# 4) shared comparison mask on the WOA23 1 degree grid
#    (build_comparison_mask with ocean_only=True, WOA23 as reference)
# ======================================================================
shared = (landsea > 1) & woa_sst.notnull()
for m in models:
    shared = shared & regridded[m].notnull()
shared = shared.fillna(False).rename("comparison_mask")

mask_count = int(shared.sum())

R_EARTH = 6371.0                                     # km
_dlat = np.deg2rad(1.0); _dlon = np.deg2rad(1.0)
_edges = np.deg2rad(np.arange(-90.0, 91.0, 1.0))
_band = (R_EARTH ** 2) * _dlon * (np.sin(_edges[1:]) - np.sin(_edges[:-1]))   # km2 per lat band
cell_area_km2 = np.repeat(_band[:, None], NLON, axis=1)
mask_area_km2 = float(np.sum(cell_area_km2[shared.values]))
print(f"shared comparison mask: {mask_count} cells, "
      f"{mask_area_km2/1e6:.2f} x10^6 km2")

# ======================================================================
# 5) area-weighted metrics on the shared cells (1 deg cell areas ~ cos(lat))
# ======================================================================
w2d = np.repeat(np.cos(np.deg2rad(woa_lat))[:, None], NLON, axis=1)
sel = shared.values
ww = w2d[sel]
ref_vals = np.asarray(woa_sst.where(shared).values, float)[sel]

sst_mean_model = {}
sst_mean_woa = {}
delta_sst = {}
rmse_sst = {}
corr_sst = {}
bias_maps = {}

for m in models:
    mod_da = regridded[m].where(shared)
    mod_vals = np.asarray(mod_da.values, float)[sel]

    sst_mean_model[m] = float(np.sum(ww * mod_vals) / np.sum(ww))
    sst_mean_woa[m] = float(np.sum(ww * ref_vals) / np.sum(ww))

    d = mod_vals - ref_vals
    delta_sst[m] = float(np.sum(ww * d) / np.sum(ww))
    rmse_sst[m] = float(np.sqrt(np.sum(ww * d * d) / np.sum(ww)))

    am = mod_vals - np.sum(ww * mod_vals) / np.sum(ww)
    bm = ref_vals - np.sum(ww * ref_vals) / np.sum(ww)
    corr_sst[m] = float(np.sum(ww * am * bm) /
                        np.sqrt(np.sum(ww * am * am) * np.sum(ww * bm * bm)))

    bias_maps[m] = (regridded[m] - woa_sst).where(shared).assign_attrs(units="degrees_celsius")

woa_sst_shared = woa_sst.where(shared)

# ======================================================================
# 6) ranking by |delta_SST| (unrounded)
# ======================================================================
ranking = sorted(models, key=lambda mm: abs(delta_sst[mm]))

# ---- summary table ---------------------------------------------------------
summary = pd.DataFrame({
    "variant": [variants[m] for m in models],
    "n_months_1971_2000": [n_months[m] for m in models],
    "native_grid_mean_SST_C": [round(native_mean[m], 2) for m in models],
    "native_ocean_area_1e12_m2": [round(native_area[m] / 1e12, 4) for m in models],
    "shared_mask_mean_model_C": [round(sst_mean_model[m], 2) for m in models],
    "shared_mask_mean_WOA_C": [round(sst_mean_woa[m], 2) for m in models],
    "delta_SST_C": [round(delta_sst[m], 2) for m in models],
    "RMSE_C": [round(rmse_sst[m], 3) for m in models],
    "spatial_r": [round(corr_sst[m], 3) for m in models],
}, index=models)
print("\nWOA23 whole-grid unweighted mean = %.2f degC" % woa_global_mean_unweighted)
print(summary.to_string())
print("ranking (smallest -> largest |delta_SST|):", " < ".join(ranking))
