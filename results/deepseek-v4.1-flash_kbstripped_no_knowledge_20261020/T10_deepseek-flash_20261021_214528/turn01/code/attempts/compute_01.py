import numpy as np
import xarray as xr

HIST = ['1971-01', '2000-12']
SSP  = ['2071-01', '2100-12']
hist_models = ['GFDL-ESM4', 'MPI-ESM1-2-LR']
ssp_models  = ['CESM2', 'GFDL-ESM4', 'MPI-ESM1-2-LR']

sig3 = lambda v: float(f"{v:.3g}")

# ------------------------------------------------------------------ #
# 1. WOA23 observed climatology, masks and vertical grid
# ------------------------------------------------------------------ #
woa_o2  = load_climatology_data(dataset='woa23', variable='oxygen',
                                period='1971-2000', time_slice='annual')
woa_pct = load_climatology_data(dataset='woa23', variable='percent_oxygen_saturation',
                                period='1971-2000', time_slice='annual')
woa_aou = load_climatology_data(dataset='woa23', variable='apparent_oxygen_utilization',
                                period='1971-2000', time_slice='annual')
woa_levels, woa_bounds = load_depth_levels(dataset='woa23', variable='oxygen',
                                           time_slice='annual')
landsea    = load_landsea_mask()
basin_mask = load_basin_mask()

# ------------------------------------------------------------------ #
# helpers
# ------------------------------------------------------------------ #
def omz_diagnostics(o2_field):
    """OMZ volume at 20 & 60 umol/kg, thickness fields and max thickness."""
    out = {}
    for thr in (20.0, 60.0):
        ex = compute_exceedance_fraction(data=o2_field, threshold=thr,
                                         comparison='less', weighting='volume', dim=None)
        out[f'volume_{int(thr)}_m3']    = float(ex['exceeding_amount'])
        out[f'volume_{int(thr)}_1e15m3'] = out[f'volume_{int(thr)}_m3'] / 1e15
    th60, m60 = identify_oxygen_minimum_zone(o2_data=o2_field, threshold=60.0)
    th20, m20 = identify_oxygen_minimum_zone(o2_data=o2_field, threshold=20.0)
    out['thickness60'] = th60
    out['mask60']      = m60
    out['thickness20'] = th20
    out['mask20']      = m20
    out['max60']       = locate_extreme(th60, op='max')
    return out


def mask_volume(mask, vol):
    """Independent cross-check: sum of native cell volumes inside the mask."""
    m = mask.astype('float64').fillna(0.0)
    return float((m * vol.fillna(0.0)).sum())


def prep_cmip(model, dataset, period):
    """Load o2/thetao/so, day-weighted long-term mean, density, o2 -> umol/kg."""
    o2 = load_regular_data(dataset=dataset, variable='o2',     model=model, time_range=period)
    th = load_regular_data(dataset=dataset, variable='thetao', model=model, time_range=period)
    so = load_regular_data(dataset=dataset, variable='so',     model=model, time_range=period)

    vol = load_volumes_for_cmip6(dataset=dataset, model=model, variable='o2')
    if 'time' in vol.dims:
        vol = vol.isel(time=0)
    areas = load_areas(dataset=dataset, model=model, variable='o2')
    if 'time' in areas.dims:
        areas = areas.isel(time=0)
    lev, bnds = load_depth_levels(dataset=dataset, model=model, variable='o2')

    o2m = temporal_aggregate(data=o2, freq='mean', op='mean')
    thm = temporal_aggregate(data=th, freq='mean', op='mean')
    som = temporal_aggregate(data=so, freq='mean', op='mean')

    sigma = compute_potential_density_anomaly(temp_data=thm, sal_data=som, ref_pressure=0)

    o2_kg = transform_unit(data=o2m, from_unit=o2m.attrs['units'],
                           to_unit='umol/kg', density=sigma)
    # keep the grid provenance the weighting tools need
    for k, v in o2m.attrs.items():
        o2_kg.attrs.setdefault(k, v)
    o2_kg.attrs['units'] = 'umol/kg'

    return dict(o2=o2_kg, sigma=sigma, vol=vol, areas=areas,
                levels=lev, bounds=bnds, src_attrs=dict(o2m.attrs))


def ocean_inventory(d):
    """Global dissolved O2 inventory in moles (umol/kg -> mol)."""
    o2_kg, sigma, vol = d['o2'], d['sigma'], d['vol']
    rho = sigma + 1000.0                       # kg m-3
    inv_manual = float((o2_kg * rho * vol).sum(skipna=True)) / 1e6
    inv_tool = float(convert_concentration_to_content(concentration=o2_kg,
                                                      density=sigma).sum(skipna=True))
    if np.isfinite(inv_tool) and inv_manual > 0 and 0.5 < inv_tool / inv_manual < 2.0:
        return inv_tool, inv_manual
    return inv_manual, inv_tool


def basin_volumes(o2_field, thr=60.0):
    out = {}
    for b in ['Pacific', 'Indian', 'Atlantic']:
        sub = select_basin(data=o2_field, basin=b, include_marginal_seas=True,
                           depth_dependent=False)
        ex = compute_exceedance_fraction(data=sub, threshold=thr, comparison='less',
                                         weighting='volume', dim=None)
        out[b] = float(ex['exceeding_amount']) / 1e15
    out['Pacific+Indian+Atlantic'] = out['Pacific'] + out['Indian'] + out['Atlantic']
    return out

# ------------------------------------------------------------------ #
# 2. WOA23 OMZ
# ------------------------------------------------------------------ #
omz = {}
omz['WOA23'] = omz_diagnostics(woa_o2)

# ------------------------------------------------------------------ #
# 3. CMIP6 historical (GFDL-ESM4, MPI-ESM1-2-LR)
# ------------------------------------------------------------------ #
cmip = {}
for m in hist_models:
    d = prep_cmip(m, 'cmip6_historical', HIST)
    cmip[f'{m}_historical'] = d
    omz[f'{m}_historical'] = omz_diagnostics(d['o2'])
    omz[f'{m}_historical']['xcheck_60_1e15m3'] = \
        mask_volume(omz[f'{m}_historical']['mask60'], d['vol']) / 1e15

# ------------------------------------------------------------------ #
# 4. CMIP6 ssp245 (CESM2, GFDL-ESM4, MPI-ESM1-2-LR)
# ------------------------------------------------------------------ #
for m in ssp_models:
    d = prep_cmip(m, 'cmip6_ssp245', SSP)
    cmip[f'{m}_ssp245'] = d
    omz[f'{m}_ssp245'] = omz_diagnostics(d['o2'])
    omz[f'{m}_ssp245']['xcheck_60_1e15m3'] = \
        mask_volume(omz[f'{m}_ssp245']['mask60'], d['vol']) / 1e15

# ------------------------------------------------------------------ #
# 5. Percentage change historical -> ssp245
# ------------------------------------------------------------------ #
pct_change = {}
for m in hist_models:
    h = omz[f'{m}_historical']
    s = omz[f'{m}_ssp245']
    for thr in (20, 60):
        pct_change[(m, thr)] = 100.0 * (s[f'volume_{thr}_m3'] - h[f'volume_{thr}_m3']) \
                               / h[f'volume_{thr}_m3']

# ------------------------------------------------------------------ #
# 6. Global dissolved oxygen inventories (moles)
# ------------------------------------------------------------------ #
inventory = {}
for m in hist_models:
    inventory[f'{m}_historical'], _ = ocean_inventory(cmip[f'{m}_historical'])
for m in ssp_models:
    inventory[f'{m}_ssp245'], _ = ocean_inventory(cmip[f'{m}_ssp245'])

# ------------------------------------------------------------------ #
# 7. Basin-resolved OMZ volume at the 60 threshold
# ------------------------------------------------------------------ #
basin_vols = {}
basin_vols['WOA23'] = basin_volumes(woa_o2)
for m in hist_models:
    basin_vols[f'{m}_historical'] = basin_volumes(cmip[f'{m}_historical']['o2'])

# ------------------------------------------------------------------ #
# 8. WOA23 AOU / O2 saturation inside the 60-threshold OMZ
# ------------------------------------------------------------------ #
m60_woa   = omz['WOA23']['mask60']
aou_in    = woa_aou.where(m60_woa)
sat_in    = woa_pct.where(m60_woa)
mean_aou  = float(compute_basic_stat(data=aou_in, op='volume_weighted_mean', dim=None))
mean_sat  = float(compute_basic_stat(data=sat_in, op='volume_weighted_mean', dim=None))

# ------------------------------------------------------------------ #
# 9. Model-vs-observation comparison
# ------------------------------------------------------------------ #
obs_vol60 = omz['WOA23']['volume_60_1e15m3']
model_err = {m: abs(omz[f'{m}_historical']['volume_60_1e15m3'] - obs_vol60)
             for m in hist_models}
better_model = min(model_err, key=model_err.get)

# ------------------------------------------------------------------ #
# 10. Common-grid thickness fields for the figures
# ------------------------------------------------------------------ #
def retag(da, src_attrs, units='m'):
    out = da.copy()
    out.attrs.update(src_attrs)
    out.attrs['units'] = units
    return out

th_gfdl_h = retag(omz['GFDL-ESM4_historical']['thickness60'],
                  cmip['GFDL-ESM4_historical']['src_attrs'])
th_mpi_h  = retag(omz['MPI-ESM1-2-LR_historical']['thickness60'],
                  cmip['MPI-ESM1-2-LR_historical']['src_attrs'])
th_woa    = retag(omz['WOA23']['thickness60'], dict(woa_o2.attrs))
th_gfdl_s = retag(omz['GFDL-ESM4_ssp245']['thickness60'],
                  cmip['GFDL-ESM4_ssp245']['src_attrs'])
th_mpi_s  = retag(omz['MPI-ESM1-2-LR_ssp245']['thickness60'],
                  cmip['MPI-ESM1-2-LR_ssp245']['src_attrs'])
th_cesm_s = retag(omz['CESM2_ssp245']['thickness60'],
                  cmip['CESM2_ssp245']['src_attrs'])

th_gfdl_h_1d = regrid_horizontal(th_gfdl_h, target='woa23', method='conservative')
th_mpi_h_1d  = regrid_horizontal(th_mpi_h,  target='woa23', method='conservative')
th_woa_1d    = regrid_horizontal(th_woa,    target='woa23', method='conservative')
th_gfdl_s_1d = regrid_horizontal(th_gfdl_s, target='woa23', method='conservative')
th_mpi_s_1d  = regrid_horizontal(th_mpi_s,  target='woa23', method='conservative')
th_cesm_s_1d = regrid_horizontal(th_cesm_s, target='woa23', method='conservative')

ch_gfdl = (th_gfdl_s_1d - th_gfdl_h_1d)
ch_mpi  = (th_mpi_s_1d  - th_mpi_h_1d)

# ------------------------------------------------------------------ #
# 11. Report
# ------------------------------------------------------------------ #
labels = ['WOA23', 'GFDL-ESM4_historical', 'MPI-ESM1-2-LR_historical',
          'CESM2_ssp245', 'GFDL-ESM4_ssp245', 'MPI-ESM1-2-LR_ssp245']

print("=" * 78)
print("OMZ volumes (1e15 m3, O2 strictly below threshold)")
print("=" * 78)
for lb in labels:
    o = omz[lb]
    print(f"{lb:24s}  <20: {sig3(o['volume_20_1e15m3']):>10}   "
          f"<60: {sig3(o['volume_60_1e15m3']):>10}   "
          f"(cross-check <60: {sig3(o['xcheck_60_1e15m3']) if 'xcheck_60_1e15m3' in o else 'n/a'})")

print("\nMax OMZ thickness at the 60 threshold")
for lb in labels:
    mx = omz[lb]['max60']
    c = mx['coordinates']
    print(f"{lb:24s}  {round(float(mx['value']),1):>8} m   "
          f"lat {c.get('latitude'):.2f}  lon {c.get('longitude'):.2f}")

print("\nPercentage change historical -> ssp245 (positive = expanding)")
for m in hist_models:
    print(f"{m:24s}  <20: {pct_change[(m,20)]:+.2f} %    <60: {pct_change[(m,60)]:+.2f} %")

print("\nGlobal dissolved oxygen inventory (mol)")
for k, v in inventory.items():
    print(f"{k:32s}  {sig3(v):.3e}")

print("\nBasin OMZ volume at the 60 threshold (1e15 m3)")
for k, bv in basin_vols.items():
    print(f"{k:24s}  Pac {sig3(bv['Pacific']):>8}  Ind {sig3(bv['Indian']):>8}  "
          f"Atl {sig3(bv['Atlantic']):>8}  Sum {sig3(bv['Pacific+Indian+Atlantic']):>8}")

print(f"\nWOA23 in-OMZ (60) volume-weighted AOU = {mean_aou:.1f} umol/kg, "
      f"O2 saturation = {mean_sat:.1f} %")
print(f"Better historical model vs WOA23: {better_model} "
      f"(|error| = {sig3(model_err[better_model])} x1e15 m3)")
