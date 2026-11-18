# -----------------------------------------------------------------------------
# Figure helpers
# -----------------------------------------------------------------------------
def regularize_for_plot(src, field):
    """Put a field on a 1-degree lon/lat grid for display only."""
    if src['grid_type'] == 'regular':
        da = field
        lon_name = src['lon_coord']
        if float(da[lon_name].min()) < 0.0:
            da = da.assign_coords({lon_name: da[lon_name] % 360.0}).sortby(lon_name)
        return da

    lat2 = np.asarray(src['lat'].values).ravel()
    lon2 = np.asarray(src['lon'].values).ravel()
    vals = np.asarray(field.values).ravel()
    valid = np.isfinite(vals)

    lon_t = np.arange(0.5, 360.0, 1.0)
    lat_t = np.arange(-89.5, 90.0, 1.0)
    LON, LAT = np.meshgrid(lon_t, lat_t)

    grid = griddata(
        (lon2[valid], lat2[valid]),
        vals[valid],
        (LON, LAT),
        method='nearest',
    )
    return xr.DataArray(
        grid,
        coords={'lon': lon_t, 'lat': lat_t},
        dims=('lat', 'lon'),
    )


# -----------------------------------------------------------------------------
# Figure 1: historical OMZ thickness at O2 < 60 umol kg-1
# -----------------------------------------------------------------------------
fig1_fields = {
    'GFDL-ESM4 historical (1971-2000)': regularize_for_plot(
        sources['GFDL-ESM4_historical'],
        omz_thickness_fields['GFDL-ESM4_historical'],
    ),
    'MPI-ESM1-2-LR historical (1971-2000)': regularize_for_plot(
        sources['MPI-ESM1-2-LR_historical'],
        omz_thickness_fields['MPI-ESM1-2-LR_historical'],
    ),
    'WOA23 1971-2000': regularize_for_plot(
        sources['WOA23'],
        omz_thickness_fields['WOA23'],
    ),
}

pool1 = np.concatenate([
    np.asarray(f.values).ravel() for f in fig1_fields.values()
])
pool1 = pool1[np.isfinite(pool1)]
vmax1 = float(np.nanpercentile(pool1, 98)) if pool1.size else 1.0
if vmax1 <= 0.0:
    vmax1 = 1.0

fig1, axes1 = plt.subplots(
    1, 3, figsize=(15, 4.2), constrained_layout=True
)

pcm1 = None
for ax, (title, da) in zip(axes1, fig1_fields.items()):
    pcm1 = ax.pcolormesh(
        da['lon'].values,
        da['lat'].values,
        da.values,
        shading='auto',
        cmap='viridis',
        vmin=0.0,
        vmax=vmax1,
    )
    ax.set_title(title, fontsize=10)
    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(-90.0, 90.0)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

fig1.colorbar(
    pcm1,
    ax=axes1,
    label='OMZ thickness (m)',
    extend='max',
)
fig1.suptitle(
    'Column thickness of water with O$_2$ < 60 $\\mu$mol kg$^{-1}$; '
    'shared sequential scale clipped at pooled 98th percentile',
    fontsize=12,
)

# -----------------------------------------------------------------------------
# Figure 2: ssp245 state and change in OMZ thickness
# -----------------------------------------------------------------------------
change_gfdl = omz_thickness_fields['GFDL-ESM4_ssp245'] - omz_thickness_fields['GFDL-ESM4_historical']
change_mpi = omz_thickness_fields['MPI-ESM1-2-LR_ssp245'] - omz_thickness_fields['MPI-ESM1-2-LR_historical']
cesm_ssp = omz_thickness_fields['CESM2_ssp245']

change_gfdl_reg = regularize_for_plot(sources['GFDL-ESM4_ssp245'], change_gfdl)
change_mpi_reg = regularize_for_plot(sources['MPI-ESM1-2-LR_ssp245'], change_mpi)
cesm_ssp_reg = regularize_for_plot(sources['CESM2_ssp245'], cesm_ssp)

pool2 = np.concatenate([
    np.asarray(change_gfdl_reg.values).ravel(),
    np.asarray(change_mpi_reg.values).ravel(),
])
pool2 = pool2[np.isfinite(pool2)]
vmax2 = float(np.nanpercentile(np.abs(pool2), 98)) if pool2.size else 1.0
if vmax2 <= 0.0:
    vmax2 = 1.0

pool3 = np.asarray(cesm_ssp_reg.values).ravel()
pool3 = pool3[np.isfinite(pool3)]
vmax3 = float(np.nanpercentile(pool3, 98)) if pool3.size else 1.0
if vmax3 <= 0.0:
    vmax3 = 1.0

fig2, axes2 = plt.subplots(
    1, 3, figsize=(15, 4.2), constrained_layout=True
)

pcm_change = None
for ax, da, title in [
    (axes2[0], change_gfdl_reg, 'GFDL-ESM4: ssp245 minus historical'),
    (axes2[1], change_mpi_reg, 'MPI-ESM1-2-LR: ssp245 minus historical'),
]:
    pcm_change = ax.pcolormesh(
        da['lon'].values,
        da['lat'].values,
        da.values,
        shading='auto',
        cmap='RdBu_r',
        vmin=-vmax2,
        vmax=vmax2,
    )
    ax.set_title(title, fontsize=10)
    ax.set_xlim(0.0, 360.0)
    ax.set_ylim(-90.0, 90.0)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

pcm_cesm = axes2[2].pcolormesh(
    cesm_ssp_reg['lon'].values,
    cesm_ssp_reg['lat'].values,
    cesm_ssp_reg.values,
    shading='auto',
    cmap='viridis',
    vmin=0.0,
    vmax=vmax3,
)
axes2[2].set_title(
    'CESM2 ssp245 state only\n(no historical member; change not formable)',
    fontsize=10,
)
axes2[2].set_xlim(0.0, 360.0)
axes2[2].set_ylim(-90.0, 90.0)
axes2[2].set_xlabel('Longitude')
axes2[2].set_ylabel('Latitude')

fig2.colorbar(
    pcm_change,
    ax=axes2[:2],
    label='$\\Delta$ OMZ thickness (m)',
    extend='both',
)
fig2.colorbar(
    pcm_cesm,
    ax=axes2[2],
    label='OMZ thickness (m)',
    extend='max',
)
fig2.suptitle(
    'Change in O$_2$ < 60 $\\mu$mol kg$^{-1}$ OMZ thickness; '
    'change panels share a zero-centred scale clipped at pooled 98th percentile',
    fontsize=12,
)
