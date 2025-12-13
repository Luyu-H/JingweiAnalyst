import numpy as np
import pandas as pd
import xarray as xr
from dask.diagnostics import ProgressBar
import zarr
import glob
import shutil
import os

def create_time_from_year_month(ds, reference_date='1965-01-01'):
    years = ds['year'].values
    months = ds['month'].values
    
    times = []
    for year in years:
        for month in months:
            date = pd.Timestamp(year=int(year), month=int(month), day=1)
            times.append(date)
    
    reference = pd.Timestamp(reference_date)
    days_since_ref = [(t - reference).days for t in times]
    
    # (year, month, depth, lat, lon) -> (time, depth, lat, lon)
    data_vars = {}
    for var_name in ds.data_vars:
        var = ds[var_name]
        
        if 'year' in var.dims and 'month' in var.dims:
            stacked = var.stack(time=('year', 'month'))
            
            dim_order = ['time']
            if 'depth' in stacked.dims:
                dim_order.append('depth')
            if 'lat' in stacked.dims:
                dim_order.append('lat')
            if 'lon' in stacked.dims:
                dim_order.append('lon')
            
            stacked = stacked.transpose(*dim_order)
            data_vars[var_name] = stacked
        else:
            data_vars[var_name] = var
    
    ds_new = xr.Dataset(data_vars)
    
    ds_new = ds_new.drop_vars(['time', 'year', 'month'], errors='ignore')
    
    ds_new = ds_new.assign_coords(
        time=xr.DataArray(
            days_since_ref,
            dims='time',
            attrs={
                'units': f'days since {reference_date} 00:00:00',
                'calendar': 'proleptic_gregorian',
                'long_name': 'time'
            }
        )
    )
    
    for coord in ['depth', 'lat', 'lon']:
        if coord in ds.coords:
            ds_new = ds_new.assign_coords({coord: ds[coord]})
    
    return ds_new

def convert_nc_to_zarr(nc_path, zarr_path, has_year_month=False):
    ds = xr.open_dataset(nc_path, chunks={})

    if has_year_month:
        print(f"  Converting year+month to time...")
        ds = create_time_from_year_month(ds)

    chunk_dict = {
        "time": 12,  # 1 year of monthly data
        "lat": 45,
        "lon": 60,
        "depth": 102  # keep depth unchunked
    }
    
    chunk_dict = {k: v for k, v in chunk_dict.items() if k in ds.dims}
    ds = ds.chunk(chunk_dict)

    if os.path.exists(zarr_path):
        shutil.rmtree(zarr_path)

    with ProgressBar():
        ds.to_zarr(zarr_path, mode="w", consolidated=True, zarr_version=2)
    
    zarr.consolidate_metadata(zarr_path)
    ds.close()

def merge_nc_files_to_zarr(nc_files, zarr_path, has_year_month=True):
    if os.path.exists(zarr_path):
        shutil.rmtree(zarr_path)
    
    print(f"[1/{len(nc_files)}] processing: {nc_files[0]}")
    convert_nc_to_zarr(nc_files[0], zarr_path, has_year_month=has_year_month)
    
    for i, nc_file in enumerate(nc_files[1:], start=2):
        print(f"[{i}/{len(nc_files)}] is added: {nc_file}")
        
        ds = xr.open_dataset(nc_file, chunks={})

        if has_year_month:
            print(f"  Converting year+month to time...")
            ds = create_time_from_year_month(ds)

        chunk_dict = {
            "time": 12,
            "lat": 45,
            "lon": 60,
            "depth": 102
        }
        
        chunk_dict = {k: v for k, v in chunk_dict.items() if k in ds.dims}
        ds = ds.chunk(chunk_dict)
        
        with ProgressBar():
            ds.to_zarr(zarr_path, mode="a", append_dim="time")
        
        ds.close()
    
    zarr.consolidate_metadata(zarr_path)

    ds_result = xr.open_zarr(zarr_path)
    print(f"\nFinal info:")
    print(f"  Dimensions: {dict(ds_result.dims)}")
    print(f"  Time range: {ds_result.time.values[0]} - {ds_result.time.values[-1]}")
    print(f"  Variables: {list(ds_result.data_vars.keys())}")
    ds_result.close()


# nc_path = "/data1/wangzhaokun-20382/Foundation_Model_results/1.1/Si_1960-2024.nc"
# zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_si_sjtu_jingwei_v1_1.zarr"
# convert_nc_to_zarr(nc_path, zarr_path)

nc_files = [
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1900-1909.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1910-1919.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1920-1929.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1930-1939.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1940-1949.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1950-1959.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1960-1969.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1970-1979.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1980-1989.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_1990-1999.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_2000-2009.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_2010-2019.nc",
    "/data1/wangzhaokun-20382/Foundation_Model_gridded_data/EN4/en4_decade_2020-2024.nc"
]
zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_temp_sjtu_jingwei_v1_1.zarr"
merge_nc_files_to_zarr(nc_files, zarr_path)