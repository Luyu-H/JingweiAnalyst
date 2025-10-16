import xarray as xr
from dask.diagnostics import ProgressBar
import zarr

def convert_nc_to_zarr(nc_path, zarr_path):
    ds = xr.open_dataset(nc_path, chunks={})
    ds = ds.chunk({
        "time": 12, # 1 year of monthly data
        "lat": 45,
        "lon": 60,
        "depth": 102 # keep depth unchunked
    })

    import shutil, os
    if os.path.exists(zarr_path):
        shutil.rmtree(zarr_path)

    with ProgressBar():
        ds.to_zarr(zarr_path, mode="w", consolidated=True, zarr_version=2)
    
    zarr.consolidate_metadata(zarr_path)


nc_path = "/data1/wangzhaokun-20382/Foundation_Model_results/1.1/Si_1960-2024.nc"
zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_si_sjtu_jingwei_v1_1.zarr"
convert_nc_to_zarr(nc_path, zarr_path)

nc_path = " /data1/ziki37-20378/jingwei/results/FMv1/train_and_inference/20250924_130043_by_direct/EN4_inference/merged.nc"
zarr_path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_dissic_sjtu_jingwei_v1_1.zarr"
convert_nc_to_zarr(nc_path, zarr_path)