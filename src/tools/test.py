from src.tools.data_loader import load_data_with_geolocator

path = "/data1/lorraine-20419/jingwei_QA/reconstruction/recon_o2_sjtu_jingwei_v1_1.zarr"
data = load_data_with_geolocator(path, time_range=("2000-01-01", "2010-12-31"), depth_range=(0, 100), 
                                geofeature="water_body", loc_name="Arctic Ocean")
print(data)