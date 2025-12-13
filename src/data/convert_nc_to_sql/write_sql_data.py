import pymysql
from netCDF4 import Dataset

import os
import math
import numpy as np
import numpy.ma as ma

conn = pymysql.connect(
    host='10.10.10.0',
    port=3306,
    database='jingwei_QA',
    user='groupleader',
    password='onlyleaders',
    charset="utf8mb4"
)
cursor = conn.cursor()

sql_insert_obs_oxygen_WOD_profile_data = """
INSERT INTO obs_oxygen_WOD_profile_data
(profile, cruise, collection, time, year, month, day, hour, lat, lon, depth, oxygen)
VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
"""

sql_insert_recon_o2_sjtu_jingwei_v1_1 = """
INSERT INTO recon_o2_sjtu_jingwei_v1_1
(time, lat, lon, depth, o2)
VALUES (%s,%s,%s,%s,%s)
"""

def clean_missing_value(value):
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    try:
        if isinstance(value, (float, int)) and np.isnan(value):
            return None
    except:
        pass
    if isinstance(value, str) and value.lower() == 'nan':
        return None
    
    return value

def insert_obs_oxygen_WOD_profile_data(nc_dir):
    for filename in os.listdir(nc_dir):
        if not filename.endswith(".nc"):
            continue
        nc_path = os.path.join(nc_dir, filename)

        try:
            nc = Dataset(nc_path, "r")

            profiles = nc.variables["profile"][:]
            cruises = nc.variables["cruise"][:]
            times = nc.variables["time"][:]
            years = nc.variables["year"][:]
            months = nc.variables["month"][:]
            days = nc.variables["day"][:]
            hours = nc.variables["hour"][:]
            lats = nc.variables["lat"][:]
            lons = nc.variables["lon"][:]
            depths = nc.variables["depth"][:]
            oxygen = nc.variables["Oxygen"][:]

            collection = 'PFL'

            batch = []
            for i in range(len(profiles)):
                for j in range(len(depths)):
                    val = oxygen[i, j]
                    if np.isnan(val) or ma.is_masked(val):
                        continue

                    row = (
                        profiles[i],
                        cruises[i],
                        collection,
                        times[i],
                        years[i],
                        months[i],
                        days[i],
                        hours[i],
                        lats[i],
                        lons[i],
                        depths[j],
                        val
                    )
                    cleaned_row = tuple(clean_missing_value(v) for v in row)
                    batch.append(cleaned_row)
            
            cursor.executemany(sql_insert_obs_oxygen_WOD_profile_data, batch)
            conn.commit()
            nc.close()
            print(f"Inserted data from {filename} successfully.")
        
        except Exception as e:
            print(f"Error processing {filename}: {e}")

    cursor.close()
    conn.close()
    print("Data insertion completed.")

def insert_recon_o2_sjtu_jingwei_v1_1(nc_path):
    try:
        nc = Dataset(nc_path, "r")

        times = nc.variables["time"][:]
        lats = nc.variables["lat"][:]
        lons = nc.variables["lon"][:]
        depths = nc.variables["depth"][:]
        o2 = nc.variables["o2"]

        print("Start processing...")
        for it, t in enumerate(times):
            data_t = o2[it, :, :, :]
            
            lat_grid, lon_grid, depth_grid = np.meshgrid(lats, lons, depths, indexing="ij")
            
            rows = np.column_stack([
                np.full(lat_grid.size, t),
                lat_grid.ravel(),
                lon_grid.ravel(),
                depth_grid.ravel(),
                data_t.ravel()
            ])

            batch = [tuple(clean_missing_value(v) for v in row) for row in rows]

            cursor.executemany(sql_insert_recon_o2_sjtu_jingwei_v1_1, batch)
            conn.commit()
            print(f"Inserted time index {it}, time={t}")

        nc.close()
        cursor.close()
        conn.close()
        print("Data insertion completed.")

    except Exception as e:
        print(f"Error processing the file: {e}")

    
if __name__ == "__main__":
    nc_file = '/data1/xinyi-20379/dataset/jingwei/results/products/o2map_v1_5500m_sjtu_jingwei_compressed.nc'
    insert_recon_o2_sjtu_jingwei_v1_1(nc_file)