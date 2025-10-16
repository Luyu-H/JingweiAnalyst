import pymysql
import pandas as pd

conn = pymysql.connect(
    host='10.10.10.0',
    port=3306,
    database='jingwei_QA',
    user='groupleader',
    password='onlyleaders',
    charset="utf8mb4"
)

cursor = conn.cursor()

sql_insert_observation_metadata = """
INSERT INTO observation_metadata (
    table_name, title, summary, comment, featureType, cdm_data_type, processing_level,
    variable_id, cf_standard_name, variable_long_name, variable_units,
    time_coverage_start, time_coverage_end,
    spatial_extent_label, geospatial_lat_min, geospatial_lat_max,
    geospatial_lon_min, geospatial_lon_max,
    geospatial_vertical_min, geospatial_vertical_max,
    geospatial_lat_units, geospatial_lon_units, geospatial_vertical_units,
    geospatial_vertical_positive, qc_definitions, missing_value,
    institution, project,
    date_created, date_modified,
    creator_name, creator_email, creator_url,
    publisher_name, publisher_url, publisher_email,
    reference, conventions, metadata_conventions, license
)
VALUES (
    %(table_name)s, %(title)s, %(summary)s, %(comment)s, %(featureType)s, %(cdm_data_type)s, %(processing_level)s,
    %(variable_id)s, %(cf_standard_name)s, %(variable_long_name)s, %(variable_units)s,
    %(time_coverage_start)s, %(time_coverage_end)s,
    %(spatial_extent_label)s, %(geospatial_lat_min)s, %(geospatial_lat_max)s,
    %(geospatial_lon_min)s, %(geospatial_lon_max)s,
    %(geospatial_vertical_min)s, %(geospatial_vertical_max)s,
    %(geospatial_lat_units)s, %(geospatial_lon_units)s, %(geospatial_vertical_units)s,
    %(geospatial_vertical_positive)s, %(qc_definitions)s, %(missing_value)s,
    %(institution)s, %(project)s,
    %(date_created)s, %(date_modified)s,
    %(creator_name)s, %(creator_email)s, %(creator_url)s,
    %(publisher_name)s, %(publisher_url)s, %(publisher_email)s,
    %(reference)s, %(conventions)s, %(metadata_conventions)s, %(license)s
)
"""

query = """
SELECT MIN(lat) as min_lat,
       MAX(lat) as max_lat,
       MIN(lon) as min_lon,
       MAX(lon) as max_lon,
       MIN(depth) as min_depth,
       MAX(depth) as max_depth
FROM obs_oxygen_WOD_profile_data
"""
query_res = pd.read_sql(query, conn).iloc[0]

geospatial_lat_min = query_res['min_lat']
geospatial_lat_max = query_res['max_lat']
geospatial_lon_min = query_res['min_lon']
geospatial_lon_max = query_res['max_lon']
geospatial_vertical_min = query_res['min_depth']
geospatial_vertical_max = query_res['max_depth']

obs_data = {
    "table_name": "obs_oxygen_WOD_profile_data",
    "title": "Oxygen profiles collected by CTD (High-resolution Conductivity-Temperature-Depth (CTD) data and high-resolution XCTD data), OSD (Bottle, low-resolution Conductivity-Temperature-Depth (CTD), low-resolution XCTD data, and plankton data) and PFL (Profiling float data, mainly from the Argo program) from 1965 to 2023 from the World Ocean Database",
    "summary": "Quality controlled Oxygen profile data collected by CTD, OSD and PFL from the World Ocean Database 2023 for Mapping Evaluation for Global Ocean Dissolved Oxygen Inventory project",
    "comment": "CTD: High-resolution Conductivity-Temperature-Depth (CTD) data and high-resolution XCTD data. OSD: Bottle, low-resolution Conductivity-Temperature-Depth (CTD), low-resolution XCTD data, and plankton data. PFL: Profiling float data, mainly from the Argo program.",
    "featureType": "Profile",
    "cdm_data_type": "Station",
    "processing_level": "raw, quality control flags applied",
    "variable_id": "o2",
    "cf_standard_name": "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water",
    "variable_long_name": "Dissolved Oxygen Concentration",
    "variable_units": "micromoles_per_kilogram",
    "time_coverage_start": "1965-01-01",
    "time_coverage_end": "2023-12-31",
    "spatial_extent_label": "global ocean",
    "geospatial_lat_min": geospatial_lat_min,
    "geospatial_lat_max": geospatial_lat_max,
    "geospatial_lon_min": geospatial_lon_min,
    "geospatial_lon_max": geospatial_lon_max,
    "geospatial_vertical_min": geospatial_vertical_min,
    "geospatial_vertical_max": geospatial_vertical_max,
    "geospatial_lat_units": "degrees_north",
    "geospatial_lon_units": "degrees_east",
    "geospatial_vertical_units": "meters",
    "geospatial_vertical_positive": "down",
    "qc_definitions": "No quality control flags in this dataset.",
    "missing_value": "Missing values in the original WOD profile NetCDF data were masked, which were excluded when storing the data in the SQL database.",
    "institution": "National Centers for Environmental Information (NCEI)",
    "project": "Mapping Evaluation for Global Ocean Dissolved Oxygen Inventory project",
    "date_created": "2024-05-30",
    "date_modified": "2024-05-30",
    "creator_name": "Ocean Climate Laboratory",
    "creator_email": "NCEI.Services@noaa.gov",
    "creator_url": "http://www.ncei.noaa.gov",
    "publisher_name": "US NOAA NATIONAL CENTERS FOR ENVIRONMENTAL INFORMATION",
    "publisher_url": "http://www.ncei.noaa.gov/",
    "publisher_email": "NCEI.Services@noaa.gov",
    "reference": "Garcia, H.E., Z. Wang, C. Bouchard, S.L. Cross, C.R. Paver, J.R. Reagan, T.P. Boyer, R.A. Locarnini, A.V. Mishonov, O.K. Baranova, D. Seidov, and D. Dukhovskoy (2023). World Ocean Atlas 2023, Volume 3: Dissolved Oxygen, Apparent Oxygen Utilization, and Oxygen Saturation. A. Mishonov Technical Ed., NOAA Atlas NESDIS 91, https://doi.org/10.25923/rb67-ns53.",
    "conventions": "CF-1.6",
    "metadata_conventions": "Unidata Dataset Discovery v1.0",
    "license": "These data are openly available to the public. Please acknowledge the use of these data with the text given in the acknowledgment attribute."
}

sql_insert_reconstruction_metadata = """
INSERT INTO observation_metadata (
    table_name, title, summary, comment, processing_level,
    date_created, date_modified,
    variable_id, cf_standard_name, variable_long_name, variable_units, missing_value,
    time_coverage_start, time_coverage_end, time_resolution,
    spatial_extent_label, 
    geospatial_lat_min, geospatial_lat_max, geospatial_lat_resolution, geospatial_lat_units,
    geospatial_lon_min, geospatial_lon_max, geospatial_lon_resolution, geospatial_lon_units,
    geospatial_vertical_min, geospatial_vertical_max, geospatial_vertical_positive, geospatial_vertical_units
)
VALUES (
    %(table_name)s, %(title)s, %(summary)s, %(comment)s, %(processing_level)s,
    %(date_created)s, %(date_modified)s,
    %(variable_id)s, %(cf_standard_name)s, %(variable_long_name)s, %(variable_units)s, %(missing_value)s,
    %(time_coverage_start)s, %(time_coverage_end)s, %(time_resolution)s,
    %(spatial_extent_label)s, 
    %(geospatial_lat_min)s, %(geospatial_lat_max)s, %(geospatial_lat_resolution)s, %(geospatial_lat_units)s,
    %(geospatial_lon_min)s, %(geospatial_lon_max)s, %(geospatial_lon_resolution)s, %(geospatial_lon_units)s,
    %(geospatial_vertical_min)s, %(geospatial_vertical_max)s, %(geospatial_vertical_positive)s, %(geospatial_vertical_units)s
)
"""

recon_data = {
    "table_name": "",
    "title": "",
    "summary": "Quality controlled Oxygen profile data collected by CTD, OSD and PFL from the World Ocean Database 2023 for Mapping Evaluation for Global Ocean Dissolved Oxygen Inventory project",
    "comment": "CTD: High-resolution Conductivity-Temperature-Depth (CTD) data and high-resolution XCTD data. OSD: Bottle, low-resolution Conductivity-Temperature-Depth (CTD), low-resolution XCTD data, and plankton data. PFL: Profiling float data, mainly from the Argo program.",
    "processing_level": "raw, quality control flags applied",
    "date_created": "2024-05-30",
    "date_modified": "2024-05-30",
    "variable_id": "o2",
    "cf_standard_name": "mole_concentration_of_dissolved_molecular_oxygen_in_sea_water",
    "variable_long_name": "Dissolved Oxygen Concentration",
    "variable_units": "micromoles_per_kilogram",
    "missing_value": "Missing values in the original WOD profile NetCDF data were masked, which were excluded when storing the data in the SQL database.",
    "time_coverage_start": "1965-01-01",
    "time_coverage_end": "2023-12-31",
    "time_resolution": "",
    "spatial_extent_label": "global ocean",
    "geospatial_lat_min": geospatial_lat_min,
    "geospatial_lat_max": geospatial_lat_max,
    "geospatial_lat_resolution": "",
    "geospatial_lat_units": "degrees_north",
    "geospatial_lon_min": geospatial_lon_min,
    "geospatial_lon_max": geospatial_lon_max,
    "geospatial_lon_resolution": "",
    "geospatial_lon_units": "degrees_east",
    "geospatial_vertical_min": geospatial_vertical_min,
    "geospatial_vertical_max": geospatial_vertical_max,
    "geospatial_vertical_positive": "down",
    "geospatial_vertical_units": "meters"
}

# cursor.execute(sql_insert_observation_metadata, obs_data)
cursor.execute(sql_insert_reconstruction_metadata, recon_data)
conn.commit()

print("Inserted one record into observation_metadata successfully.")

cursor.close()
conn.close()
