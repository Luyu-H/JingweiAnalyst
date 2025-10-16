import mysql.connector
from mysql.connector import Error

def create_connection():
    try:
        connection = mysql.connector.connect(
            host='10.10.10.0',
            port=3306,
            database='jingwei_QA',
            user='groupleader',
            password='onlyleaders'
        )
        if connection.is_connected():
            print("successfully connected to the database!")
            return connection
    except Error as e:
        print(f"Error: {e}")
        return None

def execute_sql(connection, sql):
    try:
        cursor = connection.cursor()
        cursor.execute(sql)
        connection.commit()
        print("SQL executed successfully!")
    except Error as e:
        print(f"Execution error: {e}")
    finally:
        if cursor:
            cursor.close()

sql_create_observation_metadata = """
CREATE TABLE IF NOT EXISTS observation_metadata (
    table_name VARCHAR(128) PRIMARY KEY COMMENT 'Name of the observation data table.',
    title TEXT COMMENT 'Title of the dataset.',
    summary TEXT COMMENT 'Summary of the dataset, including description of data and purpose.',
    comment TEXT COMMENT 'Additional comments about the dataset.',
    featureType VARCHAR(64) COMMENT 'The type of data model as defined by the CF Conventions. Profile indicates the data is a set of observations at a point location that vary with depth or height.',
    cdm_data_type VARCHAR(64) COMMENT 'Describes how the data is organized according to the NetCDF Common Data Model. Station means the data is grouped by measurement locations (stations).',
    processing_level VARCHAR(1000) COMMENT 'A description of the level of processing the data has undergone, such as raw, calibrated, or quality controlled.',
    variable_id VARCHAR(128) COMMENT 'A unique identifier for the variable.',
    cf_standard_name VARCHAR(500) COMMENT 'The CF standard name for the variable, following CF Conventions.',
    variable_long_name VARCHAR(500) COMMENT 'A descriptive name for the variable.',
    variable_units VARCHAR(64) COMMENT 'Units of the variable.',
    time_coverage_start DATE COMMENT 'The start date of the data coverage.',
    time_coverage_end DATE COMMENT 'The end date of the data coverage.',
    spatial_extent_label VARCHAR(128) COMMENT 'A label describing the spatial extent of the dataset.',
    geospatial_lat_min DOUBLE COMMENT 'Minimum latitude in geospatial_lat_units.',
    geospatial_lat_max DOUBLE COMMENT 'Maximum latitude in geospatial_lat_units.',
    geospatial_lon_min DOUBLE COMMENT 'Minimum longitude in geospatial_lon_units.',
    geospatial_lon_max DOUBLE COMMENT 'Maximum longitude in geospatial_lon_units.',
    geospatial_vertical_min DOUBLE COMMENT 'The minimum depth or height in geospatial_vertical_units.',
    geospatial_vertical_max DOUBLE COMMENT 'The maximum depth or height in geospatial_vertical_units.',
    geospatial_lat_units VARCHAR(128) COMMENT 'Units for latitude, typically degrees_north.',
    geospatial_lon_units VARCHAR(128) COMMENT 'Units for longitude, typically degrees_east.',
    geospatial_vertical_units VARCHAR(128) COMMENT 'Units for vertical dimension, e.g., meters.',
    geospatial_vertical_positive VARCHAR(16) COMMENT 'Direction of positive vertical values, either up or down.',
    qc_definitions TEXT COMMENT 'Definitions of quality control flags used in the dataset.',
    missing_value VARCHAR(1000) COMMENT 'Value used to represent missing data.',
    institution VARCHAR(1000) COMMENT 'Name of the institution that created the dataset.',
    project VARCHAR(128) COMMENT 'The name of the research project or program that supported the generation of this dataset.',
    date_created DATE COMMENT 'The date the dataset was created.',
    date_modified DATE COMMENT 'The date the dataset was last modified.',
    creator_name VARCHAR(128) COMMENT 'Name of the person or organization that created the dataset.',
    creator_email VARCHAR(128) COMMENT 'Email address of the creator.',
    creator_url VARCHAR(500) COMMENT 'URL for the creator or organization.',
    publisher_name VARCHAR(128) COMMENT 'Name of the person or organization that published the dataset.',
    publisher_url VARCHAR(500) COMMENT 'URL for the publisher or organization.',
    publisher_email VARCHAR(128) COMMENT 'Email address of the publisher.',
    reference TEXT COMMENT 'A citation for the scientific publication, report, or Digital Object Identifier (DOI) that describes the methodology or is the authoritative source for this data.',
    conventions VARCHAR(64) COMMENT 'Metadata or file format standards that the dataset follows.',
    metadata_conventions VARCHAR(64) COMMENT 'The conventions followed for the metadata.',
    license TEXT COMMENT 'Terms of use, copyright, and licensing under which the data is made available.'
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
"""

sql_create_observation = """
CREATE TABLE IF NOT EXISTS obs_oxygen_WOD_profile_data (
    id BIGINT AUTO_INCREMENT PRIMARY KEY COMMENT 'Unique identifier for each observation',
    profile INT COMMENT 'Unique identifier for each profile',
    cruise INT COMMENT 'Unique identifier for each cruise defined in the World Ocean Database',
    collection VARCHAR(64) COMMENT 'Connection method of the data (e.g., CTD, OSD, etc.)',
    time DOUBLE COMMENT 'Time in days since 1900-01-01 00:00:00',
    year INT COMMENT 'Year of the time',
    month INT COMMENT 'Month of the time',
    day INT COMMENT 'Day of the time',
    hour DOUBLE COMMENT 'Hour of the time',
    lat DOUBLE COMMENT 'Latitude in degrees_north',
    lon DOUBLE COMMENT 'Longitude in degrees_east',
    depth DOUBLE COMMENT 'Standard depth in meters in the World Ocean Database',
    oxygen DOUBLE COMMENT 'Values at each standard depth level for moles_of_oxygen_per_unit_mass_in_sea_water in micromoles_per_kilogram'
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
"""

sql_create_reconstruction_metadata = """
CREATE TABLE IF NOT EXISTS reconstruction_metadata (
    table_name VARCHAR(128) PRIMARY KEY COMMENT 'Name of the observation data table.',
    title TEXT COMMENT 'Title of the dataset.',
    summary TEXT COMMENT 'Summary of the dataset, including description of data and purpose.',
    comment TEXT COMMENT 'Additional comments about the dataset.',
    processing_level VARCHAR(1000) COMMENT 'A description of the level of processing the data has undergone, such as raw, calibrated, or quality controlled.',
    date_created DATE COMMENT 'The date the dataset was created.',
    date_modified DATE COMMENT 'The date the dataset was last modified.',
    variable_id VARCHAR(128) COMMENT 'A unique identifier for the variable.',
    cf_standard_name VARCHAR(500) COMMENT 'The CF standard name for the variable, following CF Conventions.',
    variable_long_name VARCHAR(500) COMMENT 'A descriptive name for the variable.',
    variable_units VARCHAR(64) COMMENT 'Units of the variable.',
    missing_value VARCHAR(1000) COMMENT 'Value used to represent missing data.',
    time_coverage_start DATE COMMENT 'The start date of the data coverage.',
    time_coverage_end DATE COMMENT 'The end date of the data coverage.',
    time_resolution VARCHAR(64) COMMENT 'The temporal resolution of the data, e.g., daily, monthly, yearly.',
    spatial_extent_label VARCHAR(128) COMMENT 'A label describing the spatial extent of the dataset.',
    geospatial_lat_min DOUBLE COMMENT 'Minimum latitude in geospatial_lat_units.',
    geospatial_lat_max DOUBLE COMMENT 'Maximum latitude in geospatial_lat_units.',
    geospatial_lat_resolution VARCHAR(64) COMMENT 'Resolution of latitude in geospatial_lat_units.',
    geospatial_lat_units VARCHAR(128) COMMENT 'Units for latitude, typically degrees_north.',
    geospatial_lon_min DOUBLE COMMENT 'Minimum longitude in geospatial_lon_units.',
    geospatial_lon_max DOUBLE COMMENT 'Maximum longitude in geospatial_lon_units.',
    geospatial_lon_resolution VARCHAR(64) COMMENT 'Resolution of longitude in geospatial_lon_units.',
    geospatial_lon_units VARCHAR(128) COMMENT 'Units for longitude, typically degrees_east.',
    geospatial_vertical_min DOUBLE COMMENT 'The minimum depth or height in geospatial_vertical_units.',
    geospatial_vertical_max DOUBLE COMMENT 'The maximum depth or height in geospatial_vertical_units.',
    geospatial_vertical_positive VARCHAR(16) COMMENT 'Direction of positive vertical values, either up or down.',
    geospatial_vertical_units VARCHAR(128) COMMENT 'Units for vertical dimension, e.g., meters.'
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
"""

sql_create_reconstruction = """
CREATE TABLE IF NOT EXISTS recon_o2_sjtu_jingwei_v1_2 (
    id BIGINT AUTO_INCREMENT PRIMARY KEY COMMENT 'Unique identifier for each observation',
    time DATETIME COMMENT 'Time in days since 1965-01-01 00:00:00 (calendar = proleptic_gregorian)',
    lat DOUBLE COMMENT 'Latitude in degrees_north',
    lon DOUBLE COMMENT 'Longitude in degrees_east',
    depth DOUBLE COMMENT 'Standard 0 - 5500 m depth levels in WOA',
    o2 DOUBLE COMMENT 'Values at each standard depth level for moles_of_oxygen_per_unit_mass_in_sea_water in μmol·kg^-1'
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
"""

if __name__ == "__main__":
    conn = create_connection()
    if conn:
        # execute_sql(conn, sql_create_observation_metadata)
        # execute_sql(conn, sql_create_observation)
        # execute_sql(conn, sql_create_reconstruction_metadata)
        execute_sql(conn, sql_create_reconstruction)
        conn.close()