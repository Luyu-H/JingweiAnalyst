import os
import logging
import tqdm
import math
from typing import List, Dict
import numpy as np
import geopandas as gpd
import shapely
from haversine import haversine_vector, Unit
from shapely.geometry import Polygon, MultiPolygon, LineString, MultiLineString, Point, MultiPoint


EXCLUDE_COUNTRIES = [] # 'Antarctica'?
EXCLUDE_STATES = [] # 'Antarctica'?
EXCLUDE_CONTINENTS = [] # 'Seven seas (open ocean)'?
EXCLUDE_WATER_BODIES = []

# The 'admin' column in the state dataframe in Natural Earth is different from the 'ADMIN' column in the country dataframe.
# This dictionary is used to map the names of countries in the state dataframe to their corresponding names in the country dataframe.
DIFF_STATE_ADMIN_TO_COUNTRY_ADMIN = {
    'Cape Verde': 'Cabo Verde',
    'Caribbean Netherlands': 'Netherlands',
    'Czech Republic': 'Czechia',
    'Gaza Strip': 'Palestine',
    'Guinea Bissau': 'Guinea-Bissau',
    'Macau S.A.R': 'Macao S.A.R',
    'Macedonia': 'North Macedonia',
    'S. Sudan': 'South Sudan',
    'Sao Tome and Principe': 'São Tomé and Principe',
    'Spratly Is.': 'Spratly Islands',
    'Swaziland': 'eSwatini',
    'West Bank': 'Palestine',
}


def get_name_col(df_type: str):
    """
    Get the name column for a given geoframe type.

    Args:
        - df_type (str): The type of geoframe. Can be 'country', 'state', 'continent', or 'water_body'.

    Returns:
        - name (str): The name of the column containing the names of the geofeatures.
    """
    if df_type == 'country':
        return 'ADMIN'
    elif df_type == 'state':
        return 'name'
    elif df_type == 'continent':
        return 'CONTINENT'
    elif df_type == 'water_body':
        return 'name'
    else:
        raise ValueError(f"Invalid geoframe type: {df_type}")

def fix_longitudes_and_split(gdf):
    """
    Fix longitudes and split any geometries that cross the 0° longitude line.
    This function takes a GeoDataFrame and returns a new GeoDataFrame with the same data but with the longitudes
    converted from -180/180 format to 0/360 format and any geometries that cross the 0° longitude line split into
    separate parts.

    Args:
        - gdf (GeoDataFrame): The GeoDataFrame to process.

    Returns:
        - (GeoDataFrame): The processed GeoDataFrame.
    """
    result_gdf = gdf.copy()
    
    def convert_coords(coords_array):
        """Convert a coordinate array from -180/180 to 0/360 format"""
        coords_array = np.array(coords_array)
        coords_array[:, 0] = np.where(coords_array[:, 0] < 0, coords_array[:, 0] + 360, coords_array[:, 0])
        return coords_array
    
    def split_geometry(geometry):
        """Split geometry into parts if it crosses the 0° longitude line."""
        if geometry.is_empty:
            return geometry
        
        if isinstance(geometry, (MultiPolygon, MultiLineString, MultiPoint)):
            list_of_parts = []
            for part in geometry.geoms:
                split_part = split_geometry(part)
                if isinstance(split_part, (Polygon, LineString, Point)):
                    list_of_parts.append(split_part)
                else:
                    list_of_parts.extend(list(split_part.geoms))
            return type(geometry)(list_of_parts)
        
        elif isinstance(geometry, Polygon):
            coords = np.array(geometry.exterior.coords)
            if (coords[:, 0] < 0).all() or (coords[:, 0] >= 0).all():
                # Convert all coordinates if they are all in one hemisphere
                return Polygon(convert_coords(coords))
            
            # Otherwise, split the polygon at the 0° longitude
            splitter = LineString([(0, -90), (0, 90)])  # Extend vertically
            split_parts = shapely.ops.split(geometry, splitter)

            # Convert coordinates in each part
            converted_parts = []
            for part in split_parts.geoms:
                part_coords = convert_coords(np.array(part.exterior.coords))
                # ugly hack to prevent nasty shapely bug
                if (part_coords[:, 0] == 0).any() and (part_coords[:, 0] > 180).any():
                    part_coords[:, 0][part_coords[:, 0] == 0] = 360
                converted_parts.append(Polygon(part_coords))
            
            return MultiPolygon(converted_parts)
        
        elif isinstance(geometry, LineString):
            return LineString(convert_coords(geometry.coords))
        
        elif isinstance(geometry, Point):
            coords = list(geometry.coords)[0]
            new_x = coords[0] + 360 if coords[0] < 0 else coords[0]
            return Point(new_x, coords[1])
        
        else:
            raise ValueError(f"Unsupported geometry type: {type(geometry)}")
    
    result_gdf['geometry'] = result_gdf['geometry'].apply(split_geometry)
    return result_gdf

def build_grid_polygons(latitudes, longitudes):
    """
    Builds a list of shapely Polygons representing grid cells.
    
    Args: 
        - latitudes (np.ndarray): 1D array of latitude values
        - longitudes (np.ndarray): 1D array of longitude values
    
    Returns:
        - grid_polygons (list[shapely.Polygon]): List of polygons representing the grid cells
        - grid_areas (list[float]): List of areas of the grid cells
    """
    grid_polygons = []
    grid_areas = []

    inter_lon_diff = longitudes[1] - longitudes[0]
    longitudes = longitudes - inter_lon_diff / 2
    longitudes = np.append(longitudes, longitudes[-1] + inter_lon_diff)
    inter_lat_diff = latitudes[1] - latitudes[0]
    latitudes = latitudes - inter_lat_diff / 2
    latitudes = np.append(latitudes, latitudes[-1] + inter_lat_diff)
    num_lat = len(latitudes)
    num_long = len(longitudes)
    for j in range(num_long-1):
        for i in range(0, num_lat-1):
            polygon = Polygon([(longitudes[j], latitudes[i]),
                               (longitudes[j], latitudes[i + 1]),
                               (longitudes[j+1], latitudes[i + 1]),
                               (longitudes[j+1], latitudes[i])])
            grid_polygons.append(polygon)
            grid_areas.append(calculate_grid_area(latitudes[i], latitudes[i + 1], longitudes[j], longitudes[j+1]))
    return grid_polygons, grid_areas

def get_intersection_areas(gdf: gpd.GeoDataFrame, grid_polygons: List[shapely.Polygon], grid_areas: List[float], 
                            name_col: str, progress_name: str, num_lats: int,num_lons: int,
                            include_admin: bool = False) -> Dict:
    """
    Calculate intersection areas between a GeoDataFrame and a set of grid polygons.

    Args:
        - gdf (gpd.GeoDataFrame): The GeoDataFrame to intersect with the grid.
        - grid_polygons (list[shapely.Polygon]): The polygons to intersect with the GeoDataFrame.
        - grid_areas (list[float]): The areas of the grid cells.
        - name_col (str): The column name of the GeoDataFrame to use for the name.
        - progress_name (str): The filename to display in the tqdm progress bar.
        - include_admin (bool, optional): Whether to include the admin country name. Default is False.

    Returns:
        - dict: A dictionary where the keys are the names and the values are numpy arrays of the areas of intersection with each grid cell.
    """
    intersection_areas = {}
    # Only cells whose bounding box meets the object's can intersect it, so an STR-tree picks those out and every other cell keeps an area of 0, which is what intersecting it would have given. This is the same computation as a loop over all 64800 cells, at a small fraction of the cost.
    cells = np.asarray(grid_polygons, dtype=object)
    tree = shapely.STRtree(cells)
    grid_areas = np.asarray(grid_areas, dtype="float64")
    for idx, obj in tqdm.tqdm(gdf.iterrows(), total=gdf.shape[0], desc=f"Building cache for {progress_name}"):
        obj_name = obj[name_col]
        if include_admin:
            country = obj['admin']
            if country:
                obj_name = f"{obj_name}, {country}"
            else:
                obj_name = f"{obj_name}"
                    
        obj_polygon = obj['geometry'].buffer(0)

        areas = np.zeros(len(cells))
        candidates = tree.query(obj_polygon)
        if candidates.size:
            overlap = shapely.area(shapely.intersection(cells[candidates], obj_polygon))
            areas[candidates] = overlap * grid_areas[candidates]
        intersection_areas[obj_name] = areas.reshape(num_lons, num_lats)
    return intersection_areas

def calculate_grid_area(lat_min, lat_max, lon_min, lon_max, radius=1):
    """
    Calculate the area of a rectangle on the Earth's surface given latitude and longitude bounds.
    
    Args:
        - lat_min, lat_max: Minimum and maximum latitude in degrees
        - lon_min, lon_max: Minimum and maximum longitude in degrees
        - radius: Radius of the Earth in kilometers 
    
    Returns:
        - Area in square kilometers
    """
    lat_min = math.radians(lat_min)
    lat_max = math.radians(lat_max)
    lon_min = math.radians(lon_min)
    lon_max = math.radians(lon_max)
    
    area = radius**2 * abs(math.sin(lat_max) - math.sin(lat_min)) * abs(lon_max - lon_min)
    
    return area

def get_nearest_grid_point(longitudes, latitudes, lon, lat):
    """
    Gets nearest grid point for any lon, lat

    Args:
        - latitudes (np.ndarray): The latitude values to use for the grid.
        - longitudes (np.ndarray): The longitude values to use for the grid.
        - lon (float): longitude
        - lat (float): latitude

    Returns:
        - (tuple): longitude_grid_index, latitude_grid_index
    """
     
    if lon < 0:
        lon += 360

    delta_lat = latitudes[1] - latitudes[0]
    delta_lon = longitudes[1] - longitudes[0]

    lon_ind = np.round((lon - longitudes[0])/delta_lon)
    lat_ind = np.round((lat - latitudes[0])/delta_lat)
    return int(lon_ind), int(lat_ind)
    
def load_and_preprocess_ne_df(natural_earth_path: str,
                              geofeature: str):
    """
    Load and preprocess a Natural Earth dataset.

    Args:
        - natural_earth_path (str): The path to the Natural Earth dataset.
        - geofeature (str): The geofeature to load, e.g. "country", "state", "continent", or "water_body".
    """
    name_col = get_name_col(geofeature)
    if geofeature == 'country':
        gdf = gpd.read_file(os.path.join(natural_earth_path, '10m_cultural', 'ne_10m_admin_0_countries.shp'))
        remove_list = EXCLUDE_COUNTRIES
    elif geofeature == 'state':
        gdf = gpd.read_file(os.path.join(natural_earth_path, '10m_cultural', 'ne_10m_admin_1_states_provinces.shp'))
        remove_list = EXCLUDE_STATES
        gdf['admin'] = gdf['admin'].replace(DIFF_STATE_ADMIN_TO_COUNTRY_ADMIN)
    elif geofeature == 'continent':
        gdf = gpd.read_file(os.path.join(natural_earth_path, '10m_cultural', 'ne_10m_admin_0_countries.shp'))
        gdf = gdf.dissolve(by='CONTINENT')
        gdf.reset_index(inplace=True)
        remove_list = EXCLUDE_CONTINENTS
    elif geofeature == 'water_body':
        gdf = gpd.read_file(os.path.join(natural_earth_path, '10m_physical', 'ne_10m_geography_marine_polys.shp'))
        remove_list = EXCLUDE_WATER_BODIES
    gdf = gdf[gdf[name_col].isnull() == False]
    gdf = gdf[~gdf[name_col].isin(remove_list)]
    return gdf

def build_caches(latitudes: np.ndarray, longitudes: np.ndarray, save_dir: str, natural_earth_path: str) -> None:

    """
    Build caches of intersection areas between grid cells and geospatial
    boundaries and saves it in the `save_dir` directory.

    Args:
        - latitudes (np.ndarray): The latitude values to use for the grid.
        - longitudes (np.ndarray): The longitude values to use for the grid.
        - save_dir (str): The directory where caches will be saved.
        - natural_earth_path (str): The path to the NaturalEarth dataset.
    """
    resolution = f"{len(longitudes)}x{len(latitudes)}"
    num_lats = len(latitudes)
    num_lons = len(longitudes)
    grid_polygons, grid_areas = build_grid_polygons(latitudes=latitudes,
                                                longitudes=longitudes)
    
    for geofeature in ['country', 'state', 'continent', 'water_body']:
        area_cache_path = os.path.join(save_dir, f"area_cache_{geofeature}_{resolution}.npz")
        if os.path.exists(area_cache_path):
            logging.info(
                f"Using existing cache for {geofeature} in {save_dir}")
            continue
        

        gdf = load_and_preprocess_ne_df(natural_earth_path, geofeature)
        name_col = get_name_col(geofeature)
        
        gdf = fix_longitudes_and_split(gdf)

        include_admin = (geofeature == 'state')

        intersection_areas = get_intersection_areas(gdf=gdf, 
                                                    grid_polygons=grid_polygons,
                                                    grid_areas=grid_areas, 
                                                    name_col=name_col, 
                                                    progress_name=geofeature,
                                                    num_lats=num_lats,
                                                    num_lons=num_lons, 
                                                    include_admin=include_admin)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        np.savez(os.path.join(save_dir, f"area_cache_{geofeature}_{resolution}.npz"), **intersection_areas)

    cm_file = os.path.join(save_dir, f'cost_matrix_{resolution}.npz')
    if os.path.exists(cm_file):
        logging.info(f"Using existing cache for {cm_file}")
    else:
        logging.info(f"Building cache for location evaluation at {cm_file}")
        all_cords = []
        for lon in longitudes:
            for lat in latitudes:
                if lon > 180:
                    all_cords.append((lat, lon - 360))
                else:
                    all_cords.append((lat, lon))
        
        cm = []
        for lat in tqdm.tqdm(latitudes):
            dists = haversine_vector(all_cords, [(lat, 0)], unit=Unit.KILOMETERS, comb=True)
            cm_grid = np.zeros((len(longitudes), len(latitudes)))
            for i, dist in enumerate(dists[0]):
                cm_grid[i//len(latitudes)][i%len(latitudes)] = dist
            cm.append(cm_grid)
        cm = np.array(cm)
        np.savez(cm_file, cost_matrix=cm)

      
        