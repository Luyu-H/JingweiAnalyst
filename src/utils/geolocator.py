import os
import numpy as np
import geopandas as gpd
from shapely.geometry import Point

from src.utils.geo_processor import build_caches, get_nearest_grid_point, get_name_col, load_and_preprocess_ne_df

class Geolocator:
    def __init__(self, latitudes, longitudes, cache_dir, natural_earth_path):
        """
        Initialize a Geolocator instance. Handles both forward and reverse geocoding

        Args:
            - latitudes (np.ndarray): 1D array of latitude values
            - longitudes (np.ndarray): 1D array of longitude values
            - cache_dir (str): Path to the directory where the caches are stored.
            - natural_earth_path (str): Path to the directory containing the Natural Earth dataset.
        """
        build_caches(latitudes, longitudes, cache_dir, natural_earth_path)

        self.country_df = gpd.read_file(os.path.join(natural_earth_path, '10m_cultural', 'ne_10m_admin_0_countries.shp'))

        # exclude the countries in the list
        self.country_df = load_and_preprocess_ne_df(natural_earth_path, 'country')
        self.state_df = load_and_preprocess_ne_df(natural_earth_path, 'state')
        self.continent_df = load_and_preprocess_ne_df(natural_earth_path, 'continent')
        self.water_body_df = load_and_preprocess_ne_df(natural_earth_path, 'water_body')

        self.countries_area_cache, self.continents_area_cache, self.states_area_cache, self.marine_area_cache, self.cost_matrix = self.read_caches(cache_dir, f"{len(longitudes)}x{len(latitudes)}")

        self.latitudes = latitudes
        self.longitudes = longitudes

    def _get_mask(self, area_cache):
        """
        Return a boolean mask indicating which grid cells intersect the given area_cache.

        Args:
            - area_cache (np.ndarray): The area cache to generate the mask for

        Returns:
            - np.ndarray: A boolean mask indicating which grid cells intersect the given area_cache
        """
        return (area_cache > 0)

    def _find_by_point(self, df, longitude, latitude, geofeature):
        """
        Find a geofeature by a given point.

        Args: 
            - df (gpd.GeoDataFrame): The geoframe to search
            - longitude (float): The longitude of the point
            - latitude (float): The latitude of the point

        Returns:
            - str: The name of the geofeature if found, otherwise None
        """
        point = Point(longitude, latitude)
        name_col = get_name_col(geofeature)
        found = df[df.contains(point)]
        return found[name_col].iloc[0] if not found.empty else None

    def read_caches(self, cache_dir, resolution):
        """
        Reads all the caches from disk and returns them as a tuple.

        Args:
            - cache_dir (str): The directory where the caches are stored.
            - resolution (str): The resolution of the caches to read, e.g. "100x100".

        Returns:
            - tuple: A tuple of 8 dictionaries, containing the area and mask caches of the countries, the continents, the states and provinces, and the water bodies.
        """
        countries_area_cache = np.load(os.path.join(
            cache_dir, f'area_cache_country_{resolution}.npz'))
        states_area_cache = np.load(os.path.join(
            cache_dir, f'area_cache_state_{resolution}.npz'))
        continents_area_cache = np.load(os.path.join(
            cache_dir, f'area_cache_continent_{resolution}.npz'))
        marine_area_cache = np.load(os.path.join(
            cache_dir, f'area_cache_water_body_{resolution}.npz'))
        cost_matrix_cache = np.load(os.path.join(cache_dir, f'cost_matrix_{resolution}.npz'))
        return countries_area_cache, continents_area_cache, states_area_cache, marine_area_cache, cost_matrix_cache['cost_matrix']

    def get_geofeature_by_point(self, longitude, latitude, geofeature):
        """
        Find the geofeature by a given point.

        Args:
            - longitude (float): The longitude of the point
            - latitude (float): The latitude of the point
            - geofeature (str): The type of geofeature to search for. Must be one of 'country', 'state', 'continent', or 'water_body'.

        Returns:
            - str: The name of the geofeature if found, otherwise None.
        
        Example:
        >>> # Assuming gl is an instance of Geolocator with proper datasets loaded
        >>> gl.get_geofeature_by_point(-74.0060, 40.7128, 'country')
        'United States of America'
        """
        df_map = {
            'country': self.country_df,
            'state': self.state_df,
            'continent': self.continent_df,
            'water_body': self.water_body_df
        }
        return self._find_by_point(df_map[geofeature], longitude, latitude, geofeature)

    def get_area_map(self, geofeature, name):
        """
        Retrieve the area map (ndarray) for a given geofeature.

        Args:
            - geofeature (str): The type of geofeature for which to retrieve the area map.
            Must be one of 'country', 'continent', 'state', or 'water_body'.
            - name (str): The name of the geofeature for which to retrieve the area map.

        Returns:
            - ndarray: 2D array with area-weighted grid cells for the given geofeature.
        
        Example:
        >>> # Assuming gl is an instance of Geolocator
        >>> area_map = gl.get_area_map('state', 'Texas, United States of America')
        """
        cache_map = {
            'country': self.countries_area_cache,
            'continent': self.continents_area_cache,
            'state': self.states_area_cache,
            'water_body': self.marine_area_cache
        }
        if geofeature not in cache_map:
            raise ValueError(f"Invalid geofeature '{geofeature}'. Supported geofeatures are: {', '.join(cache_map.keys())}.")
        
        keys_as_tuples = [(key.rsplit(',', 1)[0], key) for key in cache_map[geofeature].keys()]
        query = [k[1] for k in keys_as_tuples if k[0] == name or k[1] == name]
        if not query:
            raise ValueError(f"{name} not found in the cache for geofeature '{geofeature}'. Available keys: {list(cache_map[geofeature].keys())}")

        return cache_map[geofeature][query[0]]

    def get_mask_by_geofeature(self, geofeature, name):
        """
        Retrieve a mask for a specific geofeature and name.

        Args:
            - geofeature (str): The type of geofeature for which to retrieve the mask. 
            Must be one of 'country', 'continent', 'state', or 'water_body'.
            - name (str): The specific name of the geofeature.

        Returns:
            - np.ndarray: 2D boolean array indicating which grid cells in the area are part of the geofeature.

        Example:
        >>> # Assuming gl is an instance of Geolocator and the cache contains 'United States of America'
        >>> mask = gl.get_mask_by_geofeature('country', 'United States of America')
        """
        return self._get_mask(self.get_area_map(geofeature, name))

    def get_all_locations(self, geofeature):

        """
        Retrieve a list of all geofeatures of a given type.

        Args:
            - geofeature (str): The type of geofeature for which to retrieve the list. 
            Must be one of 'country', 'continent', 'state', or 'water_body'.

        Returns:
            - list: A list of all geofeatures of the given type.
        
        Example:
        >>> # Assuming gl is an instance of Geolocator
        >>> countries = gl.get_all_locations('country')
        >>> 'United States of America' in countries
        True    
        """
        
        cache_map = {
            'country': self.countries_area_cache,
            'continent': self.continents_area_cache,
            'state': self.states_area_cache,
            'water_body': self.marine_area_cache
        }

        return list(cache_map[geofeature].keys())

    def get_sub_geofeature(self, geofeature):
        """
        Determine the sub-geofeature for a given geofeature.

        Args: 
            - geofeature (str): The parent geofeature type (e.g., 'continent', 'country', 'water_body').

        Returns:
            - str: The sub-geofeature type (e.g., 'country' for 'continent', 'state' for 'country').
        
        Example:
        >>> # For a continent, the sub-geofeature is a country.
        >>> gl.get_sub_geofeature('continent')
        'country'
        """
        sub_geofeature_map = {
            'continent': 'country',
            'country': 'state',
            'water_body': 'water_body'
        }

        if geofeature not in sub_geofeature_map:
            raise ValueError(f"Invalid geofeature '{geofeature}'. Supported geofeatures are: {', '.join(sub_geofeature_map.keys())}.")

        return sub_geofeature_map[geofeature]

    def get_sublocations(self, geofeature, name):
        """
        Retrieve all sublocations within a given geofeature.

        Args:
            - geofeature (str): The type of geofeature for which to retrieve sublocations.
            Must be one of 'country', 'continent', or 'water_body'.
            - name (str): The specific name of the geofeature for which to retrieve sublocations.

        Returns:
            - list: A list of sublocations within the given geofeature.
            For example, if `geofeature` is "continent" and `name` is "Asia", this will return all countries in Asia.
        """
        if geofeature not in {'country', 'continent', 'water_body'}:
            raise ValueError("Invalid geofeature. Supported geofeatures are 'country', 'continent', and 'water_body'.")

        if geofeature == 'continent':
            continent_name_col = get_name_col('continent')
            if name not in self.continent_df[continent_name_col].unique():
                raise ValueError(f"Continent {name} not found.")
            sub_df = self.country_df[self.country_df[continent_name_col] == name]
            country_name_col = get_name_col('country')
            return list(sub_df[country_name_col].dropna())

        elif geofeature == 'country':
            country_name_col = get_name_col('country')
            if name not in self.country_df[country_name_col].unique():
                raise ValueError(f"Country {name} not found.")
            
            sub_df = self.state_df[self.state_df['admin'] == name]
            if sub_df.empty:
                return None
            state_name_col = get_name_col('state')
            return [f'{state}, {name}' for state in list(sub_df[state_name_col].dropna())]

        elif geofeature == 'water_body':
            valid_water_bodies = {
                'ocean', 'sea', 'reef', 'bay', 'gulf', 'strait', 'channel',
                'sound', 'fjord', 'lagoon', 'inlet', 'river'
            }
            if name.lower() not in valid_water_bodies:
                raise ValueError(f"Invalid water_body. Supported water bodies are: {', '.join(valid_water_bodies)}.")
            filtered_df = self.water_body_df[self.water_body_df['featurecla'].str.lower() == name.lower()]
            if filtered_df.empty:
                raise ValueError(f"No water bodies of type {name} found.")
            water_body_name_col = get_name_col('water_body')
            return list(filtered_df[water_body_name_col].dropna())

        return []

    def get_geodistance_for_points(self, point1, point2):
        """
        Get the shortest distance between 2 points on the globe.
        Uses precomputed cache for grid points, so do not use this if you care too much for accuracy

        Args: 
            - point1 (tuple): (longitude1, latitude1)
            - point2 (tuple): (longitude2, latitude2)

        Returns:
            - float: geodistance in KM
        """
        lon1, lat1 = point1
        lon2, lat2 = point2
        lon1, lat1 = get_nearest_grid_point(self.longitudes, self.latitudes, lon1, lat1)
        lon2, lat2 = get_nearest_grid_point(self.longitudes, self.latitudes, lon2, lat2)

        delta_lon = (lon1 - lon2)%len(self.longitudes)
        return self.cost_matrix[lat1, delta_lon, lat2]

