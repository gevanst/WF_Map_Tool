import streamlit as st
import folium
from streamlit_folium import st_folium
from folium.plugins import Draw, MeasureControl, MousePosition
import geopandas as gpd
import math
from shapely.geometry import Polygon
import requests

st.set_page_config(layout="wide")

# Set up the Streamlit interface
st.title("Wind Farm Mapping Tool")
st.write("Tool for quickly determining possible wind farm locations")
st.write("Land and Ocean class/uses/hazards are added as map layers")
st.write("Source for most map layers: https://planning.hawaii.gov/gis/download-gis-data-expanded/")

# User inputs for turbine model
turbine_choice = st.selectbox("Choose Turbine Reference Model", ["5MW", "10MW", "15MW"])
if turbine_choice == "5MW":
    hub_ht = 90
    rtr_d = 125.9
elif turbine_choice == "10MW":
    hub_ht = 119
    rtr_d = 198    
else:
    hub_ht = 150
    rtr_d = 242

clearance = st.number_input("Choose Turbine Clearance In Rotor Diameters", 1, 10)
T_Spacing = round(clearance * rtr_d)
st.write(f"Spacing between turbines (m): {T_Spacing}")

# Create a base Folium map centered on the Hawaiian Islands
hawaii_map = folium.Map(location=[20.9, -156.5], zoom_start=9)

# Add layers to the map
folium.TileLayer(tiles='https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png', name='OpenTopoMap', attr='© OpenTopoMap contributors').add_to(hawaii_map)
arcgis_ocean_basemap = 'https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}'
folium.TileLayer(tiles=arcgis_ocean_basemap, name='ArcGIS Ocean Basemap', attr='Esri, DeLorme, GEBCO, NOAA NGDC, and other contributors').add_to(hawaii_map)
folium.raster_layers.WmsTileLayer(url='https://wms.gebco.net/mapserv?', name='GEBCO Bathymetry', layers='GEBCO_LATEST', fmt='image/png', transparent=True, attr='GEBCO').add_to(hawaii_map)

# Add measurement and mouse position tools
hawaii_map.add_child(MeasureControl(primary_length_unit='meters', secondary_length_unit='miles'))
formatter = "function(num) {return L.Util.formatNum(num, 5);};"
mouse_position = MousePosition(position="topright", separator=" | ", empty_string="No Coordinates", lng_first=True, num_digits=5, prefix="Coordinates:", lat_formatter=formatter, lng_formatter=formatter)
hawaii_map.add_child(mouse_position)

# Add drawing tools for rectangle selection
draw = Draw(draw_options={"polyline": True, "polygon": True, "rectangle": True, "circle": True}, edit_options={"edit": True})
draw.add_to(hawaii_map)

# Cache loading of shapefiles and layers to improve performance
@st.cache_data
def load_layers():
    eez_file_path = "data/eez.shp"
    benthic_file_path = "data/benthic_habitat.shp"
    ocean_depth_path = "data/Ocean_Depth.shp"

    eez_gdf = gpd.read_file(eez_file_path)
    benth_gdf = gpd.read_file(benthic_file_path).to_crs("EPSG:4326")
    ocean_depth_gdf = gpd.read_file(ocean_depth_path).to_crs("EPSG:4326")

    us_eez = eez_gdf[eez_gdf['geoname'].str.contains("United States Exclusive Economic Zone", na=False)]
    us_eez_json = us_eez.to_json()
    benthic_json = benth_gdf.to_json()
    ocean_depth_json = ocean_depth_gdf.to_json()

    return us_eez_json, benthic_json, ocean_depth_json

# Load cached layers
us_eez_json, benthic_json, ocean_depth_json = load_layers()

# Define a custom style function for benthic habitats
def benthic_style(feature):
    zone = feature['properties']['ZONE']
    color_map = {
        'Back Reef': 'blue',
        'Bank/Shelf': 'purple',
        'Channel': 'orange',
        'Fore Reef': 'red',
        'Lagoon': 'cyan',
        'Reef Crest': 'yellow',
        'Reef Flat': 'pink',
        'Land': 'gray',
        'Unknown': 'black'
    }
    return {'fillColor': color_map.get(zone, 'brown'), 'color': 'black', 'weight': 1, 'fillOpacity': 0.4}

# Define a style function for the ocean depth layer
def depth_style(feature):
    depth = feature['properties'].get('DEPTH')
    depth_colors = {
        1000: '#BFEFFF',
        2000: '#7EC0EE',
        3000: '#4682B4',
        4000: '#104E8B',
        5000: '#000080'
    }
    return {'fillColor': depth_colors.get(depth, '#708090'), 'color': 'black', 'weight': 1, 'fillOpacity': 0.5}

# Create FeatureGroups for each layer
eez_layer = folium.FeatureGroup(name="US Econ. Excl. Zone")
benth_layer = folium.FeatureGroup(name="Benthic Habitats")
ocean_depth_layer = folium.FeatureGroup(name="Ocean Depths")

# Add GeoJSON layers to the respective FeatureGroups
folium.GeoJson(us_eez_json, style_function=lambda feature: {'color': 'green', 'weight': 2, 'fillOpacity': 0.1}).add_to(eez_layer)
folium.GeoJson(benthic_json, name="Benthic Habitat Zones", style_function=benthic_style).add_to(benth_layer)
folium.GeoJson(ocean_depth_json, name="Ocean Depth", style_function=depth_style).add_to(ocean_depth_layer)

# Add the FeatureGroups to the map
eez_layer.add_to(hawaii_map)
benth_layer.add_to(hawaii_map)
ocean_depth_layer.add_to(hawaii_map)

# Add LayerControl to enable toggling
folium.LayerControl().add_to(hawaii_map)

# Display the map and capture the user's drawn shapes
st_data = st_folium(hawaii_map, width=800, height=500)

#create info container
info_container = st.container()

# Caching layout function
@st.cache_data
def layout_turbines(bounds, spacing, single_turbine=False):
    min_lon, min_lat = bounds[0]
    max_lon, max_lat = bounds[2]

    # Calculate turbine positions based on spacing
    if single_turbine:
        return [((min_lon + max_lon) / 2, (min_lat + max_lat) / 2)]

    lon_spacing = spacing / 111000
    lat_spacing = spacing / 110000
    num_turbines_x = math.floor((max_lon - min_lon) / lon_spacing)
    num_turbines_y = math.floor((max_lat - min_lat) / lat_spacing)

    turbine_positions = []
    for i in range(num_turbines_x):
        for j in range(num_turbines_y):
            turbine_lon = min_lon + math.floor(i * lon_spacing * 1e5) / 1e5
            turbine_lat = min_lat + math.floor(j * lat_spacing * 1e5) / 1e5
            turbine_positions.append((turbine_lon, turbine_lat))
    
    return turbine_positions

# Check for drawn shapes
with info_container:
    if st_data and "all_drawings" in st_data and st_data["all_drawings"]:
        shape_data = st_data["all_drawings"][0]

        if shape_data["geometry"]["type"] == "Polygon":
            st.write("### Farm Space Identified")

            # Get the rectangle boundaries
            rectangle_coords = shape_data["geometry"]["coordinates"][0]
            lon_min, lat_min = rectangle_coords[0]
            lon_max, lat_max = rectangle_coords[2]

            # Calculate the center of the rectangle
            center_lat = (lat_min + lat_max) / 2
            center_lon = (lon_min + lon_max) / 2

            # Show boundaries and center coordinates
            st.write(f"**Rectangle Boundaries - SW**: ({lat_min:.4f}, {lon_min:.4f})")
            st.write(f"**Rectangle Boundaries - NE**: ({lat_max:.4f}, {lon_max:.4f})")
            st.write(f"**Center of the Rectangle**: ({center_lat:.4f}, {center_lon:.4f})")

            # Calculate turbine positions
            layout_type = st.radio("Choose Turbine Layout", ["Single Turbine", "Farm Array"])
            turbine_positions = layout_turbines(rectangle_coords, T_Spacing, single_turbine=(layout_type == "Single Turbine"))
            st.write(f"Number of Turbines: {len(turbine_positions)}")

            # Add turbines to map
            for lon, lat in turbine_positions:
                folium.Marker(location=[lat, lon], popup=f"Turbine Location\nLat: {lat}\nLon: {lon}", icon=folium.Icon(color='green', icon='wind')).add_to(hawaii_map)

            # Redisplay the updated map with turbine markers
            st_folium(hawaii_map, width=800, height=500)
