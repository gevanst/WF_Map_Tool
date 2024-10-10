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
    st.write("Rotor Diameter = 125.9m")
    st.write("Hub Height = 90m")
    st.write("cosine-loss operation model")
    st.write("Note: generator efficiency of 94.4% is assumed for the NREL 5MW turbine.")
    st.write("https://github.com/NREL/turbine-models/blob/master/Offshore/NREL_5MW_126_RWT_corrected.csv")
    hub_ht = 90
    rtr_d = 125.9
elif turbine_choice == "10MW":
    st.write("Rotor Diameter = 198m")
    st.write("Hub Height = 119m")
    st.write("cosine-loss operation model")
    st.write("Note: Generator efficiency of 94% used. Small power variations above rated removed.")
    st.write("https://github.com/NREL/turbine-models/blob/master/Offshore/IEA_10MW_198_RWT.csv")
    hub_ht = 119
    rtr_d = 198    
else:
    st.write("Rotor Diameter = 242m")
    st.write("Hub Height = 150m")
    st.write("cosine-loss operation model")
    st.write("Note: Generator efficiency of 100% assumed")
    st.write("https://github.com/IEAWindTask37/IEA-15-240-RWT/blob/master/Documentation/")
    hub_ht = 150
    rtr_d = 242

clearance = st.number_input("Choose Turbine Clearance In Rotor Diameters", 1, 10)
T_Spacing = round(clearance * rtr_d)
st.write(f"Spacing between turbines (m): {T_Spacing}")

# Create a base Folium map centered on the Hawaiian Islands
hawaii_map = folium.Map(location=[20.9, -156.5], zoom_start=9)

###############################################################################
# Add OpenTopoMap Tile Layer
folium.TileLayer(
    tiles='https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png',
    name='OpenTopoMap',
    attr='© OpenTopoMap contributors',
).add_to(hawaii_map)

# Add ArcGIS Ocean Basemap Layer
arcgis_ocean_basemap = 'https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}'
folium.TileLayer(
    tiles=arcgis_ocean_basemap,
    name='ArcGIS Ocean Basemap',
    attr='Esri, DeLorme, GEBCO, NOAA NGDC, and other contributors',
).add_to(hawaii_map)

# Add GEBCO Bathymetry Layer
folium.raster_layers.WmsTileLayer(
    url='https://wms.gebco.net/mapserv?',
    name='GEBCO Bathymetry',
    layers='GEBCO_LATEST',
    fmt='image/png',
    transparent=True,
    attr='GEBCO',
).add_to(hawaii_map)

#############################################################################

# Add measurement tool to the map
hawaii_map.add_child(MeasureControl(primary_length_unit='meters', secondary_length_unit='miles'))

# Add MousePosition plugin to show coordinates
formatter = "function(num) {return L.Util.formatNum(num, 5);};"
mouse_position = MousePosition(
    position="topright",
    separator=" | ",
    empty_string="No Coordinates",
    lng_first=True,
    num_digits=5,
    prefix="Coordinates:",
    lat_formatter=formatter,
    lng_formatter=formatter,
)
hawaii_map.add_child(mouse_position)

# Add drawing tools for rectangle selection
draw = Draw(draw_options={"polyline": True, "polygon": True, "rectangle": True, "circle": True}, edit_options={"edit": True})
draw.add_to(hawaii_map)

# Sidebar for layer selection
st.sidebar.title("Map Layers")

# Define layer information
layer_info = {
    'US Econ. Excl. Zone': {
        'file_path': 'data/eez.shp',
        'style_function': lambda feature: {'color': 'green', 'weight': 2, 'fillOpacity': 0.1},
        'tooltip_fields': ["geoname"],
        'aliases': ["Geoname"],
        'name': 'US Econ. Excl. Zone',
    },
    'Benthic Habitats': {
        'file_path': 'data/benthic_habitat.shp',
        'style_function': lambda feature: {
            'fillColor': {
                'Back Reef': 'blue',
                'Bank/Shelf': 'purple',
                'Channel': 'orange',
                'Fore Reef': 'red',
                'Lagoon': 'cyan',
                'Reef Crest': 'yellow',
                'Reef Flat': 'pink',
                'Land': 'gray',
                'Unknown': 'black'
            }.get(feature['properties']['ZONE'], 'brown'),
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.4
        },
        'tooltip_fields': ["M_STRUCT", "D_STRUCT", "M_COVER", "ZONE"],
        'aliases': ["Major Structure", "Detailed Structure", "Major Cover", "Zone"],
        'name': 'Benthic Habitats',
        'simplify_tolerance': 0.001  # Adjust as needed
    },
    # Add other layers similarly...
}

selected_layers = st.sidebar.multiselect('Select layers to display:', list(layer_info.keys()))

# Function to load and process layers
@st.cache_data
def load_layer(layer_key):
    info = layer_info[layer_key]
    gdf = gpd.read_file(info['file_path'])
    gdf = gdf.to_crs("EPSG:4326")
    if 'simplify_tolerance' in info:
        gdf['geometry'] = gdf['geometry'].simplify(info['simplify_tolerance'])
    geojson_data = gdf.to_json(drop_id=True, double_precision=3)
    return geojson_data

# Add selected layers to the map
for layer_name in selected_layers:
    info = layer_info[layer_name]
    geojson_data = load_layer(layer_name)
    feature_group = folium.FeatureGroup(name=layer_name)
    folium.GeoJson(
        geojson_data,
        name=layer_name,
        style_function=info['style_function'],
        tooltip=folium.GeoJsonTooltip(
            fields=info['tooltip_fields'],
            aliases=info['aliases'],
            localize=True
        )
    ).add_to(feature_group)
    feature_group.add_to(hawaii_map)

# Add LayerControl to enable toggling
folium.LayerControl().add_to(hawaii_map)

# Display the map and capture the user's drawn shapes
st_data = st_folium(hawaii_map, width=800, height=500)

# Create info container
info_container = st.container()

# Function to layout turbines
@st.cache_data
def layout_turbines(bounds, spacing, single_turbine=False):
    """
    Calculate turbine positions inside the selected rectangle, rounding down to ensure spacing.
    """
    min_lon, min_lat = bounds[0]
    max_lon, max_lat = bounds[2]

    # Calculate the number of turbines along each axis
    if single_turbine:
        return [((min_lon + max_lon) / 2, (min_lat + max_lat) / 2)]

    # Convert spacing from meters to approximate degrees (1 degree ≈ 111 km)
    lon_spacing = spacing / 111000
    lat_spacing = spacing / 110000

    # Calculate the number of turbines that fit, rounding down
    num_turbines_x = math.floor((max_lon - min_lon) / lon_spacing)
    num_turbines_y = math.floor((max_lat - min_lat) / lat_spacing)

    turbine_positions = []
    for i in range(num_turbines_x):
        for j in range(num_turbines_y):
            # Place turbines based on grid index and spacing, rounding down
            turbine_lon = min_lon + i * lon_spacing
            turbine_lat = min_lat + j * lat_spacing
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

            # Convert the drawn coordinates to a Polygon for area calc
            drawn_polygon = Polygon([(coord[0], coord[1]) for coord in rectangle_coords])
            # Create a GeoDataFrame and project it to a UTM zone to get the area in square meters
            gdf = gpd.GeoDataFrame(index=[0], crs="EPSG:4326", geometry=[drawn_polygon])
            gdf = gdf.to_crs(epsg=32604)  # Convert to UTM Zone 4 for accurate area calculation
            # Calculate area in square meters and convert to square kilometers
            area_sq_km = gdf.geometry.area.iloc[0] / 1e6
            # Calculate the dimensions (length and width) of the rectangle in kilometers
            width_m = gdf.geometry.iloc[0].bounds[2] - gdf.geometry.iloc[0].bounds[0]
            height_m = gdf.geometry.iloc[0].bounds[3] - gdf.geometry.iloc[0].bounds[1]

            st.write(f"**Area of Drawn Shape**: {area_sq_km:.2f} square kilometers")
            st.write(f"**Dimensions**: {width_m:.2f} m (Width) x {height_m:.2f} m (Height)")

            # Option to select single turbine or array
            layout_type = st.radio("Choose Turbine Layout", ["Single Turbine", "Farm Array"])

            # Calculate turbine positions
            turbine_positions = layout_turbines(rectangle_coords, T_Spacing, single_turbine=(layout_type == "Single Turbine"))

            st.write(f"Number of Turbines: {len(turbine_positions)}")

            # Add turbines to map
            for lon, lat in turbine_positions:
                folium.Marker(
                    location=[lat, lon],
                    popup=f"Turbine Location\nLat: {lat}\nLon: {lon}",
                    icon=folium.Icon(color='green', icon='wind')
                ).add_to(hawaii_map)

            # Redisplay the updated map with turbine markers
            st_folium(hawaii_map, width=800, height=500)
