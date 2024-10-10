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
            }.get(feature['properties'].get('ZONE', 'Unknown'), 'brown'),
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.4
        },
        'tooltip_fields': ["M_STRUCT", "D_STRUCT", "M_COVER", "ZONE"],
        'aliases': ["Major Structure", "Detailed Structure", "Major Cover", "Zone"],
        'name': 'Benthic Habitats',
        'simplify_tolerance': 0.001  # Adjust as needed
    },
    'Coral Reefs (Nautical Charts)': {
        'file_path': 'data/CoralReefs_nc.shp',
        'style_function': lambda feature: {
            'fillColor': 'teal',
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6
        },
        'tooltip_fields': ["acres", "chart1", "chart2", "notes"],
        'aliases': ["Area (Acres)", "Source Chart 1", "Source Chart 2", "Notes"],
        'name': 'Coral Reefs (Nautical Charts)',
        'simplify_tolerance': 0.001
    },
    'Whale Sanctuary': {
        'file_path': 'data/sanctuary.shp',
        'style_function': lambda feature: {
            'fillColor': 'lightblue',
            'color': 'blue',
            'weight': 1,
            'fillOpacity': 0.5
        },
        'tooltip_fields': ["sanctuary", "Shape_Leng", "Shape_Area"],
        'aliases': ["Sanctuary Name", "Perimeter (m)", "Area (sq. m)"],
        'name': 'Whale Sanctuary',
        'simplify_tolerance': 0.001
    },
    'Ocean Depth': {
        'file_path': 'data/Ocean_Depth.shp',
        'style_function': lambda feature: {
            'fillColor': {
                1000: '#BFEFFF',
                2000: '#7EC0EE',
                3000: '#4682B4',
                4000: '#104E8B',
                5000: '#000080'
            }.get(feature['properties'].get('DEPTH', 5000), '#708090'),
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.5
        },
        'tooltip_fields': ["DEPTH"],
        'aliases': ["Depth (m)"],
        'name': 'Ocean Depth',
        'simplify_tolerance': 0.001
    },
    'Dumping Areas': {
        'file_path': 'data/DumpingAreas.shp',
        'style_function': lambda feature: {
            'fillColor': 'brown',
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6
        },
        'tooltip_fields': ["id", "chart", "cfr", "notes"],
        'aliases': ["Polygon ID", "Chart Source", "Regulation", "Notes"],
        'name': 'Dumping Areas',
        'simplify_tolerance': 0.001
    },
    'Moku Ridge to Reef': {
        'file_path': 'data/moku_ridge_to_reef_dar.shp',
        'style_function': lambda feature: {
            'fillColor': {
                'Hilo': 'orange',
                'Kona': 'red',
                'Puna': 'yellow',
                'Kohala': 'green',
                'Kau': 'blue',
                'Molokai': 'purple',
                'Oahu': 'brown',
                'Maui': 'cyan'
            }.get(feature['properties'].get('Moku'), 'gray'),
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6
        },
        'tooltip_fields': ["Moku_ID", "Moku", "Mokupuni"],
        'aliases': ["Moku ID", "Moku", "Island"],
        'name': 'Moku Ridge to Reef',
        'simplify_tolerance': 0.001
    },
    'Submarine Cables': {
        'file_path': 'data/Cables.shp',
        'style_function': lambda feature: {
            'color': 'blue',
            'weight': 3,
            'opacity': 0.7,
            'dashArray': '5, 5' if 'Descript' in feature['properties'] and 'Proposed' in feature['properties']['Descript'] else None
        },
        'tooltip_fields': ["id", "length", "descriptio", "chart1", "chart2"],
        'aliases': ["Cable ID", "Length (m)", "Description", "Primary Chart", "Secondary Chart"],
        'name': 'Submarine Cables',
        'simplify_tolerance': 0.001
    },
    'Unexploded Ordnance': {
        'file_path': 'data/UnexplodedOrdnance.shp',
        'marker_color': 'red',
        'tooltip_fields': ["id", "notes", "chart"],
        'aliases': ["ID", "Notes", "Chart"],
        'name': 'Unexploded Ordnance',
    },
    'Explosive Dumping Areas': {
        'file_path': 'data/ExplosiveDumpingAreas.shp',
        'style_function': lambda feature: {
            'fillColor': 'orange',
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6
        },
        'tooltip_fields': ["chart", "notes", "acres"],
        'aliases': ["Source Chart", "Notes", "Area (Acres)"],
        'name': 'Explosive Dumping Areas',
        'simplify_tolerance': 0.001
    },
    'Sewer Lines': {
        'file_path': 'data/sewerlines.shp',
        'style_function': lambda feature: {
            'color': 'brown',
            'weight': 3,
            'opacity': 0.7,
        },
        'tooltip_fields': ["id", "notes", "island", "chart"],
        'aliases': ["Sewer ID", "Notes", "Island", "Chart"],
        'name': 'Sewer Lines',
        'simplify_tolerance': 0.001
    },
    'Submerged Buoys': {
        'file_path': 'data/SubmergedBuoys.shp',
        'marker_color': 'purple',
        'tooltip_fields': ["name", "notes"],
        'aliases': ["ID", "Description"],
        'name': 'Submerged Buoys',
    },
    'Navigation Aids': {
        'file_path': 'data/NavigationAids.shp',
        'marker_color': 'darkgreen',
        'tooltip_fields': ["id", "type", "name"],
        'aliases': ["ID", "Type", "Name"],
        'name': 'Navigation Aids',
    },
    'Obstructions': {
        'file_path': 'data/obstructions.shp',
        'marker_color': 'orange',
        'tooltip_fields': ["id", "notes"],
        'aliases': ["ID", "Notes"],
        'name': 'Obstructions',
    },
    'Fish Aggregation Devices (FADs)': {
        'file_path': 'data/fads_nmfs.shp',
        'marker_color': 'blue',
        'tooltip_fields': ["name", "island"],
        'aliases': ["ID", "Description"],
        'name': 'Fish Aggregation Devices (FADs)',
    },
    'City Zoning': {
        'file_path': 'data/cty_zoning_mau.shp',
        'style_function': lambda feature: {
            'fillColor': {
                'R-3 Residential': 'lightblue',
                'AP Airport': 'orange',
                'AG Agricultural': 'green',
                'INT Interim': 'brown',
                'PK Park': 'red',
                'M-2 Heavy Industrial': 'gray',
                'H-2 Hotel': 'purple',
                'M-1 Light Industrial': 'yellow',
                'B-2 Business-Community': 'cyan'
            }.get(feature['properties'].get('ZONE_CLASS'), 'black'),
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6,
        },
        'tooltip_fields': ['ZONE_CLASS', 'CP_AREA'],
        'aliases': ['Zoning Class:', 'Community Plan Area:'],
        'name': 'City Zoning',
        'simplify_tolerance': 0.001
    },
    'False Killer Whale Habitat': {
        'file_path': 'data/WhaleFalseKiller_MainHawaiianIslandsInsularDPS_20180724.shp',
        'style_function': lambda feature: {
            'fillColor': {
                'marine': 'blue',
                'nearshore': 'cyan',
                'offshore': 'purple',
                'unknown': 'gray'
            }.get(feature['properties'].get('HABTYPE'), 'lightgray'),
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.5,
        },
        'tooltip_fields': ['SCIENAME', 'COMNAME', 'HABTYPE', 'AREASqKm'],
        'aliases': ['Species Name', 'Common Name', 'Habitat Type', 'Area (sq. km)'],
        'name': 'False Killer Whale Habitat',
        'simplify_tolerance': 0.001
    },
    'Monk Seal Habitat': {
        'url': "https://services2.arcgis.com/C8EMgrsFcRFL6LrL/arcgis/rest/services/SealHawaiianMonk_20150821_line/FeatureServer/59/query",
        'params': {
            "where": "1=1",
            "outFields": "*",
            "f": "geojson"
        },
        'style_function': lambda feature: {
            'color': 'red',
            'weight': 3,
            'opacity': 0.7,
        },
        'tooltip_fields': ['ID', 'COMNAME', 'HABTYPE', 'AREASqKm'],
        'aliases': ['ID', 'Common Name', 'Habitat Type', 'Area (sq. km)'],
        'name': 'Monk Seal Habitat',
    },
}

selected_layers = st.sidebar.multiselect('Select layers to display:', list(layer_info.keys()), default=list(layer_info.keys()))

# Function to load and process layers
@st.cache_data
def load_layer(layer_key):
    info = layer_info[layer_key]
    if 'file_path' in info:
        gdf = gpd.read_file(info['file_path'])
        gdf = gdf.to_crs("EPSG:4326")
        if 'simplify_tolerance' in info:
            gdf['geometry'] = gdf['geometry'].simplify(info['simplify_tolerance'])
        if gdf.geometry.is_empty.any():
            gdf = gdf[~gdf.geometry.is_empty]
        geojson_data = gdf.to_json(drop_id=True, double_precision=3)
        return geojson_data, None
    elif 'url' in info:
        response = requests.get(info['url'], params=info['params'])
        geojson_data = response.json()
        return geojson_data, None
    else:
        return None, None

# Function to add point layers
def add_point_layer(gdf, layer_name, marker_color, tooltip_fields, aliases):
    feature_group = folium.FeatureGroup(name=layer_name)
    for _, row in gdf.iterrows():
        if row.geometry and row.geometry.is_valid:
            coords = [row.geometry.y, row.geometry.x]
            tooltip_content = "<br>".join(
                f"{alias}: {row[field]}" for field, alias in zip(tooltip_fields, aliases)
            )
            folium.CircleMarker(
                location=coords,
                radius=6,
                color=marker_color,
                fill=True,
                fill_color=marker_color,
                fill_opacity=0.7,
                popup=folium.Popup(tooltip_content, max_width=250)
            ).add_to(feature_group)
    feature_group.add_to(hawaii_map)

# Add selected layers to the map
for layer_name in selected_layers:
    info = layer_info[layer_name]
    geojson_data, _ = load_layer(layer_name)
    if geojson_data:
        feature_group = folium.FeatureGroup(name=layer_name)
        if 'marker_color' in info:
            # Handle point layers
            gdf = gpd.read_file(info['file_path'])
            gdf = gdf.to_crs("EPSG:4326")
            add_point_layer(
                gdf,
                layer_name,
                info['marker_color'],
                info['tooltip_fields'],
                info['aliases']
            )
        else:
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
    num_turbines_x = max(1, math.floor((max_lon - min_lon) / lon_spacing))
    num_turbines_y = max(1, math.floor((max_lat - min_lat) / lat_spacing))

    turbine_positions = []
    for i in range(num_turbines_x):
        for j in range(num_turbines_y):
            # Place turbines based on grid index and spacing
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
            # Calculate the dimensions (length and width) of the rectangle in meters
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
