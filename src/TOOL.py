import streamlit as st
import folium
from streamlit_folium import st_folium
from folium.plugins import Draw, MeasureControl, MousePosition
import geopandas as gpd
import math
from shapely.geometry import Polygon
import floris
import requests

st.set_page_config(layout="wide")

# Set up the Streamlit interface
st.title("Wind Farm Mapping Tool")
st.write("Tool for quickly determining possible wind farm locations")
st.write("Land and Ocean class/uses/hazards are added as map layers")
st.write("Source for most map layers: https://planning.hawaii.gov/gis/download-gis-data-expanded/")

# User inputs for turbine model, 
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
T_Spacing = round(clearance*rtr_d)
st.write(f"spacing between turbines (m):{T_Spacing}")

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

# Read the EEZ shapefile using GeoPandas
eez_file_path = "data/eez.shp"
benthic_file_path = "data/benthic_habitat.shp"
noaa_coral_path = "data/CoralReefs_nc.shp"
ocean_depth_path = "data/Ocean_Depth.shp"
dumping_area_path = "data/DumpingAreas.shp"
moku_file_path = "data/moku_ridge_to_reef_dar.shp"  # Moku Ridge to Reef shapefile path
bottom_type_file_path = "data/BottomType.shp"  # Bottom Type point shapefile
cable_file_path = "data/Cables.shp"  # Cable line shapefile path
ordnance_file_path = "data/UnexplodedOrdnance.shp" #point type
explosive_dumping_file_path = "data/ExplosiveDumpingAreas.shp"
offshore_sewer_path = "data/sewerlines.shp" #line type
whale_path = "data/sanctuary.shp"
buoys_sub_path = "data/SubmergedBuoys.shp" # point type
nav_aid_path = "data/NavigationAids.shp" # point type
obstruct_path = "data/obstructions.shp" # point type
fads_path = "data/fads_nmfs.shp" # point type fish aggregate devices
zoning_path = "data/cty_zoning_mau.shp" #land zoning types polygon
whale_habitat_file_path = "data/WhaleFalseKiller_MainHawaiianIslandsInsularDPS_20180724.shp"


try:
        # Construct the GeoJSON URL
    arcgis_url = "https://services2.arcgis.com/C8EMgrsFcRFL6LrL/arcgis/rest/services/SealHawaiianMonk_20150821_line/FeatureServer/59/query"
    params = {
        "where": "1=1",         # Retrieve all data
        "outFields": "*",       # Get all fields
        "f": "geojson"          # Return format as GeoJSON
    }

    # Retrieve the GeoJSON data using requests
    response = requests.get(arcgis_url, params=params)
    monk_seal_geojson = response.json()

    # Load shapefiles
    eez_gdf = gpd.read_file(eez_file_path)

    benth_gdf = gpd.read_file(benthic_file_path)
    benth_gdf = benth_gdf.to_crs("EPSG:4326") #convert to base map ref

    coral_gdf = gpd.read_file(noaa_coral_path)
    coral_gdf = coral_gdf.to_crs("EPSG:4326") #convert to base map ref

    whale_gdf = gpd.read_file(whale_path)
    whale_gdf = whale_gdf.to_crs("EPSG:4326") #convert to base map ref

    ocean_depth_gdf = gpd.read_file(ocean_depth_path)
    ocean_depth_gdf = ocean_depth_gdf.to_crs("EPSG:4326") #convert to base map ref system

    dump_area_gdf = gpd.read_file(dumping_area_path)
    dump_area_gdf = dump_area_gdf.to_crs("EPSG:4326")

    moku_gdf = gpd.read_file(moku_file_path)  # Load the Moku Ridge to Reef shapefile
    moku_gdf = moku_gdf.to_crs("EPSG:4326")

    bottom_type_gdf = gpd.read_file(bottom_type_file_path)  # Load the Bottom Type point shapefile
    bottom_type_gdf = bottom_type_gdf.to_crs("EPSG:4326")

    cables_gdf = gpd.read_file(cable_file_path)
    cables_gdf = cables_gdf.to_crs("EPSG:4326")

    ordnance_gdf = gpd.read_file(ordnance_file_path)
    ordnance_gdf = ordnance_gdf.to_crs("EPSG:4326")

    explosive_dumping_gdf = gpd.read_file(explosive_dumping_file_path)
    explosive_dumping_gdf = explosive_dumping_gdf.to_crs("EPSG:4326")

    sewer_gdf = gpd.read_file(offshore_sewer_path)
    sewer_gdf = sewer_gdf.to_crs("EPSG:4326")

    buoys_sub_gdf = gpd.read_file(buoys_sub_path)
    buoys_sub_gdf = buoys_sub_gdf.to_crs("EPSG:4326")

    nav_aid_gdf = gpd.read_file(nav_aid_path)
    nav_aid_gdf = nav_aid_gdf.to_crs("EPSG:4326")

    obstruct_gdf = gpd.read_file(obstruct_path)
    obstruct_gdf = obstruct_gdf.to_crs("EPSG:4326")

    fads_gdf = gpd.read_file(fads_path)
    fads_gdf = fads_gdf.to_crs("EPSG:4326")

    zoning_gdf = gpd.read_file(zoning_path)
    zoning_gdf = zoning_gdf.to_crs("EPSG:4326")

############
    whale_habitat_gdf = gpd.read_file(whale_habitat_file_path)
    whale_habitat_gdf = whale_habitat_gdf.to_crs("EPSG:4326")  # Convert to WGS84 (base map reference)
    whale_habitat_json = whale_habitat_gdf.to_json()

    # Create a FeatureGroup for the False Killer Whale Critical Habitat
    
###########
    # Filter for US EEZ regions around Hawaii
    us_eez = eez_gdf[eez_gdf['geoname'].str.contains("United States Exclusive Economic Zone", na=False)]
    
    # Convert to GeoJSON format
    us_eez_json = us_eez.to_json()
    benthic_json = benth_gdf.to_json()
    ocean_depth_json = ocean_depth_gdf.to_json()
    dump_area_json = dump_area_gdf.to_json()
    moku_json = moku_gdf.to_json()
    cables_json = cables_gdf.to_json()
    sewer_json = sewer_gdf.to_json()
    coral_json = coral_gdf.to_json()
    whale_json = whale_gdf.to_json()
    zoning_json = zoning_gdf.to_json()


    # Create a FeatureGroup for EEZ
    eez_layer = folium.FeatureGroup(name="US Econ. Excl. Zone")
    benth_layer = folium.FeatureGroup(name="Benthic Habitats")
    ocean_depth_layer = folium.FeatureGroup(name="Ocean Depths")
    dump_area_layer = folium.FeatureGroup(name="Dumping Grounds")
    moku_layer = folium.FeatureGroup(name="Moku Ridge to Reef")  # Feature group for Moku Ridge to Reef
    bottom_type_layer = folium.FeatureGroup(name="Marine Bottom Types")  # Feature group for bottom types
    cables_layer = folium.FeatureGroup(name="Submarine Cables")
    ordnance_layer = folium.FeatureGroup(name="Unexploded Ordnance")
    explosive_dumping_layer = folium.FeatureGroup(name="Explosive Dumping Areas")
    sewer_layer = folium.FeatureGroup(name="Sewer Lines")
    coral_layer = folium.FeatureGroup(name="NOAA Coral From Naut. Charts")
    whale_layer = folium.FeatureGroup(name="Whale Sanctuary")
    buoys_sub_layer = folium.FeatureGroup(name="Submerged Buoys")
    nav_aid_layer = folium.FeatureGroup(name="Navigation Aids")
    obstruct_layer = folium.FeatureGroup(name="Obstructions")
    fads_layer = folium.FeatureGroup(name="Fish Aggregation Devices (FADs)")
    zoning_layer = folium.FeatureGroup(name="City Zoning")
    whale_habitat_layer = folium.FeatureGroup(name="False Killer Whale Critical Habitat")
    monk_seal_layer = folium.FeatureGroup(name="Hawaiian Monk Seal Habitat")

    # Add EEZ to the FeatureGroup
    folium.GeoJson(
        us_eez_json,
        style_function=lambda feature: {'color': 'green', 'weight': 2, 'fillOpacity': 0.1}
    ).add_to(eez_layer)

    #custom style function for the benthic habitats
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
        return {
            'fillColor': color_map.get(zone, 'brown'),
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.4
        }
    
    # Define a custom style function for the Coral Reefs
    def coral_style(feature):
        return {
            'fillColor': 'teal',  # Use a distinctive color for coral reefs
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6
        }

    # Define a style function for the habitat lines
    def monk_seal_style(feature):
        return {
            'color': 'red',
            'weight': 3,
            'opacity': 0.7,
        }

    def whale_style(feature):
        return {
            'fillColor': 'lightblue',  # Use a distinctive color for the whale sanctuary
            'color': 'blue',
            'weight': 1,
            'fillOpacity': 0.5
        }

    # style funciton for ocean depths
    def depth_style(feature):
        depth = feature['properties'].get('DEPTH')
        depth_colors = {
            1000: '#BFEFFF',  # 0-1000 meters: Light Cyan
            2000: '#7EC0EE',  # 1001-2000 meters: Light Blue
            3000: '#4682B4',  # 2001-3000 meters: Steel Blue
            4000: '#104E8B',  # 3001-4000 meters: Deep Blue
            5000: '#000080'   # 4001+ meters: Navy Blue
        }
        return {
            'fillColor': depth_colors.get(depth, '#708090'),  # Default to Slate Gray for unknown depths
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.5
        }
    
    # Style function for dumping areas
    def dumping_style(feature):
        return {
            'fillColor': 'brown',  # Set a distinctive color for dumping areas
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6
        }

    def moku_style(feature):
        moku = feature['properties'].get('Moku')
        color_map = {
            'Hilo': 'orange',
            'Kona': 'red',
            'Puna': 'yellow',
            'Kohala': 'green',
            'Kau': 'blue',
            'Molokai': 'purple',
            'Oahu': 'brown',
            'Maui': 'cyan'
        }
        return {
            'fillColor': color_map.get(moku, 'gray'),  # Use distinct colors for each Moku or default to gray
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6
        }
    
    # style function for the cables based on attributes like `Descript`
    def cable_style(feature):
        return {
            'color': 'blue',
            'weight': 3,  # Adjust line thickness
            'opacity': 0.7,
            'dashArray': '5, 5' if 'Descript' in feature['properties'] and 'Proposed' in feature['properties']['Descript'] else None
        }
    
        # style function for the cables based on attributes like `Descript`
    def sewer_style(feature):
        return {
            'color': 'brown',
            'weight': 3,  # Adjust line thickness
            'opacity': 0.7,
        }
    
    # custom style function for the Explosive Dumping Areas
    def explosive_style(feature):
        return {
            'fillColor': 'orange',  # Use a distinctive color for explosive dumping areas
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6
        }
    
    def maui_zoning_style(feature):
        # Define color mapping based on ZONE_CLASS
        zone_class = feature['properties']['ZONE_CLASS']
        color_map = {
            'R-3 Residential': 'lightblue',
            'AP Airport': 'orange',
            'AG Agricultural': 'green',
            'INT Interim': 'brown',
            'PK Park': 'red',
            'M-2 Heavy Industrial': 'gray',
            'H-2 Hotel': 'purple',
            'M-1 Light Industrial': 'yellow',
            'B-2 Business-Community': 'cyan'
        }
        
        # Define default color
        fill_color = color_map.get(zone_class, 'black')
        
        return {
            'fillColor': fill_color,
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.6,
        }

    # Custom style function for habitat areas based on `HABTYPE`
    def whale_habitat_style(feature):
        hab_type = feature['properties']['HABTYPE']
        color_map = {
            'marine': 'blue',
            'nearshore': 'cyan',
            'offshore': 'purple',
            'unknown': 'gray'
        }
        return {
            'fillColor': color_map.get(hab_type, 'lightgray'),
            'color': 'black',
            'weight': 1,
            'fillOpacity': 0.5,
        }

    # Add Benthic Habitat to the FeatureGroup
    folium.GeoJson(
        benthic_json,
        name="Benthic Habitat Zones",
        style_function=benthic_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["M_STRUCT", "D_STRUCT", "M_COVER", "ZONE"],
            aliases=["Major Structure", "Detailed Structure", "Major Cover", "Zone"],
            localize=True
        )
    ).add_to(benth_layer)

    folium.GeoJson(
        zoning_json,
        name="Maui County Zoning",
        style_function=maui_zoning_style,
        tooltip=folium.GeoJsonTooltip(
            fields=['ZONE_CLASS', 'CP_AREA'],
            aliases=['Zoning Class:', 'Community Plan Area:'],
            localize=True,
            sticky=False,
            labels=True
        )
    ).add_to(zoning_layer)

        # Add the Coral Reefs layer using GeoJson with tooltips
    folium.GeoJson(
        coral_json,
        name="Coral Reefs (Nautical Charts)",
        style_function=coral_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["acres", "chart1", "chart2", "notes"],
            aliases=["Area (Acres)", "Source Chart 1", "Source Chart 2", "Notes"],
            localize=True
        )
    ).add_to(coral_layer)

    folium.GeoJson(
        whale_gdf.to_json(),
        name="Whale Sanctuary",
        style_function=whale_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["sanctuary", "Shape_Leng", "Shape_Area"],
            aliases=["Sanctuary Name", "Perimeter (m)", "Area (sq. m)"],
            localize=True
        )
    ).add_to(whale_layer)

    # Add Ocean Depth layer with depth styling
    folium.GeoJson(
        ocean_depth_json,
        name="Ocean Depth",
        style_function=depth_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["depth"],
            aliases=["Depth (m)"],
            localize=True
        )
    ).add_to(ocean_depth_layer)

    # Add Dumping Areas layer
    folium.GeoJson(
        dump_area_json,
        name="Dumping Areas",
        style_function=dumping_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["id", "chart", "cfr", "notes"],  # Display relevant fields in the tooltip
            aliases=["Polygon ID", "Chart Source", "Regulation", "Notes"],
            localize=True
        )
    ).add_to(dump_area_layer)

    # Add Moku Ridge to Reef layer
    folium.GeoJson(
        moku_json,
        name="Moku Ridge to Reef",
        style_function=moku_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["Moku_ID", "Moku", "Mokupuni"],
            aliases=["Moku ID", "Moku", "Island"],
            localize=True
        )
    ).add_to(moku_layer)

        # Add the Cables layer using GeoJson with the custom style
    folium.GeoJson(
        cables_json,
        name="Submarine Cables",
        style_function=cable_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["id", "length", "descriptio", "chart1", "chart2"],
            aliases=["Cable ID", "Length (m)", "Description", "Primary Chart", "Secondary Chart"],
            localize=True
        )
    ).add_to(cables_layer)

    #sewer layer with custom style
    folium.GeoJson(
        sewer_json,
        name="Sewer Lines",
        style_function=sewer_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["id", "notes", "island", "chart"],
            aliases=["sewer ID", "notes", "island", "Chart"],
            localize=True
        )
    ).add_to(sewer_layer)

    # Add the Explosive Dumping Areas layer using GeoJson
    folium.GeoJson(
        explosive_dumping_gdf.to_json(),
        name="Explosive Dumping Areas",
        style_function=explosive_style,
        tooltip=folium.GeoJsonTooltip(
            fields=["chart", "notes", "acres"],
            aliases=["Source Chart", "Notes", "Area (Acres)"],
            localize=True
        )
    ).add_to(explosive_dumping_layer)

    # Add False Killer Whale Habitat layer to the map with tooltips
    folium.GeoJson(
        whale_habitat_json,
        name="False Killer Whale Habitat",
        style_function=whale_habitat_style,
        tooltip=folium.GeoJsonTooltip(
            fields=['SCIENAME', 'COMNAME', 'HABTYPE', 'AREASqKm'],
            aliases=['Species Name', 'Common Name', 'Habitat Type', 'Area (sq. km)'],
            localize=True,
            sticky=False,
            labels=True
        )
    ).add_to(whale_habitat_layer)

    # Add the GeoJson layer to the map
    folium.GeoJson(
        monk_seal_geojson,
        name="Monk Seal Habitat",
        style_function=monk_seal_style,
        tooltip=folium.GeoJsonTooltip(
            fields=['ID', 'COMNAME', 'HABTYPE', 'AREASqKm'],
            aliases=['ID', 'Common Name', 'Habitat Type', 'Area (sq. km)'],
            localize=True,
            sticky=False,
            labels=True
        )
    ).add_to(monk_seal_layer)

     # Add Bottom Type Points as Markers
    for _, row in bottom_type_gdf.iterrows():
        # Get coordinates and attributes for each point
        coords = [row.geometry.y, row.geometry.x]  # Lat, Lon order for Folium
        seabed_type = row['seabed']
        code = row['code']
        point_id = row['id']
        sourcethm = row['sourcethm']

        # Create a CircleMarker for each point
        folium.CircleMarker(
            location=coords,
            radius=5,  # Marker size
            color='blue',
            fill=True,
            fill_color='blue',
            fill_opacity=0.7,
            popup=folium.Popup(f"ID: {point_id}<br>Code: {code}<br>Seabed: {seabed_type}<br>Source: {sourcethm}", max_width=250)
        ).add_to(bottom_type_layer)

    # Add Unexploded Ordnance Points as Markers
    for _, row in ordnance_gdf.iterrows():
        # Get coordinates and attributes for each point
        coords = [row.geometry.y, row.geometry.x]  # Latitude, Longitude order for Folium
        point_id = row['id']
        notes = row['notes']
        chart = row['chart']

        # Create a CircleMarker for each point
        folium.CircleMarker(
            location=coords,
            radius=6,  # Marker size
            color='red',
            fill=True,
            fill_color='red',
            fill_opacity=0.7,
            popup=folium.Popup(f"ID: {point_id}<br>Notes: {notes}<br>Source Chart: {chart}", max_width=250)
        ).add_to(ordnance_layer)

    # Add Submerged Buoys Points as Markers
    for _, row in buoys_sub_gdf.iterrows():
        if row.geometry and row.geometry.is_valid:
            coords = [row.geometry.y, row.geometry.x]
            buoy_id = row['name']
            description = row['notes']

            folium.CircleMarker(
                location=coords,
                radius=5,
                color='purple',
                fill=True,
                fill_color='purple',
                fill_opacity=0.7,
                popup=folium.Popup(f"ID: {buoy_id}<br>Description: {description}", max_width=250)
            ).add_to(buoys_sub_layer)

    # Add Navigation Aids Points as Markers
    for _, row in nav_aid_gdf.iterrows():
        if row.geometry and row.geometry.is_valid:
            coords = [row.geometry.y, row.geometry.x]
            aid_id = row['id']
            aid_type = row['type']
            name = row['name']

            folium.CircleMarker(
                location=coords,
                radius=6,
                color='darkgreen',
                fill=True,
                fill_color='darkgreen',
                fill_opacity=0.7,
                popup=folium.Popup(f"ID: {aid_id}<br>Type: {aid_type}<br>Name: {name}", max_width=250)
            ).add_to(nav_aid_layer)

    # Add Obstructions Points as Markers
    for _, row in obstruct_gdf.iterrows():
        if row.geometry and row.geometry.is_valid:
            coords = [row.geometry.y, row.geometry.x]
            obstruct_id = row['id']
            notes = row['notes']

            folium.CircleMarker(
                location=coords,
                radius=6,
                color='orange',
                fill=True,
                fill_color='orange',
                fill_opacity=0.7,
                popup=folium.Popup(f"ID: {obstruct_id}<br>Notes: {notes}", max_width=250)
            ).add_to(obstruct_layer)

    # Add FADs Points as Markers
    for _, row in fads_gdf.iterrows():
        if row.geometry and row.geometry.is_valid:
            coords = [row.geometry.y, row.geometry.x]
            fads_id = row['name']
            description = row['island']

            folium.CircleMarker(
                location=coords,
                radius=6,
                color='blue',
                fill=True,
                fill_color='blue',
                fill_opacity=0.7,
                popup=folium.Popup(f"ID: {fads_id}<br>Description: {description}", max_width=250)
            ).add_to(fads_layer)

    # Add EEZ and Benthic FeatureGroups to the map
    eez_layer.add_to(hawaii_map)
    benth_layer.add_to(hawaii_map)
    ocean_depth_layer.add_to(hawaii_map)
    dump_area_layer.add_to(hawaii_map)
    moku_layer.add_to(hawaii_map)
    bottom_type_layer.add_to(hawaii_map)  # Add bottom type layer
    cables_layer.add_to(hawaii_map)
    ordnance_layer.add_to(hawaii_map)
    explosive_dumping_layer.add_to(hawaii_map)
    sewer_layer.add_to(hawaii_map)
    buoys_sub_layer.add_to(hawaii_map)
    nav_aid_layer.add_to(hawaii_map)
    obstruct_layer.add_to(hawaii_map)
    fads_layer.add_to(hawaii_map)
    whale_layer.add_to(hawaii_map)
    zoning_layer.add_to(hawaii_map)
    coral_layer.add_to(hawaii_map)
    whale_habitat_layer.add_to(hawaii_map)
    monk_seal_layer.add_to(hawaii_map)


except Exception as e:
    st.write(f"Error loading shapefiles: {e}")

# Add LayerControl to enable toggling
folium.LayerControl().add_to(hawaii_map)

# Display the map and capture the user's drawn shapes
st_data = st_folium(hawaii_map, width=800, height=500)

#create info container
info_container = st.container()

#function to layout turbines
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

            # Convert the drawn coordinates to a Polygon for area calc
            drawn_polygon = Polygon([(coord[0], coord[1]) for coord in rectangle_coords])
            # Create a GeoDataFrame and project it to a UTM zone to get the area in square meters
            gdf = gpd.GeoDataFrame(index=[0], crs="EPSG:4326", geometry=[drawn_polygon])
            gdf = gdf.to_crs(epsg=32604)  # Convert to UTM Zone 4 for accurate area calculation
            # Calculate area in square meters and convert to square kilometers
            area_sq_km = gdf.geometry.area.iloc[0] / 1e6
            # Calculate the dimensions (length and width) of the rectangle in kilometers
            lon_min, lat_min = gdf.geometry.iloc[0].exterior.coords[0]
            lon_max, lat_max = gdf.geometry.iloc[0].exterior.coords[2]
            
            # Calculate the dimensions
            width_m = gdf.to_crs(epsg=32604).geometry.iloc[0].bounds[2] - gdf.to_crs(epsg=32604).geometry.iloc[0].bounds[0]
            height_m = gdf.to_crs(epsg=32604).geometry.iloc[0].bounds[3] - gdf.to_crs(epsg=32604).geometry.iloc[0].bounds[1]

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

