import numpy
import leafmap
import pandas
import geopandas
import pathlib
import streamlit
import streamlit_folium
import argparse
import pages.scripts.colourmaps
import datetime
import folium
import rioxarray
import odc.stac
import branca
import plotly.express
import plotly.graph_objects
import pyproj
import io


#import holoviews
#import hvplot.xarray

module_path = pathlib.Path.cwd() / 'scripts'
import sys
if str(module_path) not in sys.path:
    sys.path.append(str(module_path))
import utils
import importlib
importlib.reload(utils)

def save_geojson_with_bytesio(dataframe):
    #Function to return bytesIO of the geojson
    bindary_stream = io.BytesIO()
    dataframe.to_file(bindary_stream,  driver='GeoJSON')
    return bindary_stream


def get_map(kelp_total_extents: geopandas.GeoDataFrame, land: geopandas.GeoDataFrame, kelp_info: pandas.DataFrame):
    
    if streamlit.session_state.quarterly_index_3 != streamlit.session_state.quareterly_map_index_3:
        collection = "sentinel-2-l2a"
        rgb_bands = utils.get_band_names_from_common(['red', 'green', 'blue'])
        
        date_index = streamlit.session_state.quarterly_index_3[0]
        # Either create map or just use the existing if no change to date
        csv_file_path = pathlib.Path(kelp_info["file"].iloc[date_index])
        streamlit.subheader(f"Plot quarter {kelp_info["date"].iloc[date_index]} calculated from dates {kelp_info["dates considered"].iloc[date_index]}.")
        streamlit.caption("May take time to load...")

        folium_map = leafmap.foliumap.Map() #location=center, zoom_start=13)
        
        kelp_polygons = geopandas.read_file(csv_file_path)
        land.explore(m=folium_map, color="green", style_kwds={"fillOpacity": 0}, name="Land Outline")
        if not kelp_total_extents is None: 
            kelp_total_extents.explore(m=folium_map, color="blue", style_kwds={"fillOpacity": 0}, name="Satellite record Kelp Extents")
        kelp_polygons.explore(m=folium_map, color="magenta", style_kwds={"fillOpacity": 0}, name="Quarter Kelp Extents")

        tile_ids = kelp_info["Satellite Tile IDs"].iloc[date_index].replace(",", "").strip(" ").split(" ")
        percentiles_2 = kelp_info["Percentile 2"].iloc[date_index].replace(",", "").strip(" ").split(" ")
        percentiles_98 = kelp_info["Percentile 98"].iloc[date_index].replace(",", "").strip(" ").split(" ")
        #streamlit.text(f"Percentile 2: {percentiles_2}, Percentile 98: {percentiles_98}.")
        for index, tile_id in enumerate(tile_ids):
            folium_map.add_stac_layer(collection=collection, item=tile_id, assets=rgb_bands, name="RGB", titiler_endpoint="pc", rescale=f"{percentiles_2[index]},{percentiles_98[index]}", fit_bounds=False)

        bounds = kelp_polygons.to_crs(utils.CRS_WSG).total_bounds  # [minx, miny, maxx, maxy]
        folium_map.fit_bounds(numpy.flip(numpy.reshape(bounds, (2,2)), axis = 1).tolist())

        folium.LayerControl().add_to(folium_map)

        streamlit.session_state.quarterly_map = folium_map
        streamlit.session_state.quareterly_map_index_3 = streamlit.session_state.quarterly_index_3
    
    return streamlit.session_state.quarterly_map

def main():
    """ Create / Update the geofabric summary information and display in a dashboard.
    """
    
    streamlit.set_page_config(
        page_title="Test Sites",
        page_icon="🌏",
        layout="wide",
    )
    display_size = 700
    date_format = "%Y-%m-%d"
    
    if 'quarterly_index_3' not in streamlit.session_state:
        streamlit.session_state.quarterly_index_3 = []
        streamlit.session_state.quareterly_map_index_3 = []
        streamlit.session_state.quareterly_prev_location_3 = None
    
    data_path = pathlib.Path.cwd() / "data"
    remote_raster_path = pathlib.Path("/nesi/nobackup/niwa03660/ZBD2023_outputs/test_sites_quarterly")
    
    streamlit.button("Re-run")
    streamlit.title('Kelp Demo - click area plot to select raster display')
    
    test_sites = geopandas.read_file(data_path / "vectors" / "test_sites_offshore_3km.gpkg")
    
    location = streamlit.selectbox("Select tile to display", (test_sites["name"]), index=0,)
    
    if location != streamlit.session_state.quareterly_prev_location_3:
        streamlit.session_state.quareterly_prev_location_3 = location
        streamlit.session_state.quarterly_index_3 = streamlit.session_state.quareterly_map_index_3 = []
    
    # Define the region
    raster_path = data_path / "rasters" / "test_sites_quarterly" / f"{location}"
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    if (raster_path / "presence_absence_map.gpkg").exists():
        kelp_total_extents = geopandas.read_file(raster_path / "presence_absence_map.gpkg")
        streamlit.download_button(label=f"Download {location}'s Kelp Map 2016-2024",
                                          data=save_geojson_with_bytesio(kelp_total_extents),
                                          file_name=f"kelp_presence_absence_map_{location}.geojson",
                                          mime="application/geo+json",
                                         )
    
    streamlit.subheader("Table and Plot of areas")
    if (raster_path / "info_quarterly.csv").exists():
        kelp_info = pandas.read_csv(raster_path / "info_quarterly.csv"); kelp_info.drop(columns=["Unnamed: 0"], inplace=True, errors="ignore")
    col1, col2 = streamlit.columns([1, 5])
    with col1:
        if (raster_path / "info_quarterly.csv").exists():
            streamlit.dataframe(kelp_info[["date", "area", "dates considered", "max coverage date"]])
        else:
            streamlit.markdown("No info_quarterly.csv. Older runs displayed, but raster view of selected date is not supported")
    with col2:
        #streamlit.plotly_chart(plotly.express.line(kelp_info, x = 'date', y = ["dates considered"], markers=True, text="date"))
        figure = plotly.graph_objects.Figure()
        area_columns = []
        for older_info_file in raster_path.glob("info_quarterly*.csv"):
            name = f"{older_info_file.stem.replace('area_', '').replace("info_quarterly", "")}"
            older_kelp_info = pandas.read_csv(older_info_file)
            older_kelp_info.drop(columns=["Unnamed: 0", "file"], inplace=True, errors="ignore")
            older_kelp_info.rename(inplace=True, columns={"area": name})
            area_columns.append(f"area {older_info_file.stem.replace('area_', '')}")
            figure.add_trace(plotly.graph_objects.Scatter(x=older_kelp_info["date"], y=older_kelp_info[name], mode="lines+markers", name=name))
        
        if (raster_path / "info_quarterly.csv").exists():
            figure.add_trace(plotly.graph_objects.Scatter(x=kelp_info["date"], y=kelp_info["area"], mode="lines+markers", name="Latest run (harmonised 2-pass filtering & same thresholds)" ))
        figure.update_layout(title="Area's by date across algorithm runs", xaxis_title="date", yaxis_title="Area [m^2]", legend_title="Algorithms", template="plotly_white" )
        event = streamlit.plotly_chart(figure, on_select="rerun")
        
    selection = event["selection"]["point_indices"]
    if len(selection) > 0:
        streamlit.session_state.quarterly_index_3 = selection
    
    #streamlit.text(f"quarterly_index: {streamlit.session_state.quarterly_index_3}, quareterly_map_index: {streamlit.session_state.quareterly_map_index_3}, selection: {selection}")
    
    if len(streamlit.session_state.quarterly_index_3) and (raster_path / "info_quarterly.csv").exists():
        folium_map = get_map(kelp_total_extents, land, kelp_info)
        streamlit_folium.folium_static(folium_map, width=900)


if __name__ == '__main__':
    main()