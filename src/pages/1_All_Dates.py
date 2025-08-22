import numpy
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
import leafmap
import leafmap.foliumap


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

def get_map(kelp_total_extents: geopandas.GeoDataFrame, kelp_info: pandas.DataFrame, display_range: tuple):
    
    if streamlit.session_state.date_by_date_index != streamlit.session_state.date_by_date_map_index or streamlit.session_state.date_by_date_percentiles != display_range:
        collection = "sentinel-2-l2a"
        rgb_bands = utils.get_band_names_from_common(['red', 'green', 'blue'])
        
        date_index = streamlit.session_state.date_by_date_index[0]
        # Either create map or just use the existing if no change to date
        folium_map = leafmap.foliumap.Map() #location=center, zoom_start=13)

        #land.explore(m=folium_map, color="blue", style_kwds={"fillOpacity": 0}, name="land")
        kelp_total_extents.explore(m=folium_map, color="blue", style_kwds={"fillOpacity": 0}, name="Satellite record Kelp Extents")
        if isinstance(kelp_info["file"].iloc[date_index], str) and kelp_info["file"].iloc[date_index] != "":
            csv_file_path = pathlib.Path(kelp_info["file"].iloc[date_index])
            kelp_polygons = geopandas.read_file(csv_file_path)
            kelp_polygons.explore(m=folium_map, color="magenta", style_kwds={"fillOpacity": 0}, name="Quarter Kelp Extents")
        else:
            streamlit.text("Anomalous date flagged in manual QA/QC. Excluded")

        tile_ids = kelp_info["Satellite Tile IDs"].iloc[date_index].replace(",", "").strip(" ").split(" ")

        #streamlit.text(f"Percentile 2: {percentiles_2}, Percentile 98: {percentiles_98}.")
        for index, tile_id in enumerate(tile_ids):
            folium_map.add_stac_layer(collection=collection, item=tile_id, assets=rgb_bands, name="RGB", titiler_endpoint="pc",
                                      rescale=f"{display_range[0]},{display_range[1]}", fit_bounds=False)

        bounds = kelp_total_extents.to_crs(utils.CRS_WSG).total_bounds  # [minx, miny, maxx, maxy]
        folium_map.fit_bounds(numpy.flip(numpy.reshape(bounds, (2,2)), axis = 1).tolist())

        folium.LayerControl().add_to(folium_map)

        streamlit.session_state.date_by_date_map = folium_map
        streamlit.session_state.date_by_date_map_index = streamlit.session_state.date_by_date_index
    
    return streamlit.session_state.date_by_date_map

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
    
    if 'date_by_date_index' not in streamlit.session_state:
        streamlit.session_state.date_by_date_index = []
        streamlit.session_state.date_by_date_map_index = []
        streamlit.session_state.date_by_date_prev_location = None
        streamlit.session_state.date_by_date_percentiles = None
    
    
    data_path = pathlib.Path.cwd() / "data"
    
    streamlit.button("Re-run")
    streamlit.title('Kelp Demo - click area plot to select raster display')
    
    test_sites = geopandas.read_file(data_path / "vectors" / "test_sites_offshore_3km.gpkg")
    
    location = streamlit.selectbox("Select tile to display", (test_sites["name"]), index=0,)
    
    if location != streamlit.session_state.date_by_date_prev_location:
        streamlit.session_state.date_by_date_prev_location = location
        streamlit.session_state.date_by_date_index = streamlit.session_state.date_by_date_map_index = []
    
    # Define the region
    raster_path = data_path / "rasters" / "test_sites" / f"{location}"
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    kelp_total_extents = geopandas.read_file(raster_path / "presence_absence_map.gpkg")
    
    if 'xy' not in streamlit.session_state:
        streamlit.session_state['xy'] = []
        streamlit.session_state['prev_lat'], streamlit.session_state['prev_lon'] = numpy.nan, numpy.nan
        streamlit.session_state['spectra'] = None
    
    streamlit.subheader("Table and Plot of areas")
    if (raster_path / "info.csv").exists():
        kelp_info = pandas.read_csv(raster_path / "info.csv"); kelp_info.drop(columns=["Unnamed: 0"], inplace=True, errors="ignore")
    col1, col2 = streamlit.columns([1, 5])
    with col1:
        if (raster_path / "info.csv").exists():
            streamlit.dataframe(kelp_info[["date", "area", "ocean cloud percentage", "proportion of max coverage"]])
        else:
            streamlit.markdown("No info.csv. Older runs displayed, but raster view of selected date is not supported")
    with col2:
        #streamlit.plotly_chart(plotly.express.line(kelp_info, x = 'date', y = ["ocean cloud percentage"], markers=True, text="date"))
        figure = plotly.graph_objects.Figure()
        area_columns = []
        for older_info_file in raster_path.glob("info_*.csv"):
            name = f"{older_info_file.stem.replace('area_', '').replace('info_', '')}"
            older_kelp_info = pandas.read_csv(older_info_file)
            older_kelp_info.drop(columns=["Unnamed: 0", "file", "ocean cloud percentage"], inplace=True, errors="ignore")
            older_kelp_info.rename(inplace=True, columns={"area": name})
            area_columns.append(f"area {older_info_file.stem.replace('area_', '')}")
            figure.add_trace(plotly.graph_objects.Scatter(x=older_kelp_info["date"], y=older_kelp_info[name], mode="lines+markers", name=name))
        
        if (raster_path / "info.csv").exists():
            figure = plotly.subplots.make_subplots(specs=[[{"secondary_y": True}]])
            figure.add_trace(plotly.graph_objects.Scatter(x=kelp_info["date"], y=kelp_info["area"], mode="lines+markers", marker={'color':'blue'}, name="Area [m^2]" ), secondary_y=True)
            figure.add_trace(plotly.graph_objects.Scatter(x=kelp_info["date"], y=kelp_info["proportion of max coverage"] * 100, mode="lines+markers", marker={'color':'red'}, name="Proportion of Max Coverage [%]" ), secondary_y=False)
            
            figure.update_layout(title="Proportion of max coverage by date across algorithm runs",
                                 xaxis_title="date", yaxis2_title="Area [m^2]",
                                 yaxis_title="Max Coverage Proportion [%]", yaxis={'range': [0, 75]},
                                 legend_title="Kelp Cover", template="plotly_white" )
        event = streamlit.plotly_chart(figure, on_select="rerun")
        
    selection = event["selection"]["point_indices"]
    if len(selection) > 0:
        streamlit.session_state.date_by_date_index = selection
    
    if len(streamlit.session_state.date_by_date_index) and (raster_path / "info.csv").exists():
        date_index = streamlit.session_state.date_by_date_index[0]
        streamlit.subheader(f"Plot {kelp_info['date'].iloc[date_index]}.")
        streamlit.caption("May take time to load...")
        
        percentiles_2 = kelp_info["Percentile 2"].iloc[date_index].replace(",", "").strip(" ").split(" ")
        percentiles_98 = kelp_info["Percentile 98"].iloc[date_index].replace(",", "").strip(" ").split(" ")
        percentiles_2 = [int(percentile_2) for percentile_2 in percentiles_2]
        percentiles_98 = [int(percentile_98) for percentile_98 in percentiles_98]
        
        streamlit.session_state.date_by_date_percentiles = (int(numpy.array(percentiles_2).mean()), int(numpy.array(percentiles_98).mean()))
        
        display_range = streamlit.slider(
            f'Satellite Tile(s) default range: [{streamlit.session_state.date_by_date_percentiles[0]}, {streamlit.session_state.date_by_date_percentiles[1]}]',
            0, 10000, streamlit.session_state.date_by_date_percentiles)
        
        folium_map = get_map(kelp_total_extents, kelp_info, display_range)
        streamlit_folium.folium_static(folium_map, width=900)


if __name__ == '__main__':
    main()