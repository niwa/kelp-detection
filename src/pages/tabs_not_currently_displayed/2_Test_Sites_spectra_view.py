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
            streamlit.text(csv_file_path.name)
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
            name = f"{older_info_file.stem.replace('area_', '').replace("info_", "")}"
            older_kelp_info = pandas.read_csv(older_info_file)
            older_kelp_info.drop(columns=["Unnamed: 0", "file", "ocean cloud percentage"], inplace=True, errors="ignore")
            older_kelp_info.rename(inplace=True, columns={"area": name})
            area_columns.append(f"area {older_info_file.stem.replace('area_', '')}")
            figure.add_trace(plotly.graph_objects.Scatter(x=older_kelp_info["date"], y=older_kelp_info[name], mode="lines+markers", name=name))
        
        if (raster_path / "info.csv").exists():
            figure.add_trace(plotly.graph_objects.Scatter(x=kelp_info["date"], y=kelp_info["area"], mode="lines+markers", name="Latest run (two pass NDVI 0.03, anomaly removal using NDVI 0.213)" ))
        figure.update_layout(title="Area's by date across algorithm runs", xaxis_title="date", yaxis_title="Area [m^2]", legend_title="Algorithms", template="plotly_white" )
        event = streamlit.plotly_chart(figure, on_select="rerun")
        
    selection = event["selection"]["point_indices"]
    if len(selection) > 0:
        streamlit.session_state.date_by_date_index = selection
    
    if len(streamlit.session_state.date_by_date_index) and (raster_path / "info.csv").exists():
        date_index = streamlit.session_state.date_by_date_index[0]
        streamlit.subheader(f"Plot {kelp_info["date"].iloc[date_index]}.")
        streamlit.caption("May take time to load...")
        
        percentiles_2 = kelp_info["Percentile 2"].iloc[date_index].replace(",", "").strip(" ").split(" ")
        percentiles_98 = kelp_info["Percentile 98"].iloc[date_index].replace(",", "").strip(" ").split(" ")
        percentiles_2 = [int(percentile_2) for percentile_2 in percentiles_2]
        percentiles_98 = [int(percentile_98) for percentile_98 in percentiles_98]
        
        streamlit.session_state.date_by_date_percentiles = (int(numpy.array(percentiles_2).mean()), int(numpy.array(percentiles_98).mean()))
        
        display_range = streamlit.slider(f'Satellite Tile(s) default range: [{streamlit.session_state.date_by_date_percentiles[0]}, {streamlit.session_state.date_by_date_percentiles[1]}]',
                                                  0, 10000, streamlit.session_state.date_by_date_percentiles)
        
        folium_map = get_map(kelp_total_extents, kelp_info, display_range)
        streamlit_folium.folium_static(folium_map, width=900)
    
    ################
    
    '''if len(selection) and (raster_path / "info.csv").exists():
        csv_file_path = pathlib.Path(kelp_info["file"].iloc[selection[0]])
        data_file = csv_file_path # remote_raster_path / csv_file_path.parent.name / csv_file_path.name
        streamlit.subheader(f"Plot {kelp_info["date"].iloc[selection[0]]} with {kelp_info["ocean cloud percentage"].iloc[selection[0]]:.2f}% ocean cloud")
        streamlit.caption("May take time to load...")
        
        data = rioxarray.rioxarray.open_rasterio(data_file, chunks=True).squeeze("band", drop=True)
        transformer = pyproj.Transformer.from_crs(utils.CRS_WSG, data.rio.crs)
        
        with open(data_file, "rb") as binary_file:
            streamlit.download_button(label="Download satellite layers",
                                      data=binary_file, file_name=f"satellite_{location}_{kelp_info["date"].iloc[selection[0]]}.nc",
                                             )
        if not (data_file.parent / f'rgb_{data_file.name}').exists():
            rgb = data[["B04", "B03","B02"]].to_array("rgb", name="all images")
            utils.update_raster_defaults(rgb)
            rgb.to_netcdf(data_file.parent / f'rgb_{data_file.name}', format="NETCDF4", engine="netcdf4", encoding={"all images": {"zlib": True, "complevel": 5, "grid_mapping": rgb.encoding["grid_mapping"]}})
        with open(data_file.parent / f'rgb_{data_file.name}', "rb") as binary_file:
            streamlit.download_button(label="Download rgb layers",
                                      data=binary_file, file_name=f"satellite_rgb_{location}_{kelp_info["date"].iloc[selection[0]]}.nc",
                                             )

        col1, col2 = streamlit.columns([2, 1])
        with col1:

            folium_map = folium.Map()
            
            data["SCL"].odc.add_to(map=folium_map, name="SCL", opacity=1, cmap="viridis", vmin=0, vmax=11)
            colormap = branca.colormap.linear.viridis.scale(0, 11)
            colormap.caption = 'SCL Index'
            colormap.add_to(folium_map)
            
            data["ndvi"].odc.add_to(map=folium_map, name="ndvi", opacity=1, cmap="inferno", vmin=0, vmax=1)
            data["ndvri"].odc.add_to(map=folium_map, name="ndvri", opacity=1, cmap="inferno", vmin=0, vmax=1)
            data["ndwi"].odc.add_to(map=folium_map, name="ndwi", opacity=1, cmap="inferno", vmin=0, vmax=1)
            data["ndwi2"].odc.add_to(map=folium_map, name="ndwi2", opacity=1, cmap="inferno", vmin=0, vmax=1)
            data["ndwi3"].odc.add_to(map=folium_map, name="ndwi3", opacity=1, cmap="inferno", vmin=0, vmax=1)
            data["ndwi4"].odc.add_to(map=folium_map, name="ndwi4", opacity=1, cmap="inferno", vmin=0, vmax=1)

            data["kelp"].odc.add_to(map=folium_map, name="Kelp", opacity=1, cmap="inferno", vmin=0, vmax=1)
            colormap = branca.colormap.linear.inferno.scale(0, 1)
            colormap.caption = 'Kelp Index'
            colormap.add_to(folium_map)
            
            # In future consider https://geoviews.org or saving as a png and loading...
            rgb_dict = {
                #"False NIR color": ["nir", "green", "blue"],
                #"False rededge-2 colour": ["rededge - Band 6 - Vegetation red edge 2", "green", "blue"],
                #"False rededge-3 colour": ["rededge - Band 7 - Vegetation red edge 3", "green", "blue"],
                "Satellite RBG": ["red", "green", "blue"]
            }
            for title, names in rgb_dict.items():
                bands = utils.get_band_names_from_common(names)
                rgb = utils.normalise_rgb(data, bands)
                rgb.odc.to_rgba(bands=bands, vmin=0, vmax=.7).odc.add_to(map=folium_map, name=title)
            
            folium_map.fit_bounds(data["kelp"].odc.map_bounds())
            
            land.explore(m=folium_map, color="blue", style_kwds={"fillOpacity": 0}, name="land")
            kelp_polygons = utils.polygon_from_raster(data["kelp"])
            kelp_polygons.explore(m=folium_map, color="magenta", style_kwds={"fillOpacity": 0}, name="kelp polygon")
            
            # create callback
            folium_map.add_child(folium.LatLngPopup())
    
            folium.LayerControl().add_to(folium_map)
            st_map =  streamlit_folium.st_folium(folium_map, width=900) 
        with col2:
            
            if st_map['last_clicked'] is not None:
                lat, lon = st_map['last_clicked']['lat'], st_map['last_clicked']['lng']
                if streamlit.session_state['prev_lat'] != lat or streamlit.session_state['prev_lon'] != lon:
                    streamlit.session_state['prev_lat'], streamlit.session_state['prev_lon'] = lat, lon
                    x, y =  transformer.transform(lat, lon)
                    streamlit.session_state['xy'].append([x,y])
                    plot = utils.plot_lines(streamlit.session_state['xy'], data)
                    streamlit.pyplot(plot.get_figure())
                    streamlit.session_state['spectra'] = utils.update_spectra([x,y], data, streamlit.session_state['spectra'])
                
                streamlit.download_button(label="Export specrta from clicks",
                                          data=pandas.DataFrame(data=streamlit.session_state['spectra']).to_csv().encode("utf-8"),
                                          file_name=f"spectra_from_dashboard_{location}.csv",
                                          mime="text/csv",
                                         )
                
                if streamlit.button('Prepare specrta all dates last click'):
                    streamlit.session_state['prev_lat'], streamlit.session_state['prev_lon'] = lat, lon
                    x, y =  transformer.transform(lat, lon)
                    spectra_all_dates = utils.get_spectra_all_dates([x, y], kelp_info)
                    streamlit.download_button(label="Export specrta all dates last click",
                                              data=pandas.DataFrame(data=spectra_all_dates).to_csv().encode("utf-8"),
                                              file_name=f"spectra_all_dates_from_dashboard_{location}.csv",
                                              mime="text/csv",
                                         )'''


if __name__ == '__main__':
    main()