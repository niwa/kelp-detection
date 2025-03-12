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
import pyproj
import zipfile


#import holoviews
#import hvplot.xarray

module_path = pathlib.Path.cwd() / 'scripts'
import sys
if str(module_path) not in sys.path:
    sys.path.append(str(module_path))
import utils
import importlib
importlib.reload(utils)


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
    data_path = pathlib.Path.cwd() / "data"
    remote_raster_path = pathlib.Path("/nesi/project/niwa03660/ZBD2023_outputs")
    
    streamlit.button("Re-run")
    streamlit.title('Kelp Demo - click area plot to select raster display')
    
    test_sites = geopandas.read_file(data_path / "vectors" / "test_sites_offshore_3km.gpkg")
    
    location = streamlit.selectbox("Select tile to display", (test_sites["name"]), index=0,)
    
    # Define the region
    raster_path = data_path / "rasters" / "test_sites" / f"{location}"
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    
    if 'xy' not in streamlit.session_state:
        streamlit.session_state['xy'] = []
        streamlit.session_state['prev_lat'], streamlit.session_state['prev_lon'] = numpy.nan, numpy.nan
        streamlit.session_state['spectra'] = None
    
    streamlit.subheader("Table and Plot of areas")
    kelp_info = pandas.read_csv(raster_path / "info.csv"); kelp_info.drop(columns=["Unnamed: 0"], inplace=True, errors="ignore")
    col1, col2 = streamlit.columns([1, 5])
    with col1:
        streamlit.dataframe(kelp_info[["date", "area", "ocean cloud percentage"]])
    with col2:
        #streamlit.plotly_chart(plotly.express.line(kelp_info, x = 'date', y = ["ocean cloud percentage"], markers=True, text="date"))
        area_columns = ["area"]
        for older_info_file in raster_path.glob("info_*.csv"):
            older_kelp_info = pandas.read_csv(older_info_file)
            older_kelp_info.drop(columns=["Unnamed: 0", "file", "ocean cloud percentage"], inplace=True, errors="ignore")
            older_kelp_info.rename(inplace=True, columns={"area": f"area {older_info_file.stem.replace('area_', '')}"})
            area_columns.append(f"area {older_info_file.stem.replace('area_', '')}")
            kelp_info = pandas.merge(kelp_info, older_kelp_info, on=['date'], how="outer")
        
        plot = plotly.express.line(kelp_info, x = 'date', y = area_columns, markers=True, text="date")
        event = streamlit.plotly_chart(plot, on_select="rerun")
        
    selection = event["selection"]["point_indices"]
    
    if len(selection):
        csv_file_path = pathlib.Path(kelp_info["file"].iloc[selection[0]])
        data_file = remote_raster_path / csv_file_path.parent.name / csv_file_path.name
        streamlit.subheader(f"Plot {kelp_info["date"].iloc[selection[0]]} with {kelp_info["ocean cloud percentage"].iloc[selection[0]]:.2f}% ocean cloud")
        streamlit.caption("May take time to load...")
        
        data = rioxarray.rioxarray.open_rasterio(data_file, chunks=True).squeeze("band", drop=True, errors="ignore")
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
                                         )


if __name__ == '__main__':
    main()