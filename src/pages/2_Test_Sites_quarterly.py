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
    remote_raster_path = pathlib.Path("/nesi/nobackup/niwa03660/ZBD2023_outputs")
    
    streamlit.button("Re-run")
    streamlit.title('Kelp Demo - click area plot to select raster display')
    
    test_sites = geopandas.read_file(data_path / "vectors" / "test_sites_offshore_3km.gpkg")
    
    location = streamlit.selectbox("Select tile to display", (test_sites["name"]), index=0,)
    
    # Define the region
    raster_path = data_path / "rasters" / "test_sites_quarterly" / f"{location}"
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    
    if 'xy' not in streamlit.session_state:
        streamlit.session_state['xy'] = []
        streamlit.session_state['prev_lat'], streamlit.session_state['prev_lon'] = numpy.nan, numpy.nan
        streamlit.session_state['spectra'] = None
    
    streamlit.subheader("Table and Plot of areas")
    if (raster_path / "info_quarterly.csv").exists():
        kelp_info = pandas.read_csv(raster_path / "info_quarterly.csv"); kelp_info.drop(columns=["Unnamed: 0"], inplace=True, errors="ignore")
    col1, col2 = streamlit.columns([1, 5])
    with col1:
        if (raster_path / "info_quarterly.csv").exists():
            streamlit.dataframe(kelp_info[["date", "area", "dates considered"]])
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
            figure.add_trace(plotly.graph_objects.Scatter(x=kelp_info["date"], y=kelp_info["area"], mode="lines+markers", name="Latest run" ))
        figure.update_layout(title="Area's by date across algorithm runs", xaxis_title="date", yaxis_title="Area [m^2]", legend_title="Algorithms", template="plotly_white" )
        event = streamlit.plotly_chart(figure, on_select="rerun")
        
    selection = event["selection"]["point_indices"]
    
    if len(selection) and (raster_path / "info_quarterly.csv").exists():
        csv_file_path = pathlib.Path(kelp_info["file"].iloc[selection[0]])
        streamlit.subheader(f"Plot quarter {kelp_info["date"].iloc[selection[0]]} calculated from dates {kelp_info["dates considered"].iloc[selection[0]]}.")
        streamlit.caption("May take time to load...")
        
        #data_file = remote_raster_path / f"{csv_file_path.parent.name}_quarterly" / csv_file_path.name
        kelp_polygons = geopandas.read_file(csv_file_path)

        folium_map = folium.Map()

        land.explore(m=folium_map, color="blue", style_kwds={"fillOpacity": 0}, name="land")
        kelp_polygons.explore(m=folium_map, color="magenta", style_kwds={"fillOpacity": 0}, name="kelp polygon")


        folium.LayerControl().add_to(folium_map)
        st_map =  streamlit_folium.st_folium(folium_map, width=900) 


if __name__ == '__main__':
    main()