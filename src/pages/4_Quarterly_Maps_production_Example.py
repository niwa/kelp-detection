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
import plotly.subplots
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
    remote_raster_path = pathlib.Path("/nesi/nobackup/niwa03660/ZBD2023_outputs/test_sites_quarterly")
    
    streamlit.button("Re-run")
    streamlit.title('Kelp Demo - click area plot to select raster display')
    
    test_sites = geopandas.read_file(data_path / "vectors" / "test_sites_offshore_3km.gpkg")
    
    location = streamlit.selectbox("Select tile to display", (test_sites["name"]), index=0,)
    
    # Define the region
    raster_path = data_path / "rasters" / "test_sites_quarterly" / f"{location}"
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    kelp_total_extents = geopandas.read_file(raster_path / "presence_absence_map.gpkg")
    
    streamlit.subheader("Table and Plot of areas")
    if (raster_path / "info_quarterly.csv").exists():
        kelp_info = pandas.read_csv(raster_path / "info_quarterly.csv"); kelp_info.drop(columns=["Unnamed: 0"], inplace=True, errors="ignore")
    col1, col2 = streamlit.columns([1, 5])
    with col1:
        if (raster_path / "info_quarterly.csv").exists():
            streamlit.dataframe(kelp_info[["date", "area", "dates considered", "max coverage date", "proportion of max coverage"]])
        else:
            streamlit.markdown("No info_quarterly.csv.")
    with col2:
        #streamlit.plotly_chart(plotly.express.line(kelp_info, x = 'date', y = ["dates considered"], markers=True, text="date"))
        if (raster_path / "info_quarterly.csv").exists():
            figure = plotly.subplots.make_subplots(specs=[[{"secondary_y": True}]]) # plotly.graph_objects.Figure()
            
            figure.add_trace(plotly.graph_objects.Scatter(x=kelp_info["date"], y=kelp_info["area"] * 100, mode="lines+markers", marker={'color':'blue'}, name="Area of Kelp Coverage [m^2]" ), secondary_y=True)
            figure.add_trace(plotly.graph_objects.Scatter(x=kelp_info["date"], y=kelp_info["proportion of max coverage"] * 100, mode="lines+markers", marker={'color':'red'}, name="Proportion of Max Coverage [%]" ), secondary_y=False)
            
            figure.update_layout(title="Proportion of max coverage by date across algorithm runs",
                                 xaxis_title="date", yaxis2_title="Area [m^2]",
                                 yaxis_title="Max Coverage Proportion [%]", yaxis={'range': [0, 50]},
                                 legend_title="Algorithms", template="plotly_white" )
            event = streamlit.plotly_chart(figure, on_select="rerun")
        
    selection = event["selection"]["point_indices"]
    
    if len(selection) and (raster_path / "info_quarterly.csv").exists():
        csv_file_path = pathlib.Path(kelp_info["file"].iloc[selection[0]])
        streamlit.subheader(f"Plot quarter {kelp_info["date"].iloc[selection[0]]} calculated from dates {kelp_info["dates considered"].iloc[selection[0]]}.")
        streamlit.caption("May take time to load...")
        
        
        kelp_polygons = geopandas.read_file(csv_file_path)
        
        folium_map = folium.Map()
        #land.explore(m=folium_map, color="blue", style_kwds={"fillOpacity": 0}, name="land")
        kelp_total_extents.explore(m=folium_map, color="blue", style_kwds={"fillOpacity": 0}, name="Satellite record Kelp Extents")
        kelp_polygons.explore(m=folium_map, color="magenta", style_kwds={"fillOpacity": 0}, name="Quarter Kelp Extents")
        
        rgb_file = remote_raster_path / location / f"rgb_{kelp_info['max coverage date'].iloc[selection[0]]}.tif"
        rgb = rioxarray.rioxarray.open_rasterio(rgb_file, chunks=True)#.drop_vars("band")
        #print(rgb)
        #rgb.odc.add_to(map=folium_map, name="Satellite RBG")
        
        folium.LayerControl().add_to(folium_map)
        st_map =  streamlit_folium.st_folium(folium_map, width=900) 


if __name__ == '__main__':
    main()