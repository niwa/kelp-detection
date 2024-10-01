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



def main():
    """ Create / Update the geofabric summary information and display in a dashboard.
    """
    
    streamlit.set_page_config(
        page_title="Otago",
        page_icon="🌏",
        layout="wide",
    )
    display_size = 700
    date_format = "%Y-%m-%d"
    location = "waikouaiti"
    
    streamlit.button("Re-run")
    streamlit.title('Otago Kelp')
    
    # Define the region
    data_path = pathlib.Path.cwd() / "data"
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    #land = streamlit.session_state["land"]
    #data_path = streamlit.session_state["data_path"]
    
    streamlit.subheader("Table and Plot of areas")
    kelp_info = pandas.read_csv(data_path / "rasters" / location / "info.csv")
    col1, col2 = streamlit.columns([1, 7])
    with col1:
        streamlit.dataframe(kelp_info[["date", "area"]])
    with col2:
        plot = plotly.express.line(kelp_info, x = 'date', y = ['area'])
        streamlit.plotly_chart(plot)
    
    # Select the date to display
    files = sorted((data_path / "rasters" / location).glob(f"*.tif"))
    dates = [datetime.datetime.strptime(file.stem.strip("kelp_"), date_format) for file in files]
    date = streamlit.select_slider("Select date?", options=dates) #, format="MM/DD/YY",)
    kelp_file = data_path / "rasters" / location / f"kelp_{date.strftime(date_format)}.tif"

    streamlit.subheader("Map View")
    folium_map = folium.Map()
    land.explore(m=folium_map)

    kelp_display = rioxarray.rioxarray.open_rasterio(kelp_file, chunks=True).squeeze( "band", drop=True)
    # In future consider https://geoviews.org or saving as a png and loading...
    kelp_display.odc.add_to(folium_map, opacity=0.75, cmap="inferno", vmin=0, vmax=1) # viridis
    folium_map.fit_bounds(kelp_display.odc.map_bounds())
    
    colormap = branca.colormap.linear.inferno.scale(0, 1)
    colormap.caption = 'Kelp Index'
    colormap.add_to(folium_map)

    st_data =  streamlit_folium.st_folium(folium_map, width=900)  


if __name__ == '__main__':
    main()