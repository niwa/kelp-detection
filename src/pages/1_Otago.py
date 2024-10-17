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
import holoviews
import hvplot.xarray



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
    location = "1607" # "waikouaiti"
    
    streamlit.button("Re-run")
    streamlit.title('Kelp Demo - click area plot to select raster display')
    
    # Define the region
    data_path = pathlib.Path.cwd() / "data"
    raster_path = data_path / "rasters" / "tiles" / f"{location}_test"
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    #land = streamlit.session_state["land"]
    #data_path = streamlit.session_state["data_path"]
    
    streamlit.subheader("Table and Plot of areas")
    kelp_info = pandas.read_csv(raster_path / "info.csv")
    col1, col2 = streamlit.columns([1, 5])
    with col1:
        streamlit.dataframe(kelp_info[["date", "area", "ocean cloud percentage"]])
    with col2:
        #streamlit.plotly_chart(plotly.express.line(kelp_info, x = 'date', y = ["ocean cloud percentage"], markers=True, text="date"))
        plot = plotly.express.line(kelp_info, x = 'date', y = ["area"], markers=True, text="date")
        event = streamlit.plotly_chart(plot, on_select="rerun")
        
    selection = event["selection"]["point_indices"]
    
    if len(selection):
        kelp_file = pathlib.Path(kelp_info["file"].iloc[selection[0]])
        streamlit.subheader(f"Plot {kelp_info["date"].iloc[selection[0]]} with {kelp_info["ocean cloud percentage"].iloc[selection[0]]:.2f}% ocean cloud")
        streamlit.caption("May take time to load...")
        

        col1, col2 = streamlit.columns([1, 1])
        with col1:
            kelp_display = rioxarray.rioxarray.open_rasterio(kelp_file, chunks=True).squeeze( "band", drop=True)

            folium_map = folium.Map()
            land.explore(m=folium_map)
            # In future consider https://geoviews.org or saving as a png and loading...
            kelp_display.odc.add_to(folium_map, opacity=0.75, cmap="inferno", vmin=0, vmax=1) # viridis
            folium_map.fit_bounds(kelp_display.odc.map_bounds())

            colormap = branca.colormap.linear.inferno.scale(0, 1)
            colormap.caption = 'Kelp Index'
            colormap.add_to(folium_map)
            st_data =  streamlit_folium.st_folium(folium_map, width=900)  
        with col2:
            scl_file = kelp_file.parent / kelp_file.name.replace('kelp', 'scl')
            scl_display = rioxarray.rioxarray.open_rasterio(scl_file, chunks=True).squeeze( "band", drop=True)

            folium_map = folium.Map()
            land.explore(m=folium_map, style_kwds=dict(color="red", fill=False))
            # In future consider https://geoviews.org or saving as a png and loading...
            scl_display.odc.add_to(folium_map, opacity=0.75, vmin=0, vmax=11) # viridis
            folium_map.fit_bounds(scl_display.odc.map_bounds())

            colormap = branca.colormap.linear.viridis.scale(0, 11)
            colormap.caption = 'SCL Classification'
            colormap.add_to(folium_map)
            st_data =  streamlit_folium.st_folium(folium_map, width=900)  

        '''streamlit.subheader("hvplot Map View")
        hv_plot = kelp_display.hvplot.image(x='x', y='y', data_aspect=1, flip_yaxis=True, subplots=True, coastline="10m")
        streamlit.bokeh_chart(holoviews.render(hv_plot, backend='bokeh'))'''


if __name__ == '__main__':
    main()