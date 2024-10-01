import numpy
import streamlit
import argparse
import pathlib
import geopandas
import streamlit_folium
import pages.scripts.colourmaps
import matplotlib.colors


def geopandas_bounds_to_plot(dataframe, crs=4326):
    """ Changing bounding box representation to leaflet notation ``(lon1, lat1, lon2, lat2) -> ((lat1, lon1), (lat2, lon2))`` """
    x1, y1, x2, y2 = dataframe.to_crs(crs).total_bounds
    return ((y1, x1), (y2, x2))

def main():
    """ Create / Update the catchment summary file and display in a dashboard.
    """
    
    import streamlit as st

    streamlit.set_page_config(
        page_title="Kelp Dashboard",
        page_icon="🌊",
    )
 
    data_path = pathlib.Path.cwd() / "data"
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    regions = geopandas.read_file(data_path / "vectors" / "regions.gpkg")
    region_names = ["Otago", "Southland", "Canterbury", "Westland"]
    
    
    if "land" not in streamlit.session_state:
        streamlit.session_state["land"] = land
    if "data_path" not in streamlit.session_state:
        streamlit.session_state["data_path"] = data_path
    
     # Dashboard
    streamlit.title("Kelp Dashboard")
    
    region = streamlit.selectbox("Select Region?", region_names, index=3)
        
    streamlit.subheader("Maps - Regions")
    streamlit.caption("A caption")
    folium_map = regions.explore(column="name")#, cmap=pages.scripts.colourmaps.get_colourmap(regions, "name", region))
    folium_map.fit_bounds(geopandas_bounds_to_plot(regions))
    st_data =  streamlit_folium.st_folium(folium_map, width=900)  

    streamlit.subheader("Tables - Overall Status of Each Stage")
    streamlit.caption("Example of different regions - could give a summary statistic")
    streamlit.dataframe(regions[["name", "id"]])

if __name__ == '__main__':
    main()
