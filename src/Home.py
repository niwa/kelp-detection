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
    
    streamlit.subheader("Pelase select one of the tabs. See the left hand (cick the `>` if not visible. 'Test Sites' is currently most in use.")
    
if __name__ == '__main__':
    main()
