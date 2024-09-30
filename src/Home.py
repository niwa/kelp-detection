import numpy
import streamlit
import argparse
import pathlib
import geopandas
import streamlit_folium
import pages.scripts.colourmaps
import matplotlib.colors

def parse_args():
    """Expect a command line argument of the form:
    '--output_path path_to_output_folder'"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--data_path",
        metavar="str",
        default=r"C:\Local\repos\kelp-dashboard-demo\data", # r"/home/pearsonra/repos/kelp-dashboard-demo/data/vectors/regions.gpkg",
        action="store",
        help="The regions maps - The path to the regions geopackage file.",
    )

    return parser.parse_args()


def main(data_path: str):
    """ Create / Update the catchment summary file and display in a dashboard.
    """
    
    import streamlit as st

    streamlit.set_page_config(
        page_title="Kelp Dashboard",
        page_icon="🌊",
    )
 
    data_path = pathlib.Path(data_path)
    land = geopandas.read_file(data_path / "vectors" / "main_islands.geojson")
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
    st_data =  streamlit_folium.st_folium(folium_map, width=900)  

    streamlit.subheader("Tables - Overall Status of Each Stage")
    streamlit.caption("Example of different regions - could give a summary statistic")
    streamlit.dataframe(regions[["name", "id"]])

if __name__ == '__main__':
    args = parse_args()
    main(data_path=args.data_path)
