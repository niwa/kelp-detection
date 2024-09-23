import numpy
import pandas
import streamlit
import streamlit_folium
import argparse
import pages.scripts.colourmaps
from datetime import datetime

def parse_args():
    """Expect a command line argument of the form:
    '--tiles_subpath subpath_to_tiles_folder'"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--region_name",
        metavar="str",
        default=r"Otago",
        action="store",
        help="The region name. Used to reference the region in the regions file.",
    )

    return parser.parse_args()



def main(region_name: str):
    """ Create / Update the geofabric summary information and display in a dashboard.
    """
    
    streamlit.set_page_config(
        page_title="Otago",
        page_icon="🌏",
        layout="wide",
    )
    display_size = 700
    
    streamlit.button("Re-run")
    streamlit.title('Otago Kelp')
    
    # Define the region
    regions = streamlit.session_state["regions"]
    otago = regions[regions["name"]==region_name]
    
    # Select the date to display
    
    date = streamlit.slider("Select date?", value=datetime(2020, 1, 1), format="MM/DD/YY",)

    streamlit.subheader("Map View")
    folium_map = otago.explore(column="name")#, cmap=pages.scripts.colourmaps.get_colourmap(otago, "name"))
    st_data =  streamlit_folium.st_folium(folium_map, width=display_size, key=0)

    streamlit.subheader("Table View")
    streamlit.dataframe(regions[["name", "id"]])


if __name__ == '__main__':
    args = parse_args()
    main(region_name=args.region_name)