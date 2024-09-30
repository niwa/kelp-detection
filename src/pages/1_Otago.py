import numpy
import pandas
import streamlit
import streamlit_folium
import argparse
import pages.scripts.colourmaps
import datetime
import folium
import rioxarray
import odc.stac
import branca

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
    date_format = "%Y-%m-%d"
    location = "waikouaiti"
    
    streamlit.button("Re-run")
    streamlit.title('Otago Kelp')
    
    # Define the region
    land = streamlit.session_state["land"]
    data_path = streamlit.session_state["data_path"]
    
    # Select the date to display
    files = list((data_path / "rasters" / location).glob(f"*.tif"))
    dates = [datetime.datetime.strptime(file.stem.strip("kelp_"), date_format) for file in files]
    date = streamlit.select_slider("Select date?", options=dates) #, format="MM/DD/YY",)
    kelp_file = data_path / "rasters" / location / f"kelp_{date.strftime(date_format)}.tif"

    streamlit.subheader("Map View")
    folium_map = folium.Map()
    land.explore(m=folium_map)

    kelp_display = rioxarray.rioxarray.open_rasterio(kelp_file, chunks=True).squeeze( "band", drop=True)
    kelp_display.odc.add_to(folium_map, opacity=0.75, cmap="inferno", vmin=0, vmax=1) # viridis
    folium_map.fit_bounds(kelp_display.odc.map_bounds())
    
    colormap = branca.colormap.linear.inferno.scale(0, 1)
    colormap.caption = 'Kelp Index'
    colormap.add_to(folium_map)

    st_data =  streamlit_folium.st_folium(folium_map, width=900)  

    streamlit.subheader("Table View")
    #streamlit.dataframe(regions[["name", "id"]])


if __name__ == '__main__':
    args = parse_args()
    main(region_name=args.region_name)