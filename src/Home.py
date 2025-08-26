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
    
     # Dashboard
    streamlit.title("NIWA internal Kelp Dashboard")
    
    streamlit.subheader("This is used for data QA/QC during the data update process")

    streamlit.markdown("""
        **Dashboard Overview:** This is integrated within the kelp-detection GitHub respositry, which has several scripts for updating and analyzing kelp data. Steps for data update and QA/QC are:
        1. Create kelp detections all dates - run `create_data_NZ_wide_two_pass.py` script
        2. Review 'All Dates' tab in the Dashboard for any anomalous dates to ignore. (click the `>` if this tab is not visible.)
        3. Update `sites_dates_to_ignore.json` with dates to ignore. <span style="color: red;">IMPORTANT</span> - date format is YYYY-MM-DD
        4. Create kelp detections quarterly averaged - run `create_data_NZ_wide_two_pass_quarterly.py`
        5. Create site-wide presence-absence maps & summary info - run `summarise_NZ_wide_info.py` and `summarise_NZ_wide_info_quarterly.py`
        6. Review in the dashboard. If suitable for deployment to the external production website talk to Rose or Craig about updating in the deployment.                  
                       """)
    
if __name__ == '__main__':
    main()
