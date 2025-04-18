import pystac_client
import pystac
import odc.stac
import rioxarray
import xarray
import pathlib
import pandas
import geopandas
import shapely
import numpy
import dotenv
import datetime
import planetary_computer
import os
import copy
import geoapis.vector
import utils

def main():
    """ Create site datasets.
    """
    
    test_sites = utils.create_test_sites(distance_offshore = 3_000)
    test_sites_wsg = test_sites.to_crs(utils.CRS_WSG)
    land = geopandas.read_file(utils.DATA_PATH / "vectors" / "main_islands.gpkg")

    catalogue = {"url": "https://planetarycomputer.microsoft.com/api/stac/v1",
                 "collections": ["sentinel-2-l2a"]}

    date_format = "%Y-%m-%d"
    raster_defaults = {"resolution": 10, "nodata": 0, "dtype": "uint16"}

    bands = list(utils.SENTINEL_2B_BAND_INFO.keys()); bands.append("SCL") # bands = 
    
    filter_cloud_percentage = 30
    max_ocean_cloud_percentage = 5
    
    anomaly_detection_factor = 20
    
    # Second pass - remove beds with less than the min_pixls, then buffer outward by the specified number of pixels
    buffer = 10; min_pixels = 10 # 5

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 
    
    for site_index, row in test_sites_wsg.iterrows(): # [test_sites_wsg["name"]=="matau"]
        site_name = row['name']
        
        print(f"Test site: {site_name}") 
        raster_path = utils.DATA_PATH / "rasters" / "test_sites" / f"{site_name}"
        remote_raster_path = pathlib.Path("/nesi/nobackup/niwa03660/ZBD2023_outputs") / f"{site_name}"
        raster_path.mkdir(parents=True, exist_ok=True)
        remote_raster_path.mkdir(parents=True, exist_ok=True)
    
        # Geometry of AOI - convex hull to allow search
        site_bbox = row.geometry.bounds # shapely.box(*row.geometry.bounds) #.to_crs(utils.CRS_WSG).iloc[0].geometr
        filters = {"eo:cloud_cover":{"lt":filter_cloud_percentage}} 

        # Start from year of failure if already partial results
        if (raster_path / "info.csv").exists():
            kelp_info = pandas.read_csv(raster_path / "info.csv")
            kelp_info = kelp_info[['date', 'file', 'area', 'ocean cloud percentage']]
        else:
            print(f"No data for site {site}. Skipping RGB production")
            continue

        for index, row in kelp_info.iterrows():
            date_YYMMDD = row['date'] # datetime.datetime.strptime(row['date'], '%Y-%m')
            print(f"\tCreate RGB for date: {date_YYMMDD}")
            # run pystac client search to see available dataset
            search = client.search(
                collections=catalogue["collections"], bbox=site_bbox, datetime=date_YYMMDD, query=filters
            )

            if len(search.item_collection()) == 0:
                print(f"No valid data for for date: {date_YYMMDD}")
                continue # Nothing meeting criteria for this month

            data = odc.stac.load(search.items(), bbox=site_bbox, bands=bands,  chunks={}, groupby="solar_day", 
                                resolution = raster_defaults["resolution"], dtype=raster_defaults["dtype"], nodata=raster_defaults["nodata"])
            roi = test_sites.to_crs(data["SCL"].rio.crs).loc[[site_index]]

            data_file = pathlib.Path(row["file"])
            rgb = data[["B04", "B03","B02"]].to_array("rgb", name="all images").rio.clip(roi.geometry)
            utils.update_raster_defaults(rgb)
            rgb = rgb.squeeze('time', drop=True)
            rgb.to_netcdf(data_file.parent / f'rgb_{data_file.name}', format="NETCDF4", engine="netcdf4", encoding={"all images": {"zlib": True, "complevel": 5, "grid_mapping": rgb.encoding["grid_mapping"]}})
            
            # to display
            # rgb = rioxarray.rioxarray.open_rasterio(data_file.parent / f'rgb_{data_file.name}', chunks=True).drop_vars ("band")
            # rgb.odc.add_to(map=folium_map, name="Satellite RBG")

if __name__ == '__main__':
    main()
