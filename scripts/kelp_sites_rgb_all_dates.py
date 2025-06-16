import pystac_client
import leafmap
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
import json

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

    bands = list(utils.SENTINEL_2B_BAND_INFO.keys()); bands.append("SCL") 
    
    # Read in date to ignore file
    if (pathlib.Path.cwd() / "sites_dates_to_ignore.json").exists:
        with open(pathlib.Path.cwd() / "sites_dates_to_ignore.json", "r") as file_pointer:
            site_dates_to_ignore = json.load(file_pointer)
    else:
        site_dates_to_ignore = {}
    
    filter_cloud_percentage = 30
    rgb_bands = utils.get_band_names_from_common(["red", "green", "blue"])

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 

    for site_index, row in test_sites_wsg.iterrows(): 
        site_name = row['name']
        
        dates_to_ignore = site_dates_to_ignore[site_name] if site_name in site_dates_to_ignore.keys() else {}
        
        print(f"Test site: {site_name}") 
        raster_path = utils.DATA_PATH / "rasters" / "test_sites" / f"{site_name}"
        #remote_raster_path = pathlib.Path("/nesi/nobackup/niwa03660/ZBD2023_outputs/test_sites") / f"{site_name}"
        raster_path.mkdir(parents=True, exist_ok=True)
        #remote_raster_path.mkdir(parents=True, exist_ok=True)
    
        # Geometry of AOI
        site_bbox = row.geometry.bounds
        filters = {"eo:cloud_cover":{"lt":filter_cloud_percentage}} 
        
        # Check if any results to post process
        if (raster_path / "info.csv").exists():
            kelp_info = pandas.read_csv(raster_path / "info.csv")
            #kelp_info = kelp_info[['date','file','area','dates considered', 'max coverage date']]
        else:
            print(f"No data for site {site_name}. Skipping post processing.")
            continue

        # Save out overall presence absense map
        if not (raster_path / "presence_absence_map.gpkg").exists():
            print(f"\tWriting out presence absence for: {site_name}")
            kelp_polygons = []; file_names = []
            for file_name in kelp_info["file"].tolist():
                file_name = pathlib.Path(file_name)
                file_name = raster_path / file_name.name
                date = file_name.stem.replace('data_','')
                if date in dates_to_ignore:
                    file_names.append("")
                    continue
                
                kelp = rioxarray.rioxarray.open_rasterio(file_name, chunks=True).squeeze("band", drop=True)["kelp"]
                kelp_polygon = utils.polygon_from_raster(kelp).dissolve()
                kelp_polygon.to_file(file_name.parent / f"{date}_kelp.gpkg")
                file_names.append(file_name.parent / f"{date}_kelp.gpkg")
                kelp_polygons.append(kelp_polygon.to_crs(utils.CRS))
                '''#kelp_polygons.append(
                    geopandas.read_file(raster_path / pathlib.Path(file_name).name).to_crs(utils.CRS)
                )'''
            kelp_polygons = pandas.concat(kelp_polygons).dissolve()
            kelp_info["proportion of max coverage"] = kelp_info["area"] / kelp_polygons.area.sum()
            kelp_info["file"] = file_names
            kelp_info.to_csv(raster_path / "info.csv", index=False)
            kelp_polygons.to_file(raster_path / "presence_absence_map.gpkg")
        else:
            print(f"\tSkip presence absence for: {site_name} - already exists")
            
        # Save out RGB if not already produced
        tile_ids = []
        percentages_2 = []
        percentages_98 = []
        if "Satellite Tile IDs" not in kelp_info.columns:
            for index, row in kelp_info.iterrows():
                date_YYMMDD = row['date'] 

                # run pystac client search to see available dataset
                print(f"\tGet Tile(s) ID: {date_YYMMDD}")
                search = client.search(
                    collections=catalogue["collections"], bbox=site_bbox, datetime=date_YYMMDD, query=filters
                )

                if len(search.item_collection()) == 0:
                    print(f"No valid data for for date: {date_YYMMDD}")
                    continue # Nothing meeting criteria for this month
                
                tile_id = ""; percentage_2 = ""; percentage_98 = ""
                for item in search.items():
                    tile_id += f"{item.id}, "
                    stats=leafmap.stac_stats(collection=catalogue["collections"][0], item=item.id, titiler_endpoint="pc", assets=rgb_bands)
                    percentage_2_i = []; percentage_98_i = []
                    for key, value in stats.items():
                        percentage_2_i.append(value['percentile_2']); percentage_98_i.append(value['percentile_98'])
                    percentage_2_i = numpy.array(percentage_2_i).mean(); percentage_98_i = numpy.array(percentage_98_i).mean()
                    percentage_2 += f"{round(percentage_2_i)}, "; percentage_98 += f"{round(percentage_98_i)}, "
                tile_ids.append(tile_id); percentages_2.append(percentage_2); percentages_98.append(percentage_98); 

            # add tile id's and display ranges to the CSV
            kelp_info["Satellite Tile IDs"] = tile_ids
            kelp_info["Percentile 2"] = percentages_2
            kelp_info["Percentile 98"] = percentages_98
            kelp_info.to_csv(raster_path / "info.csv")


if __name__ == '__main__':
    main()
