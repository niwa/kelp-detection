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

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 
    
    for site_index, row in test_sites_wsg.iterrows(): # [test_sites_wsg["name"]=="matau"]
        site_name = row['name']
        
        print(f"Test site: {site_name}") 
        raster_path = utils.DATA_PATH / "rasters" / "test_sites_quarterly" / f"{site_name}"
        remote_raster_path = pathlib.Path("/nesi/nobackup/niwa03660/ZBD2023_outputs/test_sites_quarterly") / f"{site_name}"
        raster_path.mkdir(parents=True, exist_ok=True)
        remote_raster_path.mkdir(parents=True, exist_ok=True)
    
        # Geometry of AOI
        site_bbox = row.geometry.bounds
        filters = {"eo:cloud_cover":{"lt":filter_cloud_percentage}} 
        
        # Check if any results to post process
        if (raster_path / "info_quarterly.csv").exists():
            kelp_info = pandas.read_csv(raster_path / "info_quarterly.csv")
            #kelp_info = kelp_info[['date', 'file', 'max coverage date']]
        else:
            print(f"No data for site {site_name}. Skipping post processing.")
            continue

        # Save out overall presence absense map
        if not (raster_path / "presence_absence_map.gpkg").exists():
            print(f"\tWriting out presence absence for: {site_name}")
            kelp_polygons = []
            for file_name in kelp_info["file"].tolist():
                kelp_polygons.append(
                    geopandas.read_file(raster_path / pathlib.Path(file_name).name).to_crs(utils.CRS)
                )
            kelp_polygons = pandas.concat(kelp_polygons).dissolve()
            #kelp_info.assign({"proportion of max coverage": kelp_info["area"] / kelp_polygons.area.sum()})
            kelp_info["proportion of max coverage"] = kelp_info["area"] / kelp_polygons.area.sum()
            kelp_info.to_csv(raster_path / "info_quarterly.csv", index=False)
            kelp_polygons.to_file(raster_path / "presence_absence_map.gpkg")
        else:
            print(f"\tSkip presence absence for: {site_name} - already exists")
            
        # Save out RGB if not already produced
        tile_ids = []
        if len(list(remote_raster_path.glob('rgb_*.nc'))) < len(kelp_info) or "Satellite Tile IDs" not in kelp_info.columns:
            for index, row in kelp_info.iterrows():
                date_YYMMDD = row['max coverage date'] 
                filename = remote_raster_path / f'rgb_{date_YYMMDD}.nc'

                # run pystac client search to see available dataset
                print(f"\tGet Tile(s) ID: {date_YYMMDD}")
                search = client.search(
                    collections=catalogue["collections"], bbox=site_bbox, datetime=date_YYMMDD, query=filters
                )

                if len(search.item_collection()) == 0:
                    print(f"No valid data for for date: {date_YYMMDD}")
                    continue # Nothing meeting criteria for this month
                
                tile_id = ""
                for item in search.items():
                    tile_id += f"{item.id}, "
                tile_ids.append(tile_id)

                if filename.exists():
                    # Both already run - skip to the next date
                    print(f"\tSkip RGB for date: {date_YYMMDD} - already exists")
                    continue
                print(f"\tCreate RGB for date: {date_YYMMDD}")
                data = odc.stac.load(search.items(), bbox=site_bbox, bands=bands,  chunks={}, groupby="solar_day", 
                                    resolution = raster_defaults["resolution"], dtype=raster_defaults["dtype"], nodata=raster_defaults["nodata"])
                roi = test_sites.to_crs(data["SCL"].rio.crs).loc[[site_index]]

                rgb_bands = utils.get_band_names_from_common(["red", "green", "blue"])
                rgb = utils.normalise_rgb(data.isel(time=0), rgb_bands)
                utils.update_raster_defaults(rgb)
                rgb = rgb.to_array("rgb", name="Satellite RGB").rio.clip(roi.geometry)
                encoding = {"Satellite RGB": {"zlib": True, "complevel": 9, "grid_mapping": data[rgb_bands[0]].encoding["grid_mapping"]}}
                rgb.load()
                rgb.to_netcdf(filename, format="NETCDF4", engine="netcdf4", encoding=encoding)
                #rgb.odc.write_cog(filename.parent / f"{filename.stem}.tif")
            # add tile id's to the CSV
            kelp_info["Satellite Tile IDs"] = tile_ids
            kelp_info.to_csv(raster_path / "info_quarterly.csv", index=False)


if __name__ == '__main__':
    main()
