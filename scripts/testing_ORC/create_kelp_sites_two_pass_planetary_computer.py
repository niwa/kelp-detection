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
import geoapis.vector
import utils

def main():
    """ Create site datasets.
    """
    
    test_sites_andra_and_leigh_wsg_84 = utils.create_test_sites(distance_offshore = 3_000)
    #test_sites_wsg = test_sites.to_crs(utils.CRS_WSG)
    test_sites_wsg = test_sites_andra_and_leigh_wsg_84.to_crs(utils.CRS_WSG)
    land = geopandas.read_file(utils.DATA_PATH / "vectors" / "main_islands.gpkg")

    catalogue = {"url": "https://planetarycomputer.microsoft.com/api/stac/v1",
                 "collections": ["sentinel-2-l2a"]}

    date_format = "%Y-%m-%d"

    bands = list(utils.SENTINEL_2B_BAND_INFO.keys()); bands.append("SCL") # bands = ["red", "green", "blue", "nir", "SCL", "swir16", "B05", "B8A"]
    raster_defaults = {"resolution": 10, "nodata": 0, "dtype": "uint16"}
    thresholds = {"min_ndvi": 0.213, "max_ndwi": 0.1, "max_ndwi2": -0.2,} # "max_ndvi": 0.7, "min_ndvri": 0.03, "min_ndvi": 0.03
    
    filter_cloud_percentage = 30
    max_ocean_cloud_percentage = 5

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 

    for site_index, row in test_sites_wsg.iterrows(): # [test_sites_wsg["name"]=="matau"]

        site_name = row['name']
        
        print(f"Test site: {site_name}") 
        raster_path = utils.DATA_PATH / "rasters" / "test_sites" / f"{site_name}"
        #remote_raster_path = raster_path
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
            max_date = datetime.datetime.strptime(kelp_info["date"].max(), '%Y-%m-%d')
            kelp_info = kelp_info.to_dict(orient='list')
        else:
            kelp_info = {"date": [], "file": [], "area": [], "ocean cloud percentage": []}
            max_date = datetime.datetime.strptime("2015-01-31", '%Y-%m-%d')

        years = list(range(2016, 2025)) # 2016, 2025
        for year in years:
            months = [f"{year}-{str(month).zfill(2)}" for month in list(range(1, 13))]

            for month_YYMM in months:
                if max_date > datetime.datetime.strptime(month_YYMM, '%Y-%m'):
                    print(f"\tSkipping {month_YYMM} as run previously. Delete if you want a rerun")
                    continue

                print(f"\tCheck month: {month_YYMM}")
                # run pystac client search to see available dataset
                search = client.search(
                    collections=catalogue["collections"], bbox=site_bbox, datetime=month_YYMM, query=filters
                )

                if len(search.item_collection()) == 0:
                    continue # Nothing meeting criteria for this month

                data = odc.stac.load(search.items(), bbox=site_bbox, bands=bands,  chunks={}, groupby="solar_day", 
                                    resolution = raster_defaults["resolution"], dtype=raster_defaults["dtype"], nodata=raster_defaults["nodata"])

                data = utils.harmonize_post_2022(data)
                roi = test_sites_andra_and_leigh_wsg_84.to_crs(data["SCL"].rio.crs).loc[[site_index]]


                # remove if no data
                (data, ocean_cloud_percentage) = utils.screen_by_SCL_in_ROI(data, roi, max_ocean_cloud_percentage)

                if len(data.time) == 0:
                    print("\t\tNone with suitable cloud percentage.")
                    continue

                # Save out RGB
                '''rgb = data[["red", "green","blue"]].to_array("rgb", name="all images")
                utils.update_raster_defaults(rgb)
                rgb.rio.clip.to_netcdf(raster_path / f'rgb_{month_YYMM}.nc', format="NETCDF4", engine="netcdf4", encoding={"all images": {"zlib": True, "complevel": 5, "grid_mapping": rgb.encoding["grid_mapping"]}})'''

                # Convert to floats before calcualtions
                for key in data.data_vars:
                    if key == "SCL": 
                        continue
                    data[key] = data[key].astype("float32").where(data[key] != 0, numpy.nan)
                utils.update_raster_defaults(data)

                # Calculate Kelp from thresholds
                data = utils.threshold_kelp(data, thresholds, roi)

                # Save each separately
                for index in range(len(data["kelp"].time)):
                    kelp = data["kelp"].isel(time=index).load()
                    #kelp = kelp.rio.clip(roi.geometry.values)
                    filename = remote_raster_path / f'data_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.nc'

                    kelp_info["area"].append(abs(int(kelp.notnull().sum() * kelp.x.resolution * kelp.y.resolution)))
                    kelp_info["file"].append(filename)
                    kelp_info["date"].append(pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format))
                    kelp_info["ocean cloud percentage"].append(ocean_cloud_percentage[index])
                    #kelp.rio.to_raster(raster_path / f'kelp_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif', compress="deflate", driver="COG")  # missing min and max values when viewed in QGIS
                    
                    #data["SCL"].isel(time=index).rio.to_raster(raster_path / f'scl_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif', compress="deflate", driver="COG")
                    encoding = {}
                    for key in data.data_vars:
                        encoding[key] =  {"zlib": True, "complevel": 9, "grid_mapping": data[key].encoding["grid_mapping"]}
                    data.isel(time=index).to_netcdf(filename, format="NETCDF4", engine="netcdf4", encoding=encoding)
                pandas.DataFrame.from_dict(kelp_info, orient='columns').to_csv(raster_path / "info.csv", index=False)
                pandas.DataFrame.from_dict(kelp_info, orient='columns').to_csv(remote_raster_path / "info.csv", index=False)

        # Save results
        kelp_info = pandas.DataFrame.from_dict(kelp_info, orient='columns')
        kelp_info.to_csv(raster_path / "info.csv", index=False)
        kelp_info.to_csv(remote_raster_path / "info.csv", index=False)

if __name__ == '__main__':
    main()