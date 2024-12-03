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
    
    test_sites = utils.create_test_sites(distance_offshore = 3_000)
    test_sites_wsg = test_sites.to_crs(utils.CRS_WSG)
    land = geopandas.read_file(utils.DATA_PATH / "vectors" / "main_islands.gpkg")

    catalogue = {"url": "https://planetarycomputer.microsoft.com/api/stac/v1",
                 "collections": ["sentinel-2-l2a"]}

    date_format = "%Y-%m-%d"

    bands = ["red", "green", "blue", "nir", "SCL", "swir16", "B05", "B8A"]
    raster_defaults = {"resolution": 10, "nodata": 0, "dtype": "uint16"}
    scl_dict = {"no data": 0, "defective": 1, "cast shadow": 2, "cloud shadow": 3,
                "vegetation": 4, "not vegetated": 5, "water": 6, "unclassified": 7,
                "cloud medium probability": 8, "cloud high probability": 9,
                "thin cirrus": 10, "snow": 11}
    thresholds = {"min_ndvi": 0.03, "max_ndvi": 0.7, "max_ndwi": 0.1, "min_ndvri": 0.03, "max_ndwi2": -0.2,}
    
    filter_cloud_percentage = 30
    max_ocean_cloud_percentage = 10

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 
    
    for index, row in test_sites_wsg.iterrows():
        
        print(f"Test site: {row['name']}") 
        raster_path = utils.DATA_PATH / "rasters" / "test_sites" / f"{row['name']}"
        raster_path.mkdir(parents=True, exist_ok=True)
    
        # Geometry of AOI - convex hull to allow search
        geometry_query = shapely.box(*row.geometry.bounds) #row.geometry
        
        #geometry_query = tiles[tiles["name"]==name].to_crs(crs_wsg).iloc[0].geometry

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

        years = list(range(2016, 2025))
        for year in years:
            months = [f"{year}-{str(month).zfill(2)}" for month in list(range(1, 13))]

            for month_YYMM in months:

                if max_date > datetime.datetime.strptime(month_YYMM, '%Y-%m'):
                    print(f"\tSkipping {month_YYMM} as run previously. Delete if you want a rerun")
                    continue

                print(f"\tCheck month: {month_YYMM}")
                # run pystac client search to see available dataset
                search = client.search(
                    collections=catalogue["collections"], bbox=row.geometry.bounds, datetime=month_YYMM, query=filters
                ) # intersects=geometry_query, 

                if len(search.item_collection()) == 0:
                    continue # Nothing meeting criteria for this month

                data = odc.stac.load(search.items(), bbox=row.geometry.bounds, bands=bands,  chunks={}, groupby="solar_day", # geopolygon=geometry_query, 
                                    resolution = raster_defaults["resolution"], dtype=raster_defaults["dtype"], nodata=raster_defaults["nodata"])
                roi = test_sites.to_crs(data.rio.crs)[test_sites["name"] == row['name']]
                # remove if no data
                data["SCL"].load()
                #data["SCL"] = data["SCL"].rio.clip(land.to_crs(data["SCL"].rio.crs).geometry.values, invert=True)
                breakpoint()
                data["SCL"] = data["SCL"].rio.clip(roi.geometry.values)
                data["SCL"].rio.write_crs(data["SCL"].rio.crs, inplace=True);
                data = data.isel(time=(data["SCL"] != scl_dict["no data"]).any(dim=["x", "y"]))

                # remove if much cloud over ocean
                # Mask with 1s on ocean, 0 elsewhere
                ocean_mask = data["SCL"].isel(time=0).copy(deep=True)
                ocean_mask.data[:] = 1
                ocean_mask = ocean_mask.rio.clip(land.to_crs(ocean_mask.rio.crs).geometry.values, invert=True)
                # Mask by time - initially sums of cloud values then true / false by time if less than cloud threshold
                ocean_cloud_sum = (data["SCL"] == scl_dict["cloud high probability"]).sum(dim=["x", "y"]) 
                ocean_cloud_sum += (data["SCL"] == scl_dict["cloud medium probability"]).sum(dim=["x", "y"]) 
                ocean_cloud_sum += (data["SCL"] == scl_dict["cloud shadow"]).sum(dim=["x", "y"]) 
                ocean_cloud_sum += (data["SCL"] == scl_dict["cast shadow"]).sum(dim=["x", "y"]) 
                ocean_cloud_sum += (data["SCL"] == scl_dict["thin cirrus"]).sum(dim=["x", "y"])
                ocean_cloud_sum += (data["SCL"] == scl_dict["defective"]).sum(dim=["x", "y"])
                ocean_cloud_sum += (data["SCL"] == scl_dict["no data"]).sum(dim=["x", "y"]) - (ocean_mask == scl_dict["no data"]).sum(dim=["x", "y"])
                ocean_cloud_percentage = (ocean_cloud_sum / int(ocean_mask.sum())).data*100
                print(f"\t\tOcean cloud percentage {list(map('{:.2f}%'.format, ocean_cloud_percentage))}")
                cloud_mask_time = ocean_cloud_percentage < max_ocean_cloud_percentage
                data = data.isel(time=(cloud_mask_time))
                ocean_cloud_percentage = ocean_cloud_percentage[cloud_mask_time]

                if len(ocean_cloud_percentage) == 0:
                    print("\t\tNone with suitable cloud percentage.")
                    continue

                # Save out RGB
                '''rgb = data[["red", "green","blue"]].to_array("rgb", name="all images")
                utils.update_raster_defaults(rgb)
                rgb.rio.clip.to_netcdf(raster_path / f'rgb_{month_YYMM}.nc', format="NETCDF4", engine="netcdf4", encoding={"all images": {"zlib": True, "complevel": 2, "grid_mapping": rgb.encoding["grid_mapping"]}})'''

                # Convert to floats before calcualtions
                for key in data.data_vars:
                    if key == "SCL": 
                        continue
                    data[key] = data[key].astype("float32").where(data[key] != 0, numpy.nan)
                utils.update_raster_defaults(data)


                # Calculate NVDI and NVWI
                data["ndvi"] = (data.nir - data.red) / (data.nir + data.red)
                data["ndwi"] = (data.green - data.nir) / (data.green + data.nir)
                data["ndvri"] = (data.B05 - data.red) / (data.B05 + data.red);
                data["ndwi2"] = (data.swir16 + data.B05) / (data.swir16 - data.B05)
                utils.update_raster_defaults(data)

                # Calculate Kelp
                data["kelp"] = (data.nir - data.red) / (data.nir + data.red)
                data["kelp"] = data["kelp"].where(data["ndvi"].data > thresholds["min_ndvi"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["ndwi"].data < thresholds["max_ndwi"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["ndwi2"].data < thresholds["max_ndwi2"], numpy.nan)
                #data["kelp"] = data["kelp"].where(data["ndvi"].data < thresholds["max_ndvi"], numpy.nan)
                #data["kelp"] = data["kelp"].where(data["ndvri"].data > thresholds["min_ndvri"], numpy.nan)
                data["kelp"] = data["kelp"].rio.clip(land.to_crs(data["kelp"].rio.crs).geometry.values, invert=True)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cloud high probability"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["thin cirrus"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["defective"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cast shadow"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cloud shadow"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cloud medium probability"], numpy.nan)
                utils.update_raster_defaults(data)

                # Save each separately
                for index in range(len(data["kelp"].time)):
                    kelp = data["kelp"].isel(time=index).load()
                    kelp = kelp.rio.clip(roi.geometry.values)
                    filename = raster_path / f'kelp_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif'

                    kelp_info["area"].append(abs(int(kelp.notnull().sum() * kelp.x.resolution * kelp.y.resolution)))
                    kelp_info["file"].append(filename)
                    kelp_info["date"].append(pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format))
                    kelp_info["ocean cloud percentage"].append(ocean_cloud_percentage[index])

                    kelp.rio.to_raster(filename, compress="deflate", driver="COG") # missing min and max values when viewed in QGIS
                    data["SCL"].isel(time=index).rio.clip(roi.geometry.values).rio.to_raster(raster_path / f'scl_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif', compress="deflate", driver="COG")
                pandas.DataFrame.from_dict(kelp_info, orient='columns').to_csv(raster_path / "info.csv")

        # Save results
        kelp_info = pandas.DataFrame.from_dict(kelp_info, orient='columns')
        kelp_info.to_csv(raster_path / "info.csv")

def main():
    """ Create site datasets.
    """
    
    test_sites = utils.create_test_sites(distance_offshore = 3_000)
    test_sites_wsg = test_sites.to_crs(utils.CRS_WSG)
    land = geopandas.read_file(utils.DATA_PATH / "vectors" / "main_islands.gpkg")

    catalogue = {"url": "https://planetarycomputer.microsoft.com/api/stac/v1",
                 "collections": ["sentinel-2-l2a"]}

    date_format = "%Y-%m-%d"

    bands = ["red", "green", "blue", "nir", "SCL", "swir16", "B05", "B8A"]
    raster_defaults = {"resolution": 10, "nodata": 0, "dtype": "uint16"}
    scl_dict = {"no data": 0, "defective": 1, "cast shadow": 2, "cloud shadow": 3,
                "vegetation": 4, "not vegetated": 5, "water": 6, "unclassified": 7,
                "cloud medium probability": 8, "cloud high probability": 9,
                "thin cirrus": 10, "snow": 11}
    thresholds = {"min_ndvi": 0.03, "max_ndvi": 0.7, "max_ndwi": 0.1, "min_ndvri": 0.03, "max_ndwi2": -0.2,}
    
    filter_cloud_percentage = 30
    max_ocean_cloud_percentage = 10

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 
    
    for site_index, row in test_sites_wsg.iterrows():

        site_name = row['name']
        
        print(f"Test site: {site_name}") 
        raster_path = utils.DATA_PATH / "rasters" / "test_sites" / f"{site_name}"
        raster_path.mkdir(parents=True, exist_ok=True)
    
        # Geometry of AOI - convex hull to allow search
        geometry_query = shapely.box(*row.geometry.bounds) #row.geometry
        geometry = test_sites[test_sites["name"] == site_name].to_crs(utils.CRS_WSG).iloc[0].geometry
        
        #geometry_query = tiles[tiles["name"]==name].to_crs(crs_wsg).iloc[0].geometry

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

        years = list(range(2016, 2025))
        for year in years:
            months = [f"{year}-{str(month).zfill(2)}" for month in list(range(1, 13))]

            for month_YYMM in months:

                if max_date > datetime.datetime.strptime(month_YYMM, '%Y-%m'):
                    print(f"\tSkipping {month_YYMM} as run previously. Delete if you want a rerun")
                    continue

                print(f"\tCheck month: {month_YYMM}")
                # run pystac client search to see available dataset
                search = client.search(
                    collections=catalogue["collections"], intersects=geometry_query, datetime=month_YYMM, query=filters
                ) # intersects=geometry_query, 

                if len(search.item_collection()) == 0:
                    continue # Nothing meeting criteria for this month

                data = odc.stac.load(search.items(), geopolygon=geometry_query, bands=bands,  chunks={}, groupby="solar_day", # geopolygon=geometry_query, 
                                    resolution = raster_defaults["resolution"], dtype=raster_defaults["dtype"], nodata=raster_defaults["nodata"])
                roi = test_sites.to_crs(data["SCL"].rio.crs).loc[site_index]

                # remove if no data
                data = utils.screen_by_SCL_in_ROI(data, roi)
                '''data["SCL"].load()
                if not roi.geometry.intersects(shapely.box(*data.rio.bounds())):
                    print(f"\t\tWanring no intersection with ROI. Skipping.")
                    continue
                data["SCL"] = data["SCL"].rio.clip([roi.geometry])
                data["SCL"].rio.write_crs(data["SCL"].rio.crs, inplace=True);
                data = data.isel(time=(data["SCL"] != scl_dict["no data"]).any(dim=["x", "y"]))
                
                if len(data.time) == 0:
                    continue # Nothing meeting criteria for this month

                # remove if much cloud over ocean
                # Mask with 1s on ocean, 0 elsewhere
                ocean_mask = data["SCL"].isel(time=0).copy(deep=True)
                ocean_mask.data[:] = 1
                #ocean_mask = ocean_mask.rio.clip(land.to_crs(ocean_mask.rio.crs).geometry.values, invert=True)
                ocean_mask = ocean_mask.rio.clip([roi.geometry])
                # Mask by time - initially sums of cloud values then true / false by time if less than cloud threshold
                ocean_cloud_sum = (data["SCL"] == scl_dict["cloud high probability"]).sum(dim=["x", "y"]) 
                ocean_cloud_sum += (data["SCL"] == scl_dict["cloud medium probability"]).sum(dim=["x", "y"]) 
                ocean_cloud_sum += (data["SCL"] == scl_dict["cloud shadow"]).sum(dim=["x", "y"]) 
                ocean_cloud_sum += (data["SCL"] == scl_dict["cast shadow"]).sum(dim=["x", "y"]) 
                ocean_cloud_sum += (data["SCL"] == scl_dict["thin cirrus"]).sum(dim=["x", "y"])
                ocean_cloud_sum += (data["SCL"] == scl_dict["defective"]).sum(dim=["x", "y"])
                ocean_cloud_sum += (data["SCL"] == scl_dict["no data"]).sum(dim=["x", "y"]) - (ocean_mask == scl_dict["no data"]).sum(dim=["x", "y"])
                ocean_cloud_percentage = (ocean_cloud_sum / int(ocean_mask.sum())).data*100
                print(f"\t\tOcean cloud percentage {list(map('{:.2f}%'.format, ocean_cloud_percentage))}")
                cloud_mask_time = ocean_cloud_percentage < max_ocean_cloud_percentage
                data = data.isel(time=(cloud_mask_time))
                ocean_cloud_percentage = ocean_cloud_percentage[cloud_mask_time]

                if len(data.time) == 0:
                    print("\t\tNone with suitable cloud percentage.")
                    continue'''

                # Save out RGB
                '''rgb = data[["red", "green","blue"]].to_array("rgb", name="all images")
                utils.update_raster_defaults(rgb)
                rgb.rio.clip.to_netcdf(raster_path / f'rgb_{month_YYMM}.nc', format="NETCDF4", engine="netcdf4", encoding={"all images": {"zlib": True, "complevel": 2, "grid_mapping": rgb.encoding["grid_mapping"]}})'''

                # Convert to floats before calcualtions
                for key in data.data_vars:
                    if key == "SCL": 
                        continue
                    data[key] = data[key].astype("float32").where(data[key] != 0, numpy.nan)
                utils.update_raster_defaults(data)


                # Calculate NVDI and NVWI
                data["ndvi"] = (data.nir - data.red) / (data.nir + data.red)
                data["ndwi"] = (data.green - data.nir) / (data.green + data.nir)
                data["ndvri"] = (data.B05 - data.red) / (data.B05 + data.red);
                data["ndwi2"] = (data.swir16 + data.B05) / (data.swir16 - data.B05)
                utils.update_raster_defaults(data)

                # Calculate Kelp
                data["kelp"] = (data.nir - data.red) / (data.nir + data.red)
                data["kelp"] = data["kelp"].where(data["ndvi"].data > thresholds["min_ndvi"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["ndwi"].data < thresholds["max_ndwi"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["ndwi2"].data < thresholds["max_ndwi2"], numpy.nan)
                #data["kelp"] = data["kelp"].where(data["ndvi"].data < thresholds["max_ndvi"], numpy.nan)
                #data["kelp"] = data["kelp"].where(data["ndvri"].data > thresholds["min_ndvri"], numpy.nan)
                data["kelp"] = data["kelp"].rio.clip([roi.geometry]) #land.to_crs(data["kelp"].rio.crs).geometry.values, invert=True)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cloud high probability"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["thin cirrus"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["defective"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cast shadow"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cloud shadow"], numpy.nan)
                data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cloud medium probability"], numpy.nan)
                utils.update_raster_defaults(data)

                # Save each separately
                for index in range(len(data["kelp"].time)):
                    kelp = data["kelp"].isel(time=index).load()
                    #kelp = kelp.rio.clip(roi.geometry.values)
                    filename = raster_path / f'kelp_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif'

                    kelp_info["area"].append(abs(int(kelp.notnull().sum() * kelp.x.resolution * kelp.y.resolution)))
                    kelp_info["file"].append(filename)
                    kelp_info["date"].append(pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format))
                    kelp_info["ocean cloud percentage"].append(ocean_cloud_percentage[index])
                    kelp.rio.to_raster(filename, compress="deflate", driver="COG")  # missing min and max values when viewed in QGIS
                    
                    data["SCL"].isel(time=index).rio.to_raster(raster_path / f'scl_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif', compress="deflate", driver="COG")
                pandas.DataFrame.from_dict(kelp_info, orient='columns').to_csv(raster_path / "info.csv")

        # Save results
        kelp_info = pandas.DataFrame.from_dict(kelp_info, orient='columns')
        kelp_info.to_csv(raster_path / "info.csv")

if __name__ == '__main__':
    main()
