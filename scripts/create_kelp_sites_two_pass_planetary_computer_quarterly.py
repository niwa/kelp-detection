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

    bands = list(utils.SENTINEL_2B_BAND_INFO.keys()); bands.append("SCL")
    raster_defaults = {"resolution": 10, "nodata": 0, "dtype": "uint16"}
    
    # Define thresholds by year
    thresholds = {"min_ndvi": 0.03, "max_ndwi": 0.1, "max_ndwi2": -0.2,} # "max_ndwi2": -0.23 # "max_ndvi": 0.7, "min_ndvri": 0.03, 
    anomaly_thresholds = {"min_ndvi": 0.213, "max_ndwi": 0.1, "max_ndwi2": -0.2} 
    
    print(f"Thresholds by year: {thresholds}, and Anomaly thresholds by year: {anomaly_thresholds}")
    
    # Read in date to ignore file
    if (pathlib.Path.cwd() / "sites_dates_to_ignore.json").exists:
        with open(pathlib.Path.cwd() / "sites_dates_to_ignore.json", "r") as file_pointer:
            site_dates_to_ignore = json.load(file_pointer)
    else:
        site_dates_to_ignore = {}
    print(site_dates_to_ignore)
    
    filter_cloud_percentage = 30
    max_ocean_cloud_percentage = 5
    anomaly_detection_factor = 20
    debug = False
    
    # Second pass - remove beds with less than the min_pixls, then buffer outward by the specified number of pixels
    buffer = 10; min_pixels = 10 # 5

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 

    #test_sites_wsg = test_sites_wsg.iloc[20:]
    for site_index, row in test_sites_wsg.iterrows(): 
        site_name = row['name']
        
        dates_to_ignore = site_dates_to_ignore[site_name] if site_name in site_dates_to_ignore.keys() else {}
        
        print(f"Test site: {site_name}") 
        raster_path = utils.DATA_PATH / "rasters" / "test_sites_quarterly" / f"{site_name}"
        remote_raster_path = pathlib.Path("/nesi/nobackup/niwa03660/ZBD2023_outputs") / "test_sites_quarterly" / f"{site_name}_quarterly"
        raster_path.mkdir(parents=True, exist_ok=True)
        remote_raster_path.mkdir(parents=True, exist_ok=True)
    
        # Geometry of AOI - convex hull to allow search
        site_bbox = row.geometry.bounds

        filters = {"eo:cloud_cover":{"lt":filter_cloud_percentage}} 

        # Start from year of failure if already partial results
        if (raster_path / "info_quarterly.csv").exists():
            kelp_info = pandas.read_csv(raster_path / "info_quarterly.csv")
            kelp_info = kelp_info[['date', 'file', 'area', 'dates considered', "max coverage date"]]
            max_date = datetime.datetime.strptime(kelp_info["date"].max(), '%Y-%m-%d')
            kelp_info = kelp_info.to_dict(orient='list')
        else:
            kelp_info = {"date": [], "file": [], "area": [], "dates considered": [], "max coverage date": []}
            max_date = datetime.datetime.strptime("2015-01-31", '%Y-%m-%d')
        
        if (raster_path / "info_individual.csv").exists():
            kelp_info_individual = pandas.read_csv(raster_path / "info_individual.csv").to_dict()
        else:
            kelp_info_individual = {"date": [], "file": [], "area": [], "ocean cloud percentage": []}

        years = list(range(2016, 2025)) # 2016, 2025
        for year in years:
            quarters = [f"{year-1}-12/{year}-02", f"{year}-03/{year}-05", f"{year}-06/{year}-08", f"{year}-09/{year}-11"]

            for date_range in quarters:
                if max_date >= datetime.datetime.strptime(date_range.split("/")[0], '%Y-%m'):
                    print(f"\tSkipping {date_range} as run previously. Delete the info.csv if you want a rerun")
                    continue

                print(f"\tCheck month: {date_range}")
                # run pystac client search to see available dataset
                search = client.search(
                    collections=catalogue["collections"], bbox=site_bbox, datetime=date_range, query=filters
                )

                if len(search.item_collection()) == 0:
                    continue # Nothing meeting criteria for this month

                data = odc.stac.load(search.items(), bbox=site_bbox, bands=bands,  chunks={}, groupby="solar_day", 
                                    resolution = raster_defaults["resolution"], dtype=raster_defaults["dtype"], nodata=raster_defaults["nodata"])
                data = utils.harmonize_post_2022(data)
                roi = test_sites.to_crs(data["SCL"].rio.crs).loc[[site_index]]
                
                # remove any date to ignore
                mask_dates_to_ignore = data.time.isnull().data
                for date_to_ignore in dates_to_ignore:
                    mask_dates_to_ignore |= pandas.to_datetime(data.time).strftime(date_format) == date_to_ignore
                if mask_dates_to_ignore.any():
                    print(f"\t\tRemoving dates {data.time[mask_dates_to_ignore].data} of {data.time.data} as did not pass QA/QC.")
                    data = data.isel(time=(mask_dates_to_ignore==False))

                # remove if no data
                (data, ocean_cloud_percentage) = utils.screen_by_SCL_in_ROI(data, roi, max_ocean_cloud_percentage)

                if len(data.time) == 0:
                    print("\t\tNone with suitable cloud percentage.")
                    continue

                # Convert to floats before calcualtions
                for key in data.data_vars:
                    if key == "SCL": 
                        continue
                    data[key] = data[key].astype("float32").where(data[key] != 0, numpy.nan)
                utils.update_raster_defaults(data)

                # Calculate Kelp from thresholds
                data = utils.threshold_kelp(data, thresholds, roi)
                
                # Check for any big differences in kelp area using the anomaly thresholds
                data["kelp_original"] = data["kelp"].copy(deep=True)
                data = utils.threshold_kelp(data, anomaly_thresholds, roi)
                anomalous = data["kelp_original"].notnull().sum(dim=["x", "y"]) > anomaly_detection_factor * data["kelp"].notnull().sum(dim=["x", "y"])
                if anomalous.any():
                    anomalous.load()
                    print(f"\t\tKeeping only dates {data.time[anomalous==False].data} of {data.time.data} as anomalies detected in {data.time[anomalous].data}."
                          f"Area changes from {data['kelp_original'].notnull().sum(dim=['x', 'y']).load().data} to {data['kelp'].notnull().sum(dim=['x', 'y']).load().data}")
                    data = data.isel(time=(anomalous==False))  # Keep only those date that are not anomalies
                data["kelp"] = data["kelp_original"]
                data = data.drop_vars("kelp_original")
                
                if len(data.time.data) == 0:
                    print(f"\t\tNo dates remain after anomaly detection. move to next quarter.")
                    continue
                
                # Create buffered area around kelp beds
                kelp_polygons_buffered = []
                for index in range(len(data["kelp"].time)):
                    kelp_polygons_i = utils.polygon_from_raster(data["kelp"].isel(time=index)).explode(ignore_index=True, index_parts=False)
                    kelp_polygons_i = kelp_polygons_i[kelp_polygons_i.area > min_pixels*(utils.RASTER_DEFAULTS["resolution"]**2)]
                    kelp_polygons_i = kelp_polygons_i.buffer(utils.RASTER_DEFAULTS["resolution"] * buffer)
                    kelp_polygons_i = geopandas.GeoDataFrame(geometry=[kelp_polygons_i.union_all()], crs=kelp_polygons_i.crs)
                    kelp_polygons_i = kelp_polygons_i.overlay(roi, how="intersection", keep_geom_type=True)
                    kelp_polygons_buffered.append(kelp_polygons_i)
                
                # Caclulate Kelp from second thresholds - reset data first
                data = utils.threshold_kelp(data, thresholds, kelp_polygons_buffered)
                
                # Save combined polygon and area
                if debug:
                    encoding =  {"kelp": {"zlib": True, "complevel": 9, "grid_mapping": data["kelp"].encoding["grid_mapping"]}}
                    data["kelp"].to_netcdf(remote_raster_path / f'kelp_{date_range.replace("/","_")}.nc', format="NETCDF4", engine="netcdf4", encoding=encoding)
                kelp_polygons = []
                for index in range(len(data["kelp"].time)):
                    kelp_polygons.append(utils.polygon_from_raster(data["kelp"].isel(time=index)).dissolve())
                kelp_polygons = pandas.concat(kelp_polygons, ignore_index=True)
                if kelp_polygons.empty:
                    max_coverage_date = pandas.to_datetime(data.time.data[0]).strftime(date_format)
                else:
                    max_coverage_date = pandas.to_datetime(data.time.data[kelp_polygons.area.idxmax()]).strftime(date_format)
                kelp_polygons = kelp_polygons.dissolve()
                kelp_polygons.to_file(raster_path / f"{date_range.replace('/', '_')}_kelp.gpkg", index=False)

                kelp_info["area"].append(kelp_polygons.area.sum())
                kelp_info["file"].append(raster_path / f"{date_range.replace('/', '_')}_kelp.gpkg")
                kelp_info["date"].append(pandas.to_datetime(date_range.split("/")[0]).strftime(date_format))
                kelp_info["dates considered"].append(data.time.dt.strftime("%B %d, %Y").data)
                kelp_info["max coverage date"].append(max_coverage_date)
                pandas.DataFrame.from_dict(kelp_info, orient='columns').to_csv(raster_path / "info_quarterly.csv", index=False)
                if debug: 
                    pandas.DataFrame.from_dict(kelp_info, orient='columns').to_csv(remote_raster_path / "info_quarterly.csv", index=False)

                for index in range(len(data["kelp"].time)):

                    filename = remote_raster_path / f'data_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.nc'

                    data_i = data.isel(time=index)
                    kelp = data_i["kelp"].load()
                    kelp_info_individual["area"].append(abs(int(kelp.notnull().sum() * kelp.x.resolution * kelp.y.resolution)))
                    kelp_info_individual["file"].append(filename)
                    kelp_info_individual["date"].append(pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format))
                    kelp_info_individual["ocean cloud percentage"].append(ocean_cloud_percentage[index])
                    
                    if debug: # Save each date data separately
                        encoding = {}
                        for key in data_i.data_vars:
                            encoding[key] =  {"zlib": True, "complevel": 9, "grid_mapping": data[key].encoding["grid_mapping"]}
                            data_i.to_netcdf(filename, format="NETCDF4", engine="netcdf4", encoding=encoding)
                    pandas.DataFrame.from_dict(kelp_info_individual, orient='columns').to_csv(raster_path / "info.csv", index=False)
                    if debug: 
                        pandas.DataFrame.from_dict(kelp_info_individual, orient='columns').to_csv(remote_raster_path / "info.csv", index=False)

        # Save results
        kelp_info_individual = pandas.DataFrame.from_dict(kelp_info_individual, orient='columns')
        kelp_info_individual.to_csv(raster_path / "info_individual.csv", index=False)
        if debug: 
            kelp_info_individual.to_csv(remote_raster_path / "info_individual.csv", index=False)

if __name__ == '__main__':
    main()
