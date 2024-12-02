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


def update_raster_defaults(raster):
    # works on DataArrays and Datasets and for ints and floats
    if isinstance(raster, xarray.Dataset):
        for key in raster.data_vars:
            raster[key].rio.write_crs(raster[key].rio.crs, inplace=True)
            if raster[key].data.dtype == 'uint16':
                raster[key].rio.write_nodata(0, encoded=True, inplace=True)
            else: # assume float
                raster[key].rio.write_nodata(numpy.nan, encoded=True, inplace=True)
    raster.rio.write_crs(raster.rio.crs, inplace=True)
    if isinstance(raster, xarray.DataArray):
        if raster.data.dtype == 'uint16':
            raster.rio.write_nodata(0, encoded=True, inplace=True)
        else: # assume float
            raster.rio.write_nodata(numpy.nan, encoded=True, inplace=True)

def create_south_island_roi(data_path, crs):
    
    if not (data_path / "vectors" / "nz_islands.gpkg").exists():
        dotenv.load_dotenv()
        linz_key = os.environ.get("LINZ_API", None)
        fetcher = geoapis.vector.Linz(linz_key, verbose=False, crs=crs)
        islands = fetcher.run(51153)
        islands.to_file(data_path / "vectors" / "nz_islands.gpkg")
    
    import shapely
    islands = geopandas.read_file(data_path / "vectors" / "nz_islands.gpkg")
    south_islands = islands[islands["name"].isin(["South Island or Te Waipounamu", "Stewart Island/Rakiura"])]
    south_islands = islands[islands.intersects(shapely.geometry.box(*south_islands.total_bounds))]
    OFFSHORE_BUFFER = 4_000
    offshore_4km = geopandas.GeoDataFrame(geometry=south_islands.buffer(OFFSHORE_BUFFER), crs=crs).overlay(south_islands, how='difference')
    offshore_4km = offshore_4km.dissolve()
    offshore_4km.to_file(data_path / "vectors" / "offshore_4km.gpkg")

def create_tiles(data_path, crs):
    
    '''offshore_by_region = geopandas.read_file(r"https://services1.arcgis.com/3JjYDyG3oajxU6HO/arcgis/rest/services/MARINE_BioGeoRegions/FeatureServer/0/query?outFields=*&where=1%3D1&f=geojson")
    offshore_by_region = offshore_by_region.to_crs(2193)
    offshore_by_region.to_file(data_path / "vectors" / "offshore_by_region.gpkg")'''
    region_of_interest = geopandas.read_file(data_path / "vectors" / "offshore_4km.gpkg")
    
    TILE_LENGTH = 50_000
    region_of_interest = region_of_interest.to_crs(2193)
    bbox = region_of_interest.total_bounds
    n_tiles_x = int((bbox[2]-bbox[0]) / TILE_LENGTH + 1)
    n_tiles_y = int((bbox[3]-bbox[1]) / TILE_LENGTH + 1)
    origin_x = int((bbox[0] + bbox[2]) / 2) - (n_tiles_x / 2) * TILE_LENGTH
    origin_y = int((bbox[1] + bbox[3]) / 2) - (n_tiles_y / 2) * TILE_LENGTH

    tiles = {"name": [], "geometry": []}
    for i in range(n_tiles_y):
        for j in range(n_tiles_x):
            name = f"{i:02d}{j:02d}"
            tiles["name"].append(name)
            tiles["geometry"].append(
                shapely.geometry.Polygon(
                        [
                            (origin_x + j * TILE_LENGTH, origin_y + i * TILE_LENGTH),
                            (origin_x + (j + 1) * TILE_LENGTH, origin_y + i * TILE_LENGTH),
                            (origin_x + (j + 1) * TILE_LENGTH, origin_y + (i + 1) * TILE_LENGTH),
                            (origin_x + j * TILE_LENGTH, origin_y + (i + 1) * TILE_LENGTH),
                        ]
                    ))
    tiles = geopandas.GeoDataFrame(tiles, crs=crs)
    tiles = tiles[tiles.intersects(region_of_interest.dissolve().iloc[0].geometry)]
    tiles.to_file(data_path / "vectors" / f"tiles_full_size.gpkg")
    
    tiles = tiles.clip(region_of_interest)
    tiles["geometry"] = tiles["geometry"].apply(lambda geometry: shapely.geometry.box(*geometry.bounds))

    tiles.to_file(data_path / "vectors" / f"tiles.gpkg")

def main():
    """ Create Otago dataset.
    """
    
    name = "1111" # Akaroa - 0810, Waikouaiti - 0406, Kaikoura 0102, Oban 0102

    catalogue = {"url": "https://planetarycomputer.microsoft.com/api/stac/v1",
                 "collections": ["sentinel-2-l2a"]}

    crs_wsg = 4326
    crs = 2193
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

    data_path = pathlib.Path.cwd() / ".." / "data"
    raster_path = data_path / "rasters" / "tiles" / f"{name}"
    raster_path.mkdir(parents=True, exist_ok=True)
    (data_path / "vectors").mkdir(parents=True, exist_ok=True)

    # Setup vector inputs
    if not (data_path / "vectors" / "tiles.gpkg").exists() or not (data_path / "vectors" / "main_islands.gpkg").exists():
        dotenv.load_dotenv()
        linz_key = os.environ.get("LINZ_API", None)
        fetcher = geoapis.vector.Linz(linz_key, verbose=False, crs=crs)
        islands = fetcher.run(51153)
        
        main_islands = islands[islands.area > 9e8]
        main_islands.to_file(data_path / "vectors" / "main_islands.gpkg")
        
        create_tiles(data_path=data_path, crs=crs)

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 
    
    # Geometry of AOI - convex hull to allow search
    complete_roi = geopandas.read_file(data_path / "vectors" / "offshore_4km.gpkg")
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    tiles = geopandas.read_file(data_path / "vectors" / "tiles.gpkg")
    geometry_df = tiles[tiles["name"]==name]
    geometry_query = geometry_df.to_crs(crs_wsg).iloc[0].geometry

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
    
    years = list(range(2016, 2024))
    for year in years:
        months = [f"{year}-{str(month).zfill(2)}" for month in list(range(1, 13))]

        for month_YYMM in months:
            
            if max_date > datetime.datetime.strptime(month_YYMM, '%Y-%m'):
                print(f"Skipping {month_YYMM} as run previously. Delete if you want a rerun")
                continue

            print(f"Check month: {month_YYMM}")
            # run pystac client search to see available dataset
            search = client.search(
                collections=catalogue["collections"], intersects=geometry_query, datetime=month_YYMM, query=filters
            ) 

            if len(search.item_collection()) == 0:
                continue # Nothing meeting criteria for this month

            data = odc.stac.load(search.items(), geopolygon=geometry_query, bands=bands,  chunks={}, groupby="solar_day",
                                resolution = raster_defaults["resolution"], dtype=raster_defaults["dtype"], nodata=raster_defaults["nodata"])
            # remove if no data
            data["SCL"].load()
            data["SCL"] = data["SCL"].rio.clip(land.to_crs(data["SCL"].rio.crs).geometry.values, invert=True)
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
            print(f"Ocean cloud percentage {list(map('{:.2f}%'.format, ocean_cloud_percentage))}")
            cloud_mask_time = ocean_cloud_percentage < max_ocean_cloud_percentage
            data = data.isel(time=(cloud_mask_time))
            ocean_cloud_percentage = ocean_cloud_percentage[cloud_mask_time]
            
            if len(ocean_cloud_percentage) == 0:
                print("None with suitable cloud percentage.")
                continue

            # Save out RGB
            '''rgb = data[["red", "green","blue"]].to_array("rgb", name="all images")
            update_raster_defaults(rgb)
            rgb.to_netcdf(raster_path / f'rgb_{month_YYMM}.nc', format="NETCDF4", engine="netcdf4", encoding={"all images": {"zlib": True, "complevel": 2, "grid_mapping": rgb.encoding["grid_mapping"]}})'''

            # Convert to floats before calcualtions
            for key in data.data_vars:
                if key == "SCL": 
                    continue
                data[key] = data[key].astype("float32").where(data[key] != 0, numpy.nan)
            update_raster_defaults(data)


            # Calculate NVDI and NVWI
            data["ndvi"] = (data.nir - data.red) / (data.nir + data.red)
            data["ndwi"] = (data.green - data.nir) / (data.green + data.nir)
            data["ndvri"] = (data.B05 - data.red) / (data.B05 + data.red);
            data["ndwi2"] = (data.swir16 + data.B05) / (data.swir16 - data.B05)
            update_raster_defaults(data)

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
            update_raster_defaults(data)

            # Save each separately
            for index in range(len(data["kelp"].time)):
                kelp = data["kelp"].isel(time=index).load()
                filename = raster_path / f'kelp_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif'

                kelp_info["area"].append(abs(int(kelp.notnull().sum() * kelp.x.resolution * kelp.y.resolution)))
                kelp_info["file"].append(filename)
                kelp_info["date"].append(pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format))
                kelp_info["ocean cloud percentage"].append(ocean_cloud_percentage[index])

                kelp.rio.to_raster(filename, compress="deflate", driver="COG") # missing min and max values when viewed in QGIS
                data["SCL"].isel(time=index).rio.to_raster(raster_path / f'scl_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif', compress="deflate", driver="COG")
            pandas.DataFrame.from_dict(kelp_info, orient='columns').to_csv(raster_path / "info.csv")

    # Save results
    kelp_info = pandas.DataFrame.from_dict(kelp_info, orient='columns')
    kelp_info.to_csv(raster_path / "info.csv")

if __name__ == '__main__':
    main()
