import pystac_client
import pystac
import odc.stac
import rioxarray
import pathlib
import pandas
import geopandas
import numpy
import dotenv
import planetary_computer
import os
import geoapis.vector

OFFSET = -0.1
SCALE = 0.0001

def update_raster_defaults(raster):
    for key in raster.data_vars:
        if key == "SCL":
            continue
        raster[key].rio.write_crs(raster[key].rio.crs, inplace=True)
        raster[key].rio.write_nodata(numpy.nan, encoded=True, inplace=True);

def main():
    """ Create Otago dataset.
    """
    
    name = "waikouaiti"

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
    thresholds = {"min_ndvi": 0.01, "max_ndvi": 0.7, "max_ndwi": 0.2, "min_ndvri": 0.05}

    data_path = pathlib.Path.cwd() / ".." / "data"
    raster_path = data_path / "rasters" / f"{name}_planetarycomputer"
    raster_path.mkdir(parents=True, exist_ok=True)
    (data_path / "vectors").mkdir(parents=True, exist_ok=True)

    # Setup vector inputs
    if not (data_path / "vectors" / "regions.gpkg").exists() or not (data_path / "vectors" / "main_islands.gpkg").exists():
        dotenv.load_dotenv()
        linz_key = os.environ.get("LINZ_API", None)
        fetcher = geoapis.vector.Linz(linz_key, verbose=False, crs=crs)
        regions = fetcher.run(50785)
        islands = fetcher.run(51153)
        
        main_islands = islands[islands.area > 9e8]
        
        regions.to_file(data_path / "vectors" / "regions.gpkg")
        main_islands.to_file(data_path / "vectors" / "main_islands.gpkg")

    # use publically available stac link such as
    odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})
    client = pystac_client.Client.open(catalogue["url"], modifier=planetary_computer.sign_inplace) 
    
    # Geometry of AOI - convex hull to allow search
    land = geopandas.read_file(data_path / "vectors" / "main_islands.gpkg")
    geometry_df = geopandas.read_file(data_path / "vectors" / f"{name.lower()}.gpkg")
    geometry_query = geometry_df .to_crs(crs_wsg).iloc[0].geometry

    kelp_info = {"date": [], "file": [], "area": []}
    cloud_by_year = {"2016": 30, "2017": 30, "2018": 20, "2019": 20, "2020": 5, "2021": 10, "2022": 5, "2023": 5}

    years = list(range(2016, 2024))
    for year in years:
        months = [f"{year}-{str(month).zfill(2)}" for month in list(range(1, 13))]
        filters = {"eo:cloud_cover":{"lt":cloud_by_year[str(year)]}} 

        for month_YYMM in months:
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
            data["SCL"] = data["SCL"].rio.clip(land.to_crs(data["SCL"].rio.crs).geometry.values, invert=True, drop=False)
            data["SCL"].rio.write_crs(data["SCL"].rio.crs, inplace=True);
            data = data.isel(time=(data["SCL"] != scl_dict["no data"]).any(dim=["x", "y"]))

            # convert bands to actual values
            for band in bands: 
                if band == "SCL": 
                    continue
                data[band] = data[band].astype("float32").where(data[band] != 0, numpy.nan)
            update_raster_defaults(data)

            # Calculate NVDI and NVWI
            data["ndvi"] = (data.nir - data.red) / (data.nir + data.red)
            data["ndwi"] = (data.green-data.nir)/(data.green+data.nir)
            data["ndvri"] = (data.B05-data.red)/(data.B05+data.red);
            data["ndwi2"] = (data.swir16-data.B05)/(data.swir16+data.B05);
            update_raster_defaults(data)

            # Calculate Kelp
            data["kelp"] = (data.nir - data.red) / (data.nir + data.red)
            data["kelp"] = data["kelp"].where(data["ndvi"].data > thresholds["min_ndvi"], numpy.nan)
            data["kelp"] = data["kelp"].where(data["ndvi"].data < thresholds["max_ndvi"], numpy.nan)
            data["kelp"] = data["kelp"].where(data["ndwi"].data < thresholds["max_ndwi"], numpy.nan)
            data["kelp"] = data["kelp"].where(data["ndvri"].data > thresholds["min_ndvri"], numpy.nan)
            data["kelp"] = data["kelp"].rio.clip(land.to_crs(data["kelp"].rio.crs).geometry.values, invert=True)
            data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["cloud high probability"], numpy.nan)
            data["kelp"] = data["kelp"].where(data["SCL"] != scl_dict["defective"], numpy.nan)

            # Save each separately
            for index in range(len(data["kelp"].time)):
                kelp = data["kelp"].isel(time=index).load()
                filename = raster_path / f'kelp_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif'

                kelp_info["area"].append(abs(int(kelp.notnull().sum() * kelp.x.resolution * kelp.y.resolution)))
                kelp_info["file"].append(filename)
                kelp_info["date"].append(pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format))

                kelp.rio.to_raster(filename, compress="deflate", driver="COG")

    # Save results
    kelp_info = pandas.DataFrame.from_dict(kelp_info, orient='columns')
    kelp_info.to_csv(raster_path / "info.csv")

if __name__ == '__main__':
    main()
