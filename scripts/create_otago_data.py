import pystac_client
import pystac
import odc.stac
import rioxarray
import pathlib
import pandas
import geopandas
import numpy
import dotenv
import os
import geoapis.vector

OFFSET = -0.1
SCALE = 0.0001

def check_offset_scale(search, bands):
    """ Check expected offset and scale """
    for item in search.item_collection():
        for band in bands:
            if band == "scl": 
                continue
            if item.assets[band].extra_fields['raster:bands'][0]['offset'] != OFFSET:
                raise Exception("offset doesn't match expected")
            if item.assets[band].extra_fields['raster:bands'][0]['scale'] != SCALE:
                raise Exception("scale doesn't match expected")

def main():
    """ Create Otago dataset.
    """
    
    name = "waikouaiti"

    collection = "sentinel-2-c1-l2a"

    crs_wsg = 4326
    crs = 2193
    date_format = "%Y-%m-%d"

    bands = ["red", "green", "blue", "nir", "scl"]
    SCL_CIRRUS = 10
    SCL_CLOUD = 9
    SCL_DEFECTIVE = 1
    SCL_NO_DATA = 0
    thresholds = {"min_ndvi": 0.01, "max_ndvi": 0.7, "max_ndwi": 0.2}

    data_path = pathlib.Path.cwd() / ".." / "data"
    (data_path / "rasters" / name).mkdir(parents=True, exist_ok=True)
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
    client = pystac_client.Client.open("https://earth-search.aws.element84.com/v1") 
    
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
                collections=[collection], intersects=geometry_query, datetime=month_YYMM, query=filters
            ) 

            if len(search.item_collection()) == 0:
                continue # Nothing meeting criteria for this month

            check_offset_scale(search=search, bands=bands)

            data = odc.stac.load(search.items(), geopolygon=geometry_query, bands=bands,  chunks={}, groupby="solar_day")

            # remove if no data
            data["scl"].load()
            data["scl"] = data["scl"].rio.clip(land.to_crs(data["scl"].rio.crs).geometry.values, invert=True, drop=False)
            data["scl"].rio.write_crs(data["scl"].rio.crs, inplace=True);
            data = data.isel(time=(data["scl"] != SCL_NO_DATA).any(dim=["x", "y"]))

            # convert bands to actual values
            for band in bands: 
                if band == "scl": 
                    continue
                data[band].data = data[band].data * SCALE + OFFSET
                data[band].data[data[band].data <= 0] = numpy.nan
                data[band].rio.write_nodata(numpy.nan, encoded=True, inplace=True);

            # Calculate NVDI and NVWI
            data["ndvi"] = (data.nir - data.red) / (data.nir + data.red)
            data["ndwi"] = (data.green-data.nir)/(data.green+data.nir)
            data["ndvi"].rio.write_crs(data["ndvi"].rio.crs, inplace=True); data["ndvi"].rio.write_nodata(numpy.nan, encoded=True, inplace=True);
            data["ndwi"].rio.write_crs(data["ndvi"].rio.crs, inplace=True); data["ndwi"].rio.write_nodata(numpy.nan, encoded=True, inplace=True); 

            # Calculate Kelp
            data["kelp"] = (data.nir - data.red) / (data.nir + data.red)
            data["kelp"].data[(data["ndvi"].data < thresholds["min_ndvi"]) | (data["ndvi"].data > thresholds["max_ndvi"]) | (data["ndwi"].data > thresholds["max_ndwi"])] = numpy.nan
            data["kelp"].data[(data["scl"].data == SCL_CIRRUS) | (data["ndvi"].data == SCL_CLOUD) | (data["ndwi"].data == SCL_DEFECTIVE)] = numpy.nan
            data["kelp"] = data["kelp"].rio.clip(land.to_crs(data["kelp"].rio.crs).geometry.values, invert=True)
            data["kelp"].rio.write_crs(data["kelp"].rio.crs, inplace=True); data["kelp"].rio.write_nodata(numpy.nan, encoded=True, inplace=True);

            # Save each separately
            for index in range(len(data["kelp"].time)):
                kelp = data["kelp"].isel(time=index).load()
                filename = data_path / "rasters" / name / f'kelp_{pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format)}.tif'

                kelp_info["area"].append(abs(int(kelp.notnull().sum() * kelp.x.resolution * kelp.y.resolution)))
                kelp_info["file"].append(filename)
                kelp_info["date"].append(pandas.to_datetime(data["kelp"].time.data[index]).strftime(date_format))

                kelp.rio.to_raster(filename, compress="deflate")

    # Save results
    kelp_info = pandas.DataFrame.from_dict(kelp_info, orient='columns')
    kelp_info.to_csv(data_path / "rasters" / name / "info.csv")

if __name__ == '__main__':
    main()
