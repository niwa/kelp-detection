import dotenv
import os
import geoapis.vector
import shapely
import pathlib
import pandas
import geopandas
import xarray
import rioxarray
import numpy
import odc.stac
import seaborn
import matplotlib.pyplot
import planetary_computer
import rasterio.features

CRS = 2193
CRS_WSG = 4326
DATA_PATH = pathlib.Path.cwd() / ".." / "data"
DATE_FORMAT = "%Y-%m-%d"

FILTER_CLOUD_PERCENTAGE = 30
MAX_OCEAN_CLOUD_PERCENTAGE = 10

SCL_DICT = {"no data": 0, "defective": 1, "cast shadow": 2, "cloud shadow": 3,
                "vegetation": 4, "not vegetated": 5, "water": 6, "unclassified": 7,
                "cloud medium probability": 8, "cloud high probability": 9,
                "thin cirrus": 10, "snow": 11}

RASTER_DEFAULTS = {"resolution": 10, "nodata": 0, "dtype": "uint16"}




# https://sentiwiki.copernicus.eu/web/s2-processing
# https://sentinels.copernicus.eu/web/sentinel/-/copernicus-sentinel-2-major-products-upgrade-upcoming
# https://forum.step.esa.int/t/changes-in-band-data-after-25-jan-2022-baseline-04-00-harmonizevalues-sentinel-2-l2a-snappy/36270
BAND_OFFSETS_SCALING_FACTOR = 10000
BAND_OFFSET_POST_2022_01_25 = 1000
BAND_OFFSETS_POST_2022_01_25 = {
    "B01": -0.0009,
    "B02": -0.0002,
    "B03": -0.0002,
    "B04": -0.0002,
    "B05": -0.0001,
    "B06": -0.0002,
    "B07": -0.0002,
    "B08": -0.0002,
    "B8A": -0.0002,
    "B09": -0.0001,
    "B10": -0.0001,
    "B11": -0.0001,
    "B12": -0.0001,
}

# source: https://sentiwiki.copernicus.eu/web/s2-mission#S2Mission-SpectralResolutionS2-Mission-Spectral-Resolution
SENTINEL_2B_BAND_INFO = {
    "B01": {"name": "coastal", "wavelength": 442.7, "bandwidth": 20},
    "B02": {"name": "blue", "wavelength": 492.7, "bandwidth": 65},
    "B03": {"name": "green", "wavelength": 559.8, "bandwidth": 35},
    "B04": {"name": "red", "wavelength": 664.6, "bandwidth": 30},
    "B05": {"name": "rededge - Band 5 - Vegetation red edge 1", "wavelength": 704.1, "bandwidth": 14},
    "B06": {"name": "rededge - Band 6 - Vegetation red edge 2", "wavelength": 740.5, "bandwidth": 14},
    "B07": {"name": "rededge - Band 7 - Vegetation red edge 3", "wavelength": 782.8, "bandwidth": 19},
    "B08": {"name": "nir", "wavelength": 832.8, "bandwidth": 105},
    "B09": {"name": "water vapor", "wavelength": 945.1, "bandwidth": 19},
    #"B10": {"name": "", "wavelength": 1373.5, "bandwidth": 29},
    "B11": {"name": "swir16", "wavelength": 1613.7, "bandwidth": 90},
    "B12": {"name": "swir22", "wavelength": 2202.4, "bandwidth": 174},
    "B8A": {"name": "rededge - Band 8A - Vegetation red edge 4", "wavelength": 864.7, "bandwidth": 21}
}


CATALOGUE = {"planetarycomputer": {"url": "https://planetarycomputer.microsoft.com/api/stac/v1",
                                   "collections": {"sentinel": "sentinel-2-l2a", "dem": "cop-dem-glo-30"}}}

def get_band_names_from_common(common_names): 
    band_names = []
    for common_name in common_names:
        found = False
        for key, value in SENTINEL_2B_BAND_INFO.items():
            #print(f"key {key}, value {value['name']}")
            if value["name"] == common_name:
                band_names.append(key)
                found=True
                break
        if not found:
            print(f"Warning - no band name found for common name '{common_name}'")
    return band_names

def common_name_to_band(common_name): 
    band_name = None
    for key, value in SENTINEL_2B_BAND_INFO.items():
        #print(f"key {key}, value {value['name']}")
        if value["name"] == common_name:
            band_name = key
            break
    if band_name is None:
        print(f"Warning - no band name found for common name '{common_name}'")
    return band_name


def download_nz_outline():
    dotenv.load_dotenv()
    (DATA_PATH / "vectors").mkdir(parents=True, exist_ok=True)
    if not (DATA_PATH / "vectors" / "nz_islands.gpkg").exists():
        linz_key = os.environ.get("LINZ_API", None)
        fetcher = geoapis.vector.Linz(linz_key, verbose=False, crs=CRS)
        islands = fetcher.run(51153)
        islands.to_file(DATA_PATH / "vectors" / "nz_islands.gpkg")
    if not (DATA_PATH / "vectors" / "main_islands.gpkg").exists():
        islands = geopandas.read_file(DATA_PATH / "vectors" / "nz_islands.gpkg")
        big_islands = islands[islands.area > 1e8]
        main_islands = islands[islands.area > 9e8]
        main_islands.to_file(DATA_PATH / "vectors" / "main_islands.gpkg")

def create_tiles_south_island(distance_offshore = 4_000, tile_length = 10_000):

    label = "SI"
    island_names = ["South Island or Te Waipounamu", "Stewart Island/Rakiura"]
    tiles = _create_tiles(distance_offshore, tile_length, label, island_names)
    return tiles

def create_tiles_north_island(distance_offshore = 4_000, tile_length = 10_000):

    label = "NI"
    island_names = ["North Island or Te Ika-a-Māui"]

    buffer_label = (str(int(distance_offshore/1000)) if distance_offshore%1000 == 0 else f'{distance_offshore/1000:2.1f}'.replace('.', '_')) + f"km_{label}"
    tile_label = (str(int(tile_length/1000)) if tile_length%1000 == 0 else f'{tile_length/1000:2.1f}'.replace('.', '_')) + f"km2"
    tile_path = DATA_PATH / "vectors" / f"tiles_{tile_label}_buffer_{buffer_label}.gpkg"

    if not tile_path.exists():
        y0 = 5380000; x0 = 1725000
        y1 = 5480000; x1 = 1890000
        lower_NI_bbox = geopandas.GeoDataFrame(geometry=[shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])], crs=CRS)
    
        tiles = _create_tiles(distance_offshore, tile_length, label, island_names)
        tiles = tiles[tiles.intersects(lower_NI_bbox.iloc[0].geometry)]
        tiles.to_file(tile_path)
    else:
        tiles = geopandas.read_file(tile_path)
    
    return tiles

def create_tiles_chatham_island(distance_offshore = 4_000, tile_length = 10_000):

    label = "CI"
    island_names = ["Chatham Island"]
    tiles = _create_tiles(distance_offshore, tile_length, label, island_names)
    return tiles

def _create_tiles(distance_offshore, tile_length, label, island_names):

    buffer_label = (str(int(distance_offshore/1000)) if distance_offshore%1000 == 0 else f'{distance_offshore/1000:2.1f}'.replace('.', '_')) + f"km_{label}"
    tile_label = (str(int(tile_length/1000)) if tile_length%1000 == 0 else f'{tile_length/1000:2.1f}'.replace('.', '_')) + f"km2"
    tile_path = DATA_PATH / "vectors" / f"tiles_{tile_label}_buffer_{buffer_label}.gpkg"
    buffer_path = DATA_PATH / "vectors" / f"offshore_buffer_{buffer_label}.gpkg"
    
    download_nz_outline()

    if not buffer_path.exists():
        islands = geopandas.read_file(DATA_PATH / "vectors" / "nz_islands.gpkg")
        select_islands = islands[islands["name"].isin(island_names)]
        select_islands = islands[islands.intersects(shapely.geometry.box(*select_islands.total_bounds))]
        offshore_buffer = geopandas.GeoDataFrame(geometry=select_islands.buffer(distance_offshore), crs=CRS).overlay(select_islands, how='difference')
        offshore_buffer = offshore_buffer.dissolve()
        offshore_buffer.to_file(buffer_path)

    if not tile_path.exists():
        offshore_buffer = geopandas.read_file(buffer_path).to_crs(CRS)
        bbox = offshore_buffer.total_bounds
        n_tiles_x = int((bbox[2]-bbox[0]) / tile_length + 1)
        n_tiles_y = int((bbox[3]-bbox[1]) / tile_length + 1)
        origin_x = int((bbox[0] + bbox[2]) / 2) - (n_tiles_x / 2) * tile_length
        origin_y = int((bbox[1] + bbox[3]) / 2) - (n_tiles_y / 2) * tile_length
        
        tiles = {"name": [], "geometry": []}
        for i in range(n_tiles_y):
            for j in range(n_tiles_x):
                name = f"{i:02d}{j:02d}"
                tiles["name"].append(name)
                tiles["geometry"].append(
                    shapely.geometry.Polygon(
                            [
                                (origin_x + j * tile_length, origin_y + i * tile_length),
                                (origin_x + (j + 1) * tile_length, origin_y + i * tile_length),
                                (origin_x + (j + 1) * tile_length, origin_y + (i + 1) * tile_length),
                                (origin_x + j * tile_length, origin_y + (i + 1) * tile_length),
                            ]
                        ))
        tiles = geopandas.GeoDataFrame(tiles, crs=CRS)
        tiles = tiles[tiles.intersects(offshore_buffer.dissolve().iloc[0].geometry)]
        tiles = tiles.clip(offshore_buffer)
        tiles["geometry"] = tiles["geometry"].apply(lambda geometry: shapely.geometry.box(*geometry.bounds))
        tiles.to_file(tile_path)
    tiles = geopandas.read_file(tile_path)

    return tiles

def create_test_sites(distance_offshore = 4_000):

    buffer_label = (str(int(distance_offshore/1000)) if distance_offshore%1000 == 0 else f'{distance_offshore/1000:2.1f}'.replace('.', '_')) + "km"

    test_sites_path = DATA_PATH / "vectors" / f"test_sites_offshore_{buffer_label}.gpkg"
    if not test_sites_path.exists():
        # Rakiora
        y0 = 4787000; x0 = 1231000
        y1 = 4793000; x1 = 1238000
        rakiora = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        # Waikouaiti
        y0 = 4954605.549945548; x0 = 1431878.0905242788
        y1 = 4944202.117311343; x1 = 1425794.9377906742
        y2 = 4945964.990144346; x2 = 1418619.3004845034
        y3 = 4957038.81103899;  x3 = 1425819.7669855051
        waikouaiti = shapely.geometry.Polygon([[x0,y0], [x1,y1], [x2,y2], [x3,y3]])
    
        # Akaroa
        y0 = 5140000; x0 = 1592000
        y1 = 5156000; x1 = 1598000
        akaroa = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        # Matau
        y0 = 5230000; x0 = 1604000
        y1 = 5234000; x1 = 1609000
        matau = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        # Marlborogh
        y0 = 5415458; x0 = 1660774
        y1 = 5502358; x1 = 1725491
        marlborough = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        # Wellington
        y0 = 5419000; x0 = 1753000
        y1 = 5426000; x1 = 1759000
        wellington = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
    
        # Chatham
        y0 = 5076000; x0 = 2449000
        y1 = 5100000; x1 = 2470000
        chatham = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        test_sites = geopandas.GeoDataFrame({"name": ["rakiora", "waikouaiti", "akaroa", "matau", "marlborough","wellington", "chatham"], "geometry": [rakiora, waikouaiti, akaroa, matau, marlborough, wellington, chatham]}, crs=CRS)

        # clip to max distance offshore
        buffer_path = DATA_PATH / "vectors" / f"offshore_buffer_{buffer_label}_main_islands.gpkg"
        if not buffer_path.exists():
            island_names = ["North Island or Te Ika-a-Māui", "South Island or Te Waipounamu", "Stewart Island/Rakiura", "Chatham Island"]
            download_nz_outline()
            islands = geopandas.read_file(DATA_PATH / "vectors" / "nz_islands.gpkg")
            select_islands = islands[islands["name"].isin(island_names)]
            select_islands = islands[islands.intersects(shapely.geometry.box(*select_islands.total_bounds))]
            offshore_buffer = geopandas.GeoDataFrame(geometry=select_islands.buffer(distance_offshore), crs=CRS).overlay(select_islands, how='difference')
            offshore_buffer.to_file(buffer_path)
        else:
            offshore_buffer = geopandas.read_file(buffer_path)

        test_sites = test_sites.clip(offshore_buffer.dissolve().geometry.loc[0])
        test_sites.to_file(test_sites_path)
    else:
        test_sites = geopandas.read_file(test_sites_path)

    return test_sites


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

def harmonize_post_2022(raster):
    harmonize_date = "2022-01-25"
    
    for index in range(len(raster["time"])):
        if numpy.datetime64(harmonize_date) < raster["time"].isel(time=index).values:
            print(f"\t\tHarmonizing date {raster["time"].isel(time=index).values}")
            raster_i = raster.isel(time=index)

            for band in SENTINEL_2B_BAND_INFO.keys():
                #breakpoint()
                if raster_i[band].data.dtype == 'uint16':
                    raster_i[band] = raster_i[band].clip(min=1000) - BAND_OFFSET_POST_2022_01_25
                else: # Assume normalized float
                    raster_i[band] = raster_i[band].clip(min=1000) - BAND_OFFSET_POST_2022_01_25
            raster.loc[{'time': raster.time[index].values}] = raster_i
    return raster


def screen_by_SCL_in_ROI(data, roi, max_ocean_cloud_percentage):
    data["SCL"].load()
    if not roi.geometry.intersects(shapely.box(*data.rio.bounds())).any():
        print(f"\t\tWarning no intersection with ROI. Skipping.")
        return data
    data["SCL"] = data["SCL"].rio.clip(roi.geometry)
    data["SCL"].rio.write_crs(data["SCL"].rio.crs, inplace=True);
    data = data.isel(time=(data["SCL"] != SCL_DICT["no data"]).any(dim=["x", "y"]))

    if len(data.time) == 0:
        return data # Nothing meeting criteria for this month

    # remove if much cloud over ocean
    # Mask with 1s on ocean, 0 elsewhere
    ocean_mask = data["SCL"].isel(time=0).copy(deep=True)
    ocean_mask.data[:] = 1
    #ocean_mask = ocean_mask.rio.clip(land.to_crs(ocean_mask.rio.crs).geometry.values, invert=True)
    ocean_mask = ocean_mask.rio.clip(roi.geometry)
    # Mask by time - initially sums of cloud values then true / false by time if less than cloud threshold
    ocean_cloud_sum = (data["SCL"] == SCL_DICT["cloud high probability"]).sum(dim=["x", "y"]) 
    ocean_cloud_sum += (data["SCL"] == SCL_DICT["cloud medium probability"]).sum(dim=["x", "y"]) 
    #ocean_cloud_sum += (data["SCL"] == SCL_DICT["cloud shadow"]).sum(dim=["x", "y"]) 
    #ocean_cloud_sum += (data["SCL"] == SCL_DICT["cast shadow"]).sum(dim=["x", "y"]) 
    ocean_cloud_sum += (data["SCL"] == SCL_DICT["thin cirrus"]).sum(dim=["x", "y"])
    ocean_cloud_sum += (data["SCL"] == SCL_DICT["defective"]).sum(dim=["x", "y"])
    ocean_cloud_sum += (data["SCL"] == SCL_DICT["no data"]).sum(dim=["x", "y"]) - (ocean_mask == SCL_DICT["no data"]).sum(dim=["x", "y"])
    ocean_cloud_percentage = (ocean_cloud_sum / int(ocean_mask.sum())).data*100
    print(f"\t\tOcean cloud percentage {list(map('{:.2f}%'.format, ocean_cloud_percentage))}. Keep those less than {max_ocean_cloud_percentage}.")
    cloud_mask_time = ocean_cloud_percentage < max_ocean_cloud_percentage
    data = data.isel(time=(cloud_mask_time))
    ocean_cloud_percentage = ocean_cloud_percentage[cloud_mask_time]
    
    return (data, ocean_cloud_percentage)

def polygon_from_raster(data: xarray.DataArray):
    data.load()
    polygons = numpy.array([[value, shapely.geometry.shape(polygon)] for polygon, value in rasterio.features.shapes(numpy.uint8(data.notnull().values)) if value != 0])
    if len(polygons) > 0:
        polygons = geopandas.GeoDataFrame({'value': polygons[:, 0], 'geometry': polygons[:, 1]}, crs=data.rio.crs)
        transform = data.rio.transform()
        polygons["geometry"] = polygons.affine_transform([transform.a, transform.b, transform.d, transform.e, transform.xoff, transform.yoff, ]).buffer(0)
    else:
        polygons = geopandas.GeoDataFrame(geometry=[], crs=data.rio.crs)

    return polygons
    

def threshold_kelp(data, thresholds, roi):
    data["ndvi"] = (data[common_name_to_band("nir")] - data[common_name_to_band("red")]) / (data[common_name_to_band("nir")] + data[common_name_to_band("red")])
    data["ndvri"] = (data["B05"] - data[common_name_to_band("red")]) / (data["B05"] + data[common_name_to_band("red")])
    data["ndwi"] = (data[common_name_to_band("green")] - data[common_name_to_band("nir")]) / (data[common_name_to_band("green")] + data[common_name_to_band("nir")])
    data["ndwi2"] = (data[common_name_to_band("swir16")] - data["B05"]) / (data[common_name_to_band("swir16")] + data["B05"])
    data["ndwi3"] = (data[common_name_to_band("blue")] - data[common_name_to_band("coastal")]) / (data[common_name_to_band("blue")] + data[common_name_to_band("coastal")])
    data["ndwi4"] = (data[common_name_to_band("water vapor")] - data[common_name_to_band("nir")]) / (data[common_name_to_band("water vapor")] + data[common_name_to_band("nir")])
    update_raster_defaults(data)

    # Calculate Kelp
    data["kelp"] = (data[common_name_to_band("nir")] - data[common_name_to_band("red")]) / (data[common_name_to_band("nir")] + data[common_name_to_band("red")])
    for name, value in thresholds.items():
        if "min_" in name:
            data["kelp"] = data["kelp"].where(data[name.replace("min_", "")].data > value, numpy.nan)
        elif "max_" in name:
            data["kelp"] = data["kelp"].where(data[name.replace("max_", "")].data < value, numpy.nan)
        else:
            print(f"Threshold does not specify min or max: {name}")

    # Either clip all dates the same, or clip date by date
    if isinstance(roi, list):
        for index in range(len(data["kelp"].time)):
            if roi[index].area.sum() > 0:
                data.loc[{'time': data.time[index].values}] = data.isel(time=index).rio.clip(roi[index].geometry, drop=False)
            else:
                data.loc[{'time': data.time[index].values}] = numpy.nan
    elif isinstance(roi, pandas.DataFrame):
        data["kelp"] = data["kelp"].rio.clip(roi.geometry, drop=False) #land.to_crs(data["kelp"].rio.crs).geometry.values, invert=True)
    else:
        raise ValueError(f"roi should be a GeoPandas GeoDataFrame or a list of these. It is: {roi}")
    data["kelp"] = data["kelp"].where(data["SCL"] != SCL_DICT["cloud high probability"], numpy.nan)
    data["kelp"] = data["kelp"].where(data["SCL"] != SCL_DICT["thin cirrus"], numpy.nan)
    data["kelp"] = data["kelp"].where(data["SCL"] != SCL_DICT["defective"], numpy.nan)
    #data["kelp"] = data["kelp"].where(data["SCL"] != SCL_DICT["cast shadow"], numpy.nan)
    #data["kelp"] = data["kelp"].where(data["SCL"] != SCL_DICT["cloud shadow"], numpy.nan)
    data["kelp"] = data["kelp"].where(data["SCL"] != SCL_DICT["cloud medium probability"], numpy.nan)
    update_raster_defaults(data)
    
    return data

def kelp_single_day(client, date_YYMM, roi, thresholds):

    bbox = roi.to_crs(CRS_WSG).iloc[0].geometry.bounds
    
    filters = {"eo:cloud_cover":{"lt": FILTER_CLOUD_PERCENTAGE}}
    bands = list(SENTINEL_2B_BAND_INFO.keys()); bands.append("SCL") # ["red", "green", "blue", "nir", "SCL", "swir16", "B05", "B8A"] # Band 05 - Vegetation red edge 1, Band 8A - Vegetation red edge 4

    # Get data
    search_sentinel = client.search(collections=[CATALOGUE["planetarycomputer"]["collections"]["sentinel"]], bbox=bbox, datetime=date_YYMM, query=filters)
    data = odc.stac.load(search_sentinel.items(), bbox=bbox, bands=bands, chunks={}, groupby="solar_day",
                         resolution=RASTER_DEFAULTS["resolution"], dtype=RASTER_DEFAULTS["dtype"], nodata=RASTER_DEFAULTS["nodata"],
                         patch_url=planetary_computer.sign)

    # Keep only low cloud events
    roi = roi.to_crs(data["SCL"].rio.crs).iloc[[0]]
    (data, ocean_cloud_percentage) = screen_by_SCL_in_ROI(data, roi, MAX_OCEAN_CLOUD_PERCENTAGE)

    for key in data.data_vars:
        if key == "SCL": 
            continue
        data[key] = data[key].astype("float32").where(data[key] != 0, numpy.nan)

    # Kelp from thresholds
    data = threshold_kelp(data, thresholds, roi)

    return data

def plot_hists_single_day(data, date_YYMM, ocean_buffered, output_path):
    mask = data["kelp"].notnull()
    kelp_df = {}
    for key, value in SENTINEL_2B_BAND_INFO.items():
        kelp_df[f"{key}: {value['name']}"] = data[key].data[mask]
    kelp_df = pandas.DataFrame(data=kelp_df)
    kelp_df["classification"] = "kelp"
    
    mask = xarray.ones_like(data["kelp"])
    mask = mask.rio.clip(ocean_buffered.geometry, drop=False, invert=False).notnull()
    ocean_df = {}
    for key, value in SENTINEL_2B_BAND_INFO.items():
        ocean_df[f"{key}: {value['name']}"] = data[key].data[mask]
    ocean_df = pandas.DataFrame(data=ocean_df)
    ocean_df["classification"] = "ocean"
    plotting_df = pandas.concat([kelp_df, ocean_df]).melt(var_name='column', value_name='data', id_vars='classification')
    
    fig = seaborn.displot(data=kelp_df.melt(var_name='column', value_name='data', id_vars='classification'), x="data", col="column", col_wrap=4, stat='percent', bins=60, binrange=(100, 2400))
    fig.savefig(output_path / f'{date_YYMM}_kelp_hists.png')
    matplotlib.pyplot.close()
    fig = seaborn.displot(data=ocean_df.melt(var_name='column', value_name='data', id_vars='classification'), x="data", col="column", col_wrap=4, stat='percent', bins=60, binrange=(100, 2400))
    fig.savefig(output_path / f'{date_YYMM}_ocean_hists.png')
    matplotlib.pyplot.close()

def plot_lines(xys, data):
    point_df = {"name": [], "wavelength": [], "bandwidth": [], "value": [], "xy": []}
    for index, xy in enumerate(xys):
        point = data.sel(x=xy[0],y=xy[1], method="nearest")
    
        for key, value in SENTINEL_2B_BAND_INFO.items():
            point_df["name"].append(key)
            point_df["wavelength"].append(value["wavelength"])
            point_df["bandwidth"].append(value["bandwidth"])
            point_df["value"].append(float(point[key]))
            point_df["xy"].append(f"Click {index}: {xy[0]:.2f}, {xy[1]:.2f}")
    point_df = pandas.DataFrame(data=point_df)
    plot = seaborn.lineplot(data=point_df.pivot(index="wavelength", columns="xy", values="value"), markers=True)
    plot.set(title="Spectral Exploration")
    plot.set(ylabel='Reflectance')
    return plot

def update_spectra(xy, data, spectra_dict):
    if spectra_dict is None:
        spectra_dict = {"band": ["description", "wavelength", "bandwidth"], **{key: [SENTINEL_2B_BAND_INFO[key]["name"], SENTINEL_2B_BAND_INFO[key]["wavelength"], SENTINEL_2B_BAND_INFO[key]["bandwidth"]] for key in SENTINEL_2B_BAND_INFO.keys()}}
    
    point = data.sel(x=xy[0],y=xy[1], method="nearest")
    
    spectra_dict["band"].append(f"Click {len(spectra_dict['band'])-3}: {xy[0]:.2f}, {xy[1]:.2f}")
    for key, value in SENTINEL_2B_BAND_INFO.items():
        spectra_dict[key].append(float(point[key]))
    return spectra_dict
            
def get_spectra_all_dates(xy, kelp_info):
    
    data = rioxarray.rioxarray.open_rasterio(kelp_info.iloc[0]["file"], chunks=True).squeeze( "band", drop=True)
    
    indices = [key for key in data.data_vars if key not in SENTINEL_2B_BAND_INFO.keys()]
    
    spectra_dict = {"band": ["description", "wavelength", "bandwidth"], **{key: [SENTINEL_2B_BAND_INFO[key]["name"], SENTINEL_2B_BAND_INFO[key]["wavelength"], SENTINEL_2B_BAND_INFO[key]["bandwidth"]] for key in SENTINEL_2B_BAND_INFO.keys()}, **{key: ["", "", ""] for key in indices}}
    
    for index, row in kelp_info.iterrows():
        data = rioxarray.rioxarray.open_rasterio(row["file"], chunks=True).squeeze( "band", drop=True)
        point = data.sel(x=xy[0],y=xy[1], method="nearest")
    
        spectra_dict["band"].append(row["date"])
        for key in data.data_vars:
            spectra_dict[key].append(float(point[key]))
    return spectra_dict   

def normalise_rgb(data, bands):
    rgb = data[bands].copy(deep=True)
    max_val = min([rgb[band].max() for band in bands])
    for key in rgb.data_vars:
        rgb[key] = rgb[key] / max_val
    return rgb