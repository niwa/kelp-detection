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
BAND_OFFSET_POST_2022_01_25 = 1000

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
        
        # Motunau
        y0 = 5230000; x0 = 1604000
        y1 = 5234000; x1 = 1609000
        motunau = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        # Marlborough
        y0 = 5453055.733161794; x0 = 1694172.8611234897 
        y1 = 5429822.496895627; x1 = 1713207.4453088508
      #  y0 = 5415458; x0 = 1660774
        # y1 = 5502358; x1 = 1725491
        marlborough = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
    
        # Chatham
        y0 = 5076000; x0 = 2449000
        y1 = 5100000; x1 = 2470000
        chatham = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        # Pearl_Island - Polygon 1
        y0 = 4765392.547801174; x0 = 1186668.6951792038 
        y1 = 4745558.132935946; x1 = 1204618.3617616794
        pearl_island = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        # # Rakiura_Ulva - Polygon 2
        y0 = 4806923.258714756; x0 = 1223393.749804938
        y1 = 4781353.230771581; x1 = 1241770.613262257
        rakiura_ulva = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Bluff - Polygon 3
        y0 = 4830360.526196695; x0 = 1242241.2651234234 
        y1 = 4821717.456994672; x1 = 1251240.8490322009
        bluff = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Otekura - Polygon 4
        y0 = 4860347.641245389; x0 = 1346361.8815165777 
        y1 = 4845800.062247097; x1 = 1357540.5830005489
        otekura = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        #Brighton - Polygon 5
        y0 = 4908935.224387027; x0 = 1381945.325897603 
        y1 = 4882464.746047782; x1 = 1397364.623072171
        brighton = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Waikouaiti/Oamaru - Polygon 6
        y0 = 5004879.127438061; x0 = 1410502.6983560827 
        y1 = 4924655.966017282; x1 = 1444326.3632238633
        waikouaiti_oamaru = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

       # Timaru - Polygon 7
        y0 = 5087511.644697276; x0 = 1457497.571264771 
        y1 = 5060139.980299613; x1 = 1464489.8274933768
        timaru = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Akaroa/Banks Peninsula - Polygon 8
        y0 = 5152237.077727701; x0 = 1591688.587486241 
        y1 = 5137329.779619718; x1 = 1610582.2099330185
        akaroa_banks_peninsula = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Diamond Harbour - Polygon 9
        y0 = 5173971.31578453; x0 = 1577935.0334800433 
        y1 = 5163620.85565431; x1 = 1602105.9922983076
        diamond_harbour = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Black Sand Beach - Polygon 10
        y0 = 5290077.142456159; x0 = 1637051.3949440252 
        y1 = 5273017.7958546905; x1 = 1645760.0500301563
        black_sand_beach = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Wharanui_Kekerengu - Polygon 11
        y0 = 5365944.664544601; x0 = 1676008.6346885192 
        y1 = 5329490.577538376; x1 = 1696521.6617272045
        wharanui_kekerengu = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Wellington Petone - Polygon 12
        y0 = 5436416.665927516; x0 = 1747751.6254625104 
        y1 = 5408693.077236777; x1 = 1766158.9137329785
        wellington_petone = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        #Porirua - Polygon 13
        y0 = 5484483.840918069; x0 = 1746687.3415808624 
        y1 = 5444444.789643585; x1 = 1767572.332218603
        porirua = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        #Riversdale Beach - Polygon 14
        y0 = 5472042.805550492; x0 = 1849275.3141943226 
        y1 = 5427019.398448853; x1 = 1873552.6197312057
        riversdale_beach = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        #Resolution Island - Polygon 15
        y0 = 4928278.769345433; x0 = 1087321.7151045261 
        y1 = 4905362.73248211; x1 = 1106549.0823926302
        resolution_island = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        #Auckland Island - Polygon 16
        y0 = 4393640.2918102965; x0 = 1091931.9550202775 
        y1 = 4332345.107631876; x1 = 1134071.2601488307
        auckland_island = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
   
        #Campbell Island - Polygon 17
        y0 = 4183314.323317456; x0 = 1323474.7402687492
        y1 = 4158696.3337127017; x1 = 1357199.697531454
        campbell_island = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        #Long Sound - Polygon 18
        y0 = 4897091.419368935; x0 = 1089522.6079213074 
        y1 = 4865582.001213127; x1 = 1120288.4175324633
        long_sound = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        #Chatham Pitt - Polygon 19
        y0 = 5116898.692223684; x0 = 2406534.5690748165 
        y1 = 5013799.557472538; x1 = 2478598.318150406
        chatham_pitt = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])

        

        test_sites = geopandas.GeoDataFrame({"name": ["Motunau", "Marlborough", "Chatham", "Pearl Island", "Rakiura_Ulva", "Bluff", "Otekura", "Brighton", "Waikouaiti-Oamaru", "Timaru", "Akaroa-Banks Peninsula", "Diamond Harbour", "Black Sand Beach", " Wharanui Kekerengu", "Wellington-Petone", "Porirua", "Riversdale Beach", "Resolution Island", "Auckland Island", "Campbell Island", "Long Sound", "Chatham Pitt"], 
                                               "geometry": [motunau, marlborough, chatham, pearl_island, rakiura_ulva, bluff, otekura, brighton, waikouaiti_oamaru, timaru, akaroa_banks_peninsula, diamond_harbour, black_sand_beach, wharanui_kekerengu, wellington_petone, porirua, riversdale_beach, resolution_island, auckland_island, campbell_island, long_sound, chatham_pitt]}, crs=CRS)
        
        # test_sites = geopandas.GeoDataFrame({"name": ["Motunau", "Marlborough", "Chatham"], 
        #                                         "geometry": [motunau, marlborough, chatham]}, crs=CRS)

        # clip to max distance offshore
        buffer_path = DATA_PATH / "vectors" / f"offshore_buffer_{buffer_label}_main_islands.gpkg"
        if not buffer_path.exists():
            print(f"Creating buffer at: {buffer_path}")
            island_names = ["North Island or Te Ika-a-Māui", "South Island or Te Waipounamu", "Stewart Island/Rakiura", "Chatham Island", "Auckland Island", "Campbell Island/Motu Ihupuku"]
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
            print(f"\t\tHarmonizing date {raster['time'].isel(time=index).values}")
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