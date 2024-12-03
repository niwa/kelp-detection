import dotenv
import os
import geoapis.vector
import shapely
import pathlib
import geopandas
import xarray
import numpy

CRS = 2193
CRS_WSG = 4326
DATA_PATH = pathlib.Path.cwd() / ".." / "data"

def download_nz_outline():
    dotenv.load_dotenv()
    if not (DATA_PATH / "vectors" / "nz_islands.gpkg").exists():
        linz_key = os.environ.get("LINZ_API", None)
        fetcher = geoapis.vector.Linz(linz_key, verbose=False, crs=CRS)
        islands = fetcher.run(51153)
        islands.to_file(DATA_PATH / "vectors" / "nz_islands.gpkg")
    if not (DATA_PATH / "vectors" / "main_islands.gpkg"):
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
        matau = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        # Matau
        y0 = 5230000; x0 = 1604000
        y1 = 5234000; x1 = 1609000
        matau = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        # Wellington
        y0 = 5419000; x0 = 1753000
        y1 = 5426000; x1 = 1759000
        wellington = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
    
        # Chatham
        y0 = 5076000; x0 = 2449000
        y1 = 5100000; x1 = 2470000
        chatham = shapely.geometry.Polygon([[x0,y0], [x1,y0], [x1,y1], [x0,y1]])
        
        test_sites = geopandas.GeoDataFrame({"name": ["rakiora", "waikouaiti", "matau", "wellington", "chatham"], "geometry": [rakiora, waikouaiti, matau, wellington, chatham]}, crs=CRS)

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


def screen_by_SCL_in_ROI(data, roi):
    data["SCL"].load()
    if not roi.geometry.intersects(shapely.box(*data.rio.bounds())):
        print(f"\t\tWanring no intersection with ROI. Skipping.")
        return data
    data["SCL"] = data["SCL"].rio.clip([roi.geometry])
    data["SCL"].rio.write_crs(data["SCL"].rio.crs, inplace=True);
    data = data.isel(time=(data["SCL"] != scl_dict["no data"]).any(dim=["x", "y"]))

    if len(data.time) == 0:
        return data # Nothing meeting criteria for this month

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
    
    return data

