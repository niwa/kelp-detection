import pathlib
import pandas
import geopandas

import utils

def main():
    """ Combine all geometires for a site
    """
    
    site_name = "marlborough"
    quarterly_path = utils.DATA_PATH / "rasters" / "test_sites_quarterly" / site_name
    
    extent_files = list(quarterly_path.glob("*_kelp.gpkg"))
    extents = []
    for extent_file in extent_files:
        extents.append(geopandas.read_file(extent_file).to_crs(utils.CRS))
    extents = pandas.concat(extents)
    extents.dissolve().to_file(quarterly_path.parent / f"{site_name}_any_date.gpkg")

if __name__ == '__main__':
    main()
