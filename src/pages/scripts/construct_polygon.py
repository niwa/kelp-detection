# -*- coding: utf-8 -*-
"""
Update the tile instruction files. Expected to be run run in the parent repository folder (e.g. containing the python script)
"""

import json
import pathlib
import geopandas
import argparse
import time
import numpy

def parse_args():
    """Expect a command line argument of the form:
    '--output_path path_to_output_folder'"""

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--search_string_latlon",
        metavar="str",
        default=r"*y_*h_0c_latlon.nc",
        action="store",
        help="thegeneral lat lon string format.",
    )
    
    parser.add_argument(
        "--search_string_nztm",
        metavar="str",
        default=r"*y_*h_0c_nztm.nc",
        action="store",
        help="the general nztm string format.",
    )
    
    parser.add_argument(
        "--domain_path",
        metavar="str",
        required=False,
        action="store",
        help="the domain file path string - The path to the domain file defining the flood domains to simulate over.",
    )
    
    parser.add_argument(
        "--output_path",
        metavar="str",
        required=False,
        action="store",
        help="the output path for cycl workflows - should eventually be the same at the output_path.",
    )
    
    parser.add_argument(
        "--return_period",
        metavar="int",
        required=False,
        action="store",
        help="the return period.",
    )

    return parser.parse_args()


def default_search_string(domain_path: str, output_path: str, return_period: int) -> geopandas.GeoDataFrame:
    """ Create / Update the domain summary dataframe for BG-FLOOD.
    """
    args = parse_args()
    domain = main(domain_path=domain_path, output_path=output_path, search_string_latlon=args.search_string_latlon, search_string_nztm=args.search_string_nztm, return_period=return_period)
    return domain
    

def main(domain_path: str, output_path: str, search_string_latlon: str, search_string_nztm: str, return_period: int) -> geopandas.GeoDataFrame:
    """ Create / Update the domain summary dataframe for design storms.
    """

    print('Update and return the domain design storm summary file!')
    
    ## Paths
    output_path = pathlib.Path(output_path)
    domain_path = pathlib.Path(domain_path)
    
    ## Domain file
    domain = geopandas.read_file(domain_path, layer='basin')
    
    ## Update the design storm status
    def check_domain_status(domain_name):
        domain_path = output_path / domain_name
        search_model_string = f"{return_period}y_*h_0c"
        model_paths = list(domain_path.glob(search_model_string))
        
        if len(model_paths) > 1:
            print(f"Warning more than one model detected for domain, {domain_name}, and return period {return_period}")
        
        model_path = model_paths[0] / "Design_Storm" if len(model_paths) > 0 else ""

        if len(model_paths) == 0 or not model_path.exists():
            status = "not started"
        elif any(model_path.glob(search_string_latlon)) and any(model_path.glob(search_string_nztm)):
            status = "completed"
        else:
            status = "incomplete"

        return status
    domain["design storm status"] = domain.apply(lambda row: check_domain_status(row["name"]), axis = 1)
    
    ## Write out the file with a time stamp
    timestamp = time.strftime('%b-%d-%Y_%H%M', time.localtime())
    #domain.to_file(catchments_path / f"{domain_path.stem}-{timestamp}.gpkg")
    
    return domain


if __name__ == "__main__":
    """ If called as script: Read in the args and launch the main function"""
    args = parse_args()
    main(domain_path=args.domain_path, search_string_latlon=args.search_string_latlon,
         search_string_nztm=args.search_string_nztm, output_path=args.output_path, return_period=args.return_period)
