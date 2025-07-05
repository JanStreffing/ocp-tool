#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from shutil import copy2
from mpl_toolkits.basemap import Basemap

# Import the functions directly
# First, add the ocp-tool.py directory to the Python path
import sys
sys.path.append('/Users/jstreffi/software/ocp-tool')

# Then, we'll extract the functions from the ocp-tool.py file
def extract_functions_from_file():
    with open('/Users/jstreffi/software/ocp-tool/ocp-tool.py', 'r') as f:
        content = f.read()
    
    # Create a new namespace
    namespace = {}
    exec(content, namespace)
    
    # Get the functions we need
    modify_runoff_map = namespace['modify_runoff_map']
    plotting_runoff = namespace['plotting_runoff']
    return modify_runoff_map, plotting_runoff

modify_runoff_map, plotting_runoff = extract_functions_from_file()

# Parameters
res_num = 95
input_path_runoff = './input/runoff_map_default/'
output_path_runoff = './output/runoff_map_modified/'
grid_name_oce = 'CORE2'
manual_basin_removal = ['caspian-sea', 'black-sea']
verbose = True

# Create output directory if it doesn't exist
os.makedirs(output_path_runoff, exist_ok=True)

# Run the function
print("Running modify_runoff_map function...")
lons, lats = modify_runoff_map(res_num, input_path_runoff, output_path_runoff,
                              grid_name_oce, manual_basin_removal, verbose=verbose)

print("Function completed. Check the output file at:", output_path_runoff+'srunoff_maps_'+grid_name_oce+'.nc')

# Let's examine the calving point IDs to verify our changes
output_file = output_path_runoff+'srunoff_maps_'+grid_name_oce+'.nc'
if os.path.exists(output_file):
    rnffile = Dataset(output_file, 'r')
    calving = rnffile.variables[u'calving_point_id'][:]
    
    # Count occurrences of each ID in Antarctica region
    antarctica_ids = {}
    for la, lat in enumerate(lats):
        if lat < -55:  # Antarctic region
            for lo, lon in enumerate(lons):
                id_val = calving[la, lo]
                if id_val > 0:  # Ignore -1 and -2 values
                    if id_val in antarctica_ids:
                        antarctica_ids[id_val] += 1
                    else:
                        antarctica_ids[id_val] = 1
    
    print("\nAntarctica calving point IDs used and their counts:")
    for id_val, count in sorted(antarctica_ids.items()):
        print(f"ID {id_val}: {count} grid cells")
    
    rnffile.close()
else:
    print("Output file not found.")
