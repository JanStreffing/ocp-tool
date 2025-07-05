#!/usr/bin/env python3
import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from shutil import copy2
from mpl_toolkits.basemap import Basemap

# Parameters
res_num = 95
input_path_runoff = './input/runoff_map_default/'
output_path_runoff = './output/runoff_map_modified/'
grid_name_oce = 'CORE2'
manual_basin_removal = ['caspian-sea', 'black-sea']
verbose = True

# Create output directory if it doesn't exist
os.makedirs(output_path_runoff, exist_ok=True)

def modify_runoff_map(res_num, input_path_runoff, output_path_runoff,
                      grid_name_oce, manual_basin_removal, verbose=False):
    '''
    This function generates coordinate and areas fields based on
    the full and reduced gaussian gridfiles for a given truncation number.
    '''
    input_file_rnf = '%srunoff_maps.nc' % (input_path_runoff,)
    output_file_rnf = output_path_runoff+'srunoff_maps_'+grid_name_oce+'.nc'
    if os.path.exists(output_file_rnf):
        os.remove(output_file_rnf)
    copy2(input_file_rnf, output_file_rnf)

    rnffile = Dataset(output_file_rnf, 'r+')
    print(rnffile.variables.keys())

    drainage = rnffile.variables[u'drainage_basin_id'][:]
    arrival = rnffile.variables[u'arrival_point_id'][:]
    calving = rnffile.variables[u'calving_point_id'][:]

    # Set projection
    lons = rnffile.variables[u'lon'][:]
    lats = rnffile.variables[u'lat'][:]

    basin_offset_greenland = 7
    for basin in manual_basin_removal:
        if basin == 'caspian-sea':
            for lo, lon in enumerate(lons):
                if lon > 46 and lon < 56:
                    for la, lat in enumerate(lats):
                        if lat > 36 and lat < 47:
                            if drainage[la, lo] == -2:
                                drainage[la, lo] = 18+basin_offset_greenland
                                arrival[la, lo] = -1
                # adding artifical arrival points in the amazon discharge area
                # to close the global water budget
                if lon > 313 and lon < 314.5:
                    for la, lat in enumerate(lats):
                        if lat > 1 and lat < 2:
                            if arrival[la, lo] != -1:
                                arrival[la, lo] = 18+basin_offset_greenland

        if basin == 'black-sea':
            for lo, lon in enumerate(lons):
                #removing old basin
                if lon > 27 and lon < 43:
                    for la, lat in enumerate(lats):
                        if lat > 40.5 and lat < 48:
                            if drainage[la, lo] == -2:
                                drainage[la, lo] = 23+basin_offset_greenland
                                arrival[la, lo] = -1
                # adding new arrival points
                if lon > 25 and lon < 26.5:
                    for la, lat in enumerate(lats):
                        if lat > 38.5 and lat < 41:
                            if arrival[la, lo] != -1:
                                arrival[la, lo] = 23+basin_offset_greenland
                if lon > 23.5 and lon < 25:
                    for la, lat in enumerate(lats):
                        if lat > 38.5 and lat < 41:
                            if arrival[la, lo] != -1:
                                arrival[la, lo] = 28+basin_offset_greenland

    # Fix for Ob arrival
    for lo, lon in enumerate(lons):
        #removing old arrival points
        if lon > 60 and lon < 70:
            for la, lat in enumerate(lats):
                if lat > 60 and lat < 80:
                    if arrival[la, lo] == 13+basin_offset_greenland:
                        arrival[la, lo] = 6+basin_offset_greenland
        # adding new arrival points
        if lon > 72 and lon < 75:
            for la, lat in enumerate(lats):
                if lat > 65 and lat < 75:
                    if arrival[la, lo] == 6+basin_offset_greenland:
                        arrival[la, lo] = 13+basin_offset_greenland

    # Fix for Glacial calving maps
    # Antarctica
    for lo, lon in enumerate(lons):
        #removing old calving points
        for la, lat in enumerate(lats):
            if lat < -55:
                if calving[la, lo] == 66:
                    calving[la, lo] = -2

    # Adding new calving points with different IDs for different regions of Antarctica
    for lo, lon in enumerate(lons):
        # 1. Eastern Weddell Sea (ID: 73)
        if lon > 340 and lon < 20:
            for la, lat in enumerate(lats):
                if lat > -74 and lat < -65:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 73
                        
        # 2. Northwestern Weddell Sea (ID: 74)
        if lon > 300 and lon < 320:
            for la, lat in enumerate(lats):
                if lat > -70 and lat < -60:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 74
                        
        # 3. Antarctic Peninsula eastern drift route (ID: 75)
        if lon > 315 and lon < 325:
            for la, lat in enumerate(lats):
                if lat > -65 and lat < -55:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 75
                        
        # 4. Bellingshausen Sea (ID: 76)
        if lon > 260 and lon < 300:
            for la, lat in enumerate(lats):
                if lat > -75 and lat < -65:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 76
                        
        # 5. Amundsen Sea (ID: 77)
        if lon > 230 and lon < 260:
            for la, lat in enumerate(lats):
                if lat > -75 and lat < -65:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 77
                        
        # 6. Ross Sea (ID: 78)
        if lon > 170 and lon < 230:
            for la, lat in enumerate(lats):
                if lat > -75 and lat < -65:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 78
                        
        # 7. Eastern Antarctic offshore drift path (ID: 79)
        if lon > 320 and lon < 360:
            for la, lat in enumerate(lats):
                if lat > -60 and lat < -50:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 79
                        
        # 8. Wilkes Land coast (ID: 80)
        if lon > 110 and lon < 140:
            for la, lat in enumerate(lats):
                if lat > -68 and lat < -60:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 80
                        
        # 9. Prydz Bay & Davis Sea (ID: 81)
        if lon > 70 and lon < 90:
            for la, lat in enumerate(lats):
                if lat > -70 and lat < -60:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 81
                        
        # 10. Cosmonaut Sea (ID: 82)
        if lon > 40 and lon < 60:
            for la, lat in enumerate(lats):
                if lat > -70 and lat < -60:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 82
                        
        # 11. Eastern Antarctic deep south (ID: 83)
        if lon > 80 and lon < 130:
            for la, lat in enumerate(lats):
                if lat > -80 and lat < -70:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 83
                        
        # 12. Western Ross Sea (ID: 84)
        if lon > 150 and lon < 170:
            for la, lat in enumerate(lats):
                if lat > -77 and lat < -67:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 84
                        
        # 13. Eastern Ross Sea periphery (ID: 85)
        if lon > 190 and lon < 210:
            for la, lat in enumerate(lats):
                if lat > -72 and lat < -62:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 85
                        
        # 14. Marie Byrd Land coast (ID: 86)
        if lon > 220 and lon < 240:
            for la, lat in enumerate(lats):
                if lat > -78 and lat < -68:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 86
                        
        # 15. Western Antarctic deep south (ID: 87)
        if lon > 270 and lon < 330:
            for la, lat in enumerate(lats):
                if lat > -85 and lat < -75:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 87
                        
        # 16. Southern Weddell Sea (ID: 88)
        if lon > 340 and lon < 360:
            for la, lat in enumerate(lats):
                if lat > -80 and lat < -74:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 88
                        
        # 17. Southeastern Weddell Sea (ID: 89)
        if lon > 0 and lon < 30:
            for la, lat in enumerate(lats):
                if lat > -76 and lat < -70:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 89
                        
    # Greenland
    for lo, lon in enumerate(lons):
        # adding new calving points
        if lon > 300 and lon < 310:
            for la, lat in enumerate(lats):
                if lat > 50 and lat < 60:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 1

    # Saving results
    rnffile.variables[u'drainage_basin_id'][:] = drainage
    rnffile.variables[u'arrival_point_id'][:] = arrival
    rnffile.variables[u'calving_point_id'][:] = calving
    rnffile.close()

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

    return (lons, lats)

# Run the function
print("Running modify_runoff_map function to update Antarctica calving regions...")
lons, lats = modify_runoff_map(res_num, input_path_runoff, output_path_runoff,
                             grid_name_oce, manual_basin_removal, verbose=verbose)

print("\nFunction completed successfully!")
print(f"The Antarctic calving regions have been updated in the output file:")
print(f"{output_path_runoff}srunoff_maps_{grid_name_oce}.nc")
print("\nThe regions now use IDs 73-89 instead of the single ID 66 for better routing precision.")
