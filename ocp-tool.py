#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
The OpenIFS coupling preparation tool (ocp-tool) prepares input files for a 
coupled OpenIFS-FESOM2 or OpenIFS-NEMO climate simulation


To this end it performs three tasks:
1) Modifing the OpenIFS input files to fit the land sea mask and soil types
    to the ocean models land sea mask
2) Generating appropriate OASIS3-MCT input files to fit the modified land 
    sea mask
3) Modifing the runoff-mapper drainage basin and arrival point file to fit 
    the modified land sea mask


To function the scrip therefore needs the following input files:
1) Grid information txt file for the full Gaussian grid of the chosen 
    trucation number. This fileset comes with the tool.
2) Grid information txt file for the reduced Gaussian grid This fileset comes 
    with the tool.
3) OpenIFS gridpoint input file (ICMGG${EXP_NAME}INIT) containing default 
    land-sea mask that will be modified. can be requested from ECMWF: 
    openifs-support@ecmwf.int
4) Runoff-mapper input files (runoff_maps.nc) file containing default 
    drainage basin and arrival point fields. This fileset is part of the 
    EC-Earth input files and available from: 


@author: Jan Streffing (jan.streffing@awi.de), August 2019
'''

#-----------------------------------------------------------------------------
# Setup
#-----------------------------------------------------------------------------

# Loadin modules
import numpy as np
import matplotlib.pyplot as plt
from gribapi import *
from netCDF4 import Dataset

# Setting constants
earth_radius = 6371. * 1e3 #[m]
longline = ' \n ==================================================  \n'


#-----------------------------------------------------------------------------
# Function definitions
#-----------------------------------------------------------------------------

def read_grid_file(res_num, input_path_reduced_grid, input_path_full_grid):
    '''
    This function reads the reduced gaussian gridfile and returns it as a raw 
    field
    '''
    
    if 1: # linear truncation (T = NN * 2 - 1)
        NN = res_num/2 + 1
    elif 0:
        NN = res_num/4 + 1
    grid_txt = '%s/n%d_reduced.txt' % (input_path_reduced_grid,NN)
    print(' Read grid from file: %s ' % (grid_txt,) ) 
   
    print(longline)
    print(' Reading gridfiles for T%d ' % (res_num))
    print(longline)
   
    fin = open(grid_txt,'r')
    lines = fin.readlines()
    return (lines, NN)



def extract_grid_data(lines,NN):
    '''
    This function takes the raw reduced gaussian coordinate list and returns 
    coordinate and neighbour distrance lists for latitude and 
    longitude of every grindpoint, as well as the number of latitudes and 
    longitudes
    '''
    gridsize = 0
    lons_list       = []  # longitudes of each gridpoint
    lats_list       = []  # latitudes of each gridpoint
    numlons_list    = []  # number of longitude points for each latitude
    dlon_list       = []  # longitude distance in degree at each latitude
    lat_list        = []  # list of latitudes 
   
    for line in lines[3:]:
        # read latitude number, number of longitudes for red. Gaussian and regular Gauss grids
        # convert from strings to floats
        print(line)
        if NN == 320:
            [lat_num, red_points, oct_points, reg_points, lat] = [float(z) for z in line.split()]
        else:
            [lat_num, red_points, reg_points, lat] = [float(z) for z in line.split()]
         
        # longitudes for reduced Gaussian grid
        dlon = 360./red_points
        lons = np.arange(0,360,dlon)
        numlons_list.append(int(red_points))
      
        dlon_list.append(dlon)
        lat_list.append(lat)
      
        # set longitudes/latitudes for reduced Gaussian grid on this latitude
        for i in range(0,lons.shape[0]):
            lons_list.append(lons[i])
            lats_list.append(lat)
            gridsize += 1
      
    return (lons_list, lats_list, numlons_list, dlon_list, lat_list)
       

def calculate_corner_latlon(lats_list, lons_list, numlons_list, dlon_list, 
                            lat_list):
    '''
    This function calculates the latitude and longitude values at the corners
    of the gridcells based on the center values. It also saves both the corner 
    and center coordinates into a float32 arrays with oasis3-mct compatible 
    structure
    '''

    # OASIS requires grids to be 2D, but IFS grid is 1D, so we give it an 
    # extra dimension. 
    center_lons = np.array(lons_list,dtype='float32')[np.newaxis,:]
    center_lats = np.array(lats_list,dtype='float32')[np.newaxis,:]
    nx = center_lons.shape[1]
    ny = 1
   
    print(' Size of grid: nx = %d, ny = %d' % (nx,ny))
   
    # Now we calculate longitudes/latitudes of corner points for each grid cell
    crn_lons = np.zeros((4,ny,nx))
    crn_lats = np.zeros((4,ny,nx))
   
    kk = 0 # cell index
    for ii in range(0,len(numlons_list)):
        '''
        Layout of the four corners
        
        2 ---------- 1
        |            |
        |            |
        |            |
        3 -----------4
        
        ^ y
        |
        |
        |
        ----> x
        '''
        
        dlon = dlon_list[ii]
        ni   = numlons_list[ii]
        lat  = lat_list[ii]
        lons = np.arange(0,360,dlon)
      
        #     NP --- j=1 ---|--- j=2 ---|--- j=3 ---|--- j=n --- SP
        #                           <-dlat_n-> <-dlat_s->
      
        # if first latitude, the previous point was north pole
        if ii == 0:
            dlat_n = 90 - lat
            dlat_s = (lat - lat_list[ii+1]) / 2.
      
        # if last latitude, the next point is south pole
        elif ii == len(numlons_list)-1:
            dlat_n = (lat_list[ii-1] - lat) / 2.
            dlat_s = lat + 90
      
        else:
            dlat_n = (lat_list[ii-1] - lat) / 2.
            dlat_s = (lat - lat_list[ii+1]) / 2.
      
      
    for jj in range(0,ni):
        # corner 1: north-east
        crn_lons[0,0,kk] = lons[jj] + dlon/2.
        crn_lats[0,0,kk] = lat + dlat_n/2.
         
        # corner 2: north-west
        crn_lons[1,0,kk] = lons[jj] - dlon/2.
        crn_lats[1,0,kk] = lat + dlat_n/2.
         
        # corner 3: south-west
        crn_lons[2,0,kk] = lons[jj] - dlon/2.
        crn_lats[2,0,kk] = lat - dlat_s/2.
         
        # corner 4: south-east
        crn_lons[3,0,kk] = lons[jj] + dlon/2.
        crn_lats[3,0,kk] = lat - dlat_s/2.
         
        kk += 1  
        
    return (center_lats, center_lons, crn_lats, crn_lons)


def calculate_area(center_lons, numlons_list, dlon_list, lat_list):
    '''
    This function calculates the area of the gridcells based on the center 
    values and saves in into a float32 array with oasis3-mct compatible 
    structure
    '''

    # OASIS requires grids to be 2D, but IFS grid is 1D, so we give it an 
    # extra dimension. 
    nx = center_lons.shape[1]
    ny = 1
   
     # Now we calculate the cell area of each cell
    gridcell_area = np.zeros((ny,nx))
           
    kk = 0 # cell index
    for ii in range(0,len(numlons_list)):
        
        dlon = dlon_list[ii]
        ni   = numlons_list[ii]
        lat  = lat_list[ii]
        lons = np.arange(0,360,dlon)
      
        #     NP --- j=1 ---|--- j=2 ---|--- j=3 ---|--- j=n --- SP
        #                           <-dlat_n-> <-dlat_s->
      
        # if first latitude, the previous point was north pole
        if ii == 0:
            dlat_n = 90 - lat
            dlat_s = (lat - lat_list[ii+1]) / 2.
      
        # if last latitude, the next point is south pole
        elif ii == len(numlons_list)-1:
            dlat_n = (lat_list[ii-1] - lat) / 2.
            dlat_s = lat + 90
      
        else:
            dlat_n = (lat_list[ii-1] - lat) / 2.
            dlat_s = (lat - lat_list[ii+1]) / 2.
            
        # Grid cell areas in m2 width in latitude of cell is dlat_n + dlat_s
        dx = dlon * np.pi/180. * earth_radius * np.cos( np.pi/180. * lat )
        dy = (dlat_n + dlat_s) * np.pi/180. * earth_radius
        area = dx * dy
      
        for jj in range(0,ni):
            gridcell_area[0,kk] = area
            kk += 1
   
    return(gridcell_area)
    
    
    
def read_lsm(res_num, input_path_oifs, output_path_oifs, exp_name_oifs, num_fields):
    '''
    This function reads the oifs input file in grib format and save it into a 
    list of numpy arrays.
    '''
    print(' Opening Grib inpute file: %s ' % (input_path_oifs,))
    input_file_oifs = input_path_oifs + 'ICMGG' + exp_name_oifs + 'INIT'
    gid = [None] * num_fields
    gribfield = [None] * num_fields
    with open(input_file_oifs,'r+') as f:
        keys = ['N','shortName']
      
        for i in range(0,num_fields):
            gid[i] = grib_new_from_file(f)
            if gid[i] is None:
                break
 
            for key in keys:
                if not grib_is_defined(gid[i], key):
                    raise ValueError("Key '%s' was not defined" % key)
                print('%s=%s' % (key, grib_get(gid[i], key)))
         
            shortName = grib_get(gid[i],'shortName')
            
            if shortName == 'lsm':
                lsm_id = i
            if shortName == 'slt':
                slt_id = i
                
            nres = grib_get(gid[i],'N')
            gribfield[i] = grib_get_values(gid[i]) 
        f.close
    return (gribfield, lsm_id, slt_id)
         


def read_lake(res_num, input_path_lake):
    '''
    This function reads the lakemask input file in grib format and save it into a 
    a numpy array.
    '''
    print(' Opening lake inpute file: %s ' % (input_path_lake,))
    input_file_lake = input_path_lake + 'clake_' + str(res_num)
    gid = [None] * 1
    lakes = [None] * 1
    with open(input_file_lake,'r+') as f:
        keys = ['N','shortName']
      
        for i in range(0,1):
            gid[i] = grib_new_from_file(f)
            if gid[i] is None:
                break
 
            for key in keys:
                if not grib_is_defined(gid[i], key):
                    raise ValueError("Key '%s' was not defined" % key)
                print('%s=%s' % (key, grib_get(gid[i], key)))
         
            shortName = grib_get(gid[i],'shortName')
                
            nres = grib_get(gid[i],'N')
            lakes[i] = grib_get_values(gid[i]) 
        f.close
    return (lakes)



def modify_lsm(gribfield, lakes lsm_id, slt_id):
    
    lsm_binary=gribfield[lsm_id]
            
    #############################################################
    #                       Lake REMOVAL                        #
    #############################################################

        
    
def plotting():
    
        
    
def generate_coord_area(res_num, input_path_reduced_grid, input_path_full_grid):
    '''
    This function generates coordinate and areas fields based on
    the full and reduced gaussian gridfiles for a given truncation number. 
    '''
    lines, NN = read_grid_file(res_num, input_path_reduced_grid, input_path_full_grid)
    lons_list, lats_list, numlons_list, dlon_list, lat_list = extract_grid_data(lines, NN)
    center_lats, center_lons, crn_lats, crn_lons = calculate_corner_latlon(lats_list, lons_list, numlons_list, dlon_list, lat_list)
    gridcell_area = calculate_area(center_lons, numlons_list, dlon_list, lat_list)

    return (center_lats, center_lons, crn_lats, crn_lons, gridcell_area, lons_list)
    
def process_lsm(res_num, input_path_oifs, output_path_oifs, exp_name_oifs, 
                grid_name_oce, num_fields, input_path_lake):
    '''
    This function first reads, modifies and finally saves the new land 
    sea mask. Every step is mirrored for the soil type file as it has to be 
    modified in the exact same locations
    '''
    
    gribfield, lsm_id, slt_id = read_lsm(res_num, input_path_oifs, output_path_oifs, exp_name_oifs, num_fields)
    lakes = read_lake(res_num, input_path_lake)
    gribfield_mod = modify_lsm(gribfield, lsm_id, slt_id)
    #write_lsm
    
    
    

def write_oasis_files(res_num, output_path_oasis, grid_name_oce):
    '''
    This function generates coordinate and areas fields based on
    the full and reduced gaussian gridfiles for a given truncation number. 
    '''



def modify_runoff_map(res_num, input_path_runoff, output_path_runoff,
                      grid_name_oce):
    '''
    This function generates coordinate and areas fields based on
    the full and reduced gaussian gridfiles for a given truncation number. 
    '''
    
    
    
    
    
    
    
    
#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------
    
if __name__ == '__main__':
    '''
    Main program in which the tool configuration and function calls are located
    Please configure as needed
    '''
    
    # Input file directories. Modify or place files in subfolders
    input_path_reduced_grid = 'input/gaussian_grids_reduced/'
    input_path_full_grid = 'input/gaussian_grids_full/'
    input_path_oifs = 'input/openifs_input_default/'
    input_path_runoff = 'input/runoff_map_default/'
    input_path_lake = 'input/lakefiles/'
    
    # Output file directories.
    output_path_oifs = 'output/openifs_input_modified/' 
    output_path_runoff = 'output/runoff_map_modified/' 
    output_path_oasis = 'output/oasis_mct3/' 
    
    # Truncation number of desired OpenIFS grid 
    # Choose all the ones you need [63 95 159 255 319 399 511 799 1279]
    resolution_list = [159]
    
    # OpenIFS experiment name. This is a 4 digit code that is part of the 
    # name of the ICMGG????INIT file you got from EMCWF
    exp_name_oifs = 'h6mv'
    # I have not yet found a way to determine automatically the number of 
    # fields in the ICMGG????INIT file. Set it correctly!
    num_fields = 41 
    
    # Name of ocean model grid. So far supported are: 
    # FESOM2: CORE2, MR
    # NEMO:
    grid_name_oce = 'MR'
    
    for res_num in resolution_list:
        center_lats, center_lons, crn_lats, crn_lons, gridcell_area, lons_list = generate_coord_area(res_num, input_path_reduced_grid, input_path_full_grid)
        process_lsm(res_num, input_path_oifs, output_path_oifs, exp_name_oifs, grid_name_oce, num_fields, input_path_lake)
        write_oasis_files(res_num, output_path_oasis, grid_name_oce)
        modify_runoff_map(res_num, input_path_runoff, output_path_runoff, grid_name_oce)

