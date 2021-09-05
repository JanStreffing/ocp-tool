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


To function, the script therefore needs the following input files:
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


If you have trouble getting this tool to work in your python environment
you may try loading the environment.yaml with:
    conda env create -f environment.yml


@author: Jan Streffing (jan.streffing@awi.de), August 2019
'''

from __future__ import division

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from shutil import copy2
from lib import (read_grid_file, extract_grid_data, calculate_corner_latlon, 
                calculate_area, write_red_point_file, read_lsm, modify_lsm,
                select_basins, write_lsm,write_oasis_files)


#-----------------------------------------------------------------------------
# Setup
#-----------------------------------------------------------------------------




#-----------------------------------------------------------------------------
# Function definitions
#-----------------------------------------------------------------------------


def plotting_lsm(res_num, lsm_binary_l, lsm_binary_a, center_lats, center_lons):
    '''
    This function plots the final land sea mask
    '''

    fig3 = plt.figure(figsize=(12, 7))
    ax3  = fig3.add_subplot(111)
    xptsa = center_lons[np.round(lsm_binary_a[:, :])<1]
    yptsa = center_lats[np.round(lsm_binary_a[:, :])<1]
    xptsl = center_lons[np.round(lsm_binary_l[:, :])<1]
    yptsl = center_lats[np.round(lsm_binary_l[:, :])<1]
    ax3.scatter(xptsl, yptsl, s=1.5, color='red')
    ax3.scatter(xptsa, yptsa, s=2)
    figname = 'output/plots/land_points_T%d.png' % (res_num,)
    fig3.savefig(figname, format='png')


def generate_coord_area(res_num, input_path_reduced_grid, input_path_full_grid, truncation_type):
    '''
    This function generates coordinate and areas fields based on
    the full and reduced gaussian gridfiles for a given truncation number.
    '''

    lines, NN = read_grid_file.read_grid_file(res_num, input_path_reduced_grid, input_path_full_grid, truncation_type)
    lons_list, lats_list, numlons_list, dlon_list, lat_list = extract_grid_data.extract_grid_data(lines)
    center_lats, center_lons, crn_lats, crn_lons = calculate_corner_latlon.calculate_corner_latlon(lats_list, lons_list, numlons_list, dlon_list, lat_list)
    gridcell_area = calculate_area.calculate_area(center_lons, numlons_list, dlon_list, lat_list)
    write_red_point_file.write_red_point_file(lats_list, lons_list, output_path_oifs, truncation_type, NN)

    return (center_lats, center_lons, crn_lats, crn_lons, gridcell_area, lons_list, NN)


def process_lsm(res_num, input_path_oifs, output_path_oifs, exp_name_oifs,
                grid_name_oce, input_path_oce, num_fields, basin_removal, 
                coastline_addition, lons_list, center_lats, center_lons):
    '''
    This function first reads, modifies and finally saves the new land
    sea mask. Every step is mirrored for the soil type file as it has to be
    modified in the exact same locations
    '''

    gribfield, lsm_id, slt_id, cl_id, gid = read_lsm.read_lsm(res_num, input_path_oifs, 
                                                     output_path_oifs, 
                                                     exp_name_oifs, num_fields)
    lsm_ocean = read_oce.read_oce(grid_name_oce,input_path_oce) 
    lsm_binary_a, lsm_binary_l, gribfield_mod = modify_lsm.modify_lsm(gribfield, 
                                                           basin_removal, 
                                                           coastline_addition, 
                                                           lsm_id, slt_id, cl_id, 
                                                           lons_list, center_lats, 
                                                           center_lons)
    write_lsm.write_lsm(gribfield_mod, input_path_oifs, output_path_oifs, exp_name_oifs, 
              grid_name_oce, num_fields, gid)
    return (lsm_binary_a,lsm_binary_l)





def modify_runoff_map(res_num, input_path_runoff, output_path_runoff,
                      grid_name_oce, basin_removal):
    '''
    This function generates coordinate and areas fields based on
    the full and reduced gaussian gridfiles for a given truncation number.
    '''
    input_file_rnf = '%srunoff_maps.nc' % (input_path_runoff,)
    output_file_rnf = '%srunoff_maps.nc' % (output_path_runoff,)
    if os.path.exists(output_file_rnf):
        os.remove(output_file_rnf)
    copy2(input_file_rnf, output_file_rnf)

    rnffile = Dataset(output_file_rnf, 'r+')
    print (rnffile.variables.keys())

    drainage = rnffile.variables[u'drainage_basin_id'][:]
    arrival = rnffile.variables[u'arrival_point_id'][:]

    # Set projection
    lons = rnffile.variables[u'lon'][:]
    lats = rnffile.variables[u'lat'][:]

    for basin in basin_removal:

        if basin == 'caspian-sea':
            for lo, lon in enumerate(lons):
                if lon > 46 and lon < 56:
                    for la, lat in enumerate(lats):
                        if lat > 36 and lat < 47:
                            if drainage[la, lo] == -2:
                                drainage[la, lo] = 18
                                arrival[la, lo] = -1
                # adding artifical arrival points in the amazon discharge area
                # to close the global water budget
                if lon > 313 and lon < 314.5:
                    for la, lat in enumerate(lats):
                        if lat > 1 and lat < 2:
                            if arrival[la, lo] != -1:
                                arrival[la, lo] = 18

        if basin == 'black-sea':
            for lo, lon in enumerate(lons):
                #removing old basin
                if lon > 27 and lon < 43:
                    for la, lat in enumerate(lats):
                        if lat > 40.5 and lat < 48:
                            if drainage[la, lo] == -2:
                                drainage[la, lo] = 23
                                arrival[la, lo] = -1
                # adding new arrival points
                if lon > 25 and lon < 26.5:
                    for la, lat in enumerate(lats):
                        if lat > 38.5 and lat < 41:
                            if arrival[la, lo] != -1:
                                arrival[la, lo] = 23
                if lon > 23.5 and lon < 25:
                    for la, lat in enumerate(lats):
                        if lat > 38.5 and lat < 41:
                            if arrival[la, lo] != -1:
                                arrival[la, lo] = 28

    # Saving results
    rnffile.variables[u'drainage_basin_id'][:] = drainage
    rnffile.variables[u'arrival_point_id'][:] = arrival
    rnffile.close()

    plotting_runoff(drainage, arrival, lons, lats)

    return (lons, lats)


def plotting_runoff(drainage, arrival, lons, lats):

    # Split data and concatenate in reverse order to turn by 180° to Prime meridian
    ds1, ds2 = np.hsplit(np.squeeze(drainage), 2)
    drainage_cat = np.concatenate((ds2, ds1), axis=1)
    ds1, ds2 = np.hsplit(np.squeeze(arrival), 2)
    arrival_cat = np.concatenate((ds2, ds1), axis=1)

    lons = lons-180
    lon_0 = lons.mean()
    lat_0 = lats.mean()

    m = Basemap(llcrnrlon=-60., llcrnrlat=-10, urcrnrlon=-30., urcrnrlat=20., \
            resolution='l', area_thresh=1000., projection='cyl')

    #Use meshgrid to create 2D arrays from coordinates
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    fig1 = plt.figure(figsize=(12, 8))
    cmap = plt.cm.flag
    cs = m.pcolor(xi, yi, arrival_cat, cmap=cmap)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90., 120., 45.))
    m.drawmeridians(np.arange(0., 360., 90.))

    lon_0 = lons.mean()
    lat_0 = lats.mean()

    m = Basemap(llcrnrlon=20., llcrnrlat=30, urcrnrlon=80., urcrnrlat=50., \
            resolution='l', area_thresh=1000., projection='poly', \
            lat_0=0., lon_0=20.)

    #Use meshgrid to create 2D arrays from coordinates
    lon, lat = np.meshgrid(lons, lats)
    xi, yi = m(lon, lat)

    fig1 = plt.figure(figsize=(12, 8))
    cs = m.pcolor(xi, yi, drainage_cat, cmap=cmap)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90., 120., 45.))
    m.drawmeridians(np.arange(0., 360., 90.))

    fig1 = plt.figure(figsize=(12, 8))
    cs = m.pcolor(xi, yi, arrival_cat, cmap=cmap)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90., 120., 45.))
    m.drawmeridians(np.arange(0., 360., 90.))


def modify_runoff_lsm(res_num, grid_name_oce, basin_removal, lons, lats,
                      output_path_oasis):
    '''
    This function generates coordinate and areas fields based on
    the full and reduced gaussian gridfiles for a given truncation number.
    '''

    # Editing runoff mapper lsm in oasis3-mct masks file
    filename = '%smasks.nc' % (output_path_oasis,)
    oasis = Dataset(filename, 'r+')

    RnfA = oasis.variables[u'RnfA.msk'][:]
    RnfO = oasis.variables[u'RnfO.msk'][:]

    for basin in basin_removal:

        if basin == 'caspian-sea':
            for lo, lon in enumerate(lons):
                if lon > 46 and lon < 56:
                    for la, lat in enumerate(lats):
                        if lat > 36 and lat < 47:
                            RnfA[la, lo] = 0
                            RnfO[la, lo] = 1

    # Saving altered runoff mapper lsm
    oasis.variables[u'RnfA.msk'][:] = RnfA
    oasis.variables[u'RnfO.msk'][:] = RnfO
    oasis.close()


#-----------------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------------

if __name__ == '__main__':
    '''
    Main program in which the tool configuration and function calls are located
    Please configure as needed.
    '''

    # Truncation number of desired OpenIFS grid.
    # Choose the one you need e.g. [63, 95, 159, 255, 319, 399, 511, 799, 1279]
    res_num = 159

    # Choose type of trucation. linear or cubic-octahedral
    truncation_type = 'cubic-octahedral'

    # OpenIFS experiment name. This 4 digit code is part of the name of the
    # ICMGG????INIT file you got from EMCWF
    exp_name_oifs = 'h9wu' #default for linear
    #exp_name_oifs = 'h9wu'#default for cubic-octahedral
    # I have not yet found a way to determine automatically the number of
    # fields in the ICMGG????INIT file. Set it correctly or stuff will break!
    num_fields = 50

    # Name of ocean model grid. So far supported are:
    # FESOM2: CORE2, MR, HR;  NEMO:
    # Important: If you chose a supported ocean grid, manual removal of basins
    # for that grid will be applied. If you chose an new grid and wish to
    # do manual basin removal, list them in config/couplings.yaml. If you
    # want to remove or add a basin not yet supported (e.g.) for paleo 
    # simulations, add the basin to config/basins_and_coastlines.csv
    grid_name_oce = 'CORE2'

    # There is automatic removal of lakes via the lake file. To remove larger
    # features, e.g. coastal seas for low res or paleo simulations manually 
    # list them in config/grid_modifications.yaml
    # You can see available grid modifications and add your own in
    # config/basins_and_coastlines.csv
    # Besides removing basins you can also add coastlines with this routine.
    basin_removal, coastline_addition = select_basins.select_basins(grid_name_oce,
                                                                    truncation_type,
                                                                    res_num)

    # Find working directory
    dir_path = os.path.dirname(os.path.realpath(__file__))

    # IO file directories. Should work without touching it. If you want 
    # diffent input folders you can make your selections here.
    if truncation_type == 'cubic-octahedral':
        input_path_reduced_grid = 'input/gaussian_grids_octahedral_reduced/'
    elif truncation_type == 'linear':
        input_path_reduced_grid = 'input/gaussian_grids_linear_reduced/'
    else:
        sys.exit('truncation type not recognized')
    input_path_full_grid = 'input/gaussian_grids_full/'
    input_path_oifs = 'input/openifs_input_default/'
    input_path_oce = 'input/ocean_input_default/'
    input_path_runoff = 'input/runoff_map_default/'
    output_path_oifs = 'output/openifs_input_modified/'
    output_path_runoff = 'output/runoff_map_modified/'
    output_path_oasis = 'output/oasis_mct3_input/'

    
    center_lats, center_lons, \
    crn_lats, crn_lons, \
    gridcell_area, lons_list, \
    NN = generate_coord_area(res_num,
                                 input_path_reduced_grid, input_path_full_grid,
                                 truncation_type)

    lsm_binary_a,lsm_binary_l = process_lsm(res_num, input_path_oifs, output_path_oifs,
                                 exp_name_oifs, grid_name_oce, input_path_oce, num_fields,
                                 basin_removal, coastline_addition, lons_list,
                                 center_lats, center_lons)

    write_oasis_files.write_oasis_files(res_num,
                          output_path_oasis, dir_path, grid_name_oce,
                          center_lats, center_lons, crn_lats, crn_lons, gridcell_area,
                          lsm_binary_a, lsm_binary_l, NN, input_path_runoff)

    lons, lats = modify_runoff_map(res_num, input_path_runoff, output_path_runoff,
                                       grid_name_oce, basin_removal)

    modify_runoff_lsm(res_num, grid_name_oce, basin_removal, lons, lats,
                          output_path_oasis)

    plotting_lsm(res_num, lsm_binary_l, lsm_binary_a, center_lats, center_lons)



