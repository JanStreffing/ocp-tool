def modify_lsm(gribfield, manual_basin_removal, manual_coastline_addition, 
               lsm_id, slt_id, cl_id, lons_list, center_lats, center_lons):
    '''
    AUTHORS:
    Jan Streffing		2020-07-22	Split off from main tool
    Jan Streffing       2020-09-11	Reading area coordinates from file


    DESCRIPTION:
    This function firstly uses the lake mask to remove lakes from the land sea
    mask and secondly, if set, uses a preselected list of basins to manually
    alter the lsm and slt fields further. It returns both the original mask
    as well as the modified one
    

    INPUT:

    gribfield           List of all grib field arrays loaded from file
    manual_basin_removal        List of basins that will be turned to land
    manual_coastline_addition   List of coastlines that will be turned to ocean
    lsm_id              Id of field containing land sea mask
    slt_id              Id of field containing soil type
    cl_id               Id of field containing lake mask
    lons_list			List containing the longitude of each gridpoint
    center_lats			List of center latitudes
    center_lons			List of center longitudes

    
    RETURN:

    lsm_binary_a        Altered land sea mask in shape for oasis3-mct
    lsm_binary_l        Original land sea mask in shape for oasis3-mct
    gribfield_mod       OpenIFS input fields with modified lsm and soil classes 
    '''


    import copy
    import numpy as np
    import pandas as pd
    
    # Had some problems with both fields changing before
    lsm_binary_l = copy.deepcopy(gribfield[lsm_id])
    lsm_binary_l = lsm_binary_l[np.newaxis, :]

    # Automatic lake removal with lakes mask
    gribfield_mod = gribfield[:]
    # Soil class of removed lakes is set to SANDY CLAY LOAM
    for i in np.arange (0, len(gribfield_mod[slt_id])-1):
        if gribfield_mod[cl_id][i] >= 0.5:
            gribfield_mod[slt_id][i] = 6
            gribfield_mod[lsm_id][i] = 1

    # Manual lake removal for basins definded in input file
    areas = pd.read_csv("config/basins_and_coastlines.csv", delim_whitespace=True)
    if manual_basin_removal:
        print('Removing: ', manual_basin_removal)
        for basin in manual_basin_removal:
            idx = areas[areas['Area'] == basin].index[0]
            west = areas['West'][idx]
            east = areas['East'][idx]
            south = areas['South'][idx]
            north = areas['North'][idx]
            for ia in range(len(lons_list)):
                if center_lats[0, ia] > west and center_lats[0, ia] < east and center_lons[0, ia] > north and center_lons[0, ia] < south:
                    gribfield_mod[lsm_id][ia] = 1
                    gribfield_mod[slt_id][ia] = 6
                  
    # Manual addition of wet points for basins definded in input file
    if manual_coastline_addition:
        print('Adding: ', manual_coastline_addition)
        for basin in manual_coastline_addition:
            idx = areas[areas['Area'] == basin].index[0]
            west = areas['West'][idx]
            east = areas['East'][idx]
            south = areas['South'][idx]
            north = areas['North'][idx]
            for ia in range(len(lons_list)):
                if center_lats[0, ia] > west and center_lats[0, ia] < east and center_lons[0, ia] > north and center_lons[0, ia] < south:
                    gribfield_mod[lsm_id][ia] = 0
                    gribfield_mod[slt_id][ia] = 0
    
    # Mask with lakes counting as land in correct format for oasis3-mct file
    lsm_binary_a = gribfield_mod[lsm_id]
    lsm_binary_a = lsm_binary_a[np.newaxis, :]

    return (lsm_binary_a,lsm_binary_l, gribfield_mod)