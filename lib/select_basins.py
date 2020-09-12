def select_basins(grid_name_oce,truncation_type,res_num):
    '''
    AUTHORS:
    Jan Streffing		2020-09-12	Split off from main tool and reads yaml
    
    DESCRIPTION:
    This function selects the list of basins that are to be removed from the
    OpenIFS land sea mask based on a given ocean grid names. It also does the 
    same for the areas where the reverse is true and that need to be added to
    the ocean.
    
    
    INPUT:

    grid_name_oce       Name of the ocean model grid or mesh
    truncation_type     How the OpenIFS grid is pruned towards the poles
    res_num             Turncation number, minimum spherical harmonic wavelength

    
    RETURN:
        
    basin_removal       List of ocean basins that shall be turned to land
    coastline_addition  List of coastlines that shall be added to ocean
    '''

    import yaml
    
    with open(r'config/grid_modifications.yaml') as file:
        grids = yaml.full_load(file)
    basin_removal = grids["remove"][grid_name_oce]
    if (type(grids["add"][grid_name_oce])) == dict:
        coastline_addition = grids["add"][grid_name_oce][truncation_type+'-'+str(res_num)]
    else:
        coastline_addition = grids["add"][grid_name_oce]
        
    return [basin_removal,coastline_addition]