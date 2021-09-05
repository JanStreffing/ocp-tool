def read_oce(grid_name_oce,input_path_oce, lons_list, center_lats, 
             center_lons, crn_lats, crn_lons):
    '''
    AUTHORS:
    Jan Streffing		2021-09-05	Added ocean mesh reading


    DESCRIPTION:
    This function used the pyfesom2 toolbox to read the lon & lat of a fesom2 
    mesh
    

    INPUT:
    grid_name_oce       Name of the ocean model grid or mesh
    input_path_oce	    Path to folder containing ocean mesh files

    
    RETURN:
    lsm_oce			    ocean land sea mask
    '''

    import pyfesom2 as pf
    import numpy as np
    
    def sorting(l1, l2):
        # l1 and l2 has to be numpy arrays
        idx = np.argsort(l1)
        return l1[idx], l2[idx]
    
    mesh_ocean = pf.load_mesh(input_path_oce+'/'+grid_name_oce)
    
    x3, y3 = sorting(mesh_ocean.x2, mesh_ocean.y2)
    lons_list_uni = np.unique(lons_list)
    
    for idx in range(len(lons_list_uni)):
        condition1 = x3>crn_lons[0,0,idx] and x3<crn_lons[2,0,idx]
        
    return(mesh_ocean)