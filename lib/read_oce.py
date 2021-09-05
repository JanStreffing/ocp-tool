def read_oce(grid_name_oce,input_path_oce, lons_list, center_lats, center_lons, crn_lats, crn_lons):
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
    lsm_oce			    ocean land sea mask interpolared to OpenIFS grid
    '''
    
    def outer_loop(mesh_ocean,crn_lons,idx,lsm):
        lon_condition = np.logical_and(mesh_ocean.x2<crn_lons[0,0,idx] , mesh_ocean.x2>crn_lons[2,0,idx])
        reduced_idx=np.where(lon_condition)
        for ridx in reduced_idx[0]:
            lat_condition = np.logical_and(mesh_ocean.y2<crn_lats[0,0,idx] , mesh_ocean.y2>crn_lats[2,0,idx])
            if np.any(lat_condition):
                lsm[idx]=1
            else:
                lsm[idx]=0
        return lsm





    import pyfesom2 as pf
    import numpy as np
    from tqdm import tqdm
    import dask
    from dask.delayed import delayed
    from dask.diagnostics import ProgressBar
    
    lsm=[]
    ta=[]
    
    mesh_ocean = pf.load_mesh(input_path_oce+'/'+grid_name_oce)
        
    for idx in tqdm(range(len(lons_list))):
        lsm = dask.delayed(outer_loop)(mesh_ocean,crn_lons,idx,lsm)
        ta.append(lsm)
    with ProgressBar():
        lsm_oce = dask.compute(ta)
    
    return(lsm_oce)

    '''lon_condition = np.logical_and(mesh_ocean.x2<crn_lons[0,0,idx] , mesh_ocean.x2>crn_lons[2,0,idx])
        reduced_idx=np.where(lon_condition)
        for ridx in reduced_idx[0]:
            lat_condition = np.logical_and(mesh_ocean.y2<crn_lats[0,0,idx] , mesh_ocean.y2>crn_lats[2,0,idx])
            if np.any(lat_condition):
                lsm_oce[ridx]=1'''
                
    
 
