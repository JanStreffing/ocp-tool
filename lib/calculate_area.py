def calculate_area(center_lons, numlons_list, dlon_list, lat_list):
    '''
    AUTHORS:
    Jan Streffing		2020-07-22	Split off from main tool


    DESCRIPTION:
    This function calculates the area of the gridcells based on the center
    values and saves them into a float32 array with oasis3-mct compatible
    structure
    

    INPUT:
    center_lons			List of center longitudes
    numlons_list		Number of longitude points for each latitude
    dlon_list			Longitude distance in degree at each latitude
    lat_list			List of latitude rows


    RETURN:
    gridcell_area		vector containing the area of each grid cell
    '''

    import numpy as np
    
    # Setting constants
    earth_radius = 6371. * 1e3 #[m]

    # OASIS requires grids to be 2D, but IFS grid is 1D, so we give it an
    # extra dimension.
    nx = center_lons.shape[1]
    ny = 1

     # Now we calculate the cell area of each cell
    gridcell_area = np.zeros((ny, nx))

    kk = 0 # cell index
    for ii, ni in enumerate(numlons_list):

        dlon = dlon_list[ii]
        lat  = lat_list[ii]

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

        for jj in range(ni):
            gridcell_area[0, kk] = area
            kk += 1

    return(gridcell_area)