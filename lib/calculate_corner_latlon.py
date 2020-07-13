def calculate_corner_latlon(lats_list, lons_list, numlons_list, dlon_list,
                            lat_list):
    '''
    AUTHORS:
    Joakim Kjellson		2018		Original form
    Jan Streffing		2020-05-27	Split off from main tool


    DESCRIPTION:
    This function calculates the latitude and longitude values at the corners
    of the gridcells based on the center values. It also saves both the corner
    and center coordinates into a float32 arrays with oasis3-mct compatible
    structure
    

    INPUT:
    lons_list			List containing the longitude of each gridpoint
    lats_list			List containing the latitude of each gridpoint
    numlons_list		Number of longitude points for each latitude
    dlon_list			Longitude distance in degree at each latitude
    lat_list			List of latitude rows


    RETURN:
    center_lats			List of center latitudes
    center_lons			List of center longitudes
    crn_lats			List of corner latitudes
    crn_lons			List of corner longitudes
    '''

    import numpy as np

    # OASIS requires grids to be 2D, but IFS grid is 1D, so we give it an
    # extra dimension.
    center_lons = np.array(lons_list, dtype='float32')[np.newaxis, :]
    center_lats = np.array(lats_list, dtype='float32')[np.newaxis, :]
    nx = center_lons.shape[1]
    ny = 1

    print(' Size of grid: nx = %d, ny = %d' % (nx, ny))

    # Now we calculate longitudes/latitudes of corner points for each grid cell
    crn_lons = np.zeros((4, ny, nx))
    crn_lats = np.zeros((4, ny, nx))

    kk = 0 # cell index
    for ii, ni in enumerate(numlons_list):
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
        lat  = lat_list[ii]
        lons = np.arange(0, 360, dlon)

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

    for jj in range(ni):
        # corner 1: north-east
        crn_lons[0, 0, kk] = lons[jj] + dlon/2.
        crn_lats[0, 0, kk] = lat + dlat_n/2.

        # corner 2: north-west
        crn_lons[1, 0, kk] = lons[jj] - dlon/2.
        crn_lats[1, 0, kk] = lat + dlat_n/2.

        # corner 3: south-west
        crn_lons[2, 0, kk] = lons[jj] - dlon/2.
        crn_lats[2, 0, kk] = lat - dlat_s/2.

        # corner 4: south-east
        crn_lons[3, 0, kk] = lons[jj] + dlon/2.
        crn_lats[3, 0, kk] = lat - dlat_s/2.

        kk += 1

    # Make sure that longitudes are [-180, 180] and not [0, 360]
    center_lons = np.where( center_lons > 180, center_lons - 360, center_lons )
    crn_lons    = np.where( crn_lons > 180, crn_lons - 360, crn_lons )

    return (center_lats, center_lons, crn_lats, crn_lons)
