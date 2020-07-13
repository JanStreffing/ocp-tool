def extract_grid_data(lines):
    '''
    AUTHORS:
    Jan Streffing		2020-05-27	Split off from main tool


    DESCRIPTION:
    This function takes the raw reduced gaussian coordinate list and returns
    coordinate and neighbour distrance lists for latitude and
    longitude of every grindpoint, as well as the number of latitudes and
    longitudes
    

    INPUT:
    lines			List of unformated lines containing latitudes of grid


    RETURN:
    lons_list			List containing the longitude of each gridpoint
    lats_list			List containing the latitude of each gridpoint
    numlons_list		Number of longitude points for each latitude
    dlon_list			Longitude distance in degree at each latitude
    lat_list			List of latitude rows
    '''

    import numpy as np

    gridsize = 0
    lons_list       = []
    lats_list       = []
    numlons_list    = []
    dlon_list       = []
    lat_list        = []

    for line in lines[3:]:
        # read latitude number, number of longitudes for red. Gaussian and regular Gauss grids
        # convert from strings to floats
        print(line)
        _, red_points, _, lat = (float(z) for z in line.split())

        # longitudes for reduced Gaussian grid
        dlon = float(360)/red_points
        #The -0.000000001 deals with rounding errors a la 360./644*360=359.9999999999994
        lons = np.arange(0, 360-0.000000001, dlon)
        numlons_list.append(int(red_points))
        dlon_list.append(dlon)
        lat_list.append(lat)

        # set longitudes/latitudes for reduced Gaussian grid on this latitude
        lons_list.extend(lons)
        lats_list.extend([lat]*len(lons))
        gridsize += len(lons)

    return (lons_list, lats_list, numlons_list, dlon_list, lat_list)
