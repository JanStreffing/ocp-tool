def read_grid_file(res_num, truncation_type, exp_name_oifs, input_path_oifs):
    '''
    AUTHORS:
    Joakim Kjellson     2020-06-18  Original
    Jan Streffing		2021-09-05	Split off from main tool; removed txt reading


    DESCRIPTION:
    This function read the reduced gaussian from ICMGG and returns it as a raw
    field
    

    INPUT:
    res_num 			        OpenIFS truncation number
    truncation_type             OpenIFS truncation type (linear, quadratic, cubic octahedral)

    
    RETURN:
    lines			List of unformated lines containing latitudes of grid
    NN				Name of OpenIFS IO files according to ECMWF naming convention
    '''

    import os
    import numpy as np


    if truncation_type == 'linear':
        # linear truncation (T = NN * 2 - 1)
        NN = res_num/2 + 0.5
        
    elif truncation_type == 'cubic-octahedral':
        # cubic octahedral truncation (T = NN - 1)
        NN = res_num + 1

    file = '%s/ICMGG%sINIT' % (input_path_oifs, exp_name_oifs) 
    
    print(' Read grid from file: %s ' % (file,) )
    print(' Reading gridfiles for T%d ' % (res_num))
    
    latitudes = []
    nlongitudes = []
   
    # write grid description to file
    # only need to do this once
    os.system('cdo griddes %s > griddes.txt' % (file,))
   
    # read data from text file
    f = open('griddes.txt','r')
    lines = f.readlines()   
    for i in range(0,len(lines)):      
        if 'yvals' in lines[i]:
            yline = i
        elif 'rowlon' in lines[i] or 'reducedPoints' in lines[i]:
            rline = i
   
    # read from yvals until we hit rowlon   
    for i in range(yline,len(lines)):
        line = lines[i]
        print(i)
        if 'rowlon' in line or 'reducedPoints' in line:
            break
        if i == yline:
            # convert data to floats
            tmp_lat = [float(lat) for lat in line.split()[2:]]
        else:
            tmp_lat = [float(lat) for lat in line.split()]
        # append data to latitudes list
        for lat in tmp_lat: latitudes.append(lat) 
      
    for i in range(rline,len(lines)):
        line = lines[i]
        print(line)
        if 'scanningMode' in line: 
            break 
        if i == rline:
            # convert to integers
            tmp_nlon = [int(nlon) for nlon in line.split()[2:]]
        else:
            tmp_nlon = [int(nlon) for nlon in line.split()]
        # append data to nlongitudes list
        for nlon in tmp_nlon: nlongitudes.append(nlon) 
         
    f.close()
   
    print('nlon: ',nlongitudes)
    print('lat: ',latitudes)
   
    # Now construct the grid
    lons = []
    lats = []
    for ilat in range(0,len(nlongitudes)):
      
        lat  = latitudes[ilat]
        nlon = nlongitudes[ilat]
      
        lon1  = np.arange(0,360,360./nlon)
      
        for lon in lon1: 
            lons.append(lon)
            lats.append(lat)   
   
    if truncation_type == 'cubic-octahedral':
        ngrid = 'o%d' % (NN,)
        rfile = 'output/gaussian_grids_octahedral_reduced/%s_reduced.txt' % (ngrid,)
      
    elif truncation_type == 'linear':
        ngrid = 'n%d' % (NN,)
        rfile = 'output/gaussian_grids_linear_reduced/%s_reduced.txt' % (ngrid,)
   
    # Write to text file that CDO can use for interpolations
    f = open(rfile,'w')
    f.write('latitude reduced regular latitude \n')
    f.write('number points points \n')
    f.write(' ------- ------- ------- ---------- \n' )
   
    for ilat in range(0,len(nlongitudes)):
        f.write('%d %d %d %f \n' % (ilat+1, nlongitudes[ilat], len(nlongitudes)*2, latitudes[ilat]))
    f.close()
    
    fin = open(rfile, 'r')
       
    lines = fin.readlines()
    return (lines, NN)