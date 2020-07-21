def write_red_point_file(lats_list, lons_list, output_path_oifs, truncation_type, NN):
    '''
    AUTHORS:
    Jan Streffing		2020-07-22	Split off from main tool


    DESCRIPTION:
    This function writes the red_point.txt gridfile for OpenIFS remapping 
    with cdo
    

    INPUT:
    lats_list			List containing the latitude of each gridpoint
    lons_list			List containing the longitude of each gridpoint
    output_path_oifs    Path to subfolder containing finished files for OpenIFS
    truncation_type     OpenIFS grid type, linear or cubic octahedral
    NN				    Name of OpenIFS IO files according to ECMWF naming convention
    '''

    if truncation_type == "linear":
        redpoint_txt = '%s/TL%d_red_points.txt' % (output_path_oifs,NN*2-1)
    elif truncation_type == "cubic-octahedral":
        redpoint_txt = '%s/TCO%d_red_points.txt' % (output_path_oifs,NN-1)

    with open(redpoint_txt, 'w') as f:
        f.write("gridtype  = cell\n")
        f.write("gridsize  = " + str(len(lons_list)) + "\n")
        f.write("xvals     = ")
        for item in lons_list:
            f.write("%s\n" % item)
        f.write("yvals     = ")
        for item in lats_list:
            f.write("%s\n" % item)