from __future__ import division

def read_grid_file(res_num, input_path_reduced_grid, input_path_full_grid, truncation_type):
    '''
    AUTHORS:
    Jan Streffing		2020-05-27	Split off from main tool


    DESCRIPTION:
    This function reads the reduced gaussian gridfile and returns it as a raw
    field
    

    INPUT:
    res_num 			        OpenIFS truncation number
    input_path_reduced_grid	    Path to folder containing reduced Gaussian grid definitition files
    input_path_full_grid        Path to folder containing full Gaussian grid definitition files
    truncation_type             OpenIFS truncation type (linear, quadratic, cubic octahedral)

    
    RETURN:
    lines			List of unformated lines containing latitudes of grid
    NN				Name of OpenIFS IO files according to ECMWF naming convention
    '''

    if truncation_type == 'linear':
        # linear truncation (T = NN * 2 - 1)
        print(res_num)
        print(res_num/2)
        print(res_num/2 + 0.5)
        NN = res_num/2 + 0.5
        grid_txt = '%s/n%d_reduced.txt' % (input_path_reduced_grid, NN)
    elif truncation_type == 'cubic-octahedral':
        # cubic octahedral truncation (T = NN - 1)
        NN = res_num + 1
        grid_txt = '%s/o%d_reduced.txt' % (input_path_reduced_grid, NN)

    print(' Read grid from file: %s ' % (grid_txt,) )
    print(' Reading gridfiles for T%d ' % (res_num))

    fin = open(grid_txt, 'r')
    lines = fin.readlines()
    return (lines, NN)
