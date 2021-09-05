def write_lsm(gribfield_mod, input_path_oifs, output_path_oifs, exp_name_oifs,
              grid_name_oce, num_fields, gid):
    '''
    AUTHORS:
    Jan Streffing		2020-07-22	Split off from main tool
    
    DESCRIPTION:
    This function copies the input gribfile to the output folder and modifies
    it by writing the whole gribfield_mod, including the altered land sea mask
    and soil type fields into the new file
    
    
    INPUT:

    gribfield_mod       OpenIFS input fields with modified lsm and soil classes 
    input_path_oifs     Path in unmodified ICM* files
    output_path_oifs    Path where modified ICM* files will be stored
    exp_name_oifs       4 digit code is part of the file name
    grid_name_oce       Ocean grid or mesh name
    num_fields          Number of fields in GRIB file
    gid                 List of filled Ids
    '''
    
    from shutil import copy2
    #import eccdes
    import gribapi

    input_file_oifs = input_path_oifs + 'ICMGG' + exp_name_oifs + 'INIT'
    output_file_oifs = output_path_oifs + 'ICMGG' + exp_name_oifs + 'INIT_' + grid_name_oce
    copy2(input_file_oifs, output_file_oifs)

    with open(output_file_oifs, 'r+b') as f:
        for i in range(num_fields):
            gribapi.grib_set_values(gid[i], gribfield_mod[i])
            gribapi.grib_write(gid[i], f)
            gribapi.grib_release(gid[i])
