def read_lsm(res_num, input_path_oifs, output_path_oifs, exp_name_oifs, num_fields):
    '''
    AUTHORS:
    Jan Streffing		2020-07-22	Split off from main tool


    DESCRIPTION:
    This function reads the oifs input file in grib format and returns it as a
    list of numpy arrays.
    

    INPUT:
    res_num             OpenIFS truncation number
    output_path_oifs    Path to subfolder containing unmodified files for OpenIFS
    exp_name_oifs       4 digit string from ECMWF nameing conventions for OpenIFS files
    num_fields          Number of grib fields in OpenIFS inital file

    
    RETURN:
    gribfield           List of all grib field arrays loaded from file
    lsm_id              Id of field containing land sea mask
    slt_id              Id of field containing soil type
    cl_id               Id of field containing lake mask
    gid                 List of filled Ids
    '''


    import gribapi

    print(' Opening Grib inpute file: %s ' % (input_path_oifs,))
    input_file_oifs = input_path_oifs + 'ICMGG' + exp_name_oifs + 'INIT'
    gid = [None] * num_fields
    gribfield = [None] * num_fields
    with open(input_file_oifs, 'r+') as f:
        keys = ['N', 'shortName']

        for i in range(num_fields):
            gid[i] = gribapi.grib_new_from_file(f)
            if gid[i] is None:
                break

            for key in keys:
                if not gribapi.grib_is_defined(gid[i], key):
                    raise ValueError("Key '%s' was not defined" % key)
                print('%s=%s' % (key, gribapi.grib_get(gid[i], key)))

            shortName = gribapi.grib_get(gid[i], 'shortName')

            if shortName == 'lsm':
                lsm_id = i
            if shortName == 'slt':
                slt_id = i
            if shortName == 'cl':
                cl_id = i

            gribfield[i] = gribapi.grib_get_values(gid[i])

    return (gribfield, lsm_id, slt_id, cl_id, gid)

