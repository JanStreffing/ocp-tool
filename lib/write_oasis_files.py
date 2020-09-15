def write_oasis_files(res_num, output_path_oasis, dir_path, grid_name_oce, center_lats, center_lons, crn_lats, crn_lons, gridcell_area, lsm_binary_a ,lsm_binary_l , NN, input_path_runoff):
    '''
    AUTHORS:
    Jan Streffing		2020-07-22	Split off from main tool
    
    DESCRIPTION:
    This function writes the binary masks, areas and grids files for
    oasis3-mct. For OpenIFS we set up three different grids:
    All have the same lons, lats and areas
    but different land-sea masks

    atmo: Land = 0, Ocean = 1  (for atm->ocn fluxes)
    atmr: Land = 1, Ocean = 0  (for atm->runoff mapper)

    atmf is the fraction of ocean in each cell, but
    IFS always has 1 or 0, so it is the same as atmo.
    
    
    INPUT:
        
    res_num             OpenIFS truncation number
    output_path_oasis   Finished files are written here
    dir_path            Root dir of ocp-tool
    grid_name_oce       Ocean grid or mesh name
    center_lats			List of center latitudes
    center_lons			List of center longitudes
    crn_lats			List of corner latitudes
    crn_lons			List of corner longitudes
    gridcell_area		Vector containing the area of each grid cell
    lsm_binary_a        Altered land sea mask in shape for oasis3-mct
    lsm_binary_l        Original land sea mask in shape for oasis3-mct
    NN                  Truncation number for oasis grid name
    input_path_runoff   Path to copying runoff mapper fields from

    OUTPUT:
        
    
    '''

    for filebase in ['grids', 'areas', 'masks']:
        filename = '%s/%s%s.nc' % (dir_path, output_path_oasis, filebase)
        print('Writing file: %s ' % (filename,))
        nc = Dataset(filename, 'w', clobber=True)

        # For OpenIFS + NEMO + Runoffmapper we need two atmosphere grids:
        # atmo: used for atm->ocn remapping (to find ocean)
        # atmr: used for atm->runoff remapping (to find land)

        for grids_name in ('{}{:03}'.format(s, int(NN)) for s in ('A', 'L', 'R')):

            # OASIS requires certain names for the dimensions etc
            print(' Write lons, lats, corner points for grid: %s ' % (grids_name,), '(T%s)' % (res_num,))
            xname = 'x_%s' % (grids_name,)
            yname = 'y_%s' % (grids_name,)
            lonname = '%s.lon' % (grids_name,)
            latname = '%s.lat' % (grids_name,)
            nc.createDimension(xname, center_lons.shape[1])
            nc.createDimension(yname, 1)
            id_lon = nc.createVariable(lonname, 'float64', (yname, xname))
            id_lat = nc.createVariable(latname, 'float64', (yname, xname))

            # Not required, but makes grid info more readable
            id_lon.units = 'degrees_east'
            id_lon.standard_name = 'Longitude'
            id_lat.units = 'degrees_north'
            id_lat.standard_name = 'Latitude'

         # Write corner points to grids file
            if filebase == 'grids':
                crnname = 'crn_%s' % (grids_name,)
                cloname = '%s.clo' % (grids_name,)
                claname = '%s.cla' % (grids_name,)
                nc.createDimension(crnname, 4)
                id_clo = nc.createVariable(cloname, 'float64', (crnname, yname, xname))
                id_cla = nc.createVariable(claname, 'float64', (crnname, yname, xname))

         # Write land-sea masks to masks file
            elif filebase == 'masks':
                mskname = '%s.msk' % (grids_name,)
                id_msk = nc.createVariable(mskname, 'int32', (yname, xname))
                id_msk.valid_min = 0.
                id_msk.valid_max = 1

         # Write grid cell area to areas file
            elif filebase == 'areas':
                areaname = '%s.srf' % (grids_name,)
                id_area = nc.createVariable(areaname, 'float64', (yname, xname))

            id_lon[:, :] = center_lons[:, :]
            id_lat[:, :] = center_lats[:, :]
            id_lon.valid_min = center_lons.min()
            id_lon.valid_max = center_lons.max()
            id_lat.valid_min = center_lats.min()
            id_lat.valid_max = center_lats.max()

            if filebase == 'grids':
                id_clo[:, :, :] = crn_lons[:, :, :]
                id_cla[:, :, :] = crn_lats[:, :, :]
                id_clo.valid_min = crn_lons.min()
                id_clo.valid_max = crn_lons.max()
                id_cla.valid_min = crn_lats.min()
                id_cla.valid_max = crn_lats.max()

            elif filebase == 'masks':
                if grids_name.startswith('A') :
                    id_msk[:, :] = np.round(lsm_binary_a[:, :])  
                elif grids_name.startswith('L'):
                    id_msk[:, :] = np.round(lsm_binary_l[:, :])
                elif grids_name.startswith('R'):
                    id_msk[:, :] = np.abs(np.round(lsm_binary_a[:, :] - 1))
                else:
                    raise RuntimeError('Unexpected grid name: {}'.format(grids_name))

            elif filebase == 'areas':
                id_area[:, :] = gridcell_area[:, :]
                id_area.valid_min = gridcell_area.min()
                id_area.valid_max = gridcell_area.max()


        # Copying runoff mapper grids and areas into oasis3-mct files

        input_file_rnf = '%srunoff_%s.nc' % (input_path_runoff, filebase)
        rnffile = Dataset(input_file_rnf, 'r')

        nc.setncatts(rnffile.__dict__)
        for name, dimension in rnffile.dimensions.iteritems():
            nc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)

        for name, variable in rnffile.variables.iteritems():
            var_out = nc.createVariable(name, variable.datatype, variable.dimensions)
            var_out.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
            var_out[:] = variable[:]

        rnffile.close()
        nc.close()
        print(' Wrote %s ' % (filename,))