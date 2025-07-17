#!/usr/bin/env python
# coding: utf-8

import numpy as np
import eccodes
from scipy.interpolate import LinearNDInterpolator, griddata
import os
import shutil
import traceback
import time
import warnings

# Suppress some common warnings
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered')

def interpolate_co2_emissions_to_icmgg(co2_emissions_grib_file, icmgg_init_file, output_file=None, 
                                       variable_name='an', verbose=False):
    """
    Interpolate 2D CO2 emissions data from a GRIB file to the grid of a specific variable 
    in the ICMGG????INIT file.
    
    Parameters:
    -----------
    co2_emissions_grib_file : str
        Path to the CO2 emissions GRIB file to be interpolated.
    icmgg_init_file : str
        Path to the ICMGG????INIT file containing the target grid.
    output_file : str, optional
        Path to save the output GRIB file. If None, uses the icmgg_init_file path.
    variable_name : str, default 'an'
        Name of the variable in the ICMGG file to match the grid to.
    verbose : bool, default False
        Enable verbose output with detailed timing.
        
    Returns:
    --------
    str
        Path to the output file with interpolated CO2 emissions data.
    """
    print(f"Starting CO2 emissions interpolation from {co2_emissions_grib_file} to {icmgg_init_file} grid")
    
    # If no output file is specified, use the input ICMGG file
    if output_file is None:
        output_file = icmgg_init_file
        
    # Step 1: Read target grid information from ICMGG file
    target_lats, target_lons = read_icmgg_grid(icmgg_init_file, variable_name, verbose)
    
    if target_lats is None or target_lons is None:
        print(f"Error: Could not extract grid information from ICMGG file for variable {variable_name}")
        return None
    
    # Prepare target points for interpolation
    target_points = np.column_stack((target_lons.flatten(), target_lats.flatten()))
    target_shape = target_lats.shape
    
    # Step 2: Read and interpolate each CO2 emission variable
    try:
        # First pass: identify all emissions variables
        emissions_vars = []
        with open(co2_emissions_grib_file, 'rb') as f:
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                
                try:
                    short_name = eccodes.codes_get(gid, 'shortName')
                    if short_name not in emissions_vars:
                        emissions_vars.append(short_name)
                        if verbose:
                            print(f"Found CO2 emissions variable: {short_name}")
                finally:
                    eccodes.codes_release(gid)
        
        if not emissions_vars:
            print("No CO2 emissions variables found in GRIB file")
            return None
            
        print(f"Found {len(emissions_vars)} CO2 emissions variables: {', '.join(emissions_vars)}")
        
        # Second pass: read, interpolate and prepare each variable for writing
        interpolated_emissions = {}
        
        with open(co2_emissions_grib_file, 'rb') as f:
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                
                try:
                    short_name = eccodes.codes_get(gid, 'shortName')
                    
                    # Get the values and coordinates
                    values = eccodes.codes_get_array(gid, 'values')
                    lats = eccodes.codes_get_array(gid, 'latitudes')
                    lons = eccodes.codes_get_array(gid, 'longitudes')
                    
                    # Create source points for interpolation
                    source_points = np.column_stack((lons, lats))
                    
                    # Store the original GRIB message information for template
                    template_info = {
                        'short_name': short_name,
                        'param_id': eccodes.codes_get(gid, 'paramId'),
                        'level_type': eccodes.codes_get(gid, 'typeOfLevel'),
                        'level': eccodes.codes_get(gid, 'level'),
                        'grid_type': eccodes.codes_get(gid, 'gridType')
                    }
                    
                    # Linear interpolation with scipy's griddata
                    print(f"Interpolating {short_name} to target grid...")
                    interpolated_values = griddata(
                        source_points, 
                        values, 
                        target_points, 
                        method='linear', 
                        fill_value=0
                    )
                    
                    # Reshape to match target grid
                    interpolated_values = interpolated_values.reshape(target_shape)
                    
                    # Store interpolated values and template info
                    interpolated_emissions[short_name] = {
                        'values': interpolated_values,
                        'template': template_info
                    }
                    
                    print(f"Successfully interpolated {short_name}")
                    
                except Exception as e:
                    print(f"Error processing CO2 emissions variable {short_name}: {e}")
                    if verbose:
                        print(traceback.format_exc())
                finally:
                    eccodes.codes_release(gid)
        
        # Step 3: Create a copy of the original ICMGG file as backup
        backup_file = icmgg_init_file + '.bak'
        if not os.path.exists(backup_file):
            print(f"Creating backup of original ICMGG file: {backup_file}")
            shutil.copy2(icmgg_init_file, backup_file)
        
        # Step 4: Write interpolated emissions to the output file
        # First, we need to get a template message from the target file
        template_gid = None
        with open(icmgg_init_file, 'rb') as f:
            # Find a suitable template message
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                
                try:
                    # Check if this is a similar type of field (surface field)
                    level_type = eccodes.codes_get(gid, 'typeOfLevel')
                    if level_type == 'surface':
                        template_gid = eccodes.codes_clone(gid)
                        break
                finally:
                    eccodes.codes_release(gid)
        
        if template_gid is None:
            print("Could not find suitable template message in ICMGG file")
            return None
        
        # Now append the new CO2 emissions messages to the output file
        emissions_written = 0
        try:
            # Open output file for appending
            with open(output_file, 'ab') as out_file:
                for short_name, data in interpolated_emissions.items():
                    interpolated_values = data['values']
                    template_info = data['template']
                    
                    try:
                        # Create a new message based on the template
                        gid = eccodes.codes_clone(template_gid)
                        
                        # Set/override key parameters
                        eccodes.codes_set(gid, 'shortName', short_name)
                        eccodes.codes_set_long(gid, 'paramId', template_info['param_id'])
                        eccodes.codes_set(gid, 'typeOfLevel', template_info['level_type'])
                        eccodes.codes_set_long(gid, 'level', template_info['level'])
                        
                        # Try to set date to January 1, 1990
                        try:
                            eccodes.codes_set_long(gid, 'dataDate', 19900101)
                        except Exception as e:
                            if verbose:
                                print(f"Warning: Could not set dataDate: {e}")
                                
                        # Instead of trying to set individual date components which might be read-only,
                        # we'll verify the date in the final output
                        
                        # Ensure proper grid definition
                        if eccodes.codes_get(gid, 'gridType') == 'reduced_gg':
                            # For reduced Gaussian grid, we need to set number of points along parallels
                            pass  # This would require additional handling if needed
                        else:
                            # For regular lat/lon grid, set grid parameters
                            if hasattr(target_lons, 'ndim') and target_lons.ndim > 1:
                                # 2D grid
                                eccodes.codes_set_long(gid, 'Ni', target_lons.shape[1])
                                eccodes.codes_set_long(gid, 'Nj', target_lons.shape[0])
                            else:
                                # 1D grid
                                eccodes.codes_set_long(gid, 'Ni', len(target_lons))
                                eccodes.codes_set_long(gid, 'Nj', len(target_lats))
                        
                        # Set the interpolated values
                        eccodes.codes_set_values(gid, interpolated_values.flatten())
                        
                        # Write the message to the output file
                        eccodes.codes_write(gid, out_file)
                        emissions_written += 1
                        
                        print(f"Successfully wrote {short_name} to output file")
                        
                    except Exception as e:
                        print(f"Error writing CO2 emissions variable {short_name}: {e}")
                        if verbose:
                            print(traceback.format_exc())
                    finally:
                        if 'gid' in locals() and gid:
                            eccodes.codes_release(gid)
            
            print(f"Added {emissions_written} CO2 emissions variables to {output_file}")
            
        except Exception as e:
            print(f"Error writing output file: {e}")
            if verbose:
                print(traceback.format_exc())
        finally:
            # Release template
            if template_gid:
                eccodes.codes_release(template_gid)
        
        return output_file if emissions_written > 0 else None
        
    except Exception as e:
        print(f"Error during CO2 emissions interpolation: {e}")
        if verbose:
            print(traceback.format_exc())
        return None

def read_icmgg_grid(icmgg_file, variable_name='an', verbose=False):
    """Read latitude and longitude from ICMGG file for a given variable."""
    try:
        icmgg_lats = None
        icmgg_lons = None
        
        with open(icmgg_file, 'rb') as f:
            message_count = 0
            
            while True:
                # Get the next message
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break  # No more messages
                
                message_count += 1
                
                try:
                    # Check if this is the variable we're looking for
                    short_name = eccodes.codes_get(gid, 'shortName')
                    
                    if short_name == variable_name:
                        # Get lat/lon coordinates for this message
                        icmgg_lats = eccodes.codes_get_array(gid, 'latitudes')
                        icmgg_lons = eccodes.codes_get_array(gid, 'longitudes')
                        
                        # Get the shape of the grid
                        try:
                            # Try to get Ni/Nj for regular grids
                            ni = eccodes.codes_get_long(gid, 'Ni')
                            nj = eccodes.codes_get_long(gid, 'Nj')
                            
                            # Reshape lat/lon arrays if they're flattened
                            if icmgg_lats.size == ni * nj:
                                icmgg_lats = icmgg_lats.reshape(nj, ni)
                                icmgg_lons = icmgg_lons.reshape(nj, ni)
                        except:
                            # For reduced Gaussian grids, this might not work
                            pass
                        
                        if verbose:
                            print(f"Found target variable '{variable_name}', grid shape: {icmgg_lats.shape}")
                        
                        break  # We found what we need
                        
                except Exception as e:
                    if verbose:
                        print(f"Error processing message {message_count}: {e}")
                
                finally:
                    # Release the message
                    eccodes.codes_release(gid)
        
        # Check if we found the variable
        if icmgg_lats is None or icmgg_lons is None:
            print(f"Variable '{variable_name}' not found in ICMGG file")
            return None, None
            
        return icmgg_lats, icmgg_lons
        
    except Exception as e:
        print(f"Error reading ICMGG file: {str(e)}")
        return None, None

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Interpolate CO2 emissions to ICMGG grid.')
    parser.add_argument('--co2-emissions', required=True, help='Path to CO2 emissions GRIB file')
    parser.add_argument('--icmgg', required=True, help='Path to ICMGG file containing target grid')
    parser.add_argument('--output', help='Path to output file (default is to modify the ICMGG file)')
    parser.add_argument('--variable', default='lsm', help='Variable name in ICMGG to match grid (default: an)')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    interpolate_co2_emissions_to_icmgg(
        args.co2_emissions,
        args.icmgg,
        output_file=args.output,
        variable_name=args.variable,
        verbose=args.verbose
    )
