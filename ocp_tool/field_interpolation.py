#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import sys
import traceback
import shutil
import warnings
from scipy.interpolate import griddata
import eccodes
import subprocess
from pathlib import Path

# Suppress some common warnings
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered')

def interpolate_2d_fields_to_icmgg(input_grib_file, icmgg_init_file, output_file=None,
                                variable_name='lsm', field_type='co2_emissions', verbose=False):
    """
    Interpolate 2D fields (CO2 emissions or albedo) from a GRIB file to the grid of a specific variable 
    in the ICMGG????INIT file.
    
    Parameters:
    -----------
    input_grib_file : str
        Path to the input GRIB file (CO2 emissions or albedo) to be interpolated.
    icmgg_init_file : str
        Path to the ICMGG????INIT file containing the target grid.
    output_file : str, optional
        Path to save the output GRIB file. If None, uses the icmgg_init_file path.
    variable_name : str, optional
        Name of variable in ICMGG file to use for extracting target grid. Default is 'lsm'.
    field_type : str, optional
        Type of field to interpolate: 'co2_emissions' or 'albedo'. This determines the interpolation method.
        For 'co2_emissions', use python-based interpolation.
        For 'albedo', use CDO-based remapping.
    verbose : bool, optional
        Enable verbose output for debugging. Default is False.
        
    Returns:
    --------
    str or None
        Path to the output file if successful, None otherwise.
    """
    try:
        if not os.path.exists(input_grib_file):
            print(f"Error: Input GRIB file not found: {input_grib_file}")
            return None
            
        if not os.path.exists(icmgg_init_file):
            print(f"Error: ICMGG file not found: {icmgg_init_file}")
            return None
            
        # If no output file is specified, use the input ICMGG file
        if output_file is None:
            output_file = icmgg_init_file
            
        print(f"Starting {field_type} interpolation from {input_grib_file} to {icmgg_init_file} grid")
        
        # Choose the interpolation method based on field_type
        if field_type == 'albedo':
            return interpolate_albedo_fields(input_grib_file, icmgg_init_file, output_file, verbose)
        else:  # Default to CO2 emissions interpolation method
            return interpolate_co2_emissions(input_grib_file, icmgg_init_file, output_file, variable_name, verbose)
            
    except Exception as e:
        print(f"Error during {field_type} interpolation: {e}")
        if verbose:
            print(traceback.format_exc())
        return None
    
def interpolate_co2_emissions(input_grib_file, icmgg_init_file, output_file, variable_name='lsm', verbose=False):
    """Interpolate CO2 emissions data using Python-based interpolation"""
    
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
        with open(input_grib_file, 'rb') as f:
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
        
        with open(input_grib_file, 'rb') as f:
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
        # Find a message in the input file to use as a template for output
        template_gid = None
        with open(icmgg_init_file, 'rb') as f:
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                
                try:
                    # Any message will do for a template
                    template_gid = eccodes.codes_clone(gid)
                    if verbose:
                        print("Found template message in ICMGG file")
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

def interpolate_albedo_fields(input_grib_file, icmgg_init_file, output_file, verbose=False):
    """Interpolate albedo fields using Python-based interpolation (same as CO2)"""
    # Use lsm as the reference variable to extract grid information
    variable_name = 'lsm'
    
    # Step 1: Read target grid information from ICMGG file
    target_lats, target_lons = read_icmgg_grid(icmgg_init_file, variable_name, verbose)
    
    if target_lats is None or target_lons is None:
        print(f"Error: Could not extract grid information from ICMGG file for variable {variable_name}")
        return None
        
    # Prepare target points for interpolation
    target_points = np.column_stack((target_lons.flatten(), target_lats.flatten()))
    target_shape = target_lats.shape
    
    try:
        print(f"Interpolating albedo fields to match {output_file} grid using Python")
        
        # First pass: identify all albedo fields by their paramId
        albedo_params = []
        with open(input_grib_file, 'rb') as f:
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                
                try:
                    param_id = eccodes.codes_get(gid, 'paramId')
                    short_name = eccodes.codes_get(gid, 'shortName')
                    if param_id not in albedo_params:
                        albedo_params.append(param_id)
                        if verbose:
                            print(f"Found albedo field: {short_name} (paramId: {param_id})")
                finally:
                    eccodes.codes_release(gid)
        
        if not albedo_params:
            print("No albedo fields found in GRIB file")
            return None
            
        print(f"Found {len(albedo_params)} albedo fields with paramIds: {', '.join(map(str, albedo_params))}")
        
        # Second pass: read, interpolate and prepare each field for writing
        # Use paramId as the key to ensure unique identification even with duplicate shortNames
        interpolated_fields = {}
        
        with open(input_grib_file, 'rb') as f:
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                
                try:
                    short_name = eccodes.codes_get(gid, 'shortName')
                    param_id = eccodes.codes_get(gid, 'paramId')
                    
                    # Get the values and coordinates
                    values = eccodes.codes_get_array(gid, 'values')
                    lats = eccodes.codes_get_array(gid, 'latitudes')
                    lons = eccodes.codes_get_array(gid, 'longitudes')
                    
                    # Create source points for interpolation
                    source_points = np.column_stack((lons, lats))
                    
                    # Store the original GRIB message information for template
                    template_info = {
                        'short_name': short_name,
                        'param_id': param_id,
                        'level_type': eccodes.codes_get(gid, 'typeOfLevel'),
                        'grib_template': eccodes.codes_clone(gid)
                    }
                    
                    if verbose:
                        print(f"Interpolating field with paramId {param_id} (shortName: {short_name})")
                    
                    # Use nearest neighbor interpolation for albedo fields
                    interpolated_values = griddata(
                        source_points, 
                        values, 
                        target_points, 
                        method='nearest'
                    )
                    
                    # Reshape to target grid dimensions
                    interpolated_values = interpolated_values.reshape(target_shape)
                    
                    # Use paramId as the key instead of shortName to handle duplicate shortNames
                    interpolated_fields[param_id] = {
                        'values': interpolated_values,
                        'template': template_info
                    }
                        
                except Exception as e:
                    if verbose:
                        print(f"Error interpolating variable {short_name}: {e}")
                
                finally:
                    eccodes.codes_release(gid)
        
        # Make a backup of the original file
        backup_file = f"{icmgg_init_file}.bak"
        if not os.path.exists(backup_file):
            print(f"Creating backup of original ICMGG file: {backup_file}")
            shutil.copy2(icmgg_init_file, backup_file)
        
        # Create a temporary file to hold the rewritten contents
        temp_output = f"{os.path.dirname(output_file)}/temp_albedo_icmgg.grb"
        
        # Read the original file and filter out existing albedo fields (paramIds 117-120)
        print("Filtering out existing albedo fields from the original file")
        albedo_param_ids = [117, 118, 119, 120]  # The paramIds for albedo fields
        non_albedo_messages = []
        
        with open(icmgg_init_file, 'rb') as f:
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                
                try:
                    param_id = eccodes.codes_get(gid, 'paramId')
                    
                    if param_id not in albedo_param_ids:
                        # This is not an albedo field - keep it
                        non_albedo_messages.append(eccodes.codes_clone(gid))
                        if verbose:
                            print(f"Keeping non-albedo field with paramId {param_id}")
                    else:
                        if verbose:
                            print(f"Filtering out existing albedo field with paramId {param_id}")
                finally:
                    eccodes.codes_release(gid)
        
        print(f"Kept {len(non_albedo_messages)} non-albedo fields from original file")
        
        # Find a message in the input file to use as a template for output (already done via non_albedo_messages)
        if not non_albedo_messages:
            print("Error: Could not find any message in the ICMGG file to use as template")
            return None
        
        # Write all non-albedo messages to the temporary file first
        print(f"Writing non-albedo fields to temporary file {temp_output}")
        with open(temp_output, 'wb') as out_file:
            for gid in non_albedo_messages:
                eccodes.codes_write(gid, out_file)
                eccodes.codes_release(gid)
        
        # Now append the interpolated albedo fields
        print(f"Adding {len(interpolated_fields)} interpolated albedo fields to output file")
        
        with open(temp_output, 'ab') as out_file:
            for param_id, field_info in interpolated_fields.items():
                template = field_info['template']['grib_template']
                values = field_info['values']
                short_name = field_info['template']['short_name']
                
                try:
                    # Set the values and other metadata
                    eccodes.codes_set_array(template, 'values', values.flatten())
                    eccodes.codes_set(template, 'dataDate', 19900101)  # Set date to January 1, 1990
                    
                    # Write the message to the output file
                    eccodes.codes_write(template, out_file)
                    
                    if verbose:
                        print(f"Written {short_name} to output file")
                        
                except Exception as e:
                    print(f"Error writing {short_name}: {e}")
                
                finally:
                    eccodes.codes_release(template)
        
        # Move temporary file to final location
        if os.path.exists(temp_output):
            print(f"Moving temporary file to final location {output_file}")
            os.replace(temp_output, output_file)
            print(f"Successfully replaced ICMGG file with version containing new albedo fields")
            return output_file
        else:
            print("Error: Temporary file not created")
            return None
            
    except Exception as e:
        print(f"Error adding albedo fields: {e}")
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
    
    parser = argparse.ArgumentParser(description='Interpolate 2D fields (CO2 emissions or albedo) to ICMGG grid.')
    parser.add_argument('--input-grib', required=True, help='Path to input GRIB file (CO2 emissions or albedo)')
    parser.add_argument('--icmgg', required=True, help='Path to ICMGG file containing target grid')
    parser.add_argument('--output', help='Path to output file (default is to modify the ICMGG file)')
    parser.add_argument('--variable', default='lsm', help='Variable name in ICMGG to match grid (default: lsm)')
    parser.add_argument('--field-type', default='co2_emissions', choices=['co2_emissions', 'albedo'], 
                       help='Type of field to interpolate: co2_emissions or albedo')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    interpolate_2d_fields_to_icmgg(
        args.input_grib,
        args.icmgg,
        output_file=args.output,
        variable_name=args.variable,
        field_type=args.field_type,
        verbose=args.verbose
    )
