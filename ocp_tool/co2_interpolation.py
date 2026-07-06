#!/usr/bin/env python
# coding: utf-8

import numpy as np
import eccodes
from scipy.interpolate import LinearNDInterpolator, griddata, interp1d
from ocp_tool.interp_utils import ReusableGridInterp, _same_points
import argparse
import os
import shutil
import traceback
import time
import dask
import warnings

# Suppress some common warnings that might occur during parallel processing
warnings.filterwarnings('ignore', category=RuntimeWarning, message='invalid value encountered')

def get_icmgg_grid(icmgg_file, variable_name='cc', verbose=False):
    """Read latitude, longitude, and values from ICMGG file for a given variable, and collect level information"""
    icmgg_lats = None
    icmgg_lons = None
    icmgg_values = None
    level_data = {}  # Dictionary to store data at each level
    levels = set()   # Set to track all available levels
    
    try:
        with open(icmgg_file, 'rb') as f:
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                    
                # First, try to check by shortName
                try:
                    short_name = eccodes.codes_get(gid, 'shortName')
                    if verbose:
                        print(f"Found ICMGG variable: {short_name}")
                except:
                    short_name = None
                
                # If shortName doesn't match, try paramId or other identifiers
                match_found = False
                if short_name == variable_name:
                    match_found = True
                else:
                    # Check by paramId if shortName doesn't match
                    try:
                        param_id = eccodes.codes_get(gid, 'paramId')
                        if verbose:
                            print(f"Parameter ID: {param_id}")
                        # 'cc' cloud cover has paramId 248 in ECMWF GRIB files
                        if variable_name == 'cc' and param_id == 248:
                            match_found = True
                        # 'ci'/'ciwc' cloud ice typically has paramId 58/247 in ECMWF GRIB files
                        elif variable_name == 'ci' and param_id in [58, 247]:
                            match_found = True
                        # Common paramId values for temperature, specific humidity, etc.
                        elif variable_name == 'an' and param_id in [130, 131, 132, 133]:
                            match_found = True
                    except:
                        pass
                        
                    # Could add more checks here for other identification methods
                
                if match_found:
                    if verbose:
                        print(f"Found matching variable: {short_name or 'unknown'} (matching {variable_name})")
                    # We found the variable
                    try:
                        level = eccodes.codes_get(gid, 'level')
                        type_of_level = eccodes.codes_get(gid, 'typeOfLevel')
                    except:
                        level = None
                        type_of_level = None
                    
                    # Get latitude and longitude arrays
                    lats = eccodes.codes_get_array(gid, 'latitudes')
                    lons = eccodes.codes_get_array(gid, 'longitudes')
                    values = eccodes.codes_get_array(gid, 'values')
                    
                    # If level exists, store it in the level_data dictionary
                    if level is not None:
                        levels.add(level)
                        # Store grid information for each level
                        level_data[level] = {
                            'lats': lats,
                            'lons': lons,
                            'values': values,
                            'type_of_level': type_of_level
                        }
                        if verbose:
                            print(f"Target level: {level}")
                    
                    # Store the first occurrence of coordinates and values for backward compatibility
                    if icmgg_values is None:
                        icmgg_lats = lats
                        icmgg_lons = lons
                        icmgg_values = values
                        
                eccodes.codes_release(gid)
        
        # Check if we found any values
        if icmgg_values is None or len(icmgg_values) == 0:
            print(f"Variable '{variable_name}' not found in ICMGG file")
            return None, None, None, None
            
        print(f"Target grid number of horizontal points: {len(icmgg_lats)}")
        print(f"Target grid number of vertical layers: {len(levels)}")
        
        # Sort levels numerically
        sorted_levels = sorted(levels)
        
        return icmgg_lats, icmgg_lons, icmgg_values, level_data
        
    except Exception as e:
        print(f"Error reading ICMGG file: {str(e)}")
        return None, None, None, None


def interpolate_co2_to_icmgg(co2_grib_file, icmgg_iniua_file, output_file=None, variable_name='an', co2_var_name='co2', verbose=False, use_dask=False, n_workers=None):
    """
    Interpolate CO2 data from a GRIB file to the grid of a specific variable in the ICMGG????INIUA file.
    
    Parameters:
    -----------
    co2_grib_file : str
        Path to the CO2 GRIB file to be interpolated.
    icmgg_iniua_file : str
        Path to the ICMGG????INIUA file containing the target grid.
    output_file : str, optional
        Path to save the output GRIB file. If None, creates a file with '_modified' suffix.
    variable_name : str, default 'an'
        Name of the variable in the ICMGG file to match the grid to.
    co2_var_name : str, default 'co2'
        Name of the CO2 variable in the GRIB file.
    verbose : bool, default False
        Enable verbose output with detailed timing.
        
    Returns:
    --------
    str
        Path to the output file with interpolated CO2 data.
    """
    print(f"Starting CO2 interpolation from {co2_grib_file} to {icmgg_iniua_file} grid")
    
    # Step 1: Read the CO2 GRIB file using eccodes
    try:
        # We'll store information for each level
        co2_level_data = {}
        
        # Open the file and iterate through messages
        with open(co2_grib_file, 'rb') as f:
            message_count = 0
            co2_messages = 0
            
            while True:
                # Get the next message
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break  # No more messages
                
                message_count += 1
                
                try:
                    # Get the variable short name
                    short_name = eccodes.codes_get(gid, 'shortName')
                    
                    # Check if this is a CO2 message
                    if short_name == co2_var_name:
                        co2_messages += 1
                        
                        # Get the level
                        try:
                            level = eccodes.codes_get(gid, 'level')
                        except:
                            level = 0
                        
                        # Get the values
                        values = eccodes.codes_get_array(gid, 'values')
                        
                        # Get lat/lon coordinates
                        lats = eccodes.codes_get_array(gid, 'latitudes')
                        lons = eccodes.codes_get_array(gid, 'longitudes')
                        
                        # Store data for this level
                        co2_level_data[level] = {
                            'values': values,
                            'lats': lats,
                            'lons': lons
                        }
                        
                        if verbose:
                            print(f"Found CO2 data for level {level}")
                        
                        if verbose:
                            print(f"Processed CO2 message at level {level}, shape: {values.shape}")
                except Exception as e:
                    print(f"Error processing message {message_count}: {e}")
                
                # Release the message
                eccodes.codes_release(gid)
            
            if verbose:
                print(f"GRIB file summary: {message_count} total messages, {co2_messages} CO2 messages")
        
        # Check if we found any CO2 messages
        if not co2_level_data:
            print(f"Error: No {co2_var_name} data found in the GRIB file")
            return None
        
        # Use level 1 data by default or the first available level
        selected_level = 1 if 1 in co2_level_data else sorted(co2_level_data.keys())[0]
        if verbose:
            print(f"Using CO2 data from level: {selected_level}")
        
        # Get the data for the selected level
        co2_values = co2_level_data[selected_level]['values']
        co2_lats = co2_level_data[selected_level]['lats']
        co2_lons = co2_level_data[selected_level]['lons']
        
        # Create points for interpolation
        co2_points = np.column_stack((co2_lons, co2_lats))
            
        # Print CO2 grid information if verbose
        if verbose:
            print(f"CO2 data shape: {co2_values.shape}")
            if isinstance(co2_lats, np.ndarray):
                print(f"CO2 latitude range: {np.min(co2_lats)} to {np.max(co2_lats)}")
            if isinstance(co2_lons, np.ndarray):
                print(f"CO2 longitude range: {np.min(co2_lons)} to {np.max(co2_lons)}")
            
    except Exception as e:
        print(f"Error reading CO2 GRIB file: {e}")
        return None
    
    # Step 2: Extract target grid information from ICMGG file
    icmgg_lats, icmgg_lons, icmgg_values, level_data = get_icmgg_grid(icmgg_iniua_file, variable_name, verbose)
    
    if icmgg_lats is None or icmgg_lons is None or level_data is None:
        print(f"Error: Could not extract grid information from ICMGG file for variable {variable_name}")
        return None
        
    # Sort levels from the target grid
    icmgg_levels = sorted(level_data.keys())
    if verbose:
        print(f"Target grid has {len(icmgg_levels)} vertical levels")
    
    # Step 3: Interpolate CO2 data to the ICMGG grid - both horizontally and vertically
    try:
        if verbose:
            print("Interpolating CO2 data to ICMGG grid...")
            
        # Create a dictionary to store interpolated CO2 data for each level
        co2_interpolated_3d = {}
        
        # Sort CO2 levels
        co2_levels = sorted(co2_level_data.keys())
        print(f"Source CO2 data has {len(co2_levels)} vertical levels")
        if verbose:
            # Print out the CO2 level keys to check format
            print(f"CO2 level keys format: {co2_levels[:5]}...")
        
        # Step 3a: Horizontal interpolation for each level of CO2 data to the ICMGG grid
        
        # Initialize a dictionary to store the interpolated CO2 data for each level
        co2_horizontal_interpolated = {}
        
        # Debug dictionary key format
        debug_key_format = True
        
        # Get the shape of the target grid from the first level's data
        first_level = level_data[icmgg_levels[0]]
        target_shape = first_level['lats'].shape
        
        # Get the first level's lat/lon grid
        icmgg_lats = first_level['lats']
        icmgg_lons = first_level['lons']
        
        # Create target points for interpolation
        icmgg_points = np.column_stack((icmgg_lons.flatten(), icmgg_lats.flatten()))
        
        # Start timing for horizontal interpolation loop
        horizontal_start_time = time.time()
        
        # All CO2 levels share the SAME horizontal source grid (a 3D field on a
        # fixed grid), so the Delaunay triangulation griddata would rebuild per
        # level is identical every time. Build it once and reuse it across
        # levels -- bit-identical to the old per-level griddata (linear + nearest
        # fill), and fast enough that the Dask cluster is no longer needed (its
        # process/cluster startup+teardown dwarfed the ~s of actual interp).
        # The _same_points guard rebuilds only if a level's grid ever differs.
        _co2_interp = None
        _co2_ref_pts = None
        for co2_level in co2_levels:
            co2_values = co2_level_data[co2_level]['values']
            co2_points = np.column_stack((co2_level_data[co2_level]['lons'],
                                          co2_level_data[co2_level]['lats']))
            if _co2_interp is None or not _same_points(co2_points, _co2_ref_pts):
                _co2_interp = ReusableGridInterp(co2_points)
                _co2_ref_pts = co2_points
            co2_interpolated = _co2_interp.linear_with_nearest_fill(
                co2_values, icmgg_points).reshape(target_shape)
            co2_horizontal_interpolated[co2_level] = co2_interpolated
            if verbose:
                print(f"  Level {co2_level}: Horizontal interpolation complete")

        # Calculate and print horizontal interpolation timing
        horizontal_end_time = time.time()
        horizontal_duration = horizontal_end_time - horizontal_start_time
        if verbose:
            print(f"Horizontal interpolation completed in {horizontal_duration:.2f} seconds")
            
        # Step 3b: Vertical interpolation to match ICMGG levels
        if verbose:
            print("Performing vertical interpolation to match ICMGG levels...")
        
        # Start timing for vertical interpolation loop
        vertical_start_time = time.time()
        
        # Define vertical interpolation function for a single target level
        def interpolate_vertical_level(target_level, icmgg_levels, co2_levels, 
                                      co2_horizontal_interpolated, level_data, target_shape):
            # Get the pressure level for the target level
            target_pressure = None
            if target_level in level_data and 'pressure' in level_data[target_level]:
                target_pressure = level_data[target_level]['pressure']
            
            # Get array of source CO2 levels (in pressure units) and their corresponding interpolated data
            src_levels = np.array([float(level) for level in co2_levels])
            
            # Initialize the output array
            interpolated_data = np.zeros(target_shape)
            
            # If we have pressure information for the target level, use vertical interpolation
            if target_pressure is not None:
                # Reshape for easier indexing
                rows, cols = target_shape
                
                # Process the vertical interpolation for each grid point
                for i in range(rows):
                    for j in range(cols):
                        # Extract CO2 values at this grid point across all source levels
                        co2_values_at_point = []
                        valid_levels = []
                        
                        for level in co2_levels:
                            if level in co2_horizontal_interpolated:
                                value = co2_horizontal_interpolated[level][i, j]
                                if not np.isnan(value):
                                    co2_values_at_point.append(value)
                                    valid_levels.append(float(level))
                        
                        # If we have enough valid values, perform 1D interpolation
                        if len(valid_levels) >= 2:
                            try:
                                # Create interpolation function for this grid point
                                f = interp1d(valid_levels, co2_values_at_point, 
                                          kind='linear', bounds_error=False, fill_value='extrapolate')
                                
                                # Interpolate to the target level
                                interpolated_data[i, j] = f(target_level)
                            except Exception:
                                # Default value if interpolation fails
                                interpolated_data[i, j] = 0.0006
                        else:
                            # If not enough values for interpolation, use default value
                            interpolated_data[i, j] = 0.0006
            else:
                # If no pressure information, use nearest source level
                nearest_idx = np.abs(src_levels - float(target_level)).argmin()
                nearest_level = co2_levels[nearest_idx]
                interpolated_data = co2_horizontal_interpolated[nearest_level].copy()
            return interpolated_data
        
        # Initialize result dictionary
        co2_interpolated_3d = {}
        
        # Vectorised vertical interpolation (replaces the per-point interp1d
        # nested loop and the Dask cluster). The horizontal step's nearest-
        # neighbour fill leaves no NaNs, so every column shares the same full set
        # of source levels -> the per-point interp1d collapses to a single axis-0
        # interp1d over the stacked field, evaluated for all columns at once.
        # Bit-identical to the old per-point linear interpolation, no cluster.
        co2_3d = np.stack([co2_horizontal_interpolated[level] for level in co2_levels], axis=0)
        src_levels_all = np.array([float(level) for level in co2_levels])
        vinterp = interp1d(src_levels_all, co2_3d, axis=0, kind='linear',
                           bounds_error=False, fill_value='extrapolate')
        for target_level in icmgg_levels:
            has_pressure = (target_level in level_data
                            and 'pressure' in level_data[target_level])
            if has_pressure:
                co2_interpolated_3d[target_level] = vinterp(float(target_level))
            else:
                nearest_idx = int(np.abs(src_levels_all - float(target_level)).argmin())
                co2_interpolated_3d[target_level] = \
                    co2_horizontal_interpolated[co2_levels[nearest_idx]].copy()
            if verbose:
                print(f"  Level {target_level}: Vertical interpolation complete")

        # Calculate and print vertical interpolation timing
        vertical_end_time = time.time()
        vertical_duration = vertical_end_time - vertical_start_time
        # Always show total time for easy comparison
        print(f"3D interpolation complete. Generated data for {len(co2_interpolated_3d)} target levels in {vertical_duration:.2f} seconds")
        
        # Show more detailed timing if verbose
        if verbose:
            horizontal_duration = horizontal_end_time - horizontal_start_time
            print(f"  Horizontal interpolation took {horizontal_duration:.2f} seconds")
            print(f"  Vertical interpolation took {vertical_duration - horizontal_duration:.2f} seconds")
            print(f"  Using {'Dask parallelization with ' + str(n_workers) + ' workers' if use_dask else 'sequential processing'}")
            
        # Print overall summary
        total_duration = vertical_end_time - horizontal_start_time
        print(f"3D interpolation complete. Generated data for {len(co2_interpolated_3d)} target levels in {total_duration:.2f} seconds")
        
    except Exception as e:
        import traceback
        print(f"Error during interpolation: {e}")
        print("Detailed error information:")
        traceback.print_exc()
        return None
    
    # Step 4: Write the interpolated CO2 data to a new GRIB file
    if output_file is None:
        base, ext = os.path.splitext(icmgg_iniua_file)
        output_file = f"{base}_with_co2{ext}"
    
    print(f"Writing interpolated CO2 data to {output_file}...")
    
    def write_output_grib(output_file, icmgg_file, interpolated_co2_3d, level_data):
        """Write the interpolated CO2 data to a GRIB file"""
        try:
            # First, get a template message for each level from the ICMGG file
            template_gids = {}
            level_templates = {}
            original_messages = []
            
            # Read all messages from the original file
            with open(icmgg_file, 'rb') as f:
                while True:
                    gid = eccodes.codes_grib_new_from_file(f)
                    if gid is None:
                        break
                    
                    # Save this message to be written to the output file
                    original_messages.append(eccodes.codes_clone(gid))
                    
                    try:
                        short_name = eccodes.codes_get(gid, 'shortName')
                        if short_name == variable_name:
                            level = eccodes.codes_get(gid, 'level')
                            # We only need one template per level
                            if level not in template_gids:
                                template_gids[level] = eccodes.codes_clone(gid)
                                level_templates[level] = {
                                    'type_of_level': eccodes.codes_get(gid, 'typeOfLevel'),
                                    'grid_type': eccodes.codes_get(gid, 'gridType')
                                }
                    except:
                        pass
                    
                    eccodes.codes_release(gid)
            
            # If no templates were found, use the first message as fallback
            if not template_gids:
                with open(icmgg_file, 'rb') as f:
                    gid = eccodes.codes_grib_new_from_file(f)
                    if gid is not None:
                        print("No level-specific templates found, using first message as template")
                        template_gids[1] = eccodes.codes_clone(gid)
                        eccodes.codes_release(gid)
                    else:
                        print("Could not find any template message in the ICMGG file")
                        return False
            
            print(f"Found {len(template_gids)} template messages for different levels")
            
            # Write to the new output file
            with open(output_file, 'wb') as out_file:
                # First, write all original messages
                for i, gid in enumerate(original_messages):
                    if verbose and i % 50 == 0:
                        print(f"Writing original message {i+1}/{len(original_messages)}")
                    eccodes.codes_write(gid, out_file)
                
                if verbose:
                    print(f"Wrote {len(original_messages)} original messages to output file")
                
                # Now append the interpolated CO2 for each level
                # Find a CC template message to use as reference for metadata
                cc_template_gid = None
                for level, gid in template_gids.items():
                    # Use the first cc message as our reference template
                    cc_template_gid = gid
                    if verbose:
                        print(f"Using CC template from level {level} as reference for metadata")
                    break
                    
                # Create a base GRIB2 template for CO2 messages
                template_gid = None
                try:
                    # Create a new GRIB2 message from a suitable sample
                    template_gid = eccodes.codes_grib_new_from_samples("regular_ll_pl_grib2")
                    
                    if template_gid is None:
                        if verbose:
                            print("Failed to create template GRIB message")
                        return False
                        
                    if verbose:
                        print(f"Created base GRIB2 template message")
                        
                    # Extract key metadata from CC template to apply to all CO2 messages
                    if cc_template_gid:
                        # Key metadata values to copy from cc template
                        metadata_keys = [
                            "gridType", "dataType", "date", "packingType",
                            "centre", "generatingProcessIdentifier"
                        ]
                        
                        metadata_values = {}
                        for key in metadata_keys:
                            try:
                                value = eccodes.codes_get(cc_template_gid, key)
                                metadata_values[key] = value
                                if verbose:
                                    print(f"Extracted {key}={value} from CC template")
                            except Exception as e:
                                if verbose:
                                    print(f"Could not extract {key} from CC template: {e}")
                except Exception as e:
                    if verbose:
                        print(f"Error creating template GRIB message: {e}")
                        import traceback
                        print(f"Error details: {traceback.format_exc()}")
                    return False
                
                # Process each level
                co2_messages_written = 0
                for level in sorted(interpolated_co2_3d.keys()):
                    interpolated_co2 = interpolated_co2_3d[level]
                    
                    if verbose:
                        print(f"Level {level}: Creating CO2 message...")
                    
                    try:
                        # Clone the template for each level to preserve all metadata attributes
                        gid = eccodes.codes_clone(template_gid)
                        
                        # Apply extracted CC metadata values to ensure consistency
                        if 'metadata_values' in locals() and metadata_values:
                            for key, value in metadata_values.items():
                                try:
                                    if isinstance(value, (int, float)):
                                        eccodes.codes_set_long(gid, key, int(value))
                                    else:
                                        eccodes.codes_set(gid, key, str(value))
                                    if verbose:
                                        print(f"Applied {key}={value} to CO2 message")
                                except Exception as e:
                                    if verbose:
                                        print(f"Could not set {key}={value} on CO2 message: {e}")
                        
                        # Set essential identification parameters for CO2
                        eccodes.codes_set_long(gid, "discipline", 0)  # Meteorological products
                        eccodes.codes_set_long(gid, "parameterCategory", 20)  # Mass mixing ratio (kg/kg)
                        eccodes.codes_set_long(gid, "parameterNumber", 61)  # Carbon dioxide
                        eccodes.codes_set_long(gid, "paramId", 210061)  # The specific CO2 paramId
                        eccodes.codes_set(gid, "shortName", "co2")
                        
                        # Set level information
                        if level in level_templates:
                            type_of_level = level_templates[level]['type_of_level']
                        else:
                            type_of_level = 'hybrid'  # Default to hybrid level
                            
                        eccodes.codes_set(gid, "typeOfLevel", type_of_level)
                        eccodes.codes_set_long(gid, "level", level)
                        
                        # Handle grid settings based on the grid type
                        grid_type = eccodes.codes_get(gid, "gridType")
                        
                        if grid_type == "reduced_gg":
                            # For reduced Gaussian grid, we need a different approach
                            # Copy N (the number of parallels between a pole and the equator)
                            # from the cc template
                            if cc_template_gid:
                                try:
                                    n_value = eccodes.codes_get_long(cc_template_gid, "N")
                                    eccodes.codes_set_long(gid, "N", n_value)
                                    if verbose:
                                        print(f"Set N={n_value} for reduced Gaussian grid")
                                        
                                    # Copy pl array (number of points along each parallel) if available
                                    try:
                                        pl_array = eccodes.codes_get_array(cc_template_gid, "pl")
                                        eccodes.codes_set_array(gid, "pl", pl_array)
                                        if verbose:
                                            print(f"Copied pl array with length {len(pl_array)}")
                                    except Exception as e:
                                        if verbose:
                                            print(f"Could not copy pl array: {e}")
                                except Exception as e:
                                    if verbose:
                                        print(f"Could not get/set N value for Gaussian grid: {e}")
                        else:
                            # Regular lat/lon grid
                            eccodes.codes_set_long(gid, "Ni", len(lons))
                            eccodes.codes_set_long(gid, "Nj", len(lats))
                            
                            eccodes.codes_set_double(gid, "longitudeOfFirstGridPointInDegrees", lons[0])
                            eccodes.codes_set_double(gid, "longitudeOfLastGridPointInDegrees", lons[-1])
                            eccodes.codes_set_double(gid, "latitudeOfFirstGridPointInDegrees", lats[0])
                            eccodes.codes_set_double(gid, "latitudeOfLastGridPointInDegrees", lats[-1])
                            
                            # Calculate grid increments properly - ensure they are positive
                            lon_increment = abs((lons[-1] - lons[0]) / (len(lons) - 1)) if len(lons) > 1 else 0
                            lat_increment = abs((lats[-1] - lats[0]) / (len(lats) - 1)) if len(lats) > 1 else 0
                            
                            eccodes.codes_set_double(gid, "iDirectionIncrementInDegrees", lon_increment)
                            eccodes.codes_set_double(gid, "jDirectionIncrementInDegrees", lat_increment)
                        
                        # Set the data values - ensure proper flattening for 2D grid
                        eccodes.codes_set_values(gid, interpolated_co2.flatten())
                        
                        # Write the message to output file
                        eccodes.codes_write(gid, out_file)
                        co2_messages_written += 1
                        
                        if verbose:
                            print(f"Level {level}: Successfully wrote CO2 message to output file")
                            
                    except Exception as e:
                        if verbose:
                            print(f"Error creating/writing CO2 message for level {level}: {e}")
                            import traceback
                            print(f"Error details: {traceback.format_exc()}")
                    finally:
                        # Clean up the message
                        if 'gid' in locals() and gid:
                            eccodes.codes_release(gid)
                
                # Clean up the template
                if template_gid:
                    eccodes.codes_release(template_gid)
                
                # Release all template handles
                for level, gid in template_gids.items():
                    eccodes.codes_release(gid)
                
                # Release all original message handles
                for gid in original_messages:
                    eccodes.codes_release(gid)
                    
                print(f"Output GRIB file written to {output_file} with {co2_messages_written} levels of CO2 data")
            return True
        except Exception as e:
            print(f"Error writing output file: {str(e)}")
            if verbose:
                import traceback
                print(f"Error details: {traceback.format_exc()}")
            return False
    
    try:
        success = write_output_grib(output_file, icmgg_iniua_file, co2_interpolated_3d, level_data)
        
        if success:
            return output_file
        else:
            print("Data was not written successfully (this may also indicate that co2 data was already present in the output file).")
            return None
    except Exception as e:
        print(f"Error during interpolation: {e}")
        return None

def main():
    """
    Example usage of the interpolate_co2_to_icmgg function.
    """
    import argparse
    import dask
    import os
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Interpolate CO2 data from a GRIB file to an ICMGG grid')
    parser.add_argument('co2_grib', help='Path to the CO2 GRIB file')
    parser.add_argument('icmgg_iniua', help='Path to the ICMGG INIUA file')
    parser.add_argument('output_file', nargs='?', help='Output file path (default: input_name_WITH_CO2.suffix)')
    parser.add_argument('--var', default='cc', help='Variable name in the ICMGG file to use for grid reference (default: cc)')
    parser.add_argument('--co2_var', default='co2', help='Variable name in the CO2 GRIB file (default: co2)')
    parser.add_argument('--verbose', action='store_true', help="Enable verbose output with detailed timing")
    parser.add_argument('--dask', action='store_true', help="Use Dask for parallel processing")
    parser.add_argument('--workers', type=int, default=None, help="Number of Dask workers to use (default: number of CPU cores)")
    
    args = parser.parse_args()
    
    # Determine output file path
    output_file = args.output_file or f"{os.path.splitext(args.icmgg_iniua)[0]}_WITH_CO2{os.path.splitext(args.icmgg_iniua)[1]}"
    
    # Check if output file exists and delete it if it does
    if os.path.exists(output_file):
        if args.verbose:
            print(f"Output file {output_file} already exists, removing it...")
        os.remove(output_file)
        if args.verbose:
            print(f"Removed existing output file: {output_file}")
    
    if args.dask:
        dask.config.set(scheduler='processes')
        if args.workers is None:
            args.workers = os.cpu_count()
        dask.config.set(num_workers=args.workers)
    
    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        if args.verbose:
            print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)
    
    result = interpolate_co2_to_icmgg(args.co2_grib, args.icmgg_iniua, output_file, args.var, args.co2_var, args.verbose, args.dask, args.workers)

    if result:
        # Check if the CO2 shortName was set correctly
        if args.verbose:
            try:
                print("\nVerifying final metadata in output file:")
                check_cmd = f"grib_ls {result} | grep -E 'shortName|paramId' | head -n 10"
                os.system(check_cmd)
            except Exception:
                pass
    else:
        print("Data was not written successfully (this may also indicate that co2 data was already present in the output file).")

if __name__ == "__main__":
    main()
