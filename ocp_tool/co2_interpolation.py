#!/usr/bin/env python
# coding: utf-8

import numpy as np
import eccodes
from scipy.interpolate import LinearNDInterpolator, griddata, interp1d
import argparse
import os
import shutil
import traceback
import time
import dask
from dask.distributed import Client, LocalCluster
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
            
        print(f"Target grid points: {len(icmgg_lats)}")
        print(f"Found {len(levels)} vertical levels for {variable_name}")
        
        # Sort levels numerically
        sorted_levels = sorted(levels)
        print(f"Vertical levels range: {sorted_levels[0]} to {sorted_levels[-1]}")
        
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
        print("Reading CO2 GRIB file using eccodes...")
        
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
                        
                        print(f"Processed CO2 message at level {level}, shape: {values.shape}")
                except Exception as e:
                    print(f"Error processing message {message_count}: {e}")
                
                # Release the message
                eccodes.codes_release(gid)
            
            print(f"GRIB file summary: {message_count} total messages, {co2_messages} CO2 messages")
        
        # Check if we found any CO2 messages
        if not co2_level_data:
            print(f"Error: No {co2_var_name} data found in the GRIB file")
            return None
        
        # Use level 1 data by default or the first available level
        selected_level = 1 if 1 in co2_level_data else sorted(co2_level_data.keys())[0]
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
    print(f"Target grid has {len(icmgg_levels)} vertical levels")
    
    # Step 3: Interpolate CO2 data to the ICMGG grid - both horizontally and vertically
    try:
        print("Interpolating CO2 data to ICMGG grid...")
        
        # Create a dictionary to store interpolated CO2 data for each level
        co2_interpolated_3d = {}
        
        # Sort CO2 levels
        co2_levels = sorted(co2_level_data.keys())
        print(f"Source CO2 data has {len(co2_levels)} vertical levels")
        
        # Debug: Print out the CO2 level keys to check format
        debug_key_format = True
        if debug_key_format:
            print(f"CO2 level keys format: {co2_levels[:5]}...")
        
        # Step 3a: Horizontal interpolation for each level of CO2 data to the ICMGG grid
        
        # Initialize a dictionary to store the interpolated CO2 data for each level
        co2_horizontal_interpolated = {}
        
        # Debug dictionary key format
        debug_key_format = True
        
        # Get the shape of the target grid from the first level's data
        first_level = level_data[icmgg_levels[0]]
        target_shape = first_level['lats'].shape
        
        # Get the first level's lat/lon grid (assuming all levels have same horizontal grid)
        icmgg_lats = first_level['lats']
        icmgg_lons = first_level['lons']
        
        # Normalize longitude values to 0-360 range if needed
        icmgg_lons_normalized = icmgg_lons.copy()
        if np.min(icmgg_lons_normalized) < 0:
            icmgg_lons_normalized[icmgg_lons_normalized < 0] += 360
            
        # Create target points for interpolation
        icmgg_points = np.column_stack((icmgg_lons_normalized.flatten(), icmgg_lats.flatten()))
        
        # Perform horizontal interpolation for each CO2 level
        print("Performing horizontal interpolation for each CO2 level...")
        
        # Start timing for horizontal interpolation loop
        horizontal_start_time = time.time()
        
        # Define horizontal interpolation function for a single level
        def interpolate_horizontal_level(co2_level, co2_level_data, icmgg_points, target_shape):
            # Get data for this level
            co2_values = co2_level_data[co2_level]['values']
            co2_lats = co2_level_data[co2_level]['lats']
            co2_lons = co2_level_data[co2_level]['lons']
            
            # Create points for interpolation
            co2_points = np.column_stack((co2_lons, co2_lats))
            
            # Normalize CO2 longitudes
            co2_lons_normalized = co2_points[:, 0].copy()
            if np.min(co2_lons_normalized) < 0:
                co2_lons_normalized[co2_lons_normalized < 0] += 360
            co2_points_normalized = np.column_stack((co2_lons_normalized, co2_points[:, 1]))
            
            # Linear interpolation with scipy's griddata
            co2_interpolated = griddata(
                co2_points_normalized, 
                co2_values, 
                icmgg_points, 
                method='linear'
            )
            
            # Fill NaN values with nearest neighbor interpolation
            if np.isnan(co2_interpolated).any():
                nan_mask = np.isnan(co2_interpolated)
                co2_interpolated[nan_mask] = griddata(
                    co2_points_normalized, 
                    co2_values, 
                    icmgg_points[nan_mask], 
                    method='nearest'
                )
            
            # Reshape back to original grid shape
            co2_interpolated = co2_interpolated.reshape(target_shape)
            
            return co2_interpolated
        
        if use_dask:
            # Set up Dask client for parallel processing
            if verbose:
                print("Setting up Dask client for horizontal interpolation...")
            
            # Use a context manager for the Dask client to ensure proper cleanup
            n_workers = n_workers or os.cpu_count()
            with LocalCluster(n_workers=n_workers, threads_per_worker=1) as cluster, Client(cluster) as client:
                if verbose:
                    print(f"Dask cluster ready with {n_workers} workers")
                
                # Create delayed computations for each level
                delayed_results = []
                for co2_level in co2_levels:
                    delayed_result = dask.delayed(interpolate_horizontal_level)(
                        co2_level, co2_level_data, icmgg_points, target_shape
                    )
                    delayed_results.append((co2_level, delayed_result))
                
                # Compute all levels in parallel
                if verbose:
                    print(f"Computing {len(delayed_results)} horizontal interpolation tasks in parallel...")
                results = dask.compute(*[res for _, res in delayed_results])
                
                # Store results back in the dictionary
                for i, co2_level in enumerate([level for level, _ in delayed_results]):
                    co2_horizontal_interpolated[co2_level] = results[i]
                    if verbose:
                        print(f"  Level {co2_level}: Horizontal interpolation complete")
        else:
            # Sequential processing for each level
            for co2_level in co2_levels:
                # Get data for this level
                co2_values = co2_level_data[co2_level]['values']
                co2_lats = co2_level_data[co2_level]['lats']
                co2_lons = co2_level_data[co2_level]['lons']
                
                # Create points for interpolation
                co2_points = np.column_stack((co2_lons, co2_lats))
                
                # Normalize CO2 longitudes
                co2_lons_normalized = co2_points[:, 0].copy()
                if np.min(co2_lons_normalized) < 0:
                    co2_lons_normalized[co2_lons_normalized < 0] += 360
                co2_points_normalized = np.column_stack((co2_lons_normalized, co2_points[:, 1]))
                
                # Linear interpolation with scipy's griddata
                co2_interpolated = griddata(
                    co2_points_normalized, 
                    co2_values, 
                    icmgg_points, 
                    method='linear'
                )
                
                # Fill NaN values with nearest neighbor interpolation
                if np.isnan(co2_interpolated).any():
                    print(f"  Level {co2_level}: Filling NaN values with nearest neighbor interpolation...")
                    nan_mask = np.isnan(co2_interpolated)
                    co2_interpolated[nan_mask] = griddata(
                        co2_points_normalized, 
                        co2_values, 
                        icmgg_points[nan_mask], 
                        method='nearest'
                    )
                
                # Reshape back to original grid shape
                co2_interpolated = co2_interpolated.reshape(target_shape)
                
                # Store interpolated data for this level
                co2_horizontal_interpolated[co2_level] = co2_interpolated
                
                if verbose:
                    print(f"  Level {co2_level}: Horizontal interpolation complete")
            
            # Remove duplicate print for verbose output (already printed inside the loop)
            
        # Calculate and print horizontal interpolation timing
        horizontal_end_time = time.time()
        horizontal_duration = horizontal_end_time - horizontal_start_time
        if verbose:
            print(f"Horizontal interpolation completed in {horizontal_duration:.2f} seconds")
            
        # Step 3b: Vertical interpolation to match the target grid levels
        print("Performing vertical interpolation to match target grid levels...")
        
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
                if nearest_level in co2_horizontal_interpolated:
                    interpolated_data = co2_horizontal_interpolated[nearest_level].copy()
                else:
                    # If no matching level, use default value
                    interpolated_data.fill(0.0006)
            
            return interpolated_data
        
        # Initialize result dictionary
        co2_interpolated_3d = {}
        
        if use_dask:
            # Use Dask for parallel processing of vertical interpolation
            if verbose:
                print("Setting up Dask client for vertical interpolation...")
            
            # Use a context manager for the Dask client
            n_workers = n_workers or os.cpu_count()
            with LocalCluster(n_workers=n_workers, threads_per_worker=1) as cluster, Client(cluster) as client:
                if verbose:
                    print(f"Dask cluster ready with {n_workers} workers")
                
                # Create delayed computations for each target level
                delayed_results = []
                for target_level in icmgg_levels:
                    if verbose:
                        print(f"  Preparing target level {target_level} for parallel processing...")
                    
                    delayed_result = dask.delayed(interpolate_vertical_level)(
                        target_level, icmgg_levels, co2_levels, 
                        co2_horizontal_interpolated, level_data, target_shape
                    )
                    delayed_results.append((target_level, delayed_result))
                
                # Compute all levels in parallel
                if verbose:
                    print(f"Computing {len(delayed_results)} vertical interpolation tasks in parallel...")
                results = dask.compute(*[res for _, res in delayed_results])
                
                # Store results back in the dictionary
                for i, target_level in enumerate([level for level, _ in delayed_results]):
                    co2_interpolated_3d[target_level] = results[i]
                    if verbose:
                        print(f"  Level {target_level}: Vertical interpolation complete")
                        print(f"  Level {target_level}: CO2 range: {np.min(co2_interpolated_3d[target_level])} to {np.max(co2_interpolated_3d[target_level])}")
        else:
            # Sequential processing for each target level
            for target_level in icmgg_levels:
                print(f"  Processing target level {target_level}...")
                
                # Initialize array for this level with the same shape as the target grid
                co2_interpolated_3d[target_level] = np.zeros(target_shape)
                
                # Loop through each grid point and perform 1D vertical interpolation
                total_points = target_shape[0] * target_shape[1] if len(target_shape) > 1 else target_shape[0]
                
                # For efficiency, reshape arrays to 2D if they're not already
                flat_shape = (total_points,)
                
                # Check if array is 1D or 2D
                is_2d = len(target_shape) > 1
                
                # Create a vertical profile for each grid point
                for i in range(target_shape[0]):
                    if is_2d:
                        for j in range(target_shape[1]):
                            # Extract CO2 values for this grid point at all source levels
                            source_levels = np.array(co2_levels, dtype=float)
                            
                            # Extract CO2 values for this grid point at all source levels
                            source_values = np.array([co2_horizontal_interpolated[level][i, j] for level in co2_levels])
                            
                            # Create 1D interpolator for the vertical profile
                            try:
                                if len(source_levels) > 1:  # Need at least 2 points for interpolation
                                    vertical_interp = interp1d(source_levels, source_values, 
                                                            kind='linear', bounds_error=False, fill_value="extrapolate")
                                    # Interpolate to the target level
                                    co2_interpolated_3d[target_level][i, j] = vertical_interp(float(target_level))
                                else:
                                    # If only one source level, use that value for all target levels
                                    co2_interpolated_3d[target_level][i, j] = source_values[0]
                            except Exception as e:
                                # In case of interpolation errors, use nearest source level
                                if len(source_levels) > 0:
                                    nearest_idx = np.abs(source_levels - float(target_level)).argmin()
                                    co2_interpolated_3d[target_level][i, j] = source_values[nearest_idx]
                                else:
                                    # If no source levels, set to default value
                                    co2_interpolated_3d[target_level][i, j] = 0.0006  # Default CO2 value
                    else:
                        # For 1D arrays
                        source_levels = np.array(co2_levels, dtype=float)
                        # Extract CO2 values for this grid point at all source levels
                        source_values = np.array([co2_horizontal_interpolated[level][i] for level in co2_levels])
                        
                        # Create 1D interpolator for the vertical profile
                        try:
                            if len(source_levels) > 1:  # Need at least 2 points for interpolation
                                vertical_interp = interp1d(source_levels, source_values, 
                                                        kind='linear', bounds_error=False, fill_value="extrapolate")
                                # Interpolate to the target level
                                co2_interpolated_3d[target_level][i] = vertical_interp(float(target_level))
                            else:
                                # If only one source level, use that value for all target levels
                                co2_interpolated_3d[target_level][i] = source_values[0]
                        except Exception as e:
                            # In case of interpolation errors, use nearest source level
                            if len(source_levels) > 0:
                                nearest_idx = np.abs(source_levels - float(target_level)).argmin()
                                co2_interpolated_3d[target_level][i] = source_values[nearest_idx]
                            else:
                                # If no source levels, set to default value
                                co2_interpolated_3d[target_level][i] = 0.0006  # Default CO2 value
                
                if verbose:
                    print(f"  Level {target_level}: Vertical interpolation complete")
                    print(f"  Level {target_level}: CO2 range: {np.min(co2_interpolated_3d[target_level])} to {np.max(co2_interpolated_3d[target_level])}")
        
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
    
    try:
        print(f"Writing interpolated CO2 data to {output_file}...")
        
        def write_output_grib(output_file, icmgg_file, interpolated_co2_3d, level_data):
            """Write the interpolated CO2 data to a GRIB file"""
            try:
                # Copy the ICMGG file to the output file
                shutil.copy2(icmgg_file, output_file)
                
                # First, get a template message for each level from the ICMGG file
                template_gids = {}
                level_templates = {}
                
                # Read template messages for each level
                with open(icmgg_file, 'rb') as f:
                    while True:
                        gid = eccodes.codes_grib_new_from_file(f)
                        if gid is None:
                            break
                        
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
                
                # Append the interpolated CO2 to the output file for each level
                with open(output_file, 'ab') as out_file:
                    # Process each level
                    for level in sorted(interpolated_co2_3d.keys()):
                        interpolated_co2 = interpolated_co2_3d[level]
                        
                        # Choose the best template for this level
                        template_level = level if level in template_gids else list(template_gids.keys())[0]
                        template_gid = template_gids[template_level]
                        
                        # Clone the template
                        gid = eccodes.codes_clone(template_gid)
                            
                        # Create a completely new message for CO2 to avoid eccodes parameter database mapping
                        if verbose:
                            print(f"Level {level}: Creating new CO2 message from scratch...")
                        
                        try:
                            # Get the sample code for a new grib message (similar to grib_copy example)
                            # Get original data values and key attributes we need to preserve
                            original_values = eccodes.codes_get_values(gid)  # Store original data values
                            
                            # Store essential keys for grid definition
                            keys_to_preserve = [
                                'gridType', 'gridDefinitionTemplateNumber', 'Ni', 'Nj',
                                'iDirectionIncrementInDegrees', 'jDirectionIncrementInDegrees',
                                'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees',
                                'latitudeOfLastGridPointInDegrees', 'longitudeOfLastGridPointInDegrees',
                                'numberOfPointsAlongAParallel', 'numberOfPointsAlongAMeridian',
                                'binaryScaleFactor', 'decimalScaleFactor', 'packingType',
                                'referenceValue', 'bitsPerValue',
                                'typeOfLevel', 'gridDefinitionDescription', 'resolutionAndComponentFlags',
                                'dataRepresentationType', 'bitmapPresent'
                            ]
                            
                            preserved_values = {}
                            for key in keys_to_preserve:
                                try:
                                    preserved_values[key] = eccodes.codes_get(gid, key)
                                except Exception:
                                    # Not all keys may exist in all templates
                                    pass
                            
                            # Don't try to unset paramId as it's causing errors
                            # Instead, directly set the CO2 parameters without trying to clear first
                                
                            # Now set the CO2-specific metadata
                            # Set with low-level keys rather than high-level ones
                            eccodes.codes_set_long(gid, 'table2Version', 253)  # A table that works with co2
                            eccodes.codes_set_long(gid, 'indicatorOfParameter', 162)  # Direct parameter code
                            eccodes.codes_set_string(gid, 'shortName', 'co2')
                            eccodes.codes_set_string(gid, 'name', 'Carbon Dioxide')
                            eccodes.codes_set_string(gid, 'units', 'kg kg**-1')
                            
                            # Set level
                            eccodes.codes_set_long(gid, 'level', level)
                            
                            # Read back to verify
                            try:
                                current_shortname = eccodes.codes_get(gid, 'shortName')
                                if verbose:
                                    print(f"Level {level}: Final shortName is '{current_shortname}'")
                            except:
                                if verbose:
                                    print(f"Level {level}: Could not read back shortName")
                                
                        except Exception as e:
                            if verbose:
                                print(f"Level {level}: Error with custom CO2 message: {str(e)}")
                            else:
                                print(f"Error with CO2 metadata: {str(e)}")
                            
                        # Set level information
                        if level in level_templates:
                            type_of_level = level_templates[level]['type_of_level']
                        else:
                            type_of_level = 'hybrid'  # Default to hybrid level
                            
                        eccodes.codes_set(gid, 'typeOfLevel', type_of_level)
                        eccodes.codes_set(gid, 'level', level)
                        
                        # Print what grid type we're using if verbose
                        try:
                            grid_type = eccodes.codes_get(gid, 'gridType')
                            if verbose:
                                print(f"Level {level}: Using grid type: {grid_type}")
                        except:
                            if verbose:
                                print(f"Level {level}: Could not get grid type from template")
                        
                        # Set the data values
                        eccodes.codes_set_values(gid, interpolated_co2)
                        
                        # Write the message to the file
                        eccodes.codes_write(gid, out_file)
                        
                        # Release this GRIB handle
                        eccodes.codes_release(gid)
                        
                        if verbose:
                            print(f"Added CO2 data for level {level}")
                    
                    # Release all template handles
                    for gid in template_gids.values():
                        eccodes.codes_release(gid)
                    
                print(f"Output GRIB file written to {output_file} with {len(interpolated_co2_3d)} levels of CO2 data")
                return True
            except Exception as e:
                print(f"Error writing output file: {str(e)}")
                print(traceback.format_exc())
                return False
        
        write_output_grib(output_file, icmgg_iniua_file, co2_interpolated_3d, level_data)
        
        print(f"Successfully wrote interpolated CO2 data to {output_file}")
        return output_file
        
    except Exception as e:
        print(f"Error writing output file: {e}")
        return None

def main():
    """
    Example usage of the interpolate_co2_to_icmgg function.
    """
    import argparse
    import dask
    
    parser = argparse.ArgumentParser(description="Interpolate CO2 data from a GRIB file to an ICMGG grid")
    parser.add_argument('co2_grib', help='Path to CO2 GRIB file')
    parser.add_argument('icmgg_iniua', help='Path to ICMGG INIUA file')
    parser.add_argument('--output', help="Path to save the output GRIB file")
    parser.add_argument('--var', default='cc', help="Name of the variable in the ICMGG file to match the grid to")
    parser.add_argument('--co2_var', default='co2', help="Name of the CO2 variable in the GRIB file")
    parser.add_argument('--verbose', action='store_true', help="Enable verbose output with detailed timing")
    parser.add_argument('--dask', action='store_true', help="Use Dask for parallel processing")
    parser.add_argument('--workers', type=int, default=None, help="Number of Dask workers to use (default: number of CPU cores)")
    
    args = parser.parse_args()
    
    output_file = args.output or f"{os.path.splitext(args.icmgg_iniua)[0]}_WITH_CO2{os.path.splitext(args.icmgg_iniua)[1]}"
    
    if args.dask:
        dask.config.set(scheduler='processes')
        if args.workers is None:
            args.workers = os.cpu_count()
        dask.config.set(num_workers=args.workers)
    
    result = interpolate_co2_to_icmgg(args.co2_grib, args.icmgg_iniua, output_file, args.var, args.co2_var, args.verbose, args.dask, args.workers)

    if result:
        print(f"Successfully interpolated CO2 data to {result}")
        
        # Check if the CO2 shortName was set correctly
        if args.verbose:
            try:
                print("\nVerifying final metadata in output file:")
                check_cmd = f"grib_ls {result} | grep -E 'shortName|paramId' | head -n 10"
                os.system(check_cmd)
            except Exception:
                pass
    else:
        print("Interpolation failed.")

if __name__ == "__main__":
    main()
