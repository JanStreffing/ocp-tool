"""
Land-Sea Mask processing module for OCP-Tool.
Handles reading, modifying, and writing GRIB land-sea mask files.
"""

import copy
from pathlib import Path
from typing import Tuple, List, Optional
from dataclasses import dataclass

import numpy as np
import gribapi
from shutil import copy2

from .config import OCPConfig
from .gaussian_grids import GaussianGrid


@dataclass
class LSMData:
    """Container for land-sea mask processing results."""
    lsm_binary_atm: np.ndarray      # LSM for atmosphere (lakes as land)
    lsm_binary_land: np.ndarray     # LSM for land points only
    lsm_binary_runoff: np.ndarray   # LSM for runoff mapping
    gribfield_modified: List[np.ndarray]  # Modified GRIB fields
    land_lats: List[float]          # Latitudes of land points
    land_lons: List[float]          # Longitudes of land points


def count_grib_fields(grib_file: Path) -> int:
    """
    Auto-detect the number of fields in a GRIB file.
    
    Args:
        grib_file: Path to GRIB file
        
    Returns:
        Number of GRIB messages in file
    """
    count = 0
    with open(grib_file, 'rb') as f:
        while True:
            gid = gribapi.grib_new_from_file(f)
            if gid is None:
                break
            gribapi.grib_release(gid)
            count += 1
    return count


def read_grib_fields(
    input_file: Path,
    num_fields: Optional[int] = None,
    verbose: bool = False
) -> Tuple[List[np.ndarray], int, int, int, List, int]:
    """
    Read GRIB fields from OpenIFS ICMGG file.
    
    Args:
        input_file: Path to ICMGG GRIB file
        num_fields: Number of fields to read (auto-detected if None)
        verbose: Print debug info
        
    Returns:
        Tuple of (gribfield, lsm_id, slt_id, cl_id, gid_list, num_fields)
    """
    print(f"\n {'='*50} \n")
    print(f" Opening Grib input file: {input_file}")
    
    # Auto-detect number of fields if not specified
    if num_fields is None:
        num_fields = count_grib_fields(input_file)
        print(f" Auto-detected {num_fields} fields in GRIB file")
    
    gid = [None] * num_fields
    gribfield = [None] * num_fields
    lsm_id = slt_id = cl_id = -1
    lsm_saved = False
    
    with open(input_file, 'rb') as f:
        keys = ['N', 'shortName']
        
        if verbose:
            print(f'num_fields: {num_fields}')
        
        for i in range(num_fields):
            gid[i] = gribapi.grib_new_from_file(f)
            if gid[i] is None:
                break
            
            for key in keys:
                if not gribapi.grib_is_defined(gid[i], key):
                    raise ValueError(f"Key '{key}' was not defined")
                if verbose:
                    print(f'{key}={gribapi.grib_get(gid[i], key)}')
            
            short_name = gribapi.grib_get(gid[i], 'shortName')
            
            if short_name == 'lsm':
                lsm_id = i
                num_timesteps = gribapi.grib_get(gid[i], 'numberOfDataPoints')
                print(f'number of lsm timesteps: {num_timesteps}')
            elif short_name == 'slt':
                slt_id = i
            elif short_name == 'cl':
                cl_id = i
            
            if not lsm_saved or short_name != 'lsm':
                values = gribapi.grib_get_values(gid[i])
                if verbose:
                    print(f"Shape of '{short_name}' values: {np.shape(values)}")
                gribfield[i] = values
    
    # Filter out None values
    gribfield = [field for field in gribfield if field is not None]
    print(f"\n {'='*50} \n")
    print(f'Shape of Gribfield: {np.shape(gribfield)}')
    
    return gribfield, lsm_id, slt_id, cl_id, gid, num_fields


def modify_lsm(
    gribfield: List[np.ndarray],
    ocean_lsm: np.ndarray,
    ocean_grid_name: str,
    lsm_id: int,
    slt_id: int,
    cl_id: int,
    grid: GaussianGrid,
    verbose: bool = False
) -> LSMData:
    """
    Modify land-sea mask based on ocean model grid.
    
    - Removes lakes that are ocean in the ocean model
    - Sets soil type to SANDY CLAY LOAM (6) for removed lakes
    - Creates masks for different purposes (atmosphere, land, runoff)
    
    Args:
        gribfield: List of GRIB field arrays
        ocean_lsm: Ocean model land-sea mask
        ocean_grid_name: Name of ocean grid
        lsm_id: Index of LSM field in gribfield
        slt_id: Index of soil type field
        cl_id: Index of lake cover field
        grid: Gaussian grid data
        verbose: Print debug info
        
    Returns:
        LSMData with modified masks and land point coordinates
    """
    # Create copies for different mask types
    lsm_binary_land = copy.deepcopy(gribfield[lsm_id])
    lsm_binary_land = lsm_binary_land[np.newaxis, :]
    lsm_binary_runoff = lsm_binary_land.copy()
    gribfield_mod = gribfield[:]
    
    if ocean_grid_name != 'AMIP':
        # Automatic lake removal based on ocean mask
        # Polygon method: ocean_lsm = 1 means land, 0 means ocean
        n_points = len(gribfield_mod[slt_id])
        
        for i in range(n_points - 1):
            # Point is land in IFS but ocean in FESOM -> make it ocean
            if gribfield_mod[lsm_id][i] >= 0.5 and ocean_lsm[i] < 0.5:
                gribfield_mod[slt_id][i] = 0
                gribfield_mod[lsm_id][i] = 0
            # Point is ocean in IFS but land in FESOM -> make it land
            elif gribfield_mod[lsm_id][i] <= 0.5 and ocean_lsm[i] >= 0.5:
                gribfield_mod[slt_id][i] = 6  # SANDY CLAY LOAM
                gribfield_mod[lsm_id][i] = 1
    else:
        print(" Skipped modifying OpenIFS grid, because we are in AMIP mode")
    
    # Mask with lakes counting as land (for atmosphere coupling)
    lsm_binary_atm = gribfield_mod[lsm_id]
    lsm_binary_atm = lsm_binary_atm[np.newaxis, :]
    
    # Generate list of land point coordinates for LPJ-GUESS
    land_lats = []
    land_lons = []
    n_points = len(gribfield_mod[slt_id])
    
    for i in range(n_points - 1):
        if gribfield_mod[lsm_id][i] >= 0.5:
            land_lats.append(grid.center_lats[0, i])
            land_lons.append(grid.center_lons[0, i])
    
    print(f'Number of land points: {len(land_lats)}')
    
    return LSMData(
        lsm_binary_atm=lsm_binary_atm,
        lsm_binary_land=lsm_binary_land,
        lsm_binary_runoff=lsm_binary_runoff,
        gribfield_modified=gribfield_mod,
        land_lats=land_lats,
        land_lons=land_lons
    )


def write_modified_grib(
    gribfield_mod: List[np.ndarray],
    input_file: Path,
    output_file: Path,
    num_fields: int,
    gid_list: List,
    lsm_id: int,
    verbose: bool = False
) -> None:
    """
    Write modified GRIB fields to output file.
    
    Args:
        gribfield_mod: Modified GRIB field arrays
        input_file: Original input GRIB file
        output_file: Output file path
        num_fields: Number of fields
        gid_list: GRIB handles from reading
        lsm_id: Index of LSM field
        verbose: Print debug info
    """
    # Copy original file to output location
    copy2(input_file, output_file)
    
    # Read LSM datadates for multi-timestep handling
    with open(output_file, 'rb') as f:
        lsm_datadates = []
        while True:
            gidi = gribapi.grib_new_from_file(f)
            if gidi is None:
                break
            short_name = gribapi.grib_get(gidi, 'shortName')
            if short_name == 'lsm':
                data_date = gribapi.grib_get(gidi, 'dataDate')
                lsm_datadates.append((gidi, data_date))
        
        print(f'lsm_datadates: {lsm_datadates}')
    
    # Write modified fields
    with open(output_file, 'r+b') as f:
        lsm_index = 0
        
        # Write LSM fields with correct dates
        for gidi, data_date in lsm_datadates:
            print(f'Overwriting lsm for date: {data_date}')
            
            if lsm_index < len(gribfield_mod):
                gribapi.grib_set_values(gidi, gribfield_mod[lsm_id].flatten())
                gribapi.grib_set(gidi, 'dataDate', data_date)
                gribapi.grib_write(gidi, f)
                lsm_index += 1
            else:
                print(f"Error: Index {lsm_index} out of range for gribfield_mod.")
            
            gribapi.grib_release(gidi)
        
        # Write non-LSM fields
        for i in range(num_fields):
            short_name = gribapi.grib_get(gid_list[i], 'shortName')
            if short_name != 'lsm':
                gribapi.grib_set_values(gid_list[i], gribfield_mod[i])
                gribapi.grib_write(gid_list[i], f)
            gribapi.grib_release(gid_list[i])


def process_land_sea_mask(
    config: OCPConfig,
    grid: GaussianGrid,
    ocean_lsm: np.ndarray,
    resolution: int
) -> LSMData:
    """
    Complete land-sea mask processing pipeline.
    
    Args:
        config: OCP configuration
        grid: Gaussian grid data
        ocean_lsm: Ocean model land-sea mask
        resolution: Truncation number
        
    Returns:
        LSMData with all mask data
    """
    input_file = config.get_icmgg_input_file()
    output_file = config.get_icmgg_output_file()
    
    # Read GRIB fields (auto-detect num_fields)
    gribfield, lsm_id, slt_id, cl_id, gid_list, num_fields = read_grib_fields(
        input_file,
        num_fields=None,  # Auto-detect
        verbose=config.options.verbose
    )
    
    # Modify LSM based on ocean grid
    lsm_data = modify_lsm(
        gribfield=gribfield,
        ocean_lsm=ocean_lsm,
        ocean_grid_name=config.ocean.grid_name,
        lsm_id=lsm_id,
        slt_id=slt_id,
        cl_id=cl_id,
        grid=grid,
        verbose=config.options.verbose
    )
    
    # Write modified GRIB
    write_modified_grib(
        gribfield_mod=lsm_data.gribfield_modified,
        input_file=input_file,
        output_file=output_file,
        num_fields=num_fields,
        gid_list=gid_list,
        lsm_id=lsm_id,
        verbose=config.options.verbose
    )
    
    return lsm_data


def create_slt_output_for_lpjg(
    config: OCPConfig,
    resolution: int
) -> bool:
    """
    Create separate SLT (soil type) output file for LPJ-GUESS.
    
    Uses eccodes (grib_copy) and CDO to extract soil type field
    and convert to NetCDF format.
    
    Args:
        config: OCP configuration
        resolution: Truncation number
        
    Returns:
        True if successful, False otherwise
    """
    import os
    
    print(f"\n {'='*50} \n")
    print(" Creating SLT output file for LPJG")
    print(f"\n {'='*50} \n")
    
    input_file = config.get_icmgg_output_file()
    
    # Generate output filename with ocean grid name
    if config.atmosphere.truncation_type == 'cubic-octahedral':
        slt_output_name = f'slt_TCO{resolution}_{config.ocean.grid_name}.nc'
    elif config.atmosphere.truncation_type == 'linear':
        slt_output_name = f'slt_TL{resolution}_{config.ocean.grid_name}.nc'
    else:
        raise ValueError(f"Unknown truncation type: {config.atmosphere.truncation_type}")
    
    slt_output_path = config.output_paths.lpj_guess / slt_output_name
    temp_grib_file = config.output_paths.lpj_guess / 'temp_slt_var43.grb'
    
    if not input_file.exists():
        print(f"Error: Input ICMGG file not found: {input_file}")
        return False
    
    try:
        print(f"Extracting SLT field (variable 43) from {input_file}")
        
        # Extract SLT field using grib_copy
        cmd_extract = f"grib_copy -w shortName=slt {input_file} {temp_grib_file}"
        if config.options.verbose:
            print(f"Running: {cmd_extract}")
        
        result = os.system(cmd_extract)
        if result != 0:
            print(f"Error: grib_copy command failed with exit code {result}")
            return False
        
        if not temp_grib_file.exists():
            print(f"Error: Temporary GRIB file not created: {temp_grib_file}")
            return False
        
        print(f"Converting GRIB to NetCDF: {slt_output_path}")
        
        # Convert to NetCDF
        cmd_convert = f"cdo -f nc copy {temp_grib_file} {slt_output_path}"
        if config.options.verbose:
            print(f"Running: {cmd_convert}")
        
        result = os.system(cmd_convert)
        if result != 0:
            print(f"Error: cdo command failed with exit code {result}")
            if temp_grib_file.exists():
                temp_grib_file.unlink()
            return False
        
        # Rename variable from 'slt' to 'var43' (LPJ-GUESS expects var43)
        cmd_rename = f"ncrename -v slt,var43 {slt_output_path}"
        if config.options.verbose:
            print(f"Running: {cmd_rename}")
        
        result = os.system(cmd_rename)
        if result != 0:
            print(f"Warning: ncrename failed, variable may still be named 'slt'")
        
        # Cleanup
        if temp_grib_file.exists():
            temp_grib_file.unlink()
            if config.options.verbose:
                print(f"Removed temporary file: {temp_grib_file}")
        
        if slt_output_path.exists():
            print(f"Successfully created SLT file: {slt_output_path}")
            if config.options.verbose:
                print(f"File size: {slt_output_path.stat().st_size} bytes")
                os.system(f"ncdump -h {slt_output_path}")
            return True
        else:
            print(f"Error: Output NetCDF file was not created: {slt_output_path}")
            return False
            
    except Exception as e:
        print(f"Unexpected error creating SLT file: {e}")
        import traceback
        traceback.print_exc()
        
        if temp_grib_file.exists():
            temp_grib_file.unlink()
        
        return False
