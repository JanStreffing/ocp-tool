"""
Runoff map processing module for OCP-Tool.
Handles drainage basin and arrival point modifications.
"""

from pathlib import Path
from typing import List, Tuple

import numpy as np
from netCDF4 import Dataset
from shutil import copy2

from .config import OCPConfig


def modify_runoff_map(
    config: OCPConfig,
    resolution: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Modify runoff map drainage basins and arrival points.
    
    Applies basin-specific corrections for closed seas (Caspian, Black Sea)
    and fixes for specific rivers (Ob).
    
    Args:
        config: OCP configuration
        resolution: Truncation number
        
    Returns:
        Tuple of (lons, lats) arrays from runoff file
    """
    input_file = config.input_paths.runoff_default / 'runoff_maps.nc'
    output_file = config.output_paths.runoff_modified / f'srunoff_maps_{config.ocean.grid_name}.nc'
    
    # Remove existing output file if present
    if output_file.exists():
        output_file.unlink()
    
    # Copy input to output
    copy2(input_file, output_file)
    
    # Open for modification
    rnffile = Dataset(str(output_file), 'r+')
    print(rnffile.variables.keys())
    
    drainage = rnffile.variables['drainage_basin_id'][:]
    arrival = rnffile.variables['arrival_point_id'][:]
    calving = rnffile.variables['calving_point_id'][:]
    lons = rnffile.variables['lon'][:]
    lats = rnffile.variables['lat'][:]
    
    # Apply basin modifications
    for basin in config.runoff.manual_basin_removal:
        if basin == 'caspian-sea':
            _modify_caspian_sea(drainage, arrival, lons, lats)
        elif basin == 'black-sea':
            _modify_black_sea(drainage, arrival, lons, lats)
    
    # Fix Ob river arrival points
    _fix_ob_river(drainage, arrival, lons, lats)
    
    # Fix glacial calving maps
    _fix_calving_maps(calving, lons, lats)
    
    # Save modifications
    rnffile.variables['drainage_basin_id'][:] = drainage
    rnffile.variables['arrival_point_id'][:] = arrival
    rnffile.variables['calving_point_id'][:] = calving
    rnffile.close()
    
    return lons, lats


def _modify_caspian_sea(
    drainage: np.ndarray,
    arrival: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray
) -> None:
    """Modify Caspian Sea basin - redirect to basin 18 with Amazon discharge."""
    for lo, lon in enumerate(lons):
        if 46 < lon < 56:
            for la, lat in enumerate(lats):
                if 36 < lat < 47:
                    if drainage[la, lo] == -2:
                        drainage[la, lo] = 18
                        arrival[la, lo] = -1
        # Add artificial arrival points in Amazon discharge area
        if 313 < lon < 314.5:
            for la, lat in enumerate(lats):
                if 1 < lat < 2:
                    if arrival[la, lo] != -1:
                        arrival[la, lo] = 18


def _modify_black_sea(
    drainage: np.ndarray,
    arrival: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray
) -> None:
    """Modify Black Sea basin - redirect to basin 23/28."""
    for lo, lon in enumerate(lons):
        # Remove old basin
        if 27 < lon < 43:
            for la, lat in enumerate(lats):
                if 40.5 < lat < 48:
                    if drainage[la, lo] == -2:
                        drainage[la, lo] = 23
                        arrival[la, lo] = -1
        # Add new arrival points
        if 25 < lon < 26.5:
            for la, lat in enumerate(lats):
                if 38.5 < lat < 41:
                    if arrival[la, lo] != -1:
                        arrival[la, lo] = 23
        if 23.5 < lon < 25:
            for la, lat in enumerate(lats):
                if 38.5 < lat < 41:
                    if arrival[la, lo] != -1:
                        arrival[la, lo] = 28


def _fix_ob_river(
    drainage: np.ndarray,
    arrival: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray
) -> None:
    """Fix Ob river arrival points."""
    for lo, lon in enumerate(lons):
        # Remove old arrival points
        if 60 < lon < 70:
            for la, lat in enumerate(lats):
                if 60 < lat < 80:
                    if arrival[la, lo] == 13:
                        arrival[la, lo] = 6
        # Add new arrival points
        if 72 < lon < 75:
            for la, lat in enumerate(lats):
                if 65 < lat < 75:
                    if arrival[la, lo] == 6:
                        arrival[la, lo] = 13


def _fix_calving_maps(
    calving: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray
) -> None:
    """Fix glacial calving maps for Antarctica and Greenland."""
    # Antarctica - remove old calving points
    for lo, lon in enumerate(lons):
        for la, lat in enumerate(lats):
            if lat < -55:
                if calving[la, lo] == 66:
                    calving[la, lo] = -2
    
    # Antarctica - add new calving points
    for lo, lon in enumerate(lons):
        if 300 < lon < 320:
            for la, lat in enumerate(lats):
                if -70 < lat < -60:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 66
        if 315 < lon < 325:
            for la, lat in enumerate(lats):
                if -65 < lat < -55:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 66
        if 320 < lon < 360:
            for la, lat in enumerate(lats):
                if -60 < lat < -50:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 66
        if 170 < lon < 180:
            for la, lat in enumerate(lats):
                if -75 < lat < -65:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 66
    
    # Greenland
    for lo, lon in enumerate(lons):
        if 300 < lon < 310:
            for la, lat in enumerate(lats):
                if 50 < lat < 60:
                    if calving[la, lo] == -2:
                        calving[la, lo] = 1


def modify_runoff_lsm(
    config: OCPConfig,
    lons: np.ndarray,
    lats: np.ndarray
) -> None:
    """
    Modify runoff mapper land-sea mask in OASIS masks file.
    
    Args:
        config: OCP configuration
        lons: Longitude array from runoff file
        lats: Latitude array from runoff file
    """
    masks_file = config.output_paths.oasis / 'masks.nc'
    oasis = Dataset(str(masks_file), 'r+')
    
    RnfA = oasis.variables['RnfA.msk'][:]
    RnfO = oasis.variables['RnfO.msk'][:]
    
    for basin in config.runoff.manual_basin_removal:
        if basin == 'caspian-sea':
            for lo, lon in enumerate(lons):
                if 46 < lon < 56:
                    for la, lat in enumerate(lats):
                        if 36 < lat < 47:
                            RnfA[la, lo] = 0
                            RnfO[la, lo] = 1
    
    oasis.variables['RnfA.msk'][:] = RnfA
    oasis.variables['RnfO.msk'][:] = RnfO
    oasis.close()
