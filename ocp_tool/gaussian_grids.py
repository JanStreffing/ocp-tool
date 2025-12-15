"""
Grid processing module for OCP-Tool.
Handles Gaussian grid reading, coordinate calculation, and FESOM grid processing.
"""

import os
import subprocess
from pathlib import Path
from typing import Tuple, List, Optional
from dataclasses import dataclass

import numpy as np
from netCDF4 import Dataset

from .config import OCPConfig, EARTH_RADIUS_M


@dataclass
class GaussianGrid:
    """Container for Gaussian grid data."""
    center_lons: np.ndarray  # Shape (1, n_points)
    center_lats: np.ndarray  # Shape (1, n_points)
    corner_lons: np.ndarray  # Shape (4, 1, n_points)
    corner_lats: np.ndarray  # Shape (4, 1, n_points)
    cell_areas: np.ndarray   # Shape (1, n_points)
    lons_list: List[float]   # Flat list of longitudes
    nn: int                  # Grid parameter N


def read_grid_file(
    resolution: int,
    reduced_grid_path: Path,
    full_grid_path: Path,
    truncation_type: str,
    experiment_name: str = 'hagw',
    openifs_path: Optional[Path] = None,
    verbose: bool = False
) -> Tuple[List[str], int]:
    """
    Read reduced Gaussian grid file and return raw lines.
    
    Args:
        resolution: Truncation number (e.g., 95 for TCO95)
        reduced_grid_path: Path to reduced grid files
        full_grid_path: Path to full grid files
        truncation_type: 'linear' or 'cubic-octahedral'
        experiment_name: OpenIFS experiment name
        openifs_path: Path to OpenIFS input files
        verbose: Print debug info
        
    Returns:
        Tuple of (lines from grid file, NN parameter)
    """
    if truncation_type == 'linear':
        nn = int(resolution / 2 + 0.5)
        grid_file = reduced_grid_path / f'n{nn}_reduced.txt'
    elif truncation_type == 'cubic-octahedral':
        nn = resolution + 1
        grid_file = reduced_grid_path / f'o{nn}_reduced.txt'
    else:
        raise ValueError(f"Unknown truncation type: {truncation_type}")
    
    print(f"\n {'='*50} \n")
    print(f" Reading OpenIFS gridfile for T{resolution}")
    print(f"\n {'='*50} \n")
    
    if grid_file.exists():
        print(f" Read grid from file: {grid_file}")
        with open(grid_file, 'r') as f:
            lines = f.readlines()
    else:
        # Fall back to extracting from ICMGG file
        if openifs_path is None:
            raise FileNotFoundError(f"Grid file not found: {grid_file}")
        icmgg_file = openifs_path / f'ICMGG{experiment_name}INIT'
        print(f" Read grid from file: {icmgg_file}")
        grid_file = _read_grid_from_icmgg(icmgg_file, nn, truncation_type, reduced_grid_path)
        with open(grid_file, 'r') as f:
            lines = f.readlines()
    
    return lines, nn


def _read_grid_from_icmgg(
    icmgg_file: Path,
    nn: int,
    truncation_type: str,
    output_path: Path
) -> Path:
    """
    Extract grid information from ICMGG GRIB file using CDO.
    
    Returns path to generated grid file.
    """
    # Try GRIB2 first, then GRIB1
    try:
        os.system(f'grib_copy -w edition=2 {icmgg_file} {icmgg_file}.grib2')
        os.system(f'cdo griddes {icmgg_file}.grib2 > griddes.txt')
    except:
        os.system(f'grib_copy -w edition=1 {icmgg_file} {icmgg_file}.grib1')
        os.system(f'cdo griddes {icmgg_file}.grib1 > griddes.txt')
    
    # Parse grid description
    with open('griddes.txt', 'r') as f:
        lines = f.readlines()
    
    latitudes = []
    nlongitudes = []
    yline = rline = 0
    
    for i, line in enumerate(lines):
        if 'yvals' in line:
            yline = i
        elif 'rowlon' in line or 'reducedPoints' in line:
            rline = i
    
    # Read latitudes
    for i in range(yline, len(lines)):
        line = lines[i]
        if 'rowlon' in line or 'reducedPoints' in line:
            break
        if i == yline:
            tmp_lat = [float(lat) for lat in line.split()[2:]]
        else:
            tmp_lat = [float(lat) for lat in line.split()]
        latitudes.extend(tmp_lat)
    
    # Read number of longitudes per latitude
    for i in range(rline, len(lines)):
        line = lines[i]
        if 'scanningMode' in line:
            break
        if i == rline:
            tmp_nlon = [int(nlon) for nlon in line.split()[2:]]
        else:
            tmp_nlon = [int(nlon) for nlon in line.split()]
        nlongitudes.extend(tmp_nlon)
    
    # Generate output file
    if truncation_type == 'cubic-octahedral':
        rfile = output_path / f'o{nn}_reduced.txt'
    else:
        rfile = output_path / f'n{nn}_reduced.txt'
    
    with open(rfile, 'w') as f:
        f.write('latitude reduced regular latitude \n')
        f.write('number points points \n')
        f.write(' ------- ------- ------- ---------- \n')
        for ilat, (nlon, lat) in enumerate(zip(nlongitudes, latitudes)):
            f.write(f'{ilat+1} {nlon} {len(nlongitudes)*2} {lat} \n')
    
    return rfile


def extract_grid_data(
    lines: List[str],
    verbose: bool = False
) -> Tuple[List[float], List[float], List[int], List[float], List[float]]:
    """
    Parse reduced Gaussian grid file and extract coordinate data.
    
    Returns:
        Tuple of (lons_list, lats_list, numlons_list, dlon_list, lat_list)
    """
    lons_list = []
    lats_list = []
    numlons_list = []
    dlon_list = []
    lat_list = []
    
    for line in lines[3:]:  # Skip header
        if verbose:
            print(line)
        _, red_points, _, lat = (float(z) for z in line.split())
        
        dlon = 360.0 / red_points
        lons = np.arange(0, 360 - 1e-9, dlon)
        
        numlons_list.append(int(red_points))
        dlon_list.append(dlon)
        lat_list.append(lat)
        lons_list.extend(lons)
        lats_list.extend([lat] * len(lons))
    
    return lons_list, lats_list, numlons_list, dlon_list, lat_list


def calculate_corner_coordinates(
    lats_list: List[float],
    lons_list: List[float],
    numlons_list: List[int],
    dlon_list: List[float],
    lat_list: List[float],
    verbose: bool = False
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate corner coordinates for each grid cell.
    
    Corner layout:
        2 ---------- 1
        |            |
        |            |
        3 ---------- 4
    
    Returns:
        Tuple of (center_lats, center_lons, corner_lats, corner_lons)
    """
    center_lons = np.array(lons_list, dtype='float32')[np.newaxis, :]
    center_lats = np.array(lats_list, dtype='float32')[np.newaxis, :]
    nx = center_lons.shape[1]
    ny = 1
    
    print(f" Size of grid: nx = {nx}, ny = {ny}")
    
    corner_lons = np.zeros((4, ny, nx))
    corner_lats = np.zeros((4, ny, nx))
    
    kk = 0
    for ii, ni in enumerate(numlons_list):
        dlon = dlon_list[ii]
        lat = lat_list[ii]
        lons = np.arange(0, 360, dlon)
        
        # Calculate latitude deltas
        if ii == 0:  # First latitude - north pole
            dlat_n = 90 - lat
            dlat_s = (lat - lat_list[ii + 1]) / 2.0
        elif ii == len(numlons_list) - 1:  # Last latitude - south pole
            dlat_n = (lat_list[ii - 1] - lat) / 2.0
            dlat_s = lat + 90
        else:
            dlat_n = (lat_list[ii - 1] - lat) / 2.0
            dlat_s = (lat - lat_list[ii + 1]) / 2.0
        
        for jj in range(ni):
            # Corner 1: north-east
            corner_lons[0, 0, kk] = lons[jj] + dlon / 2.0
            corner_lats[0, 0, kk] = lat + dlat_n
            # Corner 2: north-west
            corner_lons[1, 0, kk] = lons[jj] - dlon / 2.0
            corner_lats[1, 0, kk] = lat + dlat_n
            # Corner 3: south-west
            corner_lons[2, 0, kk] = lons[jj] - dlon / 2.0
            corner_lats[2, 0, kk] = lat - dlat_s
            # Corner 4: south-east
            corner_lons[3, 0, kk] = lons[jj] + dlon / 2.0
            corner_lats[3, 0, kk] = lat - dlat_s
            kk += 1
    
    # Convert to [-180, 180] range
    center_lons = np.where(center_lons > 180, center_lons - 360, center_lons)
    corner_lons = np.where(corner_lons > 180, corner_lons - 360, corner_lons)
    
    return center_lats, center_lons, corner_lats, corner_lons


def calculate_cell_areas(
    center_lons: np.ndarray,
    numlons_list: List[int],
    dlon_list: List[float],
    lat_list: List[float],
    verbose: bool = False
) -> np.ndarray:
    """
    Calculate grid cell areas in square meters.
    
    Returns:
        Array of cell areas with shape (1, n_points)
    """
    nx = center_lons.shape[1]
    cell_areas = np.zeros((1, nx))
    
    kk = 0
    for ii, ni in enumerate(numlons_list):
        dlon = dlon_list[ii]
        lat = lat_list[ii]
        
        # Calculate latitude deltas
        if ii == 0:
            dlat_n = 90 - lat
            dlat_s = (lat - lat_list[ii + 1]) / 2.0
        elif ii == len(numlons_list) - 1:
            dlat_n = (lat_list[ii - 1] - lat) / 2.0
            dlat_s = lat + 90
        else:
            dlat_n = (lat_list[ii - 1] - lat) / 2.0
            dlat_s = (lat - lat_list[ii + 1]) / 2.0
        
        # Cell area in m²
        dx = dlon * np.pi / 180.0 * EARTH_RADIUS_M * np.cos(np.pi / 180.0 * lat)
        dy = (dlat_n + dlat_s) * np.pi / 180.0 * EARTH_RADIUS_M
        area = dx * dy
        
        for jj in range(ni):
            cell_areas[0, kk] = area
            kk += 1
    
    return cell_areas


def generate_gaussian_grid(
    config: OCPConfig,
    resolution: int
) -> GaussianGrid:
    """
    Generate complete Gaussian grid data for a given resolution.
    
    Args:
        config: OCP configuration
        resolution: Truncation number
        
    Returns:
        GaussianGrid dataclass with all grid data
    """
    lines, nn = read_grid_file(
        resolution=resolution,
        reduced_grid_path=config.input_paths.gaussian_grids_reduced,
        full_grid_path=config.input_paths.gaussian_grids_full,
        truncation_type=config.atmosphere.truncation_type,
        experiment_name=config.atmosphere.experiment_name,
        openifs_path=config.input_paths.openifs_default,
        verbose=config.options.verbose
    )
    
    lons_list, lats_list, numlons_list, dlon_list, lat_list = extract_grid_data(
        lines, verbose=config.options.verbose
    )
    
    center_lats, center_lons, corner_lats, corner_lons = calculate_corner_coordinates(
        lats_list, lons_list, numlons_list, dlon_list, lat_list,
        verbose=config.options.verbose
    )
    
    cell_areas = calculate_cell_areas(
        center_lons, numlons_list, dlon_list, lat_list,
        verbose=config.options.verbose
    )
    
    return GaussianGrid(
        center_lons=center_lons,
        center_lats=center_lats,
        corner_lons=corner_lons,
        corner_lats=corner_lats,
        cell_areas=cell_areas,
        lons_list=lons_list,
        nn=nn
    )


def read_fesom_grid(
    config: OCPConfig,
    verbose: bool = False
) -> np.ndarray:
    """
    Read FESOM ocean grid and generate OpenIFS-compatible land-sea mask.
    
    Uses prep_fesom.sh script to process the mesh file.
    
    Returns:
        Array of land-sea mask values on the OpenIFS grid
    """
    print(f"\n {'='*50} \n")
    print(" Using cdo and nco to interpolate from ocean to OpenIFS land sea mask")
    
    fesom_mesh_dir = config.input_paths.fesom_mesh
    ocean_grid = config.ocean.grid_name
    mesh_file = config.ocean.mesh_file
    interp_res = config.ocean.intermediate_resolution
    cavity = str(config.ocean.has_ice_cavities)
    exp_name = config.atmosphere.experiment_name
    
    # Save current directory and change to fesom_mesh directory
    original_dir = os.getcwd()
    os.chdir(fesom_mesh_dir)
    
    try:
        print(f" Does FESOM grid path exist? {mesh_file.exists()}")
        print(f" Is the overwrite active? {config.ocean.force_overwrite_griddes}")
        
        if mesh_file.exists() and not config.ocean.force_overwrite_griddes:
            print(f" Using existing grid description file '{mesh_file}'")
            cmd = f'./prep_fesom.sh {mesh_file} {ocean_grid} {interp_res} ../openifs_input_default/ICMGG{exp_name}INIT {cavity}'
        else:
            # Would need pyfesom2 to create mesh - for now, require existing file
            raise FileNotFoundError(f"Mesh file not found and force_overwrite_griddes=False: {mesh_file}")
        
        print(f"\n {'='*50} \n")
        print(" Using the following command to generate OpenIFS lsm based on FESOM mesh description file:")
        print(cmd)
        print(f" Reading ocean based land sea mask: {ocean_grid}")
        
        os.system(cmd)
        
        # Read the generated file
        oifs_file = fesom_mesh_dir / f'{ocean_grid}_oifs.nc'
        mesh = Dataset(str(oifs_file), 'r')
        
        if verbose:
            print(mesh.variables.keys())
        
        fesom_lsm = mesh.variables['cell_area'][:]
        mesh.close()
        
        return fesom_lsm
        
    finally:
        os.chdir(original_dir)
