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
from collections import defaultdict
from scipy.spatial import ConvexHull

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


def read_fesom_grid_polygon(
    config: OCPConfig,
    grid: 'GaussianGrid',
    verbose: bool = False
) -> np.ndarray:
    """
    Read FESOM ocean grid and generate OpenIFS-compatible land-sea mask using polygon method.
    
    This method uses exact triangulation to determine if each OpenIFS grid point
    falls inside the FESOM ocean mesh. It scales well for high-resolution grids
    as it doesn't require an intermediate regular grid.
    
    Algorithm:
    1. Build ocean boundary from coastal edges + calving front edges
    2. Create convex hulls for small disconnected cavity components
    3. Point is LAND if: outside ocean polygon OR inside cavity hull
    
    Args:
        config: OCP configuration
        grid: Gaussian grid data with center coordinates
        verbose: Print debug info
        
    Returns:
        Array of land-sea mask values (1=land, 0=ocean) on the OpenIFS grid
    """
    from shapely.geometry import Point, Polygon
    from shapely.ops import unary_union
    from shapely.prepared import prep
    
    print(f"\n {'='*50}")
    print(" Using POLYGON method for land-sea mask (no intermediate grid)")
    print(f" {'='*50}\n")
    
    mesh_file = config.ocean.mesh_file
    has_cavities = config.ocean.has_ice_cavities
    
    print(f" Loading FESOM mesh: {mesh_file}")
    mesh = Dataset(str(mesh_file), 'r')
    lon = mesh.variables['lon'][:]
    lat = mesh.variables['lat'][:]
    triag = mesh.variables['triag_nodes'][:] - 1  # Convert to 0-indexed
    
    if has_cavities and 'cav_nod_mask' in mesh.variables:
        cav_mask = mesh.variables['cav_nod_mask'][:]
        print(f" Cavity nodes: {int(np.sum(cav_mask == 1))}")
    else:
        cav_mask = np.zeros(len(lon))
        print(" No cavity mask found or cavities disabled")
    
    mesh.close()
    
    # Build edge information
    print(" Building edge information...")
    edge_triangles = defaultdict(int)
    for tri in triag:
        for e in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]:
            edge_triangles[tuple(sorted(e))] += 1
    
    # Classify edges
    cavity_ocean_edges = []
    coastal_edges = []
    
    for (n0, n1), count in edge_triangles.items():
        c0, c1 = cav_mask[n0], cav_mask[n1]
        if count == 2 and c0 != c1:
            # Internal edge between cavity and ocean
            cavity_ocean_edges.append((n0, n1))
        elif count == 1 and c0 == 0 and c1 == 0:
            # Boundary edge with no cavity nodes (coastline)
            coastal_edges.append((n0, n1))
    
    print(f" Calving front edges: {len(cavity_ocean_edges)}")
    print(f" Coastal edges: {len(coastal_edges)}")
    
    # Build ocean boundary and trace polygons
    new_boundary = coastal_edges + cavity_ocean_edges
    adj = defaultdict(set)
    for n0, n1 in new_boundary:
        adj[n0].add(n1)
        adj[n1].add(n0)
    
    visited = set()
    polygon_rings = []
    
    def trace_ring(start):
        ring = [start]
        current, prev = start, None
        while True:
            neighbors = adj[current] - {prev} if prev else adj[current]
            if not neighbors:
                break
            next_n = min(neighbors)
            edge = tuple(sorted([current, next_n]))
            if edge in visited:
                neighbors = neighbors - {next_n}
                if not neighbors:
                    break
                next_n = min(neighbors)
                edge = tuple(sorted([current, next_n]))
                if edge in visited:
                    break
            visited.add(edge)
            if next_n == start:
                break
            ring.append(next_n)
            prev, current = current, next_n
            if len(ring) > 100000:
                break
        return ring
    
    print(" Tracing polygon boundaries...")
    for node in sorted(adj.keys()):
        if any(tuple(sorted([node, n])) not in visited for n in adj[node]):
            ring = trace_ring(node)
            if len(ring) >= 3:
                polygon_rings.append(ring)
    
    print(f" Found {len(polygon_rings)} polygon rings")
    
    # Create shapely polygons
    shapely_polys = []
    for ring in polygon_rings:
        coords = [(lon[n], lat[n]) for n in ring]
        lons_ring = [c[0] for c in coords]
        if max(lons_ring) - min(lons_ring) > 180:
            continue  # Skip dateline-crossing polygons
        try:
            p = Polygon(coords)
            if p.is_valid and p.area > 0.0001:
                shapely_polys.append(p)
        except:
            pass
    
    print(f" Valid shapely polygons: {len(shapely_polys)}")
    
    ocean_poly = unary_union(shapely_polys)
    ocean_prep = prep(ocean_poly)
    
    # Build cavity hulls for small disconnected components
    cavity_prep = None
    if has_cavities and np.sum(cav_mask == 1) > 0:
        print(" Building cavity hulls for small components...")
        
        # Find connected components of cavity nodes
        node_adj = defaultdict(set)
        for tri in triag:
            for i in range(3):
                for j in range(i+1, 3):
                    node_adj[tri[i]].add(tri[j])
                    node_adj[tri[j]].add(tri[i])
        
        cavity_nodes = set(np.where(cav_mask == 1)[0])
        visited_nodes = set()
        components = []
        
        def bfs(start):
            comp = set()
            queue = [start]
            while queue:
                node = queue.pop(0)
                if node in visited_nodes:
                    continue
                visited_nodes.add(node)
                comp.add(node)
                for neighbor in node_adj[node]:
                    if neighbor in cavity_nodes and neighbor not in visited_nodes:
                        queue.append(neighbor)
            return comp
        
        for node in cavity_nodes:
            if node not in visited_nodes:
                components.append(bfs(node))
        
        print(f" Cavity components: {len(components)}")
        
        # Create convex hulls for small components
        cavity_hulls = []
        for comp in components:
            if len(comp) < 50 and len(comp) >= 3:
                nodes = list(comp)
                coords = np.column_stack([lon[nodes], lat[nodes]])
                try:
                    hull = ConvexHull(coords)
                    p = Polygon(coords[hull.vertices])
                    if p.is_valid:
                        cavity_hulls.append(p)
                except:
                    pass
        
        if cavity_hulls:
            cavity_poly = unary_union(cavity_hulls)
            cavity_prep = prep(cavity_poly)
            print(f" Small cavity hulls: {len(cavity_hulls)}")
    
    # Get OpenIFS grid coordinates
    oifs_lons = np.array(grid.lons_list)
    oifs_lats = grid.center_lats.flatten()
    
    # Normalize longitudes to -180 to 180
    oifs_lons_std = np.where(oifs_lons > 180, oifs_lons - 360, oifs_lons)
    
    print(f" Testing {len(oifs_lons)} OpenIFS grid points...")
    
    # Use matplotlib triangulation for EXACT point location (no search radius)
    from matplotlib.tri import Triangulation
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    
    print(" Building exact triangulation lookup (parallelized)...")
    
    # Check if triangle has ANY cavity node
    tri_has_cavity = np.any(cav_mask[triag] == 1, axis=1)
    
    # Split triangles by dateline
    tri_lons_all = lon[triag]
    lon_spans = tri_lons_all.max(axis=1) - tri_lons_all.min(axis=1)
    valid_std = lon_spans <= 180  # Non-dateline triangles
    valid_dl = lon_spans > 180    # Dateline-crossing triangles
    
    fesom_lsm = np.ones(len(oifs_lons))  # Default: land (1)
    
    # Pass 1: Standard coordinates for non-dateline triangles
    print(f" Pass 1: {np.sum(valid_std)} triangles (standard coords)...")
    triag_std = triag[valid_std]
    cavity_std = tri_has_cavity[valid_std]
    
    try:
        tri_std = Triangulation(lon, lat, triangles=triag_std)
        finder_std = tri_std.get_trifinder()
        tri_indices = finder_std(oifs_lons_std, oifs_lats)
        
        for i in range(len(oifs_lons)):
            tri_idx = tri_indices[i]
            if tri_idx >= 0:
                if cavity_std[tri_idx]:
                    fesom_lsm[i] = 1  # Cavity -> land
                else:
                    fesom_lsm[i] = 0  # Ocean
        
        print(f"   Found {np.sum(tri_indices >= 0)} points in standard triangles")
    except Exception as e:
        print(f" Warning: Standard triangulation failed: {e}")
    
    # Pass 2: Shifted coordinates for dateline triangles
    n_dl = np.sum(valid_dl)
    if n_dl > 0:
        print(f" Pass 2: {n_dl} triangles (dateline, shifted coords)...")
        
        # Shift longitudes to 0-360 for dateline handling
        lon_shifted = np.where(lon < 0, lon + 360, lon)
        oifs_lons_shifted = np.where(oifs_lons_std < 0, oifs_lons_std + 360, oifs_lons_std)
        
        triag_dl = triag[valid_dl]
        cavity_dl = tri_has_cavity[valid_dl]
        
        try:
            tri_dl = Triangulation(lon_shifted, lat, triangles=triag_dl)
            finder_dl = tri_dl.get_trifinder()
            
            # Only check points not yet classified as ocean
            unclassified = fesom_lsm == 1
            tri_indices_dl = finder_dl(oifs_lons_shifted[unclassified], oifs_lats[unclassified])
            
            found_dl = tri_indices_dl >= 0
            unclass_indices = np.where(unclassified)[0]
            
            for j, idx in enumerate(unclass_indices):
                if found_dl[j]:
                    if cavity_dl[tri_indices_dl[j]]:
                        fesom_lsm[idx] = 1
                    else:
                        fesom_lsm[idx] = 0
            
            print(f"   Found {np.sum(found_dl)} points in dateline triangles")
        except Exception as e:
            print(f" Warning: Dateline triangulation failed: {e}")
    
    land_count = int(np.sum(fesom_lsm == 1))
    ocean_count = int(np.sum(fesom_lsm == 0))
    print(f" Result: Land={land_count}, Ocean={ocean_count}")
    
    return fesom_lsm
