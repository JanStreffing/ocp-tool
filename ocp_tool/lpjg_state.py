#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LPJ-GUESS State File Interpolation

This module handles reading, interpolating, and writing LPJ-GUESS binary
state files. It treats gridcell state blobs as opaque binary data and
uses nearest-neighbor interpolation to remap between resolutions.

State file format:
- meta.bin: num_processes, vegmode, npft, pft_names
- N.state: gridcell blobs + index (coords + file positions)
"""

import os
import struct
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from scipy.spatial import cKDTree


@dataclass
class LPJGMetadata:
    """Metadata from meta.bin file."""
    num_processes: int
    vegmode: int
    npft: int
    pft_names: List[str]


@dataclass
class GridcellBlob:
    """A single gridcell's state as opaque binary data."""
    lon: float
    lat: float
    data: bytes


def read_meta_bin(state_dir: Path) -> LPJGMetadata:
    """
    Read the meta.bin file from an LPJ-GUESS state directory.
    
    Args:
        state_dir: Path to the state directory containing meta.bin
        
    Returns:
        LPJGMetadata with process count, vegetation mode, and PFT info
    """
    meta_path = state_dir / "meta.bin"
    
    with open(meta_path, 'rb') as f:
        # Read header
        num_processes = struct.unpack('i', f.read(4))[0]
        vegmode = struct.unpack('i', f.read(4))[0]
        npft = struct.unpack('i', f.read(4))[0]
        
        # Read PFT names
        pft_names = []
        for _ in range(npft):
            name_len = struct.unpack('Q', f.read(8))[0]  # unsigned long
            name = f.read(name_len).decode('ascii')
            pft_names.append(name)
    
    return LPJGMetadata(
        num_processes=num_processes,
        vegmode=vegmode,
        npft=npft,
        pft_names=pft_names
    )


def write_meta_bin(state_dir: Path, metadata: LPJGMetadata, num_processes: int) -> None:
    """
    Write a meta.bin file with updated process count.
    
    Args:
        state_dir: Path to output state directory
        metadata: Original metadata (vegmode, npft, pft_names preserved)
        num_processes: New number of processes for output
    """
    meta_path = state_dir / "meta.bin"
    
    with open(meta_path, 'wb') as f:
        f.write(struct.pack('i', num_processes))
        f.write(struct.pack('i', metadata.vegmode))
        f.write(struct.pack('i', metadata.npft))
        
        for name in metadata.pft_names:
            name_bytes = name.encode('ascii')
            f.write(struct.pack('Q', len(name_bytes)))
            f.write(name_bytes)


def read_state_file(state_path: Path) -> List[GridcellBlob]:
    """
    Read a single .state file and extract all gridcell blobs.
    
    Args:
        state_path: Path to a .state file
        
    Returns:
        List of GridcellBlob objects with coordinates and binary data
    """
    with open(state_path, 'rb') as f:
        # Read entire file
        f.seek(0, 2)  # Seek to end
        file_size = f.tell()
        
        # Read number of elements from end of file
        f.seek(file_size - 8)  # size_t is 8 bytes
        num_elements = struct.unpack('Q', f.read(8))[0]
        
        if num_elements == 0:
            return []
        
        # Calculate index start position
        # Each index entry: 2 doubles (lon, lat) + 1 streamsize (8 bytes)
        index_entry_size = 8 + 8 + 8  # lon, lat, file_position
        index_size = num_elements * index_entry_size + 8  # +8 for num_elements
        index_start = file_size - index_size
        
        # Read index
        f.seek(index_start)
        index = []
        for _ in range(num_elements):
            lon = struct.unpack('d', f.read(8))[0]
            lat = struct.unpack('d', f.read(8))[0]
            file_pos = struct.unpack('q', f.read(8))[0]  # streamsize is signed
            index.append((lon, lat, file_pos))
        
        # Sort index by file position to read blobs sequentially
        index_sorted = sorted(enumerate(index), key=lambda x: x[1][2])
        
        # Read blobs
        blobs = [None] * num_elements
        for i, (orig_idx, (lon, lat, file_pos)) in enumerate(index_sorted):
            # Determine blob size
            if i + 1 < len(index_sorted):
                next_pos = index_sorted[i + 1][1][2]
            else:
                next_pos = index_start
            
            blob_size = next_pos - file_pos
            
            f.seek(file_pos)
            data = f.read(blob_size)
            
            blobs[orig_idx] = GridcellBlob(lon=lon, lat=lat, data=data)
        
        return blobs


def read_all_state_files(state_dir: Path) -> Tuple[LPJGMetadata, List[GridcellBlob]]:
    """
    Read all state files from a directory.
    
    Args:
        state_dir: Path to state directory
        
    Returns:
        Tuple of (metadata, list of all gridcell blobs)
    """
    metadata = read_meta_bin(state_dir)
    
    all_blobs = []
    for rank in range(metadata.num_processes):
        state_path = state_dir / f"{rank}.state"
        if state_path.exists():
            blobs = read_state_file(state_path)
            all_blobs.extend(blobs)
            print(f"  Read {len(blobs)} gridcells from {rank}.state")
    
    print(f"  Total: {len(all_blobs)} gridcells from {metadata.num_processes} state files")
    return metadata, all_blobs


def write_state_file(state_path: Path, blobs: List[GridcellBlob]) -> None:
    """
    Write a single .state file with the given gridcell blobs.
    
    Args:
        state_path: Output path for the .state file
        blobs: List of GridcellBlob objects to write
    """
    with open(state_path, 'wb') as f:
        # Write blobs and track positions
        index = []
        for blob in blobs:
            file_pos = f.tell()
            f.write(blob.data)
            index.append((blob.lon, blob.lat, file_pos))
        
        # Sort index by coordinates (as LPJ-GUESS does)
        index_sorted = sorted(index)
        
        # Write index
        for lon, lat, file_pos in index_sorted:
            f.write(struct.pack('d', lon))
            f.write(struct.pack('d', lat))
            f.write(struct.pack('q', file_pos))
        
        # Write number of elements
        f.write(struct.pack('Q', len(blobs)))


def interpolate_lpjg_state(
    source_dir: Path,
    target_dir: Path,
    target_lons: np.ndarray,
    target_lats: np.ndarray,
    target_land_mask: np.ndarray,
    num_output_files: int = 1,
    verbose: bool = True
) -> Dict[str, any]:
    """
    Interpolate LPJ-GUESS state from source to target grid using nearest-neighbor.
    
    Uses the same KDTree-based NN approach as ocp-tool's other interpolations.
    
    Args:
        source_dir: Path to source state directory
        target_dir: Path to output state directory (will be created)
        target_lons: 1D array of target grid longitudes
        target_lats: 1D array of target grid latitudes
        target_land_mask: Boolean array, True = land (needs vegetation)
        num_output_files: Number of output .state files (MPI ranks)
        verbose: Print progress information
        
    Returns:
        Dictionary with statistics about the interpolation
    """
    if verbose:
        print(f"\n{'='*60}")
        print(" LPJ-GUESS State Interpolation")
        print(f"{'='*60}")
        print(f"  Source: {source_dir}")
        print(f"  Target: {target_dir}")
    
    # Create output directory
    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    
    # Read source state
    if verbose:
        print("\nStep 1: Reading source state files...")
    metadata, source_blobs = read_all_state_files(source_dir)
    
    if len(source_blobs) == 0:
        raise ValueError("No gridcells found in source state directory")
    
    # Extract source coordinates
    source_lons = np.array([b.lon for b in source_blobs])
    source_lats = np.array([b.lat for b in source_blobs])
    
    # Get target land points
    target_lons = np.asarray(target_lons).flatten()
    target_lats = np.asarray(target_lats).flatten()
    target_land_mask = np.asarray(target_land_mask).flatten().astype(bool)
    
    land_lons = target_lons[target_land_mask]
    land_lats = target_lats[target_land_mask]
    n_target_land = len(land_lons)
    
    if verbose:
        print(f"\nStep 2: Building KDTree for nearest-neighbor interpolation...")
        print(f"  Source land cells: {len(source_blobs)}")
        print(f"  Target land cells: {n_target_land}")
    
    # Build KDTree (reusing ocp-tool pattern)
    src_points = np.column_stack((source_lons, source_lats))
    target_points = np.column_stack((land_lons, land_lats))
    tree = cKDTree(src_points)
    distances, indices = tree.query(target_points)
    
    if verbose:
        print(f"  Max NN distance: {distances.max():.4f} degrees")
        print(f"  Mean NN distance: {distances.mean():.4f} degrees")
    
    # Create target blobs by copying source blobs with new coordinates
    if verbose:
        print(f"\nStep 3: Creating interpolated gridcells...")
    
    target_blobs = []
    for i, (lon, lat, src_idx) in enumerate(zip(land_lons, land_lats, indices)):
        # Copy blob data from nearest source, but use target coordinates
        src_blob = source_blobs[src_idx]
        target_blobs.append(GridcellBlob(
            lon=lon,
            lat=lat,
            data=src_blob.data  # Binary data copied unchanged
        ))
    
    # Distribute blobs across output files
    if verbose:
        print(f"\nStep 4: Writing {num_output_files} output state files...")
    
    blobs_per_file = len(target_blobs) // num_output_files
    remainder = len(target_blobs) % num_output_files
    
    start_idx = 0
    for rank in range(num_output_files):
        # Distribute remainder evenly
        n_blobs = blobs_per_file + (1 if rank < remainder else 0)
        end_idx = start_idx + n_blobs
        
        file_blobs = target_blobs[start_idx:end_idx]
        state_path = target_dir / f"{rank}.state"
        write_state_file(state_path, file_blobs)
        
        if verbose:
            print(f"  Wrote {len(file_blobs)} gridcells to {rank}.state")
        
        start_idx = end_idx
    
    # Write metadata
    write_meta_bin(target_dir, metadata, num_output_files)
    if verbose:
        print(f"  Wrote meta.bin (vegmode={metadata.vegmode}, npft={metadata.npft})")
    
    # Compute statistics
    stats = {
        'source_gridcells': len(source_blobs),
        'target_gridcells': n_target_land,
        'target_total_points': len(target_lons),
        'land_fraction': n_target_land / len(target_lons),
        'max_nn_distance': float(distances.max()),
        'mean_nn_distance': float(distances.mean()),
        'output_files': num_output_files,
        'pft_names': metadata.pft_names,
    }
    
    if verbose:
        print(f"\n{'='*60}")
        print(" Interpolation Complete")
        print(f"{'='*60}")
        print(f"  Source → Target: {stats['source_gridcells']} → {stats['target_gridcells']} gridcells")
        print(f"  Land fraction: {stats['land_fraction']:.2%}")
        print(f"  Max NN distance: {stats['max_nn_distance']:.4f}°")
    
    return stats


def interpolate_lpjg_state_with_config(
    source_dir: Path,
    target_dir: Path,
    grid,  # GaussianGrid from ocp-tool
    lsm_data,  # LSMData from ocp-tool
    num_output_files: int = 1,
    verbose: bool = True
) -> Dict[str, any]:
    """
    Convenience wrapper that uses ocp-tool's grid and LSM data.
    
    This integrates with the existing ocp-tool workflow, using the
    polygon-method land mask already computed.
    
    Args:
        source_dir: Path to source LPJ-GUESS state directory
        target_dir: Path to output state directory
        grid: GaussianGrid object with target grid coordinates
        lsm_data: LSMData object with land-sea mask (1=land, 0=ocean)
        num_output_files: Number of output .state files
        verbose: Print progress
        
    Returns:
        Dictionary with interpolation statistics
    """
    # Extract coordinates and land mask from ocp-tool objects
    target_lons = np.squeeze(grid.center_lons)
    target_lats = np.squeeze(grid.center_lats)
    
    # lsm_binary_land: 1 = land, 0 = ocean
    land_mask = np.squeeze(lsm_data.lsm_binary_land) == 1
    
    return interpolate_lpjg_state(
        source_dir=Path(source_dir),
        target_dir=Path(target_dir),
        target_lons=target_lons,
        target_lats=target_lats,
        target_land_mask=land_mask,
        num_output_files=num_output_files,
        verbose=verbose
    )
