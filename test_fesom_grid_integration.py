#!/usr/bin/env python3
"""Test FESOM grid integration into ocp-tool OASIS files"""

from pathlib import Path
import sys

# Add ocp-tool to path
sys.path.insert(0, '/work/ab0246/a270092/software/ocp-tool')

from ocp_tool.config import load_config

# Load TCO95_CORE2 config
config = load_config('/work/ab0246/a270092/software/ocp-tool/configs/TCO95_CORE2.yaml')

print("Configuration loaded:")
print(f"  Ocean grid: {config.ocean.grid_name}")
print(f"  Mesh file: {config.ocean.mesh_file}")
print(f"  OASIS output dir: {config.output_paths.oasis}")

# Check if OASIS files exist from previous run
oasis_dir = config.output_paths.oasis
if not oasis_dir.exists():
    print(f"\nError: OASIS directory does not exist: {oasis_dir}")
    print("Run ocp-tool first to generate OpenIFS OASIS grids")
    sys.exit(1)

grids_file = oasis_dir / 'grids.nc'
masks_file = oasis_dir / 'masks.nc'
areas_file = oasis_dir / 'areas.nc'

if not all(f.exists() for f in [grids_file, masks_file, areas_file]):
    print("\nError: OASIS files not found. Run ocp-tool first.")
    sys.exit(1)

print(f"\nOASIS files found:")
print(f"  grids.nc: {grids_file.stat().st_size / (1024*1024):.2f} MB")
print(f"  masks.nc: {masks_file.stat().st_size / (1024*1024):.2f} MB")
print(f"  areas.nc: {areas_file.stat().st_size / (1024*1024):.2f} MB")

# Test adding FESOM grid
print("\nAdding FESOM grid to existing OASIS files...")

try:
    import pyfesom2 as pf
    
    mesh_path = Path(config.ocean.mesh_file).parent
    print(f"Loading FESOM mesh from: {mesh_path}")
    mesh = pf.load_mesh(str(mesh_path))
    print(f"  Loaded: {mesh.n2d} nodes, {mesh.e2d} elements")
    
    # Add FESOM grid
    pf.write_fesom_oasis_files(
        mesh=mesh,
        output_dir=str(oasis_dir),
        prefix='feom',
        overwrite=True
    )
    
    print("\n✓ FESOM grid added successfully!")
    
    # Verify files were updated
    print(f"\nUpdated OASIS files:")
    print(f"  grids.nc: {grids_file.stat().st_size / (1024*1024):.2f} MB")
    print(f"  masks.nc: {masks_file.stat().st_size / (1024*1024):.2f} MB")
    print(f"  areas.nc: {areas_file.stat().st_size / (1024*1024):.2f} MB")
    
    # Check grid contents
    from netCDF4 import Dataset
    with Dataset(str(grids_file), 'r') as nc:
        print(f"\nGrids in grids.nc:")
        grids = [v for v in nc.variables.keys() if '.lon' in v]
        for grid in sorted(grids):
            grid_name = grid.replace('.lon', '')
            print(f"  - {grid_name}")
    
except Exception as e:
    print(f"\nError: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
