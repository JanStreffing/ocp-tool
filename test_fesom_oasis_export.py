#!/usr/bin/env python
"""Test pyfesom2 OASIS export functionality"""

import pyfesom2 as pf
import os
from pathlib import Path

# Get FESOM mesh path from config (directory containing mesh.nc)
mesh_path = Path('/work/ab0246/a270092/input/fesom2/CORE2/')

print(f"Loading FESOM mesh from: {mesh_path}")
mesh = pf.load_mesh(str(mesh_path))

print(f"\nMesh info:")
print(f"  Number of 2D nodes: {mesh.n2d}")
print(f"  Number of 2D elements: {mesh.e2d}")

# Create test output directory
output_dir = Path('/work/ab0246/a270092/software/ocp-tool/test_oasis_output')
output_dir.mkdir(exist_ok=True)

print(f"\nWriting OASIS files to: {output_dir}")
result = pf.write_fesom_oasis_files(
    mesh=mesh,
    output_dir=str(output_dir),
    prefix='feom',
    overwrite=True
)

print("\nOutput files created:")
for key, path in result.items():
    size = os.path.getsize(path) / (1024*1024)  # MB
    print(f"  {key}: {path} ({size:.2f} MB)")

print("\n✓ FESOM→OASIS export successful!")
