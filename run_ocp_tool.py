#!/usr/bin/env python3
"""
OCP-Tool: OpenIFS Coupling Preparation Tool (Refactored Version)

Prepares input files for coupled OpenIFS-FESOM2 or OpenIFS-NEMO climate simulations.

Usage:
    python run_ocp_tool.py [config.yaml]
    
If no config file is specified, uses config.yaml in the current directory.
"""

import sys
import time
from os import makedirs
from pathlib import Path

from ocp_tool.config import load_config, OCPConfig
from ocp_tool.gaussian_grids import generate_gaussian_grid, read_fesom_grid_polygon
from ocp_tool.lsm import process_land_sea_mask, create_slt_output_for_lpjg
from ocp_tool.oasis_writer import write_oasis_grid_files, interpolate_vegin_data
from ocp_tool.runoff import modify_runoff_map, modify_runoff_lsm
from ocp_tool.plotting import plot_land_sea_mask, plot_runoff_maps
from ocp_tool.co2_interpolation import interpolate_co2_to_icmgg
from ocp_tool.field_interpolation import interpolate_2d_fields_to_icmgg
from ocp_tool.create_outputdirs import create_outputdirs


def run_ocp_tool(config: OCPConfig) -> None:
    """
    Main OCP-Tool processing pipeline.
    
    Args:
        config: Loaded OCP configuration
    """
    print("=" * 60)
    print(" OCP-Tool: OpenIFS Coupling Preparation Tool")
    print("=" * 60)
    print(f"\nConfiguration:")
    print(f"  Atmosphere: {config.atmosphere.truncation_type} T{config.atmosphere.resolution_list}")
    print(f"  Ocean grid: {config.ocean.grid_name}")
    print(f"  Ice cavities: {config.ocean.has_ice_cavities}")
    print(f"  Experiment: {config.atmosphere.experiment_name}")
    
    
    # Process each resolution
    for resolution in config.atmosphere.resolution_list:
        print(f"\n{'='*60}")
        print(f" Processing resolution T{resolution}")
        print(f" Output: ./output/TCO{resolution}_{config.ocean.grid_name}")
        print(f"{'='*60}\n")
        
        # Step 0: Create Output directories
        print("Step 0: Creating output directories...")
        create_outputdirs(config, resolution)

        # Step 1: Generate Gaussian grid
        print("Step 1: Generating Gaussian grid coordinates...")
        grid = generate_gaussian_grid(config, resolution)
        
        # Step 2: Read ocean grid (if not AMIP)
        if config.ocean.grid_name != 'AMIP':
            print("\nStep 2: Reading FESOM ocean grid (polygon method)...")
            ocean_lsm = read_fesom_grid_polygon(config, grid, verbose=config.options.verbose)
        else:
            print("\nStep 2: Skipped reading FESOM mesh (AMIP mode)")
            ocean_lsm = []
        
        # Step 3: Process land-sea mask
        print("\nStep 3: Processing land-sea mask...")
        lsm_data = process_land_sea_mask(config, grid, ocean_lsm, resolution)
        
        # Step 4: Write OASIS files
        print("\nStep 4: Writing OASIS grid/mask/area files...")
        write_oasis_grid_files(
            config, grid, lsm_data, resolution,
            parallel=config.options.parallel_workers > 1
        )
        
        # Step 5: Interpolate vegetation and CO2 data
        print("\nStep 5: Interpolating vegetation and CO2 restart data...")
        interpolate_vegin_data(config, grid)
        
        # Step 6: Plot land-sea mask
        print("\nStep 6: Creating land-sea mask plots...")
        plot_land_sea_mask(config, grid, lsm_data, resolution)
        
        # Step 7: Interpolate 3D CO2 to INIUA file
        print("\nStep 7: Interpolating 3D CO2 concentrations...")
        icmgg_iniua_file = config.get_icmgg_iniua_file()
        interpolate_co2_to_icmgg(
            str(config.co2_grib_file),
            str(icmgg_iniua_file),
            output_file=str(icmgg_iniua_file),
            use_dask=config.options.use_dask,
            n_workers=config.options.parallel_workers
        )
        
        # Step 8: Interpolate 2D CO2 emissions to INIT file
        print("\nStep 8: Processing CO2 emissions data...")
        icmgg_init_file = config.get_icmgg_output_file()
        interpolate_2d_fields_to_icmgg(
            str(config.co2_emissions_grib_file),
            str(icmgg_init_file),
            output_file=str(icmgg_init_file),
            variable_name='lsm',
            field_type='co2_emissions',
            verbose=config.options.verbose
        )
        
        # Step 9: Add albedo fields
        print("\nStep 9: Adding bare soil albedo fields...")
        interpolate_2d_fields_to_icmgg(
            str(config.albedo_file),
            str(icmgg_init_file),
            field_type='albedo',
            verbose=config.options.verbose
        )
        
        # Step 10: Modify runoff maps
        print("\nStep 10: Modifying runoff maps...")
        lons, lats = modify_runoff_map(config, resolution)
        
        # Step 11: Modify runoff LSM
        print("\nStep 11: Modifying runoff land-sea mask...")
        modify_runoff_lsm(config, lons, lats)
        
        # Step 12: Create SLT output for LPJ-GUESS
        print("\nStep 12: Creating SLT output for LPJ-GUESS...")
        create_slt_output_for_lpjg(config, resolution)
    
    print("\n" + "=" * 60)
    print(" OCP-Tool processing complete!")
    print("=" * 60)


def main():
    """Main entry point."""
    # Determine config file path
    if len(sys.argv) > 1:
        config_path = Path(sys.argv[1])
    else:
        config_path = Path('config.yaml')
    
    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        print("Usage: python run_ocp_tool.py [config.yaml]")
        sys.exit(1)
    
    print(f"Loading configuration from: {config_path}")
    
    # Load configuration
    config = load_config(config_path)
    
    # Run the tool with timing
    start_time = time.time()
    run_ocp_tool(config)
    elapsed_time = time.time() - start_time
    
    print(f"\nTotal execution time: {elapsed_time:.2f} seconds")


if __name__ == '__main__':
    main()
