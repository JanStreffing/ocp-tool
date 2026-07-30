"""Entry point for running ocp-tool as a module: python3 -m ocp_tool [config.yaml]"""

import sys
import time
from pathlib import Path
from shutil import copy2

from ocp_tool.config import load_config
from ocp_tool.gaussian_grids import generate_gaussian_grid, read_fesom_grid_polygon
from ocp_tool.lsm import process_land_sea_mask, create_slt_output_for_lpjg
from ocp_tool.oasis_writer import write_oasis_grid_files, interpolate_vegin_data
from ocp_tool.runoff import modify_runoff_map, modify_runoff_lsm
from ocp_tool.plotting import plot_land_sea_mask, plot_runoff_maps
from ocp_tool.co2_interpolation import interpolate_co2_to_icmgg
from ocp_tool.field_interpolation import interpolate_2d_fields_to_icmgg
from ocp_tool.create_outputdirs import create_outputdirs


def run_ocp_tool(config):
    """Main OCP-Tool processing pipeline."""
    print("=" * 60)
    print(" OCP-Tool: OpenIFS Coupling Preparation Tool")
    print("=" * 60)
    print(f"\nConfiguration:")
    print(f"  Atmosphere: {config.atmosphere.truncation_type} T{config.atmosphere.resolution_list}")
    print(f"  Ocean grid: {config.ocean.grid_name}")
    print(f"  Ice cavities: {config.ocean.has_ice_cavities}")
    print(f"  Experiment: {config.atmosphere.experiment_name}")

    for resolution in config.atmosphere.resolution_list:
        print(f"\n{'='*60}")
        print(f" Processing resolution T{resolution}")
        print(f"{'='*60}\n")

        print("Step 0: Creating output directories...")
        create_outputdirs(config, resolution)

        print("Step 1: Generating Gaussian grid coordinates...")
        grid = generate_gaussian_grid(config, resolution)

        if config.ocean.grid_name != 'AMIP':
            print("\nStep 2: Reading FESOM ocean grid (polygon method)...")
            ocean_lsm = read_fesom_grid_polygon(config, grid, verbose=config.options.verbose)
        else:
            print("\nStep 2: Skipped reading FESOM mesh (AMIP mode)")
            ocean_lsm = []

        print("\nStep 3: Processing land-sea mask...")
        lsm_data = process_land_sea_mask(config, grid, ocean_lsm, resolution)

        print("\nStep 4: Writing OASIS grid/mask/area files...")
        write_oasis_grid_files(config, grid, lsm_data, resolution, parallel=False)

        print("\nStep 5: Interpolating vegetation and CO2 restart data...")
        interpolate_vegin_data(config, grid)

        print("\nStep 6: Creating land-sea mask plots...")
        plot_land_sea_mask(config, grid, lsm_data, resolution)

        print("\nStep 7: Interpolating 3D CO2 concentrations...")
        icmgg_iniua_file = config.get_icmgg_iniua_file()
        # Stage the source INIUA into the output dir (mirrors how the INIT
        # file is copied in process_land_sea_mask); CO2 is then added in place.
        icmgg_iniua_input = config.get_icmgg_iniua_input_file()
        copy2(icmgg_iniua_input, icmgg_iniua_file)
        interpolate_co2_to_icmgg(
            str(config.co2_grib_file),
            str(icmgg_iniua_file),
            output_file=str(icmgg_iniua_file),
            use_dask=config.options.use_dask,
            n_workers=config.options.parallel_workers
        )

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

        print("\nStep 9: Adding bare soil albedo fields...")
        interpolate_2d_fields_to_icmgg(
            str(config.albedo_file),
            str(icmgg_init_file),
            field_type='albedo',
            verbose=config.options.verbose
        )

        if config.ocean.grid_name != 'AMIP':
            print("\nStep 10: Modifying runoff maps...")
            lons, lats = modify_runoff_map(config, resolution)

            print("\nStep 11: Modifying runoff land-sea mask...")
            modify_runoff_lsm(config, lons, lats)
        else:
            print("\nStep 10-11: Skipped runoff processing (AMIP mode)")

        print("\nStep 12: Creating SLT output for LPJ-GUESS...")
        create_slt_output_for_lpjg(config, resolution)

    print("\n" + "=" * 60)
    print(" OCP-Tool processing complete!")
    print("=" * 60)


def main():
    """Main entry point."""
    if len(sys.argv) > 1:
        config_path = Path(sys.argv[1])
    else:
        config_path = Path('config.yaml')

    if not config_path.exists():
        print(f"Error: Configuration file not found: {config_path}")
        print("Usage: python3 -m ocp_tool [config.yaml]")
        sys.exit(1)

    print(f"Loading configuration from: {config_path}")
    config = load_config(config_path)

    start_time = time.time()
    run_ocp_tool(config)
    elapsed_time = time.time() - start_time
    print(f"\nTotal execution time: {elapsed_time:.2f} seconds")


if __name__ == '__main__':
    main()
