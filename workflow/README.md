# OCP-tool Snakemake Workflow

This directory contains the Snakemake workflow for automated climate model preparation.

## Usage

Configure your settings in `config/config.yaml`, then run:

```bash
# Local execution  
snakemake --cores 4

# With pixi
pixi run workflow-run

# Dry run to check workflow
pixi run workflow-dry
```

## Configuration

Edit `config/config.yaml` to specify:
- Grid resolution (`res_num`)  
- Experiment name (`exp_name_oifs`)
- Ocean grid name (`grid_name_oce`)
- Input/output paths

## Rules

- `prepare_gaussian_grids`: Process OpenIFS grid files
- `process_fesom_grid`: Handle ocean model grids
- `modify_land_sea_mask`: Core LSM modification  
- `generate_oasis_files`: Create OASIS3-MCT files
- `modify_runoff_maps`: Adjust runoff routing
- `generate_plots`: Create visualizations