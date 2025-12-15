"""
Configuration module for OCP-Tool.
Handles loading YAML config and provides typed configuration dataclasses.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional
import yaml


@dataclass
class AtmosphereConfig:
    """Atmosphere grid configuration."""
    resolution_list: List[int]
    truncation_type: str  # 'linear' or 'cubic-octahedral'
    experiment_name: str  # 4-digit ECMWF experiment code


@dataclass
class OceanConfig:
    """Ocean grid configuration."""
    grid_name: str
    has_ice_cavities: bool
    mesh_file: Path
    force_overwrite_griddes: bool


@dataclass
class RunoffConfig:
    """Runoff basin modification configuration."""
    manual_basin_removal: List[str]


@dataclass
class InputPaths:
    """Input directory paths."""
    fesom_mesh: Path
    gaussian_grids_full: Path
    gaussian_grids_reduced: Path  # Selected based on truncation_type
    openifs_default: Path
    runoff_default: Path
    lpj_guess: Path


@dataclass
class OutputPaths:
    """Output directory paths."""
    openifs_modified: Path
    runoff_modified: Path
    oasis: Path
    lpj_guess: Path
    plots: Path


@dataclass
class ProcessingOptions:
    """Processing options."""
    verbose: bool
    parallel_workers: int
    use_dask: bool


@dataclass
class OCPConfig:
    """Main configuration container."""
    atmosphere: AtmosphereConfig
    ocean: OceanConfig
    runoff: RunoffConfig
    input_paths: InputPaths
    output_paths: OutputPaths
    options: ProcessingOptions
    root_dir: Path
    
    @property
    def co2_grib_file(self) -> Path:
        return self.input_paths.openifs_default / 'cams_co2_initial.grib'
    
    @property
    def co2_emissions_grib_file(self) -> Path:
        return self.input_paths.openifs_default / 'cams_co2_emissions.grib'
    
    @property
    def albedo_file(self) -> Path:
        return self.input_paths.openifs_default / 'bare_soil_albedos.grb'
    
    def get_icmgg_input_file(self) -> Path:
        """Get path to input ICMGG file."""
        return self.input_paths.openifs_default / f'ICMGG{self.atmosphere.experiment_name}INIT'
    
    def get_icmgg_output_file(self) -> Path:
        """Get path to output ICMGG INIT file."""
        return self.output_paths.openifs_modified / f'ICMGG{self.atmosphere.experiment_name}INIT_{self.ocean.grid_name}'
    
    def get_icmgg_iniua_file(self) -> Path:
        """Get path to output ICMGG INIUA file."""
        return self.output_paths.openifs_modified / f'ICMGG{self.atmosphere.experiment_name}INIUA'


def load_config(config_path: str | Path) -> OCPConfig:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to config.yaml file
        
    Returns:
        OCPConfig dataclass with all configuration
    """
    config_path = Path(config_path)
    
    with open(config_path, 'r') as f:
        raw = yaml.safe_load(f)
    
    root_dir = Path(raw['paths']['root_dir'])
    
    # Determine reduced grid path based on truncation type
    truncation_type = raw['atmosphere']['truncation_type']
    if truncation_type == 'cubic-octahedral':
        reduced_grid_key = 'gaussian_grids_octahedral_reduced'
    elif truncation_type == 'linear':
        reduced_grid_key = 'gaussian_grids_linear_reduced'
    else:
        raise ValueError(f"Unknown truncation type: {truncation_type}")
    
    # Build paths - resolve relative paths against root_dir
    def resolve_path(path_str: str) -> Path:
        p = Path(path_str)
        if p.is_absolute():
            return p
        return root_dir / p
    
    input_paths = InputPaths(
        fesom_mesh=resolve_path(raw['paths']['input']['fesom_mesh']),
        gaussian_grids_full=resolve_path(raw['paths']['input']['gaussian_grids_full']),
        gaussian_grids_reduced=resolve_path(raw['paths']['input'][reduced_grid_key]),
        openifs_default=resolve_path(raw['paths']['input']['openifs_default']),
        runoff_default=resolve_path(raw['paths']['input']['runoff_default']),
        lpj_guess=resolve_path(raw['paths']['input']['lpj_guess']),
    )
    
    output_paths = OutputPaths(
        openifs_modified=resolve_path(raw['paths']['output']['openifs_modified']),
        runoff_modified=resolve_path(raw['paths']['output']['runoff_modified']),
        oasis=resolve_path(raw['paths']['output']['oasis']),
        lpj_guess=resolve_path(raw['paths']['output']['lpj_guess']),
        plots=resolve_path(raw['paths']['output']['plots']),
    )
    
    return OCPConfig(
        atmosphere=AtmosphereConfig(
            resolution_list=raw['atmosphere']['resolution_list'],
            truncation_type=truncation_type,
            experiment_name=raw['atmosphere']['experiment_name'],
        ),
        ocean=OceanConfig(
            grid_name=raw['ocean']['grid_name'],
            has_ice_cavities=raw['ocean']['has_ice_cavities'],
            mesh_file=Path(raw['ocean']['mesh_file']),
            force_overwrite_griddes=raw['ocean']['force_overwrite_griddes'],
        ),
        runoff=RunoffConfig(
            manual_basin_removal=raw['runoff']['manual_basin_removal'],
        ),
        input_paths=input_paths,
        output_paths=output_paths,
        options=ProcessingOptions(
            verbose=raw['options']['verbose'],
            parallel_workers=raw['options']['parallel_workers'],
            use_dask=raw['options']['use_dask'],
        ),
        root_dir=root_dir,
    )


# Constants
EARTH_RADIUS_M = 6371.0 * 1e3  # Earth radius in meters
