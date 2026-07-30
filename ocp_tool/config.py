"""
Configuration module for OCP-Tool.
Handles loading YAML config and provides typed configuration dataclasses.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Union
import yaml

from .cycles import AUTO, CycleSpec, resolve_cycle


@dataclass
class AtmosphereConfig:
    """Atmosphere grid configuration."""
    resolution_list: List[int]
    truncation_type: str  # 'linear' or 'cubic-octahedral'
    experiment_name: str  # 4-digit ECMWF experiment code
    # OpenIFS cycle of the ICMGG/ICMSH input files: '43r3', '48r1' or 'auto'.
    # See ocp_tool/cycles.py — governs the snow-field layout.
    model_cycle: str = AUTO


@dataclass
class OceanConfig:
    """Ocean grid configuration."""
    grid_name: str
    has_ice_cavities: bool
    mesh_file: Optional[Path]
    intermediate_resolution: str = "r360x181"
    force_overwrite_griddes: bool = False


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
class PaleoConfig:
    """Paleo reconstruction configuration."""
    enabled: bool
    experiment_id: str  # e.g. "EP", "LP" — used in filenames
    reconstruction_dir: Path  # directory containing reconstruction NetCDF files
    modern_reference_dir: Path  # directory containing modern reference files
    icmsh_input_file: Path  # path to ICMSHaackINIT_CORE2 (spectral topo)
    calnoro_binary: Optional[Path]  # path to compiled calnoro binary
    # Reconstruction file names (relative to reconstruction_dir)
    ice_mask_file: str = "{exp_id}_icemask_v1.0.nc"
    topography_file: str = "{exp_id}_topo_v1.0.nc"
    lsm_file: str = "{exp_id}_LSM_v1.0.nc"
    lake_file: str = "{exp_id}_lake_v1.0.nc"
    soil_file: str = "{exp_id}_soil_v1.0.nc"
    biome_file: str = "{exp_id}_mbiome_v1.0.nc"
    # Modern reference file names (relative to modern_reference_dir)
    modern_topo_file: str = "Modern_std_topo_v1.0.nc"
    modern_lake_file: str = "Modern_std_soil_lake_v1.0.nc"

    def get_reconstruction_file(self, file_attr: str) -> Path:
        """Get full path to a reconstruction file, substituting experiment_id."""
        template = getattr(self, file_attr)
        filename = template.format(exp_id=self.experiment_id)
        return self.reconstruction_dir / filename

    def get_modern_file(self, file_attr: str) -> Path:
        """Get full path to a modern reference file."""
        filename = getattr(self, file_attr)
        return self.modern_reference_dir / filename


@dataclass
class ProcessingOptions:
    """Processing options."""
    verbose: bool
    parallel_workers: int
    use_dask: bool
    generate_rmp: bool = True  # Generate OASIS remapping weight files


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
    paleo: Optional[PaleoConfig] = None
    
    @property
    def co2_grib_file(self) -> Path:
        return self.input_paths.openifs_default / 'cams_co2_initial.grib'
    
    @property
    def co2_emissions_grib_file(self) -> Path:
        return self.input_paths.openifs_default / 'cams_co2_emissions.grib'
    
    @property
    def albedo_file(self) -> Path:
        return self.input_paths.openifs_default / 'bare_soil_albedos.grb'
    
    @property
    def cycle(self) -> CycleSpec:
        """
        OpenIFS cycle of the input files, resolved once and cached.

        Either taken from ``atmosphere.model_cycle`` or auto-detected from the
        input ICMGG file. An explicit setting that contradicts the file raises.
        """
        cached = getattr(self, '_cycle_spec', None)
        if cached is None:
            cached = resolve_cycle(
                self.atmosphere.model_cycle,
                self.get_icmgg_input_file(),
                verbose=self.options.verbose,
            )
            self._cycle_spec = cached
        return cached

    def get_icmgg_input_file(self) -> Path:
        """Get path to input ICMGG file."""
        return self.input_paths.openifs_default / f'ICMGG{self.atmosphere.experiment_name}INIT'

    def get_icmgg_output_file(self) -> Path:
        """Get path to output ICMGG INIT file."""
        return self.output_paths.openifs_modified / f'ICMGG{self.atmosphere.experiment_name}INIT_{self.ocean.grid_name}'
    
    def get_icmgg_iniua_input_file(self) -> Path:
        """Get path to input ICMGG INIUA file."""
        return self.input_paths.openifs_default / f'ICMGG{self.atmosphere.experiment_name}INIUA'

    def get_icmgg_iniua_file(self) -> Path:
        """Get path to output ICMGG INIUA file."""
        return self.output_paths.openifs_modified / f'ICMGG{self.atmosphere.experiment_name}INIUA'
    
    def get_icmsh_input_file(self) -> Path:
        """Get path to input ICMSH file (spectral topography)."""
        if self.paleo and self.paleo.icmsh_input_file:
            return self.paleo.icmsh_input_file
        return self.input_paths.openifs_default / f'ICMSH{self.atmosphere.experiment_name}INIT'
    
    def get_icmsh_output_file(self) -> Path:
        """Get path to output ICMSH file."""
        suffix = f'_{self.paleo.experiment_id}' if self.paleo else f'_{self.ocean.grid_name}'
        return self.output_paths.openifs_modified / f'ICMSH{self.atmosphere.experiment_name}INIT{suffix}'


def load_config(config_path: Union[str, Path]) -> OCPConfig:
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to config.yaml file
        
    Returns:
        OCPConfig dataclass with all configuration
    """
    config_path = Path(config_path).resolve()
    
    with open(config_path, 'r') as f:
        raw = yaml.safe_load(f)
    
    # Auto-detect root_dir: use config file's parent's parent (configs/ -> root)
    # or explicit path if provided
    if 'root_dir' in raw.get('paths', {}) and raw['paths']['root_dir']:
        root_dir = Path(raw['paths']['root_dir'])
    else:
        # Assume config is in <root>/configs/
        root_dir = config_path.parent.parent
    
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
    
    # Auto-compute output directory name from atmosphere and ocean settings
    # Format: TCO{resolution}_{ocean_grid} (e.g., TCO95_CORE2, TCO319_CORE3)
    resolution = raw['atmosphere']['resolution_list'][0]  # Use first resolution
    ocean_grid = raw['ocean']['grid_name']
    trunc_prefix = "TCO" if truncation_type == "cubic-octahedral" else "TL"
    output_subdir = f"{trunc_prefix}{resolution}_{ocean_grid}"
    
    output_base = root_dir / "output" / output_subdir
    output_paths = OutputPaths(
        openifs_modified=output_base / "openifs_input_modified",
        runoff_modified=output_base / "runoff_map_modified",
        oasis=output_base / "oasis_mct3_input",
        lpj_guess=output_base / "lpj-guess",
        plots=output_base / "plots",
    )
    
    return OCPConfig(
        atmosphere=AtmosphereConfig(
            resolution_list=raw['atmosphere']['resolution_list'],
            truncation_type=truncation_type,
            experiment_name=raw['atmosphere']['experiment_name'],
            model_cycle=raw['atmosphere'].get('model_cycle', AUTO),
        ),
        ocean=OceanConfig(
            grid_name=raw['ocean']['grid_name'],
            has_ice_cavities=raw['ocean']['has_ice_cavities'],
            mesh_file=Path(raw['ocean']['mesh_file']) if raw['ocean']['mesh_file'] is not None else None,
            intermediate_resolution=raw['ocean'].get('intermediate_resolution', 'r360x181'),
            force_overwrite_griddes=raw['ocean'].get('force_overwrite_griddes', False),
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
            generate_rmp=raw['options'].get('generate_rmp', True),
        ),
        root_dir=root_dir,
        paleo=_load_paleo_config(raw, resolve_path) if 'paleo' in raw else None,
    )


def _load_paleo_config(raw: dict, resolve_path) -> Optional[PaleoConfig]:
    """Parse paleo section from raw YAML config."""
    paleo_raw = raw.get('paleo', {})
    if not paleo_raw.get('enabled', False):
        return None
    return PaleoConfig(
        enabled=True,
        experiment_id=paleo_raw['experiment_id'],
        reconstruction_dir=resolve_path(paleo_raw['reconstruction_dir']),
        modern_reference_dir=resolve_path(paleo_raw['modern_reference_dir']),
        icmsh_input_file=resolve_path(paleo_raw['icmsh_input_file']),
        calnoro_binary=resolve_path(paleo_raw['calnoro_binary']) if paleo_raw.get('calnoro_binary') else None,
        ice_mask_file=paleo_raw.get('ice_mask_file', '{exp_id}_icemask_v1.0.nc'),
        topography_file=paleo_raw.get('topography_file', '{exp_id}_topo_v1.0.nc'),
        lsm_file=paleo_raw.get('lsm_file', '{exp_id}_LSM_v1.0.nc'),
        lake_file=paleo_raw.get('lake_file', '{exp_id}_lake_v1.0.nc'),
        soil_file=paleo_raw.get('soil_file', '{exp_id}_soil_v1.0.nc'),
        biome_file=paleo_raw.get('biome_file', '{exp_id}_mbiome_v1.0.nc'),
        modern_topo_file=paleo_raw.get('modern_topo_file', 'Modern_std_topo_v1.0.nc'),
        modern_lake_file=paleo_raw.get('modern_lake_file', 'Modern_std_soil_lake_v1.0.nc'),
    )


# Constants
EARTH_RADIUS_M = 6371.0 * 1e3  # Earth radius in meters
