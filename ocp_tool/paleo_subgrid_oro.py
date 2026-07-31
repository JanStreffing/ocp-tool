"""
Paleo subgrid-scale orography module for OCP-Tool.

Replaces convert_topodata.sh + calnoro + create_subgrid_oro.sh + mod_oro_oifs.sh:
1. Prepare topography input for calnoro (remap to r720x360, rename to OROMEA)
2. Run compiled calnoro binary via subprocess
3. Post-process calnoro output (SERVICE format → NetCDF, rename variables)
4. Apply SSO parameters (sdor, isor, anor, slor) to ICMGGaackINIT

All interpolation done in Python — no CDO dependency.
Calnoro itself remains a Fortran binary called via subprocess.
"""

import os
import struct
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
import eccodes

from .config import OCPConfig, PaleoConfig
from .gaussian_grids import GaussianGrid
from .paleo_input import (
    _read_netcdf_field,
    _interp_to_gaussian,
    _replace_grib_fields,
)


# SSO variable mapping: calnoro output code → (OpenIFS GRIB code, name, scale_factor)
SSO_MAPPING = {
    52: (160, 'OROSTD', 1.0),       # standard deviation of orography
    53: (163, 'OROSIG', 10.0),      # slope of sub-gridscale orography (×10)
    54: (161, 'OROGAM', 1.0),       # anisotropy of sub-gridscale orography
    55: (162, 'OROTHE', np.pi/180), # angle of sub-gridscale orography (deg→rad)
}


def _create_topodata_netcdf(
    topo_data: np.ndarray,
    topo_lons: np.ndarray,
    topo_lats: np.ndarray,
    output_path: Path,
) -> None:
    """
    Create topodata.nc for calnoro on a regular 0.5° grid (r720x360).

    Replaces convert_topodata.sh: remap to r720x360 and rename variable to OROMEA.
    """
    # Define target regular grid: 0.5° resolution, 720x360
    tgt_lons_1d = np.linspace(0.25, 359.75, 720)
    tgt_lats_1d = np.linspace(89.75, -89.75, 360)
    tgt_lons_2d, tgt_lats_2d = np.meshgrid(tgt_lons_1d, tgt_lats_1d)

    # Interpolate topography to regular grid
    src_pts = np.column_stack([topo_lons.ravel(), topo_lats.ravel()])
    src_vals = topo_data.ravel()

    # Remove NaN / masked values
    valid = np.isfinite(src_vals)
    if hasattr(src_vals, 'mask'):
        valid &= ~src_vals.mask
    src_pts = src_pts[valid]
    src_vals = src_vals[valid]

    tgt_pts = np.column_stack([tgt_lons_2d.ravel(), tgt_lats_2d.ravel()])
    topo_regular = griddata(src_pts, src_vals, tgt_pts, method='nearest', fill_value=0)
    topo_regular = topo_regular.reshape(360, 720)

    # Write NetCDF
    ds = Dataset(str(output_path), 'w', format='NETCDF4')
    ds.createDimension('lon', 720)
    ds.createDimension('lat', 360)

    lon_var = ds.createVariable('lon', 'f8', ('lon',))
    lon_var[:] = tgt_lons_1d
    lon_var.units = 'degrees_east'

    lat_var = ds.createVariable('lat', 'f8', ('lat',))
    lat_var[:] = tgt_lats_1d
    lat_var.units = 'degrees_north'

    oromea = ds.createVariable('OROMEA', 'f8', ('lat', 'lon'))
    oromea[:] = topo_regular
    oromea.long_name = 'mean orography'
    oromea.units = 'm'

    ds.close()
    print(f"  Created topodata.nc: {output_path}")


def _read_service_file(filepath: Path, nlat: int = 64, nlon: int = 128) -> Dict[int, np.ndarray]:
    """
    Read a Fortran SERVICE format file (as produced by calnoro).

    SERVICE format: for each field:
      - 8-int header (code, level, date, time, nlon, nlat, ?, ?)
      - nlon*nlat float32 values

    Returns dict mapping variable code → 2D array.

    Parameters
    ----------
    filepath : path to .srv file
    nlat, nlon : expected grid dimensions (T63 Gaussian = 64 lats × 128 lons)
    """
    fields = {}
    filesize = os.path.getsize(filepath)

    with open(filepath, 'rb') as f:
        while f.tell() < filesize:
            # Fortran record: 4-byte length prefix
            rec_len_bytes = f.read(4)
            if len(rec_len_bytes) < 4:
                break
            rec_len = struct.unpack('>i', rec_len_bytes)[0]

            # Read header (8 ints)
            header = struct.unpack(f'>{rec_len // 4}i', f.read(rec_len))
            f.read(4)  # trailing record length

            code = header[0]
            grid_nlon = header[4] if len(header) > 4 else nlon
            grid_nlat = header[5] if len(header) > 5 else nlat

            # Read data record
            rec_len_bytes = f.read(4)
            if len(rec_len_bytes) < 4:
                break
            rec_len = struct.unpack('>i', rec_len_bytes)[0]
            n_vals = rec_len // 4
            data = np.array(struct.unpack(f'>{n_vals}f', f.read(rec_len)))
            f.read(4)  # trailing record length

            fields[code] = data.reshape(grid_nlat, grid_nlon)

    return fields


def prepare_calnoro_input(
    config: OCPConfig,
    work_dir: Path,
) -> Path:
    """
    Prepare the topodata.nc file that calnoro expects as input.

    Returns path to topodata.nc.
    """
    paleo = config.paleo

    print("\n  Preparing calnoro input...")

    topo_file = paleo.get_reconstruction_file('topography_file')
    if not topo_file.exists():
        raise FileNotFoundError(f"Paleo topography not found: {topo_file}")

    topo_data, topo_lons, topo_lats = _read_netcdf_field(topo_file)

    topodata_path = work_dir / 'topodata.nc'
    _create_topodata_netcdf(topo_data, topo_lons, topo_lats, topodata_path)

    return topodata_path


def run_calnoro(
    config: OCPConfig,
    work_dir: Path,
    truncation: int = 63,
) -> Path:
    """
    Run the calnoro Fortran binary.

    Parameters
    ----------
    config : with paleo.calnoro_binary set
    work_dir : directory containing topodata.nc, output goes here
    truncation : spectral truncation (63 for T63)

    Returns path to the output sso_par_fil.srv file.
    """
    paleo = config.paleo

    if paleo.calnoro_binary is None or not paleo.calnoro_binary.exists():
        raise FileNotFoundError(
            f"Calnoro binary not found: {paleo.calnoro_binary}\n"
            "Compile calnoro from esm_tools: ./install_calnoro.sh levante"
        )

    print(f"\n  Running calnoro (truncation T{truncation})...")
    print(f"  Binary: {paleo.calnoro_binary}")
    print(f"  Working dir: {work_dir}")

    result = subprocess.run(
        [str(paleo.calnoro_binary)],
        input=f"{truncation}\n",
        capture_output=True,
        text=True,
        cwd=str(work_dir),
        timeout=600,
    )

    if result.returncode != 0:
        print(f"  calnoro stderr:\n{result.stderr}")
        raise RuntimeError(f"calnoro failed with exit code {result.returncode}")

    sso_file = work_dir / 'sso_par_fil.srv'
    if not sso_file.exists():
        raise FileNotFoundError(f"calnoro did not produce {sso_file}")

    print(f"  calnoro output: {sso_file}")
    return sso_file


def postprocess_calnoro_output(
    sso_file: Path,
    truncation: int = 63,
    verbose: bool = False,
) -> Dict[int, np.ndarray]:
    """
    Read and post-process calnoro SERVICE output.

    Returns dict mapping calnoro variable code → 2D array (lat-flipped).
    Calnoro for T63 uses a 64-latitude × 128-longitude Gaussian grid.
    """
    print("\n  Post-processing calnoro output...")

    # T63 Gaussian grid dimensions
    nlat = truncation + 1  # approximate
    nlon = nlat * 2

    fields = _read_service_file(sso_file, nlat=nlat, nlon=nlon)

    # Flip latitudes (calnoro outputs S→N, OpenIFS expects N→S)
    for code in fields:
        fields[code] = fields[code][::-1, :]
        if verbose:
            print(f"    Code {code}: shape {fields[code].shape}, "
                  f"range [{fields[code].min():.4f}, {fields[code].max():.4f}]")

    print(f"  Read {len(fields)} SSO fields from calnoro output")
    return fields


def apply_subgrid_orography(
    config: OCPConfig,
    grid: GaussianGrid,
    masks: Dict[str, np.ndarray],
    sso_fields: Dict[int, np.ndarray],
    truncation: int = 63,
) -> Dict[int, np.ndarray]:
    """
    Map calnoro SSO fields to OpenIFS GRIB codes and interpolate to Gaussian grid.

    Applies unit conversions:
    - OROSIG (slope): ×10
    - OROTHE (angle): degrees → radians

    Returns dict mapping OpenIFS GRIB paramId → values array on Gaussian grid.
    """
    verbose = config.options.verbose

    print("\n  Applying subgrid-scale orography to Gaussian grid...")

    tgt_lons = np.array(grid.lons_list)
    tgt_lats = grid.center_lats.flatten()
    land_mask = masks['land']

    # Build T63 Gaussian grid coordinates for the calnoro output
    nlat = truncation + 1
    nlon = nlat * 2
    # Simple regular grid approximation for the T63 output
    sso_lons_1d = np.linspace(0, 360 - 360/nlon, nlon)
    sso_lats_1d = np.linspace(90, -90, nlat)
    sso_lons_2d, sso_lats_2d = np.meshgrid(sso_lons_1d, sso_lats_1d)

    results = {}
    for calnoro_code, (oifs_code, name, scale) in SSO_MAPPING.items():
        if calnoro_code not in sso_fields:
            if verbose:
                print(f"    Skipping {name} (code {calnoro_code} not in calnoro output)")
            continue

        field_2d = sso_fields[calnoro_code]

        # Interpolate to Gaussian grid
        field_gauss = _interp_to_gaussian(
            field_2d, sso_lons_2d, sso_lats_2d,
            tgt_lons, tgt_lats,
            method='nearest', fill_value=0,
        )

        # Apply scale factor
        field_gauss *= scale

        # Mask ocean to 0
        field_gauss[~land_mask] = 0

        results[oifs_code] = field_gauss
        print(f"    {name} (code {oifs_code}): range [{field_gauss.min():.4f}, {field_gauss.max():.4f}]")

    return results


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def modify_paleo_subgrid_orography(
    config: OCPConfig,
    grid: GaussianGrid,
    masks: Dict[str, np.ndarray],
) -> Path:
    """
    Full subgrid-scale orography pipeline.

    1. Prepare calnoro input (topodata.nc)
    2. Run calnoro
    3. Post-process output
    4. Apply SSO fields to ICMGG

    Parameters
    ----------
    config : OCPConfig with paleo section
    grid : Gaussian grid
    masks : derivative masks from create_paleo_masks()

    Returns
    -------
    Path to the modified ICMGG file
    """
    paleo = config.paleo
    verbose = config.options.verbose
    icmgg_file = config.get_icmgg_output_file()

    print("\n" + "=" * 60)
    print(f" Paleo Subgrid-Scale Orography: {paleo.experiment_id}")
    print("=" * 60)

    # Use a temporary working directory for calnoro
    work_dir = config.output_paths.openifs_modified / f'calnoro_{paleo.experiment_id}'
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Prepare input
    prepare_calnoro_input(config, work_dir)

    # 2. Run calnoro
    sso_file = run_calnoro(config, work_dir)

    # 3. Post-process
    sso_fields = postprocess_calnoro_output(sso_file, verbose=verbose)

    # 4. Apply to Gaussian grid
    sso_replacements = apply_subgrid_orography(config, grid, masks, sso_fields)

    # 5. Write to ICMGG
    if sso_replacements:
        print(f"\n  Writing {len(sso_replacements)} SSO fields to {icmgg_file}...")
        _replace_grib_fields(icmgg_file, icmgg_file, sso_replacements, verbose=verbose)
        print(f"  Subgrid-scale orography applied: {icmgg_file}")
    else:
        print("  Warning: No SSO fields to write")

    return icmgg_file
