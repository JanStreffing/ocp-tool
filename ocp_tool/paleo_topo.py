"""
Paleo topography modification module for OCP-Tool.

Replaces mod_topo_oifs.sh: modifies the ICMSH spectral topography file
using the anomaly method (paleo topography − modern topography) applied
to the existing OpenIFS spectral orography.

Uses CDO (via Python bindings) for spectral ↔ gridpoint transforms
(sp2gpl / gp2spl), eccodes for GRIB I/O, and scipy for interpolation.
"""

import os
import tempfile
from pathlib import Path
from shutil import copy2
from typing import Dict, Tuple

import numpy as np
from scipy.interpolate import griddata
import eccodes

from .config import OCPConfig
from .gaussian_grids import GaussianGrid
from .paleo_input import (
    _read_netcdf_field,
    _interp_to_gaussian,
    _distance_weighted_fill,
)

# GRIB code of the spectral orography (geopotential at the surface)
Z_PARAM_ID = 129


# ── CDO helpers for spectral ↔ gridpoint ────────────────────────

def _sp2gp(spectral_grib: str, gridpoint_nc: str) -> None:
    """Convert spectral GRIB → Gaussian gridpoint NetCDF using CDO sp2gpl."""
    from cdo import Cdo
    cdo = Cdo()
    cdo.sp2gpl(input=spectral_grib, output=gridpoint_nc, options='-f nc4 --eccodes')


def _gp2sp(gridpoint_nc: str, spectral_grib: str, grib_opts: str = '-f grb2 --eccodes') -> None:
    """Convert Gaussian gridpoint NetCDF → spectral GRIB using CDO gp2spl."""
    from cdo import Cdo
    cdo = Cdo()
    cdo.gp2spl(input=gridpoint_nc, output=spectral_grib, options=grib_opts)


def _extract_z_grib(icmsh_file: Path, output_file: str) -> None:
    """Extract the 'z' (geopotential / orography) message from an ICMSH file."""
    import subprocess
    subprocess.run(
        ['grib_copy', '-w', 'shortName=z', str(icmsh_file), output_file],
        check=True,
    )


def _replace_z_in_icmsh(icmsh_input: Path, z_spectral_grib: str, icmsh_output: Path) -> None:
    """
    Write the modified orography back into ICMSH, reusing the original message.

    CDO's gp2spl emits a fresh GRIB2 message carrying no parameter identity —
    ecCodes reports paramId 0. Merging that in place of the original leaves the
    file without code 129, and OpenIFS aborts at the very first read of CNT3::

        SPECTRAL 2D FIELD MISSING:  129
        ABOR1 : IOSTREAM_MIX:SPEC_IN - MISSING FIELD

    Copying only the spectral coefficients into the original message keeps its
    edition, paramId, truncation and packing, so the file stays readable.
    """
    with open(z_spectral_grib, 'rb') as f:
        gid = eccodes.codes_grib_new_from_file(f)
        if gid is None:
            raise RuntimeError(f"No GRIB message in {z_spectral_grib}")
        try:
            new_values = eccodes.codes_get_array(gid, 'values')
        finally:
            eccodes.codes_release(gid)

    replaced = 0
    with open(icmsh_input, 'rb') as fin, open(icmsh_output, 'wb') as fout:
        while True:
            gid = eccodes.codes_grib_new_from_file(fin)
            if gid is None:
                break
            try:
                if int(eccodes.codes_get(gid, 'paramId')) == Z_PARAM_ID:
                    n_orig = int(eccodes.codes_get(gid, 'numberOfValues'))
                    if n_orig != new_values.size:
                        raise RuntimeError(
                            f"Spectral truncation mismatch writing orography: "
                            f"original message holds {n_orig} coefficients, "
                            f"gp2spl produced {new_values.size}"
                        )
                    eccodes.codes_set_array(gid, 'values', new_values)
                    replaced += 1
                eccodes.codes_write(gid, fout)
            finally:
                eccodes.codes_release(gid)

    if replaced != 1:
        raise RuntimeError(
            f"Expected exactly one orography message (paramId {Z_PARAM_ID}) in "
            f"{icmsh_input}, found {replaced}"
        )


# ── Topography reading helpers ──────────────────────────────────

def _read_z_gridpoint(icmsh_file: Path, verbose: bool = False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read spectral 'z' from ICMSH, convert to gridpoint via CDO sp2gpl.

    Returns (topo_m, lats_1d, lons_1d) on a Gaussian grid.
    """
    import xarray as xr

    with tempfile.TemporaryDirectory(prefix='ocp_topo_') as tmpdir:
        z_grib = os.path.join(tmpdir, 'z_spectral.grib')
        z_nc = os.path.join(tmpdir, 'z_gridpoint.nc')

        _extract_z_grib(icmsh_file, z_grib)
        _sp2gp(z_grib, z_nc)

        ds = xr.open_dataset(z_nc)
        # Variable may be called 'z' or 'orog' depending on CDO version
        for vname in ['z', 'orog', 'Z']:
            if vname in ds.data_vars:
                data = ds[vname].values.squeeze()
                break
        else:
            # Fallback: take the first data variable
            vname = list(ds.data_vars)[0]
            data = ds[vname].values.squeeze()

        lats = ds['lat'].values
        lons = ds['lon'].values

        ds.close()

    topo_m = data / 9.80665  # geopotential → metres
    if verbose:
        print(f"  sp2gp: {vname} → {data.shape}, "
              f"topo range [{topo_m.min():.1f}, {topo_m.max():.1f}] m")
    return topo_m, lats, lons


# ── Main entry point ────────────────────────────────────────────

def modify_topography(
    config: OCPConfig,
    grid: GaussianGrid,
    masks: Dict[str, np.ndarray],
) -> Path:
    """
    Modify ICMSH topography using the anomaly method.

    Steps:
    1. Read modern + paleo topography from reconstruction NetCDF files
    2. Clip below sea level to 0
    3. Distance-weighted fill for land mask gaps
    4. Compute anomaly (paleo − modern) on Gaussian grid
    5. Read existing OpenIFS spectral topo via CDO sp2gpl → gridpoint
    6. Add anomaly on the gridpoint grid, mask by land, clip below 0
    7. Convert back to spectral via CDO gp2spl, replace 'z' in ICMSH

    Returns path to the modified ICMSH file.
    """
    import xarray as xr

    paleo = config.paleo
    verbose = config.options.verbose

    print("\n" + "=" * 60)
    print(f" Paleo Topography Modification: {paleo.experiment_id}")
    print("=" * 60)

    tgt_lons = np.array(grid.lons_list)
    tgt_lats = grid.center_lats.flatten()
    land_mask = masks['land']

    # --- 1. Read reconstruction topographies ---
    modern_topo_file = paleo.get_modern_file('modern_topo_file')
    paleo_topo_file = paleo.get_reconstruction_file('topography_file')

    print(f"  Modern topo: {modern_topo_file}")
    print(f"  Paleo topo:  {paleo_topo_file}")

    modern_data, m_lons, m_lats = _read_netcdf_field(modern_topo_file)
    paleo_data, p_lons, p_lats = _read_netcdf_field(paleo_topo_file)

    # --- 2. Clip below sea level ---
    modern_data = np.where(modern_data < 0, 0, modern_data)
    paleo_data = np.where(paleo_data < 0, 0, paleo_data)

    # --- 3. Interpolate to Gaussian grid ---
    modern_gauss = _interp_to_gaussian(modern_data, m_lons, m_lats,
                                       tgt_lons, tgt_lats, method='nearest', fill_value=0)
    paleo_gauss = _interp_to_gaussian(paleo_data, p_lons, p_lats,
                                      tgt_lons, tgt_lats, method='nearest', fill_value=0)

    # Distance-weighted fill for land points missing paleo topo
    paleo_gauss[~np.isfinite(paleo_gauss)] = 0
    paleo_on_land = np.where(land_mask, paleo_gauss, 0)
    paleo_on_land[land_mask & (paleo_on_land == 0)] = 200  # default elevation
    paleo_filled = _distance_weighted_fill(paleo_on_land, land_mask & (paleo_on_land > 0),
                                           tgt_lons, tgt_lats)
    paleo_filled[~land_mask] = 0

    # --- 4. Compute anomaly ---
    topo_anom = paleo_filled - modern_gauss
    topo_anom[~land_mask] = 0

    print(f"  Topo anomaly range: [{topo_anom.min():.1f}, {topo_anom.max():.1f}] m")

    # --- 5. Read existing spectral topo → gridpoint via CDO ---
    icmsh_input = config.get_icmsh_input_file()
    print(f"  Reading spectral topo from: {icmsh_input}")

    oifs_topo, gp_lats, gp_lons = _read_z_gridpoint(icmsh_input, verbose=verbose)
    # oifs_topo is 2D (lat, lon)

    gp_lons2d, gp_lats2d = np.meshgrid(gp_lons, gp_lats)

    # Interpolate anomaly from Gaussian grid to the gridpoint grid
    anom_pts = np.column_stack([tgt_lons, tgt_lats])
    gp_pts = np.column_stack([gp_lons2d.ravel(), gp_lats2d.ravel()])
    anom_on_gp = griddata(anom_pts, topo_anom, gp_pts, method='nearest', fill_value=0)
    anom_on_gp = anom_on_gp.reshape(oifs_topo.shape)

    # --- 6. Add anomaly, clip ---
    new_topo = oifs_topo + anom_on_gp
    new_topo = np.where(new_topo < 0, 0, new_topo)

    # Mask ocean (where the paleo LSM says ocean)
    land_on_gp = griddata(anom_pts, land_mask.astype(float), gp_pts,
                          method='nearest', fill_value=0).reshape(oifs_topo.shape)
    new_topo[land_on_gp < 0.5] = 0

    print(f"  New topo range: [{new_topo.min():.1f}, {new_topo.max():.1f}] m")

    # --- 7. Convert modified gridpoint → spectral, write ICMSH ---
    icmsh_output = config.get_icmsh_output_file()
    print(f"  Writing modified ICMSH: {icmsh_output}")

    new_geopotential = new_topo * 9.80665

    with tempfile.TemporaryDirectory(prefix='ocp_topo_') as tmpdir:
        # Write modified gridpoint field as NetCDF
        mod_nc = os.path.join(tmpdir, 'z_modified.nc')

        # Read the sp2gp output to get metadata, then overwrite data
        z_grib_tmp = os.path.join(tmpdir, 'z_orig.grib')
        z_nc_tmp = os.path.join(tmpdir, 'z_orig.nc')
        _extract_z_grib(icmsh_input, z_grib_tmp)
        _sp2gp(z_grib_tmp, z_nc_tmp)

        ds = xr.open_dataset(z_nc_tmp)
        # Find the variable name
        for vname in ['z', 'orog', 'Z']:
            if vname in ds.data_vars:
                break
        else:
            vname = list(ds.data_vars)[0]

        ds_mod = ds.copy(deep=True)
        ds_mod[vname].values[:] = new_geopotential.reshape(ds_mod[vname].shape)
        ds_mod.to_netcdf(mod_nc)
        ds.close()
        ds_mod.close()

        # gp2sp: gridpoint NetCDF → spectral GRIB
        z_new_grib = os.path.join(tmpdir, 'z_new_spectral.grib')
        _gp2sp(mod_nc, z_new_grib)

        # Replace 'z' in the original ICMSH
        _replace_z_in_icmsh(icmsh_input, z_new_grib, icmsh_output)

    print(f"  Topography modification complete: {icmsh_output}")
    return icmsh_output
