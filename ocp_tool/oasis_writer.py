"""
OASIS file generation module for OCP-Tool.
Handles writing grids.nc, masks.nc, and areas.nc files.
"""

from pathlib import Path
from typing import Optional
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import xarray as xr
from netCDF4 import Dataset
from scipy.spatial import cKDTree

from .config import OCPConfig
from .gaussian_grids import GaussianGrid
from .lsm import LSMData


def write_oasis_grid_files(
    config: OCPConfig,
    grid: GaussianGrid,
    lsm_data: LSMData,
    resolution: int,
    parallel: bool = True
) -> None:
    """
    Write OASIS grids, masks, and areas files.
    
    Args:
        config: OCP configuration
        grid: Gaussian grid data
        lsm_data: Land-sea mask data
        resolution: Truncation number
        parallel: Use parallel file writing
    """
    nn = grid.nn
    if len(str(nn)) > 4:
        nn = int(str(nn)[:-1])
    
    # Generate LPJG OASIS grid name
    if config.atmosphere.truncation_type == 'cubic-octahedral':
        lpjg_oasis_name = f'TCO{nn-1}-land'
    else:
        lpjg_oasis_name = f'TL{nn*2-1}-land'
    
    # Grid names to process
    grid_names = [
        f'{s}{nn:03d}' if s in ('A', 'L', 'R') else s
        for s in ('A', 'L', 'R', lpjg_oasis_name)
    ]
    
    file_types = ['grids', 'areas', 'masks']
    
    if parallel:
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = [
                executor.submit(
                    _write_single_oasis_file,
                    config, grid, lsm_data, resolution, nn, grid_names, file_type
                )
                for file_type in file_types
            ]
            for future in futures:
                future.result()
    else:
        for file_type in file_types:
            _write_single_oasis_file(
                config, grid, lsm_data, resolution, nn, grid_names, file_type
            )
    
    # Copy runoff files into OASIS files (skip for AMIP mode - no ocean coupling)
    if config.ocean.grid_name != 'AMIP':
        _append_runoff_to_oasis_files(config, file_types)
        if config.ocean.write_oasis_grid:
            _append_fesom_grid_to_oasis_files(config)
        else:
            print(' Skipping FESOM grid export (FESOM writes feom at runtime)')
    else:
        # Add AMIP forcing-reader grid for SST/SIC forcing
        _append_amip_grid_to_oasis_files(config, file_types)

    # Add ice-sheet grids, off unless ice.enabled is set in the config
    if config.ice is not None:
        _append_pism_grid_to_oasis_files(config, file_types)
        _write_ismm_restart(config)


def _write_single_oasis_file(
    config: OCPConfig,
    grid: GaussianGrid,
    lsm_data: LSMData,
    resolution: int,
    nn: int,
    grid_names: list,
    file_type: str
) -> None:
    """Write a single OASIS file (grids, masks, or areas)."""
    output_file = config.output_paths.oasis / f'{file_type}.nc'
    print(f'Writing file: {output_file}')
    
    nc = Dataset(str(output_file), 'w', clobber=True)
    
    for grids_name in grid_names:
        print(f' Write lons, lats, corner points for grid: {grids_name} (T{resolution})')
        
        xname = f'x_{grids_name}'
        yname = f'y_{grids_name}'
        lonname = f'{grids_name}.lon'
        latname = f'{grids_name}.lat'
        
        # Only create dimensions if they don't exist
        if xname not in nc.dimensions:
            nc.createDimension(xname, grid.center_lons.shape[1])
        if yname not in nc.dimensions:
            nc.createDimension(yname, 1)
        
        # Skip if variables already exist (avoid duplicates)
        if lonname in nc.variables:
            continue
        
        id_lon = nc.createVariable(lonname, 'float64', (yname, xname))
        id_lat = nc.createVariable(latname, 'float64', (yname, xname))
        
        id_lon.units = 'degrees_east'
        id_lon.standard_name = 'Longitude'
        id_lat.units = 'degrees_north'
        id_lat.standard_name = 'Latitude'
        
        if file_type == 'grids':
            crnname = f'crn_{grids_name}'
            cloname = f'{grids_name}.clo'
            claname = f'{grids_name}.cla'
            if crnname not in nc.dimensions:
                nc.createDimension(crnname, 4)
            id_clo = nc.createVariable(cloname, 'float64', (crnname, yname, xname))
            id_cla = nc.createVariable(claname, 'float64', (crnname, yname, xname))
        
        elif file_type == 'masks':
            mskname = f'{grids_name}.msk'
            id_msk = nc.createVariable(mskname, 'int32', (yname, xname))
            id_msk.coordinates = f'{grids_name}.lat {grids_name}.lon'
            id_msk.valid_min = 0.
            id_msk.valid_max = 1
        
        elif file_type == 'areas':
            areaname = f'{grids_name}.srf'
            id_area = nc.createVariable(areaname, 'float64', (yname, xname))
            id_area.coordinates = f'{grids_name}.lat {grids_name}.lon'
        
        # Write coordinate data
        id_lon[:, :] = grid.center_lons[:, :]
        id_lat[:, :] = grid.center_lats[:, :]
        id_lon.valid_min = float(grid.center_lons.min())
        id_lon.valid_max = float(grid.center_lons.max())
        id_lat.valid_min = float(grid.center_lats.min())
        id_lat.valid_max = float(grid.center_lats.max())
        
        if file_type == 'grids':
            id_clo[:, :, :] = grid.corner_lons[:, :, :]
            id_cla[:, :, :] = grid.corner_lats[:, :, :]
            id_clo.valid_min = float(grid.corner_lons.min())
            id_clo.valid_max = float(grid.corner_lons.max())
            id_cla.valid_min = float(grid.corner_lats.min())
            id_cla.valid_max = float(grid.corner_lats.max())
        
        elif file_type == 'masks':
            if grids_name.startswith('A'):
                id_msk[:, :] = np.round(lsm_data.lsm_binary_atm[:, :])
            elif grids_name.startswith('L'):
                id_msk[:, :] = np.round(lsm_data.lsm_binary_land[:, :])
            elif grids_name.startswith('R'):
                if config.ocean.grid_name == 'ORCA05':
                    id_msk[:, :] = np.abs(np.round(lsm_data.lsm_binary_runoff[:, :] - 1))
                else:
                    id_msk[:, :] = np.abs(np.round(lsm_data.lsm_binary_atm[:, :] - 1))
            elif '-land' in grids_name:
                id_msk[:, :] = np.abs(np.round(lsm_data.lsm_binary_atm[:, :] - 1))
            else:
                raise RuntimeError(f'Unexpected grid name: {grids_name}')
        
        elif file_type == 'areas':
            id_area[:, :] = grid.cell_areas[:, :]
            id_area.valid_min = float(grid.cell_areas.min())
            id_area.valid_max = float(grid.cell_areas.max())
    
    nc.close()
    print(f' Wrote {output_file}')


def _append_runoff_to_oasis_files(config: OCPConfig, file_types: list) -> None:
    """Append runoff mapper grids/areas/masks to OASIS files."""
    for file_type in file_types:
        oasis_file = config.output_paths.oasis / f'{file_type}.nc'
        runoff_file = config.input_paths.runoff_default / f'runoff_{file_type}.nc'
        
        print(f' Adding {runoff_file} to oasis files')
        
        nc = Dataset(str(oasis_file), 'a')
        rnffile = Dataset(str(runoff_file), 'r')
        
        nc.setncatts(rnffile.__dict__)
        
        for name, dimension in rnffile.dimensions.items():
            if name not in nc.dimensions:
                nc.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
        
        for name, variable in rnffile.variables.items():
            if name not in nc.variables:
                var_out = nc.createVariable(name, variable.datatype, variable.dimensions)
                var_out.setncatts({k: variable.getncattr(k) for k in variable.ncattrs()})
                var_out[:] = variable[:]
        
        rnffile.close()
        nc.close()
        
        print(f"\n {'='*50} \n")


def _append_fesom_grid_to_oasis_files(config: OCPConfig) -> None:
    """Append FESOM ocean grid to OASIS files using pyfesom2."""
    try:
        import pyfesom2 as pf
    except ImportError:
        print("Warning: pyfesom2 not available, skipping FESOM grid export")
        return
    
    if config.ocean.mesh_file is None:
        print("Warning: No FESOM mesh file specified, skipping FESOM grid export")
        return
    
    mesh_path = Path(config.ocean.mesh_file).parent
    print(f" Loading FESOM mesh from: {mesh_path}")
    
    try:
        mesh = pf.load_mesh(str(mesh_path))
        print(f"  FESOM mesh loaded: {mesh.n2d} nodes, {mesh.e2d} elements")
        
        # Write FESOM grid to OASIS files with 'feom' prefix
        output_dir = config.output_paths.oasis
        print(f" Adding FESOM grid to OASIS files in: {output_dir}")
        
        pf.write_fesom_oasis_files(
            mesh=mesh,
            output_dir=str(output_dir),
            prefix='feom',
            overwrite=True
        )
        
        print(" FESOM grid successfully added to OASIS files")
        print(f"\n {'='*50} \n")
        
    except Exception as e:
        print(f"Warning: Failed to add FESOM grid to OASIS files: {e}")
        import traceback
        traceback.print_exc()


def _write_ismm_restart(config: OCPConfig, restart_name: str = "rstism.nc") -> None:
    """Seed the ISM-mapper OASIS restart with the ice mask.

    OASIS restarts hold the SOURCE field on the SOURCE grid and remap at read
    time, so only <prefix>_plit is needed here; the atmosphere fluxes are safely
    zero-filled by NNOREST. A mask is not: with LAG the receiver takes its first
    value from this file, and an empty mask strips the glacier tiles from the
    atmosphere for one coupling interval.
    """
    out = config.output_paths.oasis / restart_name
    print(f' Writing ISM-mapper restart {out}')

    wrote = []
    for region in config.ice.regions:
        try:
            with Dataset(str(region.grid_file)) as src:
                if "thk" not in src.variables:
                    print(f'Warning: no thk in {region.grid_file}, skipping {region.name}')
                    continue
                thk = np.array(src.variables["thk"][:])
        except Exception as e:
            print(f'Warning: could not read {region.grid_file}: {e}')
            continue

        if thk.ndim == 3:
            thk = thk[-1]
        mask = np.where(thk > 0.0, 1.0, 0.0)

        ny, nx = mask.shape
        vn = f"{region.prefix}_plit"
        mode = "a" if out.exists() else "w"
        nc = Dataset(str(out), mode, format="NETCDF3_CLASSIC")
        try:
            yname, xname = f"ny_{region.name}", f"nx_{region.name}"
            if yname not in nc.dimensions:
                nc.createDimension(yname, ny)
            if xname not in nc.dimensions:
                nc.createDimension(xname, nx)
            if vn not in nc.variables:
                v = nc.createVariable(vn, "f8", (yname, xname))
                v.units = "1"
                v.long_name = "ice sheet mask, 1 = ice"
                v[:] = mask
                wrote.append(f"{vn} ({nx}x{ny}, {int(mask.sum())} ice cells)")
        finally:
            nc.close()

    print(f"  {', '.join(wrote) if wrote else '(nothing written)'}")
    print(f"\n {'='*50} \n")


def _append_pism_grid_to_oasis_files(config: OCPConfig, file_types: list) -> None:
    """Append each configured PISM ice-sheet domain to the OASIS files."""
    from .grids import PISM

    for region in config.ice.regions:
        print(f' Adding PISM ice grid {region.name} from {region.grid_file}')
        try:
            grid = PISM(str(region.grid_file), name=region.name)
            lats = grid.cell_latitudes()
            lons = grid.cell_longitudes()
            corners = grid.cell_corners()
            areas = grid.cell_areas()
            masks = grid.cell_masks()
        except Exception as e:
            print(f'Warning: failed to build PISM grid {region.name}: {e}')
            import traceback
            traceback.print_exc()
            continue

        ny, nx = lats.shape
        name = region.name
        xname, yname, crnname = f'x_{name}', f'y_{name}', f'crn_{name}'

        for file_type in file_types:
            oasis_file = config.output_paths.oasis / f'{file_type}.nc'
            nc = Dataset(str(oasis_file), 'a')

            if xname not in nc.dimensions:
                nc.createDimension(xname, nx)
            if yname not in nc.dimensions:
                nc.createDimension(yname, ny)

            if file_type == 'grids':
                if crnname not in nc.dimensions:
                    nc.createDimension(crnname, 4)
                for vn, data, units, sname in (
                    (f'{name}.lon', lons, 'degrees_east', 'Longitude'),
                    (f'{name}.lat', lats, 'degrees_north', 'Latitude'),
                ):
                    if vn not in nc.variables:
                        v = nc.createVariable(vn, 'float64', (yname, xname))
                        v.units = units
                        v.standard_name = sname
                        v[:] = data
                for vn, data, units, sname in (
                    (f'{name}.clo', corners[1], 'degrees_east', 'Longitude corners'),
                    (f'{name}.cla', corners[0], 'degrees_north', 'Latitude corners'),
                ):
                    if vn not in nc.variables:
                        v = nc.createVariable(vn, 'float64', (crnname, yname, xname))
                        v.units = units
                        v.standard_name = sname
                        v[:] = data

            elif file_type == 'areas':
                vn = f'{name}.srf'
                if vn not in nc.variables:
                    v = nc.createVariable(vn, 'float64', (yname, xname))
                    v.units = 'm2'
                    v.standard_name = 'Grid cell area'
                    v[:] = areas

            elif file_type == 'masks':
                vn = f'{name}.msk'
                if vn not in nc.variables:
                    v = nc.createVariable(vn, 'int32', (yname, xname))
                    v.units = '1'
                    v.standard_name = 'Grid mask (0=active, 1=masked)'
                    v[:] = masks

            nc.close()

        print(f'  {name}: {nx}x{ny}, total area {areas.sum()/1e12:.2f} x10^6 km2')

    print(' PISM ice grid(s) added successfully')
    print(f"\n {'='*50} \n")


def _append_amip_grid_to_oasis_files(config: OCPConfig, file_types: list) -> None:
    """Append AMIP forcing-reader grid (360x180 regular) to OASIS files for SST/SIC forcing."""
    print(' Adding AMIP forcing-reader grid (360x180) to OASIS files')
    
    # Create regular lat-lon grid: 180 lats x 360 lons
    nlats = 180
    nlons = 360
    
    # Cell centers: lon = 0.5 to 359.5, lat = -89.5 to 89.5
    lons_1d = np.arange(0.5, 360.0, 1.0)  # 0.5, 1.5, 2.5, ..., 359.5
    lats_1d = np.arange(-89.5, 90.0, 1.0)  # -89.5, -88.5, ..., 89.5
    lons, lats = np.meshgrid(lons_1d, lats_1d)
    
    # Cell corners (4 corners per cell)
    lons_corners = np.arange(0.0, 360.0, 1.0)  # 0, 1, 2, ..., 359
    lats_corners = np.arange(-90.0, 91.0, 1.0)  # -90, -89, ..., 90
    
    # Build corner arrays: shape (nlats, nlons, 4)
    corners_lon = np.zeros((nlats, nlons, 4))
    corners_lat = np.zeros((nlats, nlons, 4))
    
    for i in range(nlats):
        for j in range(nlons):
            # Corners in counter-clockwise order
            # Use 360° for the right edge of the last cell, not 0°
            lon_right = lons_corners[j+1] if j < nlons-1 else 360.0
            corners_lon[i, j, 0] = lons_corners[j]       # bottom-left
            corners_lon[i, j, 1] = lon_right             # bottom-right
            corners_lon[i, j, 2] = lon_right             # top-right
            corners_lon[i, j, 3] = lons_corners[j]       # top-left
            
            corners_lat[i, j, 0] = lats_corners[i]       # bottom-left
            corners_lat[i, j, 1] = lats_corners[i]       # bottom-right
            corners_lat[i, j, 2] = lats_corners[i+1]     # top-right
            corners_lat[i, j, 3] = lats_corners[i+1]     # top-left
    
    # Calculate cell areas (in m^2)
    earth_radius = 6371000.0  # meters
    cell_area = np.zeros((nlats, nlons))
    for i in range(nlats):
        lat_min = lats_corners[i] * np.pi / 180.0
        lat_max = lats_corners[i+1] * np.pi / 180.0
        dlon_rad = 1.0 * np.pi / 180.0
        # Area = R^2 * dlon * (sin(lat_max) - sin(lat_min))
        cell_area[i, :] = earth_radius**2 * dlon_rad * (np.sin(lat_max) - np.sin(lat_min))
    
    # Write to each file type
    for file_type in file_types:
        oasis_file = config.output_paths.oasis / f'{file_type}.nc'
        nc = Dataset(str(oasis_file), 'a')
        
        grid_name = 'AMIP'
        xname = f'x_{grid_name}'
        yname = f'y_{grid_name}'
        
        # Create dimensions
        if xname not in nc.dimensions:
            nc.createDimension(xname, nlons)
        if yname not in nc.dimensions:
            nc.createDimension(yname, nlats)
        
        if file_type == 'grids':
            # Write lon/lat centers
            lonname = f'{grid_name}.lon'
            latname = f'{grid_name}.lat'
            
            if lonname not in nc.variables:
                var_lon = nc.createVariable(lonname, 'float64', (yname, xname))
                var_lon.units = 'degrees_east'
                var_lon.standard_name = 'Longitude'
                var_lon[:] = lons
            
            if latname not in nc.variables:
                var_lat = nc.createVariable(latname, 'float64', (yname, xname))
                var_lat.units = 'degrees_north'
                var_lat.standard_name = 'Latitude'
                var_lat[:] = lats
            
            # Write corners
            crnname = f'crn_{grid_name}'
            cloname = f'{grid_name}.clo'
            claname = f'{grid_name}.cla'
            
            if crnname not in nc.dimensions:
                nc.createDimension(crnname, 4)
            
            if cloname not in nc.variables:
                var_clo = nc.createVariable(cloname, 'float64', (crnname, yname, xname))
                var_clo.units = 'degrees_east'
                var_clo.standard_name = 'Longitude corners'
                var_clo[:] = np.transpose(corners_lon, (2, 0, 1))  # (nlats, nlons, 4) -> (4, nlats, nlons)
            
            if claname not in nc.variables:
                var_cla = nc.createVariable(claname, 'float64', (crnname, yname, xname))
                var_cla.units = 'degrees_north'
                var_cla.standard_name = 'Latitude corners'
                var_cla[:] = np.transpose(corners_lat, (2, 0, 1))  # (nlats, nlons, 4) -> (4, nlats, nlons)
        
        elif file_type == 'areas':
            areaname = f'{grid_name}.srf'
            if areaname not in nc.variables:
                var_area = nc.createVariable(areaname, 'float64', (yname, xname))
                var_area.units = 'm2'
                var_area.standard_name = 'Grid cell area'
                var_area[:] = cell_area
        
        elif file_type == 'masks':
            maskname = f'{grid_name}.msk'
            if maskname not in nc.variables:
                var_mask = nc.createVariable(maskname, 'int32', (yname, xname))
                var_mask.units = '1'
                var_mask.standard_name = 'Grid mask (0=land, 1=ocean)'
                # All land (0) for AMIP forcing reader
                var_mask[:] = np.zeros((nlats, nlons), dtype=np.int32)
        
        nc.close()
    
    print(' AMIP forcing-reader grid added successfully')
    print(f"\n {'='*50} \n")


def interpolate_vegin_data(
    config: OCPConfig,
    grid: GaussianGrid
) -> None:
    """
    Interpolate vegetation input data to target grid using KD-tree nearest neighbor.
    
    Also interpolates CO2 restart files and nitrogen deposition files.
    
    Args:
        config: OCP configuration
        grid: Target Gaussian grid
    """
    print("Read input files")
    vegin_grid = xr.open_dataset(config.input_paths.lpj_guess / 'vegin_grid.nc')
    vegin = xr.open_dataset(config.input_paths.lpj_guess / 'vegin.nc')
    
    # Extract source grid
    print("Extract source grid and data")
    src_lon = np.squeeze(vegin_grid["A128.lon"].values)
    src_lat = np.squeeze(vegin_grid["A128.lat"].values)
    
    # Target grid
    print("Prepare target grid")
    target_lon = np.squeeze(grid.center_lons)
    target_lat = np.squeeze(grid.center_lats)
    
    # Build KD-tree
    print("Create a KDTree for nearest-neighbor interpolation")
    src_points = np.column_stack((src_lon, src_lat))
    target_points = np.column_stack((target_lon, target_lat))
    tree = cKDTree(src_points)
    _, idx = tree.query(target_points)
    
    # Interpolate vegin
    print("Creating interpolated dataset")
    new_ds = _interpolate_dataset(vegin, src_lon, target_lon, idx)
    
    # Write output
    vegin_out_path = config.output_paths.oasis / 'vegin.nc'
    _write_netcdf3_file(new_ds, vegin_out_path, "vegin.nc")
    print(f"Interpolated data successfully saved to {vegin_out_path}")
    
    vegin.close()
    vegin_grid.close()
    
    # Interpolate CO2 restart files
    _interpolate_co2_restarts(config, src_lon, target_lon, idx)
    
    # Interpolate nitrogen deposition files
    _interpolate_ndep_files(config, src_lon, target_lon, target_lat, idx)


def _interpolate_dataset(
    src_ds: xr.Dataset,
    src_lon: np.ndarray,
    target_lon: np.ndarray,
    idx: np.ndarray
) -> xr.Dataset:
    """Interpolate xarray dataset using pre-computed KD-tree indices."""
    new_ds = xr.Dataset()
    new_ds.attrs.update(src_ds.attrs)
    
    for var in src_ds.data_vars:
        try:
            src_data = src_ds[var].values.squeeze()
            original_dims = src_ds[var].dims
            
            if len(original_dims) == 1 and original_dims[0].endswith("_ncnt"):
                continue
            
            if src_data.size == 1:
                interpolated = np.full(target_lon.shape, src_data.item())
                new_ds[var] = (["y", "x"], interpolated[np.newaxis, :])
                continue
            
            if src_data.size != len(src_lon):
                continue
            
            interpolated = src_data.ravel()[idx].reshape(target_lon.shape)
            new_ds[var] = (["y", "x"], interpolated[np.newaxis, :])
            
            if hasattr(src_ds[var], 'attrs'):
                new_ds[var].attrs.update(src_ds[var].attrs)
            
            print(f"Updated variable: {var}")
            
        except Exception as e:
            print(f"Error interpolating variable {var}: {e}")
    
    return new_ds


def _interpolate_co2_restarts(
    config: OCPConfig,
    src_lon: np.ndarray,
    target_lon: np.ndarray,
    idx: np.ndarray
) -> None:
    """Interpolate CO2 restart files."""
    import os
    
    for fname in ['rst_co2_ao.nc', 'rst_co2_av.nc']:
        input_file = config.input_paths.lpj_guess / fname
        if not input_file.exists():
            print(f"\nFile {input_file} not found, skipping.")
            continue
        
        print(f"\nInterpolating {fname}...")
        src_ds = xr.open_dataset(input_file)
        
        # Special handling for rst_co2_av.nc
        if fname == 'rst_co2_av.nc':
            if 'GUE_CFIR' in src_ds.data_vars and 'GUE_CNAT' in src_ds.data_vars:
                print("  Adding fire flux (GUE_CFIR) to natural flux (GUE_CNAT)...")
                gue_cnat_with_fire = src_ds['GUE_CNAT'].values + src_ds['GUE_CFIR'].values
                src_ds['GUE_CNAT'] = (src_ds['GUE_CNAT'].dims, gue_cnat_with_fire)
                src_ds = src_ds.drop_vars('GUE_CFIR')
            
            if 'GUE_CNPP' not in src_ds.data_vars and 'GUE_CNAT' in src_ds.data_vars:
                print("  Adding NPP field (GUE_CNPP)...")
                gue_cnpp_data = np.zeros_like(src_ds['GUE_CNAT'].values)
                src_ds['GUE_CNPP'] = (src_ds['GUE_CNAT'].dims, gue_cnpp_data)
                src_ds['GUE_CNPP'].attrs['long_name'] = 'Net Primary Production flux'
                src_ds['GUE_CNPP'].attrs['units'] = 'kg m-2 s-1'
        
        # Interpolate
        co2_ds = xr.Dataset()
        co2_ds.attrs.update(src_ds.attrs)
        
        for var in src_ds.data_vars:
            src_data = src_ds[var].values.squeeze()
            if src_data.size == len(src_lon):
                interpolated = src_data.ravel()[idx].reshape(target_lon.shape)
                co2_ds[var] = (["y", "x"], interpolated[np.newaxis, :])
                co2_ds[var].attrs.update(src_ds[var].attrs)
        
        out_path = config.output_paths.oasis / fname
        _write_netcdf3_file(co2_ds, out_path, fname)
        print(f"Saved {fname} to {out_path}")
        
        src_ds.close()


def _interpolate_ndep_files(
    config: OCPConfig,
    src_lon: np.ndarray,
    target_lon: np.ndarray,
    target_lat: np.ndarray,
    idx: np.ndarray
) -> None:
    """Interpolate nitrogen deposition files."""
    truncation_type = config.atmosphere.truncation_type
    resolution = config.atmosphere.resolution_list[0]  # Use first resolution
    ifs_grid = f"TCO{resolution}" if truncation_type == "cubic-octahedral" else f"TL{resolution}"
    
    ndep_files = [
        'drynhx_TL255_hist_d1.nc', 'wetnhx_TL255_hist_d1.nc',
        'drynoy_TL255_hist_d1.nc', 'wetnoy_TL255_hist_d1.nc'
    ]
    
    for ndep_fname in ndep_files:
        ndep_input = config.input_paths.lpj_guess / ndep_fname
        if not ndep_input.exists():
            print(f"\nNitrogen deposition file {ndep_input} not found, skipping.")
            continue
        
        print(f"\nInterpolating nitrogen deposition file {ndep_fname}...")
        ndep_ds = xr.open_dataset(ndep_input)
        
        data_vars = [v for v in ndep_ds.data_vars if v not in ['time_bnds', 'lon_bnds', 'lat_bnds']]
        if not data_vars:
            print(f"  No data variables found in {ndep_fname}, skipping.")
            ndep_ds.close()
            continue
        
        var_name = data_vars[0]
        print(f"  Interpolating variable: {var_name}")
        
        src_data = ndep_ds[var_name].values
        n_times = src_data.shape[0]
        print(f"  Number of timesteps: {n_times}")
        
        interpolated_data = src_data[:, idx]
        
        ndep_out = xr.Dataset()
        ndep_out['time'] = ndep_ds['time']
        if 'time_bnds' in ndep_ds:
            ndep_out['time_bnds'] = ndep_ds['time_bnds']
        
        ndep_out[var_name] = (['time', 'ncells'], interpolated_data)
        ndep_out[var_name].attrs.update(ndep_ds[var_name].attrs)
        
        ndep_out['lon'] = (['ncells'], target_lon)
        ndep_out['lon'].attrs = {'standard_name': 'longitude', 'long_name': 'longitude', 'units': 'degrees_east'}
        ndep_out['lat'] = (['ncells'], target_lat)
        ndep_out['lat'].attrs = {'standard_name': 'latitude', 'long_name': 'latitude', 'units': 'degrees_north'}
        
        ndep_out.attrs.update(ndep_ds.attrs)
        ndep_out.attrs['history'] = f"Interpolated from TL255 to {ifs_grid} using nearest-neighbor; " + ndep_out.attrs.get('history', '')
        
        out_fname = ndep_fname.replace('TL255', ifs_grid)
        out_path = config.output_paths.lpj_guess / out_fname
        
        print(f"  Writing to {out_path}...")
        ndep_out.to_netcdf(out_path)
        print(f"  Saved {out_fname}")
        
        ndep_ds.close()
        ndep_out.close()


def _write_netcdf3_file(dataset: xr.Dataset, output_path: Path, file_label: str) -> bool:
    """Write xarray dataset to NetCDF3_64BIT_OFFSET format."""
    print(f"Writing {file_label} to {output_path} in NetCDF3_64BIT_OFFSET format")
    try:
        nc_out = Dataset(str(output_path), "w", format="NETCDF3_64BIT_OFFSET")
        
        for dim_name, dim_size in dataset.dims.items():
            nc_out.createDimension(dim_name, size=dim_size)
        
        for attr_name, attr_value in dataset.attrs.items():
            try:
                nc_out.setncattr(attr_name, attr_value)
            except Exception as e:
                print(f"Warning: Could not copy global attribute {attr_name}: {e}")
        
        for var_name, var in dataset.variables.items():
            var_data = var.values
            var_dims = var.dims
            var_dtype = var_data.dtype
            
            if var_dtype == np.float64:
                out_dtype = 'f8'
            elif var_dtype == np.float32:
                out_dtype = 'f4'
            elif var_dtype == np.int64:
                out_dtype = 'i4'
                var_data = var_data.astype(np.int32)
            elif var_dtype == np.int32:
                out_dtype = 'i4'
            elif var_dtype == np.int16:
                out_dtype = 'i2'
            elif var_dtype in (np.int8, np.uint8):
                out_dtype = 'i1'
            elif var_dtype == np.bool_:
                out_dtype = 'i1'
                var_data = var_data.astype(np.int8)
            else:
                out_dtype = 'f8'
                var_data = var_data.astype(np.float64)
            
            fill_value = None
            if '_FillValue' in var.attrs:
                fill_value = var.attrs['_FillValue']
                if out_dtype.startswith('i') and isinstance(fill_value, float):
                    fill_value = np.int32(fill_value)
                elif out_dtype.startswith('f') and isinstance(fill_value, int):
                    fill_value = float(fill_value)
            
            if fill_value is not None:
                var_out = nc_out.createVariable(var_name, out_dtype, var_dims, fill_value=fill_value)
            else:
                var_out = nc_out.createVariable(var_name, out_dtype, var_dims)
            
            for attr_name, attr_value in var.attrs.items():
                if attr_name != '_FillValue':
                    try:
                        var_out.setncattr(attr_name, attr_value)
                    except Exception as e:
                        print(f"Warning: Could not copy attribute {attr_name}: {e}")
            
            var_out[:] = var_data
        
        nc_out.close()
        print(f"NetCDF3 file creation complete for {file_label}")
        return True
        
    except Exception as e:
        print(f"Error writing NetCDF3 file for {file_label}: {e}")
        print("Falling back to standard xarray output")
        dataset.to_netcdf(output_path)
        return False
