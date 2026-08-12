"""
Land-Sea Mask processing module for OCP-Tool.
Handles reading, modifying, and writing GRIB land-sea mask files.
"""

import copy
from pathlib import Path
from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass

import numpy as np
import gribapi
from shutil import copy2

from .config import OCPConfig
from .cycles import check_lake_fields
from .gaussian_grids import GaussianGrid


# FLake does not integrate the lake depth it is given. It clips it:
#
#   flakeene_mod.F90:211  ZDEPTH_W_MAX = 50.0   ! Maximum lake depth simulated
#   flakeene_mod.F90:212  ZDEPTH_W_MIN =  2.0   ! Minimum lake depth simulated
#   flakeene_mod.F90:238  ZDEPTH_W(JL) = MIN(ZDEPTH_W_MAX, MAX(ZDEPTH_W_MIN, PDEPTH_W(JL)))
#
# and independently again when seeding the mixed layer,
#
#   surftstp_ctl_mod.F90:882  ZHLML(JL) = MIN(PLDEPTH(JL), 50.0_JPRB)
#
# so 2 m and 50 m are the only depths at the ends of the range that mean
# anything, and a deeper number is indistinguishable from 50 m.
FLAKE_DEPTH_MIN_M = 2.0
FLAKE_DEPTH_MAX_M = 50.0

# FLake prognostics, restored alongside the cover so a recovered lake starts
# from its own state. Ordering is irrelevant; presence is checked per field.
FLAKE_PROGNOSTICS = ('lmlt', 'lblt', 'lmld', 'lict', 'licd', 'ltlt', 'lshf')


@dataclass
class LSMData:
    """Container for land-sea mask processing results."""
    lsm_binary_atm: np.ndarray      # LSM for atmosphere (lakes as land)
    lsm_binary_land: np.ndarray     # LSM for land points only
    lsm_binary_runoff: np.ndarray   # LSM for runoff mapping
    gribfield_modified: List[np.ndarray]  # Modified GRIB fields
    land_lats: List[float]          # Latitudes of land points
    land_lons: List[float]          # Longitudes of land points


def count_grib_fields(grib_file: Path) -> int:
    """
    Auto-detect the number of fields in a GRIB file.
    
    Args:
        grib_file: Path to GRIB file
        
    Returns:
        Number of GRIB messages in file
    """
    count = 0
    with open(grib_file, 'rb') as f:
        while True:
            gid = gribapi.grib_new_from_file(f)
            if gid is None:
                break
            gribapi.grib_release(gid)
            count += 1
    return count


def read_grib_fields(
    input_file: Path,
    num_fields: Optional[int] = None,
    verbose: bool = False
) -> Tuple[List[np.ndarray], int, int, int, List, int]:
    """
    Read GRIB fields from OpenIFS ICMGG file.
    
    Args:
        input_file: Path to ICMGG GRIB file
        num_fields: Number of fields to read (auto-detected if None)
        verbose: Print debug info
        
    Returns:
        Tuple of (gribfield, lsm_id, slt_id, cl_id, gid_list, num_fields)
    """
    print(f"\n {'='*50} \n")
    print(f" Opening Grib input file: {input_file}")
    
    # Auto-detect number of fields if not specified
    if num_fields is None:
        num_fields = count_grib_fields(input_file)
        print(f" Auto-detected {num_fields} fields in GRIB file")
    
    gid = [None] * num_fields
    gribfield = [None] * num_fields
    lsm_id = slt_id = cl_id = -1
    lsm_saved = False
    
    with open(input_file, 'rb') as f:
        keys = ['N', 'shortName']
        
        if verbose:
            print(f'num_fields: {num_fields}')
        
        for i in range(num_fields):
            gid[i] = gribapi.grib_new_from_file(f)
            if gid[i] is None:
                break
            
            for key in keys:
                if not gribapi.grib_is_defined(gid[i], key):
                    raise ValueError(f"Key '{key}' was not defined")
                if verbose:
                    print(f'{key}={gribapi.grib_get(gid[i], key)}')
            
            short_name = gribapi.grib_get(gid[i], 'shortName')
            
            if short_name == 'lsm':
                lsm_id = i
                num_timesteps = gribapi.grib_get(gid[i], 'numberOfDataPoints')
                print(f'number of lsm timesteps: {num_timesteps}')
            elif short_name == 'slt':
                slt_id = i
            elif short_name == 'cl':
                cl_id = i
            
            if not lsm_saved or short_name != 'lsm':
                values = gribapi.grib_get_values(gid[i])
                if verbose:
                    print(f"Shape of '{short_name}' values: {np.shape(values)}")
                gribfield[i] = values
    
    # Filter out None values
    gribfield = [field for field in gribfield if field is not None]
    print(f"\n {'='*50} \n")
    print(f'Shape of Gribfield: {np.shape(gribfield)}')
    
    return gribfield, lsm_id, slt_id, cl_id, gid, num_fields


def fill_flipped_from_nearest_neighbour(
    gribfield_mod: List[np.ndarray],
    lsm_id: int,
    grid: GaussianGrid,
    changed_to_land: List[int],
    changed_to_ocean: List[int],
    verbose: bool = False,
) -> None:
    """
    Overwrite every surface field of each land-sea-flipped cell with the value
    from the nearest stable neighbour of its NEW type (great-circle nearest).

    The land-sea mask itself (``lsm_id``) is left untouched so the freshly set
    coastline is preserved. Donor cells exclude all flipped cells, so a flipped
    cell is never seeded from another flipped cell.
    """
    changed = set(changed_to_land) | set(changed_to_ocean)
    if not changed:
        return

    lats = np.deg2rad(np.asarray(grid.center_lats[0], dtype=float))
    lons = np.deg2rad(np.asarray(grid.center_lons[0], dtype=float))
    # unit vectors on the sphere -> nearest = largest dot product
    vx = np.cos(lats) * np.cos(lons)
    vy = np.cos(lats) * np.sin(lons)
    vz = np.sin(lats)

    lsm = np.asarray(gribfield_mod[lsm_id], dtype=float)
    npts = lsm.shape[0]
    changed_mask = np.zeros(npts, dtype=bool)
    changed_mask[list(changed)] = True
    is_land = lsm >= 0.5
    donor_land = np.where(is_land & ~changed_mask)[0]
    donor_ocean = np.where(~is_land & ~changed_mask)[0]

    def _nearest(i, donors):
        dots = vx[donors] * vx[i] + vy[donors] * vy[i] + vz[donors] * vz[i]
        return int(donors[int(np.argmax(dots))])

    def _seed(cells, donors, label):
        if len(cells) and donors.size == 0:
            print(f'  WARNING: no donor cells of type {label}; '
                  f'{len(cells)} flipped cells left unchanged')
            return
        for i in cells:
            j = _nearest(i, donors)
            for f in range(len(gribfield_mod)):
                if f == lsm_id:
                    continue
                if len(gribfield_mod[f]) == npts:
                    gribfield_mod[f][i] = gribfield_mod[f][j]

    _seed(changed_to_land, donor_land, 'land')
    _seed(changed_to_ocean, donor_ocean, 'ocean')

    if verbose:
        print(f'  NN-filled {len(changed_to_land)} new-land and '
              f'{len(changed_to_ocean)} new-ocean cells from nearest neighbours')


def read_field_ids(input_file: Path) -> Dict[str, int]:
    """
    shortName -> position in the gribfield list, in file order.

    ``read_grib_fields`` returns indices for only the three fields it needs.
    Anything else that has to be addressed by name -- the FLake set, here --
    needs the same mapping for the rest of the file. Built by a second scan
    rather than threaded through the existing return tuple, which several
    callers unpack positionally.

    On a duplicated shortName the FIRST occurrence wins, matching
    ``read_grib_fields``.
    """
    ids: Dict[str, int] = {}
    with open(input_file, 'rb') as f:
        i = 0
        while True:
            gid = gribapi.grib_new_from_file(f)
            if gid is None:
                break
            try:
                ids.setdefault(gribapi.grib_get(gid, 'shortName'), i)
            finally:
                gribapi.grib_release(gid)
            i += 1
    return ids


def _clip_flake_depth(depth: float) -> float:
    """Clip to the band FLake actually integrates over."""
    return float(min(FLAKE_DEPTH_MAX_M, max(FLAKE_DEPTH_MIN_M, depth)))


def restore_lake_cover(
    gribfield_pristine: Dict[int, np.ndarray],
    gribfield_mod: List[np.ndarray],
    lake_ids: Dict[str, int],
    changed_to_land: List[int],
    verbose: bool = False,
) -> int:
    """
    Put back the lake fields that the nearest-neighbour fill overwrote.

    Must run AFTER ``fill_flipped_from_nearest_neighbour``, because that is what
    destroyed them. The pristine snapshot is needed as well as the modified
    list: the point is to recover what the input file said before the fill
    replaced it with a dry neighbour's value.

    ONE RULE, NO THRESHOLD. Every flipped cell gets its own pristine ``cl``,
    ``dl`` and FLake prognostics back. Nothing is invented and nothing is
    reclassified, because the input file has already done the classification
    and it is more informative than any threshold applied here:

    ``cl`` is the lake share of the grid box, capped by the water share -- on
    TCO95 ``cl <= 1 - lsm`` holds at every one of the 40 320 cells. So at the
    cells that were genuinely lake (the Caspian, the Great Lakes) the file says
    ``cl = 1 - lsm`` exactly, mean 0.799; and at the coastal cells whose water
    is sea rather than lake it says ``cl = 0.019`` against a water share of
    0.663. Restoring the pristine value therefore gives each population what
    ECMWF classified it as, in one operation.

    An earlier version of this thresholded on ``cl >= 0.5`` and wrote ``1.0``.
    That overstated lake area at those cells by 20 % on average, since 27 of
    the 70 sit between 0.5 and 0.7, and it needed a companion option to decide
    what to do with everything else. Both problems are artefacts of the
    threshold, not of the data.

    ``dl`` is clipped to the band FLake actually integrates,
    ``MIN(50, MAX(2, dl))`` (``flakeene_mod.F90:238`` and
    ``surftstp_ctl_mod.F90:882``), so the file says what the model will do.

    Returns
    -------
    int : how many cells were restored.
    """
    if not changed_to_land:
        return 0

    cl_id, dl_id = lake_ids['cl'], lake_ids['dl']
    idx = np.asarray(changed_to_land, dtype=int)

    # Everything FLake needs, so a recovered lake starts from its own state
    # rather than from whatever dry neighbour the fill happened to pick.
    fields = [lake_ids[n] for n in ('cl', 'dl') + FLAKE_PROGNOSTICS
              if n in lake_ids and lake_ids[n] in gribfield_pristine]

    for i in idx:
        for f in fields:
            if len(gribfield_mod[f]) > i and len(gribfield_pristine[f]) > i:
                gribfield_mod[f][i] = gribfield_pristine[f][i]
        gribfield_mod[dl_id][i] = _clip_flake_depth(gribfield_pristine[dl_id][i])

    if verbose:
        pris = np.asarray(gribfield_pristine[cl_id], float)[idx]
        print(f'  Lake fields restored at {len(idx)} flipped cells '
              f'(mean cl {pris.mean():.3f}; '
              f'{int((pris >= 0.5).sum())} of them are majority lake)')
    return len(idx)


def modify_lsm(
    gribfield: List[np.ndarray],
    ocean_lsm: np.ndarray,
    ocean_grid_name: str,
    lsm_id: int,
    slt_id: int,
    cl_id: int,
    grid: GaussianGrid,
    verbose: bool = False,
    lake_ids: Optional[Dict[str, int]] = None,
    lakes: Optional['LakeConfig'] = None,
) -> LSMData:
    """
    Modify land-sea mask based on ocean model grid.

    - Flips the land-sea mask to agree with the ocean model grid
    - Rebuilds each flipped cell's full surface column from the nearest stable
      neighbour of its new type (so it is physically self-consistent)
    - Optionally gives back the lake cover that rebuild removed
    - Creates masks for different purposes (atmosphere, land, runoff)

    Args:
        gribfield: List of GRIB field arrays
        ocean_lsm: Ocean model land-sea mask
        ocean_grid_name: Name of ocean grid
        lsm_id: Index of LSM field in gribfield
        slt_id: Index of soil type field
        cl_id: Index of lake cover field
        grid: Gaussian grid data
        verbose: Print debug info
        lake_ids: shortName -> gribfield index, needed only for the lake step
        lakes: LakeConfig; the lake step is skipped when it is None or off

    Returns:
        LSMData with modified masks and land point coordinates
    """
    # Create copies for different mask types
    lsm_binary_land = copy.deepcopy(gribfield[lsm_id])
    lsm_binary_land = lsm_binary_land[np.newaxis, :]
    lsm_binary_runoff = lsm_binary_land.copy()
    gribfield_mod = gribfield[:]
    
    # ocean_lsm is the ocean model grid when coupled, the reconstruction mask
    # for a paleo AMIP run, and empty otherwise.
    if len(ocean_lsm):
        # Automatic lake removal based on ocean mask
        # Polygon method: ocean_lsm = 1 means land, 0 means ocean
        n_points = len(gribfield_mod[slt_id])

        # gribfield_mod is a SHALLOW copy, so the NN fill below mutates these
        # arrays in place and the input values are gone once it has run. The
        # lake step needs to know what the file said before that, so snapshot
        # the fields it reads while they are still pristine.
        do_lakes = lakes is not None and lakes.restore_flipped_lakes
        pristine_lake = {}
        if do_lakes:
            if lake_ids is None or 'cl' not in lake_ids or 'dl' not in lake_ids:
                raise ValueError(
                    'lakes.restore_flipped_lakes is on but the cl/dl field '
                    'indices were not supplied to modify_lsm().'
                )
            for name in ('cl', 'dl') + FLAKE_PROGNOSTICS:
                f = lake_ids.get(name)
                if f is not None and f < len(gribfield_mod):
                    pristine_lake[f] = np.array(gribfield_mod[f], copy=True)

        changed_to_ocean = []
        changed_to_land = []
        for i in range(n_points - 1):
            # Point is land in IFS but ocean in FESOM -> make it ocean
            if gribfield_mod[lsm_id][i] >= 0.5 and ocean_lsm[i] < 0.5:
                gribfield_mod[lsm_id][i] = 0
                changed_to_ocean.append(i)
            # Point is ocean in IFS but land in FESOM -> make it land
            elif gribfield_mod[lsm_id][i] <= 0.5 and ocean_lsm[i] >= 0.5:
                gribfield_mod[lsm_id][i] = 1
                changed_to_land.append(i)

        # A flipped cell otherwise keeps the whole surface column of its OLD
        # type: a new land point retains ocean's zero soil moisture, leftover
        # sea-ice fraction and SST; a new ocean point retains land soil/skin
        # state. That mix is physically inconsistent and makes the atmosphere
        # surface/moist physics produce NaNs. Rebuild each flipped cell from the
        # nearest stable neighbour of its NEW type so every surface field looks
        # like its surroundings.
        fill_flipped_from_nearest_neighbour(
            gribfield_mod, lsm_id, grid,
            changed_to_land, changed_to_ocean, verbose=verbose,
        )

        # ...and the fill has just overwritten the lake cover of every cell that
        # used to be one, because cl is only another surface field to it.
        if do_lakes:
            restore_lake_cover(
                gribfield_pristine=pristine_lake,
                gribfield_mod=gribfield_mod,
                lake_ids=lake_ids,
                changed_to_land=changed_to_land,
                verbose=verbose,
            )
    else:
        print(" Skipped modifying OpenIFS grid, because we are in AMIP mode")
    
    # Mask with lakes counting as land (for atmosphere coupling)
    lsm_binary_atm = gribfield_mod[lsm_id]
    lsm_binary_atm = lsm_binary_atm[np.newaxis, :]
    
    # Generate list of land point coordinates for LPJ-GUESS
    land_lats = []
    land_lons = []
    n_points = len(gribfield_mod[slt_id])
    
    for i in range(n_points - 1):
        if gribfield_mod[lsm_id][i] >= 0.5:
            land_lats.append(grid.center_lats[0, i])
            land_lons.append(grid.center_lons[0, i])
    
    print(f'Number of land points: {len(land_lats)}')
    
    return LSMData(
        lsm_binary_atm=lsm_binary_atm,
        lsm_binary_land=lsm_binary_land,
        lsm_binary_runoff=lsm_binary_runoff,
        gribfield_modified=gribfield_mod,
        land_lats=land_lats,
        land_lons=land_lons
    )


def write_modified_grib(
    gribfield_mod: List[np.ndarray],
    input_file: Path,
    output_file: Path,
    num_fields: int,
    gid_list: List,
    lsm_id: int,
    verbose: bool = False
) -> None:
    """
    Write modified GRIB fields to output file.
    
    Args:
        gribfield_mod: Modified GRIB field arrays
        input_file: Original input GRIB file
        output_file: Output file path
        num_fields: Number of fields
        gid_list: GRIB handles from reading
        lsm_id: Index of LSM field
        verbose: Print debug info
    """
    # Copy original file to output location
    copy2(input_file, output_file)
    
    # Read LSM datadates for multi-timestep handling
    with open(output_file, 'rb') as f:
        lsm_datadates = []
        while True:
            gidi = gribapi.grib_new_from_file(f)
            if gidi is None:
                break
            short_name = gribapi.grib_get(gidi, 'shortName')
            if short_name == 'lsm':
                data_date = gribapi.grib_get(gidi, 'dataDate')
                lsm_datadates.append((gidi, data_date))
        
        print(f'lsm_datadates: {lsm_datadates}')
    
    # Write modified fields
    with open(output_file, 'r+b') as f:
        lsm_index = 0
        
        # Write LSM fields with correct dates
        for gidi, data_date in lsm_datadates:
            print(f'Overwriting lsm for date: {data_date}')
            
            if lsm_index < len(gribfield_mod):
                gribapi.grib_set_values(gidi, gribfield_mod[lsm_id].flatten())
                gribapi.grib_set(gidi, 'dataDate', data_date)
                gribapi.grib_write(gidi, f)
                lsm_index += 1
            else:
                print(f"Error: Index {lsm_index} out of range for gribfield_mod.")
            
            gribapi.grib_release(gidi)
        
        # Write non-LSM fields
        for i in range(num_fields):
            short_name = gribapi.grib_get(gid_list[i], 'shortName')
            if short_name != 'lsm':
                gribapi.grib_set_values(gid_list[i], gribfield_mod[i])
                gribapi.grib_write(gid_list[i], f)
            gribapi.grib_release(gid_list[i])


def process_land_sea_mask(
    config: OCPConfig,
    grid: GaussianGrid,
    ocean_lsm: np.ndarray,
    resolution: int
) -> LSMData:
    """
    Complete land-sea mask processing pipeline.
    
    Args:
        config: OCP configuration
        grid: Gaussian grid data
        ocean_lsm: Ocean model land-sea mask
        resolution: Truncation number
        
    Returns:
        LSMData with all mask data
    """
    input_file = config.get_icmgg_input_file()
    output_file = config.get_icmgg_output_file()
    
    # Read GRIB fields (auto-detect num_fields)
    gribfield, lsm_id, slt_id, cl_id, gid_list, num_fields = read_grib_fields(
        input_file,
        num_fields=None,  # Auto-detect
        verbose=config.options.verbose
    )
    
    # The lake step addresses fields by name, which read_grib_fields does not
    # provide beyond lsm/slt/cl, and it must not run against a file that cannot
    # support it -- so resolve the indices and verify the file up front, before
    # anything has been modified.
    lake_ids = None
    if config.lakes.restore_flipped_lakes:
        check_lake_fields(input_file, config.cycle)
        lake_ids = read_field_ids(input_file)

    # Modify LSM based on ocean grid
    lsm_data = modify_lsm(
        gribfield=gribfield,
        ocean_lsm=ocean_lsm,
        ocean_grid_name=config.ocean.grid_name,
        lsm_id=lsm_id,
        slt_id=slt_id,
        cl_id=cl_id,
        grid=grid,
        verbose=config.options.verbose,
        lake_ids=lake_ids,
        lakes=config.lakes,
    )
    
    # Write modified GRIB
    write_modified_grib(
        gribfield_mod=lsm_data.gribfield_modified,
        input_file=input_file,
        output_file=output_file,
        num_fields=num_fields,
        gid_list=gid_list,
        lsm_id=lsm_id,
        verbose=config.options.verbose
    )
    
    return lsm_data


def create_slt_output_for_lpjg(
    config: OCPConfig,
    resolution: int
) -> bool:
    """
    Create separate SLT (soil type) output file for LPJ-GUESS.
    
    Uses eccodes (grib_copy) and CDO to extract soil type field
    and convert to NetCDF format.
    
    Args:
        config: OCP configuration
        resolution: Truncation number
        
    Returns:
        True if successful, False otherwise
    """
    import os
    
    print(f"\n {'='*50} \n")
    print(" Creating SLT output file for LPJG")
    print(f"\n {'='*50} \n")
    
    input_file = config.get_icmgg_output_file()
    
    # Generate output filename with ocean grid name
    suffix = config.options.output_suffix
    if config.atmosphere.truncation_type == 'cubic-octahedral':
        slt_output_name = f'slt_TCO{resolution}_{config.ocean.grid_name}{suffix}.nc'
    elif config.atmosphere.truncation_type == 'linear':
        slt_output_name = f'slt_TL{resolution}_{config.ocean.grid_name}{suffix}.nc'
    else:
        raise ValueError(f"Unknown truncation type: {config.atmosphere.truncation_type}")
    
    slt_output_path = config.output_paths.lpj_guess / slt_output_name
    temp_grib_file = config.output_paths.lpj_guess / 'temp_slt_var43.grb'
    
    if not input_file.exists():
        print(f"Error: Input ICMGG file not found: {input_file}")
        return False
    
    try:
        print(f"Extracting SLT field (variable 43) from {input_file}")
        
        # Extract SLT field using grib_copy
        cmd_extract = f"grib_copy -w shortName=slt {input_file} {temp_grib_file}"
        if config.options.verbose:
            print(f"Running: {cmd_extract}")
        
        result = os.system(cmd_extract)
        if result != 0:
            print(f"Error: grib_copy command failed with exit code {result}")
            return False
        
        if not temp_grib_file.exists():
            print(f"Error: Temporary GRIB file not created: {temp_grib_file}")
            return False
        
        print(f"Converting GRIB to NetCDF: {slt_output_path}")
        
        # Convert to NetCDF
        cmd_convert = f"cdo -f nc copy {temp_grib_file} {slt_output_path}"
        if config.options.verbose:
            print(f"Running: {cmd_convert}")
        
        result = os.system(cmd_convert)
        if result != 0:
            print(f"Error: cdo command failed with exit code {result}")
            if temp_grib_file.exists():
                temp_grib_file.unlink()
            return False
        
        # Rename variable from 'slt' to 'var43' (LPJ-GUESS expects var43)
        cmd_rename = f"ncrename -v slt,var43 {slt_output_path}"
        if config.options.verbose:
            print(f"Running: {cmd_rename}")
        
        result = os.system(cmd_rename)
        if result != 0:
            print(f"Error: ncrename failed (exit {result}). Is NCO on PATH?")
            return False

        # LPJ-GUESS reads var43 and falls back to slt, so a file left with the
        # wrong name still runs and the mistake surfaces much later. Check.
        import netCDF4
        with netCDF4.Dataset(slt_output_path) as ds:
            if 'var43' not in ds.variables:
                print(f"Error: {slt_output_path} has no var43, only "
                      f"{list(ds.variables)}")
                return False


        # Cleanup
        if temp_grib_file.exists():
            temp_grib_file.unlink()
            if config.options.verbose:
                print(f"Removed temporary file: {temp_grib_file}")
        
        if slt_output_path.exists():
            print(f"Successfully created SLT file: {slt_output_path}")
            if config.options.verbose:
                print(f"File size: {slt_output_path.stat().st_size} bytes")
                os.system(f"ncdump -h {slt_output_path}")
            return True
        else:
            print(f"Error: Output NetCDF file was not created: {slt_output_path}")
            return False
            
    except Exception as e:
        print(f"Unexpected error creating SLT file: {e}")
        import traceback
        traceback.print_exc()
        
        if temp_grib_file.exists():
            temp_grib_file.unlink()
        
        return False
