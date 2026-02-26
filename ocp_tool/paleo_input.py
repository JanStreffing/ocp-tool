"""
Paleo land surface modification module for OCP-Tool.

Replaces mod_input_oifs.sh: modifies ICMGG land surface variables
(ice, lakes, soils, soil water, vegetation type/cover/LAI, and remaining
fields) based on paleoclimate reconstruction data.

All interpolation done in Python (scipy) — no CDO dependency.
Operates directly on the OpenIFS Gaussian grid.
"""

import copy
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata, NearestNDInterpolator
import eccodes

from .config import OCPConfig, PaleoConfig
from .gaussian_grids import GaussianGrid
from .lsm import LSMData, read_grib_fields

# ---------------------------------------------------------------------------
# Lookup tables: PlioMIP3 reconstruction → OpenIFS categories
# Derived from mod_input_oifs.sh
# ---------------------------------------------------------------------------

# PlioMIP3 soil code → OpenIFS soil type (code 43)
SOIL_LOOKUP = {
    28: 1, 29: 1,   # → coarse
    30: 2, 31: 2,   # → medium
    32: 4,           # → organic
    33: 2,
    34: 3,           # → fine
    35: 4,
    36: 1,
    37: 3,
    38: 2,
    39: 2,
    40: 1,
    41: 3,
    42: 1,
    43: 3,
}

# OpenIFS soil type → volumetric soil water content per layer
SWVL_LOOKUP = {
    # layer: {soil_type: water_content}
    39: {1: 0.210, 2: 0.283, 3: 0.373, 4: 0.209},  # Layer 1, level 0
    40: {1: 0.212, 2: 0.279, 3: 0.374, 4: 0.248},  # Layer 2, level 7
    41: {1: 0.208, 2: 0.270, 3: 0.375, 4: 0.262},  # Layer 3, level 28
    42: {1: 0.188, 2: 0.300, 3: 0.399, 4: 0.255},  # Layer 4, level 100
}

# PlioMIP3 biome → OpenIFS high vegetation type (code 30)
HIGH_VEG_TYPE_LOOKUP = {1: 6, 2: 3, 6: 5, 7: 4}

# PlioMIP3 biome → OpenIFS low vegetation type (code 29)
LOW_VEG_TYPE_LOOKUP = {3: 7, 4: 2, 5: 11, 8: 9}

# High vegetation type → cover fraction (code 28)
HIGH_VEG_COVER_LOOKUP = {3: 0.95, 4: 0.925, 5: 0.95, 6: 0.95}

# Low vegetation type → cover fraction (code 27)
LOW_VEG_COVER_LOOKUP = {2: 0.87, 7: 0.90, 9: 0.40, 11: 0.10}

# High vegetation type → LAI (code 67)
HIGH_VEG_LAI_LOOKUP = {3: 5.9, 4: 4.89, 5: 5.97, 6: 6.44}

# Low vegetation type → LAI (code 66)
LOW_VEG_LAI_LOOKUP = {2: 2.92, 7: 3.79, 9: 2.95, 11: 3.11}

# GRIB codes for remaining fields that need drowned/added interpolation
REMAINING_FIELD_CODES = [
    139, 170, 183, 236,  # soil temperatures
    198,                  # skin reservoir content
    235,                  # skin temperature
    31,                   # sea ice fraction
    8, 9, 10, 11, 12, 13, 14,  # lake variables
    238,                  # snow layer temperature
    32, 33, 34,           # snow albedo, density, SST
    35, 36, 37, 38,       # ice temperatures
    148,                  # Charnock parameter
    174,                  # albedo
    15, 16, 17, 18,       # UV/IR albedos
    74,                   # sdfor
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _read_netcdf_field(filepath: Path, var_name: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read a 2D field from a NetCDF file.

    Returns (values, lons, lats) on the file's native grid.
    If var_name is None, reads the first non-coordinate variable.
    """
    ds = Dataset(str(filepath), 'r')
    coord_names = {'lon', 'lat', 'longitude', 'latitude', 'x', 'y', 'time',
                   'lon_bnds', 'lat_bnds', 'time_bnds', 'bounds_lon', 'bounds_lat'}

    if var_name is None:
        for name in ds.variables:
            if name.lower() not in coord_names:
                var_name = name
                break
        if var_name is None:
            ds.close()
            raise ValueError(f"No data variable found in {filepath}")

    data = np.squeeze(np.array(ds.variables[var_name][:]))

    # Find lon/lat variables
    lon_name = lat_name = None
    for candidate in ['lon', 'longitude', 'x']:
        if candidate in ds.variables:
            lon_name = candidate
            break
    for candidate in ['lat', 'latitude', 'y']:
        if candidate in ds.variables:
            lat_name = candidate
            break

    lons_1d = np.array(ds.variables[lon_name][:])
    lats_1d = np.array(ds.variables[lat_name][:])
    ds.close()

    # Build 2D coordinate arrays
    if data.ndim == 2:
        lons_2d, lats_2d = np.meshgrid(lons_1d, lats_1d)
    else:
        lons_2d = lons_1d
        lats_2d = lats_1d

    return data, lons_2d, lats_2d


def _interp_to_gaussian(
    src_data: np.ndarray,
    src_lons: np.ndarray,
    src_lats: np.ndarray,
    tgt_lons: np.ndarray,
    tgt_lats: np.ndarray,
    method: str = 'nearest',
    fill_value: float = 0.0,
) -> np.ndarray:
    """
    Interpolate source data to the OpenIFS Gaussian grid points.

    Parameters
    ----------
    src_data : 2D source field
    src_lons, src_lats : 2D coordinate arrays matching src_data
    tgt_lons, tgt_lats : 1D target coordinate arrays (Gaussian grid)
    method : 'nearest', 'linear', or 'cubic'
    fill_value : value for points outside convex hull (linear/cubic only)
    """
    src_pts = np.column_stack([src_lons.ravel(), src_lats.ravel()])
    src_vals = src_data.ravel()

    # Remove NaN / masked values
    valid = np.isfinite(src_vals)
    if hasattr(src_vals, 'mask'):
        valid &= ~src_vals.mask
    src_pts = src_pts[valid]
    src_vals = src_vals[valid]

    tgt_pts = np.column_stack([tgt_lons, tgt_lats])

    return griddata(src_pts, src_vals, tgt_pts, method=method, fill_value=fill_value)


def _distance_weighted_fill(
    field: np.ndarray,
    mask: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray,
) -> np.ndarray:
    """
    Fill missing/masked points using nearest-neighbour from valid points.

    Replaces CDO's ``setmisstodis``.

    Parameters
    ----------
    field : 1D array on Gaussian grid (values at valid points, anything at others)
    mask : boolean, True where values are valid (known)
    lons, lats : 1D coordinate arrays
    """
    if np.all(mask) or np.sum(mask) == 0:
        return field.copy()

    known_pts = np.column_stack([lons[mask], lats[mask]])
    known_vals = field[mask]
    interp = NearestNDInterpolator(known_pts, known_vals)

    result = field.copy()
    missing = ~mask
    query_pts = np.column_stack([lons[missing], lats[missing]])
    result[missing] = interp(query_pts)
    return result


def _apply_lookup(field: np.ndarray, lookup: Dict, default: float = 0.0) -> np.ndarray:
    """Map integer categories through a lookup table."""
    result = np.full_like(field, default, dtype=np.float64)
    for key, val in lookup.items():
        result[np.isclose(field, key, atol=0.5)] = val
    return result


def _get_grib_field_by_code(grib_file: Path, param_code: int) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Read a single GRIB field by parameter code, returning (values, lats, lons).
    """
    with open(grib_file, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                code = eccodes.codes_get(gid, 'paramId')
                if code == param_code:
                    vals = eccodes.codes_get_array(gid, 'values')
                    lats = eccodes.codes_get_array(gid, 'latitudes')
                    lons = eccodes.codes_get_array(gid, 'longitudes')
                    return vals, lats, lons
            finally:
                eccodes.codes_release(gid)
    return None, None, None


def _read_all_grib_fields(grib_file: Path) -> Dict[int, np.ndarray]:
    """Read all fields from a GRIB file, keyed by paramId."""
    fields = {}
    with open(grib_file, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                code = eccodes.codes_get(gid, 'paramId')
                vals = eccodes.codes_get_array(gid, 'values')
                fields[code] = vals
            finally:
                eccodes.codes_release(gid)
    return fields


def _replace_grib_fields(
    input_file: Path,
    output_file: Path,
    replacements: Dict[int, np.ndarray],
    verbose: bool = False,
) -> None:
    """
    Copy a GRIB file, replacing specific fields identified by paramId.

    Parameters
    ----------
    input_file : source GRIB
    output_file : destination GRIB
    replacements : dict mapping paramId → new values array
    """
    # Read all messages into memory first
    messages = []
    with open(input_file, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            messages.append(gid)

    with open(output_file, 'wb') as f:
        for gid in messages:
            code = eccodes.codes_get(gid, 'paramId')
            if code in replacements:
                eccodes.codes_set_array(gid, 'values', replacements[code].astype(np.float64))
                if verbose:
                    print(f"  Replaced paramId {code}")
            eccodes.codes_write(gid, f)
            eccodes.codes_release(gid)


# ---------------------------------------------------------------------------
# Main paleo input modification pipeline
# ---------------------------------------------------------------------------

def create_paleo_masks(
    config: OCPConfig,
    grid: GaussianGrid,
    lsm_data: LSMData,
) -> Dict[str, np.ndarray]:
    """
    Create derivative masks on the Gaussian grid.

    Uses the pre-modification (CORE2) and post-modification (paleo) LSM
    already available from the ocp-tool pipeline — no intermediate grid.

    Returns dict with keys: 'land', 'ocean', 'added', 'drowned', 'reduced_ice'
    """
    paleo = config.paleo
    verbose = config.options.verbose

    print("\n" + "=" * 60)
    print(" Creating paleo derivative masks")
    print("=" * 60)

    tgt_lons = np.array(grid.lons_list)
    tgt_lats = grid.center_lats.flatten()

    # Paleo LSM (post-modification) — already in lsm_data from step 3
    paleo_lsm = lsm_data.lsm_binary_atm.flatten()

    # Read CORE2 (pre-modification) LSM from original ICMGG file
    core2_lsm, _, _ = _get_grib_field_by_code(config.get_icmgg_input_file(), 172)
    if core2_lsm is None:
        raise RuntimeError("Could not read LSM (code 172) from input ICMGG file")

    # Binary masks on Gaussian grid
    land_mask = paleo_lsm >= 0.5
    ocean_mask = paleo_lsm < 0.5
    added_land = (core2_lsm < 0.5) & (paleo_lsm >= 0.5)   # ocean → land
    drowned = (core2_lsm >= 0.5) & (paleo_lsm < 0.5)       # land → ocean

    print(f"  Land points:  {np.sum(land_mask)}")
    print(f"  Ocean points: {np.sum(ocean_mask)}")
    print(f"  Added land:   {np.sum(added_land)}")
    print(f"  Drowned:      {np.sum(drowned)}")

    # Reduced ice sheet mask
    reduced_ice = np.zeros_like(paleo_lsm, dtype=bool)
    ice_file = paleo.get_reconstruction_file('ice_mask_file')
    if ice_file.exists():
        paleo_ice_data, ice_lons, ice_lats = _read_netcdf_field(ice_file)
        paleo_ice = _interp_to_gaussian(paleo_ice_data, ice_lons, ice_lats,
                                        tgt_lons, tgt_lats, method='nearest')

        # Read modern ice (snow depth, code 141) from ICMGG
        modern_ice, _, _ = _get_grib_field_by_code(config.get_icmgg_input_file(), 141)
        if modern_ice is not None:
            # Modern ice > 0 but paleo ice absent → reduced ice sheet
            modern_has_ice = modern_ice > 0.0
            paleo_has_ice = paleo_ice >= 2.0  # value 2 = ice sheet in PlioMIP3
            reduced_ice = modern_has_ice & ~paleo_has_ice & land_mask
            print(f"  Reduced ice:  {np.sum(reduced_ice)}")
    else:
        print(f"  Warning: ice mask file not found: {ice_file}")

    masks = {
        'land': land_mask,
        'ocean': ocean_mask,
        'added': added_land,
        'drowned': drowned,
        'reduced_ice': reduced_ice,
        'paleo_lsm': paleo_lsm,
        'core2_lsm': core2_lsm,
    }
    return masks


def modify_ice_snow_depth(
    config: OCPConfig,
    grid: GaussianGrid,
    masks: Dict[str, np.ndarray],
) -> np.ndarray:
    """
    Modify snow depth (code 141) from reconstruction ice mask.

    PlioMIP3 convention: ice_mask value 1 → no ice (snow=0),
    value 2 → ice sheet (snow depth = 10 m w.e.)
    """
    paleo = config.paleo

    print("\n  Modifying ice / snow depth (code 141)...")

    tgt_lons = np.array(grid.lons_list)
    tgt_lats = grid.center_lats.flatten()
    n_pts = len(tgt_lons)

    ice_file = paleo.get_reconstruction_file('ice_mask_file')
    if not ice_file.exists():
        print(f"  Warning: ice mask not found: {ice_file}, skipping")
        return None

    ice_data, ice_lons, ice_lats = _read_netcdf_field(ice_file)
    ice_on_gauss = _interp_to_gaussian(ice_data, ice_lons, ice_lats,
                                       tgt_lons, tgt_lats, method='nearest')

    # Map: 1 → 0, 2 → 10
    snow_depth = np.zeros(n_pts, dtype=np.float64)
    snow_depth[np.isclose(ice_on_gauss, 2, atol=0.5)] = 10.0

    print(f"    Ice sheet points (snow=10): {np.sum(snow_depth > 0)}")
    return snow_depth


def modify_lakes(
    config: OCPConfig,
    grid: GaussianGrid,
    masks: Dict[str, np.ndarray],
    existing_lake: np.ndarray,
) -> np.ndarray:
    """
    Modify lake cover (code 26) using anomaly method.

    PlioMIP3 convention: lake fraction indicated by value 1100.
    """
    paleo = config.paleo

    print("\n  Modifying lake cover (code 26)...")

    tgt_lons = np.array(grid.lons_list)
    tgt_lats = grid.center_lats.flatten()

    lake_file = paleo.get_reconstruction_file('lake_file')
    modern_lake_file = paleo.get_modern_file('modern_lake_file')

    if not lake_file.exists() or not modern_lake_file.exists():
        print(f"  Warning: lake files not found, skipping")
        return None

    # Read and create binary lake masks (value 1100 = lake)
    paleo_lake_data, plk_lons, plk_lats = _read_netcdf_field(lake_file)
    modern_lake_data, mlk_lons, mlk_lats = _read_netcdf_field(modern_lake_file)

    paleo_lake_binary = (np.isclose(paleo_lake_data, 1100, atol=0.5)).astype(float)
    modern_lake_binary = (np.isclose(modern_lake_data, 1100, atol=0.5)).astype(float)

    # Interpolate both to Gaussian grid
    paleo_lake_gauss = _interp_to_gaussian(paleo_lake_binary, plk_lons, plk_lats,
                                           tgt_lons, tgt_lats, method='nearest')
    modern_lake_gauss = _interp_to_gaussian(modern_lake_binary, mlk_lons, mlk_lats,
                                            tgt_lons, tgt_lats, method='nearest')

    # Anomaly method
    lake_anom = paleo_lake_gauss - modern_lake_gauss
    new_lake = existing_lake + lake_anom

    # Clamp to 0 or 1
    new_lake = np.where(new_lake >= 0.5, 1.0, 0.0)

    print(f"    Lake points: {np.sum(new_lake > 0)}")
    return new_lake


def modify_soils(
    config: OCPConfig,
    grid: GaussianGrid,
    masks: Dict[str, np.ndarray],
) -> Optional[np.ndarray]:
    """
    Modify soil type (code 43) from reconstruction data with lookup table.
    """
    paleo = config.paleo

    print("\n  Modifying soil type (code 43)...")

    tgt_lons = np.array(grid.lons_list)
    tgt_lats = grid.center_lats.flatten()
    land_mask = masks['land']

    soil_file = paleo.get_reconstruction_file('soil_file')
    if not soil_file.exists():
        print(f"  Warning: soil file not found: {soil_file}, skipping")
        return None

    soil_data, s_lons, s_lats = _read_netcdf_field(soil_file)

    # Interpolate to Gaussian grid
    soil_gauss = _interp_to_gaussian(soil_data, s_lons, s_lats,
                                     tgt_lons, tgt_lats, method='nearest')

    # Distance-weighted fill over land (replace missing with code 31 = default)
    soil_gauss[~np.isfinite(soil_gauss)] = 0
    soil_on_land = np.where(land_mask, soil_gauss, 0)
    soil_on_land[land_mask & (soil_on_land == 0)] = 31  # default fill
    soil_filled = _distance_weighted_fill(soil_on_land, land_mask & (soil_on_land > 0),
                                          tgt_lons, tgt_lats)

    # Apply PlioMIP3 → OpenIFS lookup
    oifs_soil = _apply_lookup(soil_filled, SOIL_LOOKUP, default=0.0)

    # Ensure integers and mask ocean
    oifs_soil = np.round(oifs_soil).astype(np.float64)
    oifs_soil[~land_mask] = 0
    # Clamp to valid range 1–4 on land
    oifs_soil[land_mask & (oifs_soil < 1)] = 1
    oifs_soil[land_mask & (oifs_soil > 4)] = 1

    print(f"    Soil types on land — 1:{np.sum(oifs_soil==1)}, 2:{np.sum(oifs_soil==2)}, "
          f"3:{np.sum(oifs_soil==3)}, 4:{np.sum(oifs_soil==4)}")
    return oifs_soil


def modify_soil_water(
    soil_type: np.ndarray,
    masks: Dict[str, np.ndarray],
) -> Dict[int, np.ndarray]:
    """
    Compute volumetric soil water layers (codes 39–42) from soil type.
    """
    print("\n  Computing volumetric soil water layers (codes 39–42)...")

    results = {}
    for code, lut in SWVL_LOOKUP.items():
        swvl = np.zeros_like(soil_type)
        for stype, wval in lut.items():
            swvl[np.isclose(soil_type, stype, atol=0.5)] = wval
        swvl[~masks['land']] = 0
        results[code] = swvl

    return results


def modify_vegetation(
    config: OCPConfig,
    grid: GaussianGrid,
    masks: Dict[str, np.ndarray],
) -> Dict[int, np.ndarray]:
    """
    Modify vegetation type, cover, and LAI from reconstruction biome data.

    Returns dict mapping GRIB code → values array:
      30 (tvh), 29 (tvl), 28 (cvh), 27 (cvl), 67 (lai_hv), 66 (lai_lv)
    """
    paleo = config.paleo

    print("\n  Modifying vegetation (codes 27–30, 66–67)...")

    tgt_lons = np.array(grid.lons_list)
    tgt_lats = grid.center_lats.flatten()
    land_mask = masks['land']

    biome_file = paleo.get_reconstruction_file('biome_file')
    if not biome_file.exists():
        print(f"  Warning: biome file not found: {biome_file}, skipping")
        return {}

    biome_data, b_lons, b_lats = _read_netcdf_field(biome_file)

    # Interpolate to Gaussian grid
    biome_gauss = _interp_to_gaussian(biome_data, b_lons, b_lats,
                                      tgt_lons, tgt_lats, method='nearest')

    # Distance-weighted fill over land
    biome_gauss[~np.isfinite(biome_gauss)] = 0
    biome_on_land = np.where(land_mask, biome_gauss, 0)
    biome_on_land[land_mask & (biome_on_land == 0)] = 8  # default biome
    biome_filled = _distance_weighted_fill(biome_on_land, land_mask & (biome_on_land > 0),
                                           tgt_lons, tgt_lats)
    biome_filled[~land_mask] = 0

    # High vegetation type (code 30)
    thv = _apply_lookup(biome_filled, HIGH_VEG_TYPE_LOOKUP, default=0.0)
    thv[~land_mask] = 0

    # Low vegetation type (code 29)
    tlv = _apply_lookup(biome_filled, LOW_VEG_TYPE_LOOKUP, default=0.0)
    tlv[~land_mask] = 0

    # High vegetation cover (code 28)
    cvh = _apply_lookup(thv, HIGH_VEG_COVER_LOOKUP, default=0.0)
    cvh[~land_mask] = 0

    # Low vegetation cover (code 27)
    cvl = _apply_lookup(tlv, LOW_VEG_COVER_LOOKUP, default=0.0)
    cvl[~land_mask] = 0

    # High vegetation LAI (code 67)
    lai_hv = _apply_lookup(thv, HIGH_VEG_LAI_LOOKUP, default=0.0)
    lai_hv[~land_mask] = 0

    # Low vegetation LAI (code 66)
    lai_lv = _apply_lookup(tlv, LOW_VEG_LAI_LOOKUP, default=0.0)
    lai_lv[~land_mask] = 0

    results = {30: thv, 29: tlv, 28: cvh, 27: cvl, 67: lai_hv, 66: lai_lv}

    for code, arr in results.items():
        nonzero = np.sum(arr > 0)
        print(f"    Code {code}: {nonzero} non-zero points")

    return results


def interpolate_remaining_fields(
    config: OCPConfig,
    grid: GaussianGrid,
    masks: Dict[str, np.ndarray],
) -> Dict[int, np.ndarray]:
    """
    Interpolate remaining GRIB fields for added/drowned grid points.

    For each field:
    - Drowned points: fill from surrounding ocean values
    - Added land + reduced ice sheet points: fill from surrounding land values
    - Combine both contributions
    """
    verbose = config.options.verbose
    icmgg_file = config.get_icmgg_output_file()

    print("\n  Interpolating remaining fields for added/drowned points...")

    tgt_lons = np.array(grid.lons_list)
    tgt_lats = grid.center_lats.flatten()
    land_mask = masks['land']
    ocean_mask = masks['ocean']
    added = masks['added']
    drowned = masks['drowned']
    reduced_ice = masks['reduced_ice']

    # Points that need filling
    needs_land_fill = added | reduced_ice
    needs_ocean_fill = drowned

    if not np.any(needs_land_fill) and not np.any(needs_ocean_fill):
        print("    No added/drowned points — nothing to interpolate")
        return {}

    # Read all current fields from the (already LSM-adapted) ICMGG
    all_fields = _read_all_grib_fields(icmgg_file)

    results = {}
    for code in REMAINING_FIELD_CODES:
        if code not in all_fields:
            if verbose:
                print(f"    Code {code} not found in ICMGG, skipping")
            continue

        field = all_fields[code].copy()

        # --- Ocean fill for drowned points ---
        if np.any(needs_ocean_fill):
            ocean_known = ocean_mask & ~drowned
            if np.any(ocean_known):
                filled_ocean = _distance_weighted_fill(field, ocean_known, tgt_lons, tgt_lats)
                field[needs_ocean_fill] = filled_ocean[needs_ocean_fill]

        # --- Land fill for added / reduced-ice points ---
        if np.any(needs_land_fill):
            land_known = land_mask & ~added & ~reduced_ice
            if np.any(land_known):
                filled_land = _distance_weighted_fill(field, land_known, tgt_lons, tgt_lats)
                field[needs_land_fill] = filled_land[needs_land_fill]

        results[code] = field
        if verbose:
            print(f"    Interpolated code {code}")

    print(f"    Processed {len(results)} remaining fields")
    return results


# ---------------------------------------------------------------------------
# Top-level entry point
# ---------------------------------------------------------------------------

def modify_paleo_input(
    config: OCPConfig,
    grid: GaussianGrid,
    lsm_data: LSMData,
) -> Path:
    """
    Full paleo land-surface modification pipeline for ICMGG.

    This runs after the standard ocp-tool LSM adaptation (steps 1–3).
    It modifies the ICMGGaackINIT_<grid> file in-place, producing the
    final file with all reconstruction-based changes.

    Parameters
    ----------
    config : OCPConfig with paleo section enabled
    grid : Gaussian grid from step 1
    lsm_data : LSMData from step 3 (land-sea mask processing)

    Returns
    -------
    Path to the modified ICMGG file
    """
    paleo = config.paleo
    if not paleo or not paleo.enabled:
        raise ValueError("Paleo configuration is not enabled")

    icmgg_file = config.get_icmgg_output_file()
    verbose = config.options.verbose

    print("\n" + "=" * 60)
    print(f" Paleo Input Modification: {paleo.experiment_id}")
    print(f" Target: {icmgg_file}")
    print("=" * 60)

    # 1. Create derivative masks
    masks = create_paleo_masks(config, grid, lsm_data)

    # 2. Collect all field replacements
    replacements: Dict[int, np.ndarray] = {}

    # 3. Ice / snow depth (code 141)
    snow = modify_ice_snow_depth(config, grid, masks)
    if snow is not None:
        replacements[141] = snow

    # 4. Lakes (code 26)
    existing_lake, _, _ = _get_grib_field_by_code(icmgg_file, 26)
    if existing_lake is not None:
        new_lake = modify_lakes(config, grid, masks, existing_lake)
        if new_lake is not None:
            replacements[26] = new_lake

    # 5. Soils (code 43)
    soil_type = modify_soils(config, grid, masks)
    if soil_type is not None:
        replacements[43] = soil_type

        # 6. Volumetric soil water (codes 39–42)
        swvl = modify_soil_water(soil_type, masks)
        replacements.update(swvl)

    # 7. Vegetation type, cover, LAI (codes 27–30, 66–67)
    veg = modify_vegetation(config, grid, masks)
    replacements.update(veg)

    # 8. Remaining fields (drowned / added interpolation)
    remaining = interpolate_remaining_fields(config, grid, masks)
    replacements.update(remaining)

    # 9. Write all replacements into the ICMGG file
    print(f"\n  Writing {len(replacements)} modified fields to {icmgg_file}...")
    _replace_grib_fields(icmgg_file, icmgg_file, replacements, verbose=verbose)

    print(f"\n  Paleo input modification complete: {icmgg_file}")
    return icmgg_file
