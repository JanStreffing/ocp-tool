"""
Plotting module for OCP-Tool.
Handles visualization of land-sea masks and runoff maps.
Uses polygon-based plotting with cell corners for accurate representation.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from .config import OCPConfig
from .gaussian_grids import GaussianGrid
from .lsm import LSMData


def _build_cell_vertices(grid: GaussianGrid) -> np.ndarray:
    """
    Build cell vertices array from grid corner coordinates.
    
    Args:
        grid: GaussianGrid with corner_lons and corner_lats
        
    Returns:
        Array of shape (n_cells, 4, 2) with [lon, lat] for each corner
    """
    # corner_lons/lats have shape (4, 1, n_points), squeeze and transpose
    clo = grid.corner_lons.squeeze().T  # (n_points, 4)
    cla = grid.corner_lats.squeeze().T  # (n_points, 4)
    cell_verts = np.stack([clo, cla], axis=-1)  # (n_points, 4, 2)
    return cell_verts


def _get_valid_cell_mask(grid: GaussianGrid) -> np.ndarray:
    """
    Create mask for cells that don't wrap around the date line.
    
    Args:
        grid: GaussianGrid with corner coordinates
        
    Returns:
        Boolean mask array
    """
    clo = grid.corner_lons.squeeze().T  # (n_points, 4)
    lon_range = np.max(clo, axis=1) - np.min(clo, axis=1)
    return lon_range < 180


def plot_land_sea_mask(
    config: OCPConfig,
    grid: GaussianGrid,
    lsm_data: LSMData,
    resolution: int
) -> None:
    """
    Plot the final land-sea mask using polygon cells.
    
    Args:
        config: OCP configuration
        grid: Gaussian grid data with corners
        lsm_data: Land-sea mask data
        resolution: Truncation number
    """
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    
    # Build cell vertices from corners
    cell_verts = _build_cell_vertices(grid)
    valid_mask = _get_valid_cell_mask(grid)
    
    # Get LSM values (flatten to 1D)
    lsm_atm = lsm_data.lsm_binary_atm.flatten()
    lsm_land = lsm_data.lsm_binary_land.flatten()
    
    # Create color array: 0=wet (blue), 1=land (tan), new dry points (red)
    # wet points: lsm_atm < 1
    # new dry points: lsm_land < 1 but lsm_atm >= 1 (points that became dry)
    colors = np.zeros(len(lsm_atm))
    colors[lsm_atm >= 0.5] = 1  # Land
    colors[lsm_atm < 0.5] = 0   # Wet
    
    # Identify new dry points (were wet in atm, now dry in land)
    new_dry = (np.round(lsm_land) < 1) & (np.round(lsm_atm) >= 0.5)
    colors[new_dry] = 2  # New dry points
    
    # Apply valid mask
    valid_verts = cell_verts[valid_mask]
    valid_colors = colors[valid_mask]
    
    # Create custom colormap: blue=wet, tan=land, red=new dry
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(['#4169E1', '#D2B48C', '#FF4444'])  # Blue, Tan, Red
    
    # Create PolyCollection
    collection = PolyCollection(
        valid_verts,
        array=valid_colors,
        cmap=cmap,
        edgecolors='none',
        linewidths=0,
        transform=ccrs.PlateCarree(),
    )
    collection.set_clim(0, 2)
    ax.add_collection(collection)
    
    # Set global extent
    ax.set_global()
    
    # Add coastlines for reference
    ax.coastlines(linewidth=0.5, color='black', zorder=5)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#4169E1', label='Wet points'),
        Patch(facecolor='#D2B48C', label='Land'),
        Patch(facecolor='#FF4444', label='New dry points'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')
    
    ax.set_title(f'Land-Sea Mask T{resolution}', fontsize=14)
    
    output_file = config.output_paths.plots / f'land_points_T{resolution}.png'
    fig.savefig(str(output_file), format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Saved LSM plot to {output_file}")


def plot_runoff_maps(
    config: OCPConfig,
    drainage: np.ndarray,
    arrival: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray
) -> None:
    """
    Plot runoff drainage basin and arrival point maps.
    
    Args:
        config: OCP configuration
        drainage: Drainage basin ID array
        arrival: Arrival point ID array
        lons: Longitude array
        lats: Latitude array
    """
    # Shift data to center on Prime meridian
    ds1, ds2 = np.hsplit(np.squeeze(drainage), 2)
    drainage_cat = np.concatenate((ds2, ds1), axis=1)
    ds1, ds2 = np.hsplit(np.squeeze(arrival), 2)
    arrival_cat = np.concatenate((ds2, ds1), axis=1)
    
    lons_shifted = lons - 180
    lon, lat = np.meshgrid(lons_shifted, lats)
    
    cmap = plt.cm.flag
    
    # Amazon region - arrival
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.set_extent([-60, -30, -10, 20], crs=ccrs.PlateCarree())
    ax.pcolormesh(lon, lat, arrival_cat, cmap=cmap, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.gridlines(draw_labels=True)
    output_file = config.output_paths.plots / 'runoff_amazon_arrival.png'
    fig.savefig(str(output_file), format='png')
    plt.close(fig)
    
    # Ob region - drainage
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.set_extent([50, 110, 40, 80], crs=ccrs.PlateCarree())
    ax.pcolormesh(lon, lat, drainage_cat, cmap=cmap, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.gridlines(draw_labels=True)
    output_file = config.output_paths.plots / 'runoff_ob_drainage.png'
    fig.savefig(str(output_file), format='png')
    plt.close(fig)
    
    # Ob region - arrival
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())
    ax.set_extent([50, 110, 40, 80], crs=ccrs.PlateCarree())
    ax.pcolormesh(lon, lat, arrival_cat, cmap=cmap, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.gridlines(draw_labels=True)
    output_file = config.output_paths.plots / 'runoff_ob_arrival.png'
    fig.savefig(str(output_file), format='png')
    plt.close(fig)
    
    print(f"Saved runoff plots to {config.output_paths.plots}")
