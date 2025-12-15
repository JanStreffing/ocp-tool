"""
Plotting module for OCP-Tool.
Handles visualization of land-sea masks and runoff maps.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

from .config import OCPConfig
from .gaussian_grids import GaussianGrid
from .lsm import LSMData


def plot_land_sea_mask(
    config: OCPConfig,
    grid: GaussianGrid,
    lsm_data: LSMData,
    resolution: int
) -> None:
    """
    Plot the final land-sea mask showing wet and dry points.
    
    Args:
        config: OCP configuration
        grid: Gaussian grid data
        lsm_data: Land-sea mask data
        resolution: Truncation number
    """
    fig = plt.figure(figsize=(24, 14))
    ax = fig.add_subplot(111)
    
    # Extract wet and dry points
    lsm_atm = lsm_data.lsm_binary_atm
    lsm_land = lsm_data.lsm_binary_land
    
    xpts_atm = grid.center_lons[np.round(lsm_atm[:, :]) < 1]
    ypts_atm = grid.center_lats[np.round(lsm_atm[:, :]) < 1]
    xpts_land = grid.center_lons[np.round(lsm_land[:, :]) < 1]
    ypts_land = grid.center_lats[np.round(lsm_land[:, :]) < 1]
    
    ax.scatter(xpts_land, ypts_land, s=100/resolution, color='red', marker='.', label='New dry points')
    ax.scatter(xpts_atm, ypts_atm, s=200/resolution, marker='.', label='Wet points')
    ax.legend(loc="lower right")
    
    output_file = config.output_paths.plots / f'land_points_T{resolution}.png'
    fig.savefig(str(output_file), format='png', dpi=600)
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
    
    # Amazon region
    m = Basemap(
        llcrnrlon=-60., llcrnrlat=-10,
        urcrnrlon=-30., urcrnrlat=20.,
        resolution='l', area_thresh=1000., projection='cyl'
    )
    xi, yi = m(lon, lat)
    
    fig = plt.figure(figsize=(12, 8))
    m.pcolor(xi, yi, arrival_cat, cmap=cmap)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90., 120., 45.))
    m.drawmeridians(np.arange(0., 360., 90.))
    output_file = config.output_paths.plots / 'runoff_amazon_arrival.png'
    fig.savefig(str(output_file), format='png')
    plt.close(fig)
    
    # Ob region - drainage
    m = Basemap(
        llcrnrlon=50., llcrnrlat=40,
        urcrnrlon=110., urcrnrlat=80.,
        resolution='l', area_thresh=1000., projection='cyl'
    )
    
    fig = plt.figure(figsize=(12, 8))
    m.pcolor(xi, yi, drainage_cat, cmap=cmap)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90., 120., 45.))
    m.drawmeridians(np.arange(0., 360., 90.))
    output_file = config.output_paths.plots / 'runoff_ob_drainage.png'
    fig.savefig(str(output_file), format='png')
    plt.close(fig)
    
    # Ob region - arrival
    fig = plt.figure(figsize=(12, 8))
    m.pcolor(xi, yi, arrival_cat, cmap=cmap)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90., 120., 45.))
    m.drawmeridians(np.arange(0., 360., 90.))
    output_file = config.output_paths.plots / 'runoff_ob_arrival.png'
    fig.savefig(str(output_file), format='png')
    plt.close(fig)
    
    print(f"Saved runoff plots to {config.output_paths.plots}")
