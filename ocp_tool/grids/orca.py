import numpy as np
from netCDF4 import Dataset


def _valid_subgrid(subgrid):
    return subgrid in ('t', 'u', 'v')


class ORCA:

    _orca_names = {
        (362, 292, 75): 'ORCA1L75',
    }

    def __init__(self, domain_cfg, masks=None):

        self.domain_cfg = domain_cfg
        self.masks = masks

        with Dataset(domain_cfg) as nc:
            if not {'x', 'y'}.issubset(nc.dimensions):
                raise RuntimeError(
                    'Missing dimensions in NEMO domain config'
                )
            if not {
                'glamt', 'glamu', 'glamv', 'glamf',
                'gphit', 'gphiu', 'gphiv', 'gphif',
                'e1t', 'e1u', 'e1v', 'e1f',
                'e2t', 'e2u', 'e2v', 'e2f',
                'top_level',
            }.issubset(nc.variables):
                raise RuntimeError(
                    'Missing variables in NEMO domain config'
                )
            try:
                self.name = self._orca_names[
                    (
                        nc.dimensions['x'].size,
                        nc.dimensions['y'].size,
                        nc.dimensions['z'].size,
                    )
                ]
            except KeyError:
                raise RuntimeError(
                    'Unknown dimensions in NEMO domain config'
                )
        if self.masks is not None:
            with Dataset(self.masks) as nc:
                if not {
                    'tmaskutil', 'umaskutil', 'vmaskutil'
                }.issubset(nc.variables):
                    raise RuntimeError(
                        'Missing variables in NEMO masks file'
                    )

    def cell_latitudes(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.domain_cfg) as nc:
            return nc.variables[f'gphi{subgrid}'][0, ...].data

    def cell_longitudes(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.domain_cfg) as nc:
            return nc.variables[f'glam{subgrid}'][0, ...].data

    def cell_areas(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')
        with Dataset(self.domain_cfg) as nc:
            return \
                nc.variables[f'e1{subgrid}'][0, ...].data \
                * nc.variables[f'e2{subgrid}'][0, ...].data

    def cell_masks(self, subgrid='t'):
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')

        def mask_borders(mask):
            mask[-1, :] = 1  # mask north-fold line
            mask[:, (0, -1)] = 1  # mask east+west borders

        # If a NEMO mask file is provided, just read T, U, V masks
        if self.masks is not None:
            with Dataset(self.masks) as nc:
                mask = np.where(
                    nc.variables[f'{subgrid}maskutil'][0, ...].data > 0, 0, 1
                )
                mask_borders(mask)
                return mask

        # Without a NEMO mask file, compute masks from top_level in domain_cfg
        with Dataset(self.domain_cfg) as nc:
            tmask = np.where(
                nc.variables['top_level'][0, ...].data == 0, 1, 0
            )
            if subgrid == 't':
                mask_borders(tmask)
                return tmask
            elif subgrid == 'u':
                umask = tmask \
                        * tmask.take(
                            range(1, tmask.shape[1]+1), axis=1, mode='wrap'
                          )
                mask_borders(umask)
                return umask
            elif subgrid == 'v':
                vmask = tmask \
                        * tmask.take(
                            range(1, tmask.shape[0]+1), axis=0, mode='clip'
                          )
                mask_borders(vmask)
                return vmask

    def cell_corners(self, subgrid='t'):
        """For the ORCA grid and staggered subgrids, see
        NEMO reference manual, section 4 'Space Domain (DOM)'

        Corner numbering used here:
        j
        ^  1 ------- 0
        |  |         |
        |  |         |
        |  |         |
        |  2 --------3
        +------------> i
        """
        if not _valid_subgrid(subgrid):
            raise ValueError(f'Invalid NEMO subgrid: {subgrid}')

        with Dataset(self.domain_cfg) as nc:
            if subgrid == 't':
                lats = nc.variables['gphif'][0, ...].data
                lons = nc.variables['glamf'][0, ...].data
            elif subgrid == 'u':
                lats = nc.variables['gphiv'][0, ...].data
                lons = nc.variables['glamv'][0, ...].data
            elif subgrid == 'v':
                lats = nc.variables['gphiu'][0, ...].data
                lons = nc.variables['glamu'][0, ...].data

        if lats.shape != lons.shape:
            raise ValueError(
                f'Incompatible lat/lon arrays in {self.domain_cfg}'
            )

        # Note that some corner lats/lons will be left undefined (set to an
        # invalid initial value), because we do not handle the north-fold
        # (v-grid) or the southern end of the grid over the Antarctic (t,
        # u-grids).
        corners = np.full((2, 4, *lats.shape), -99999.0)

        # give indices in the corners array sensible names
        lat, lon = 0, 1
        nw, ne, se, sw = 0, 1, 2, 3

        if subgrid == 't':
            corners[lat, ne, :, :] = np.roll(lats, 1, axis=1)
            corners[lon, ne, :, :] = np.roll(lons, 1, axis=1)
            corners[lat, se, 1:, :] = np.roll(lats[:-1, :], 1, axis=1)
            corners[lon, se, 1:, :] = np.roll(lons[:-1, :], 1, axis=1)

        elif subgrid == 'u':
            corners[lat, ne, :, :] = lats
            corners[lon, ne, :, :] = lons
            corners[lat, se, 1:, :] = lats[:-1, :]
            corners[lon, se, 1:, :] = lons[:-1, :]

        elif subgrid == 'v':
            corners[lat, ne, :-1, :] = np.roll(lats[1:, :], 1, axis=1)
            corners[lon, ne, :-1, :] = np.roll(lons[1:, :], 1, axis=1)
            corners[lat, se, :, :] = np.roll(lats, 1, axis=1)
            corners[lon, se, :, :] = np.roll(lons, 1, axis=1)

        # the westly corners are just copies of the eastly corners of the
        # "next" cell
        corners[:, nw, :, :] = np.roll(corners[:, ne, :, :], -1, axis=2)
        corners[:, sw, :, :] = np.roll(corners[:, se, :, :], -1, axis=2)

        return corners
