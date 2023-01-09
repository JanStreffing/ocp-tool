import numpy as np

from .earth import RADIUS as EARTH_RADIUS


def _longitudes(N, *, loc='c'):
    if loc in ('center', 'c'):
        return np.linspace(0, 360, N+1)[:-1]
    bounds = np.linspace(0, 360, 2*N+1)[1::2]
    if loc in ('west', 'w', 'left', 'l'):
        return np.roll(bounds, 1)
    elif loc in ('east', 'e', 'right', 'r'):
        # Note: east bound of first cell is defined >0, not negative!
        return bounds


def _latitude_bounds(lats, *, loc):
    centers = 0.5*(lats[:-1]+lats[1:])
    if loc in ('south', 's', 'lower', 'l'):
        return np.array((*centers, -90))
    elif loc in ('north', 'n', 'upper', 'u'):
        return np.array((90, *centers))


class ReducedGaussianGrid:

    def __init__(self, lats, nlons):
        self.lats = np.array(lats)
        self.nlons = nlons

    def _repeat(self, values):
        """Repeats the values of a list (one value per latitude row), according
        to the number of cells (longitudes) in each row. The resulting vector
        has one value per grid cell."""
        return np.block(
            [np.repeat(v, n) for v, n in zip(values, self.nlons)]
        )

    def _tile(self, func, *args, **kwargs):
        """Creates, for each latitude row, a vector by calling func(n, ...),
        where n is the number of cells (longitudes) in that respective row. The
        results are concatenated to a vector that has one value for each grid
        cell."""
        return np.block(
            [func(n, *args, **kwargs) for n in self.nlons]
        )

    def cell_latitudes(self):
        return self._repeat(self.lats)

    def cell_longitudes(self):
        return self._tile(_longitudes)

    def _cell_corner_latitudes(self):
        return np.array(
            [
                self._repeat(_latitude_bounds(self.lats, loc='n')),
                self._repeat(_latitude_bounds(self.lats, loc='n')),
                self._repeat(_latitude_bounds(self.lats, loc='s')),
                self._repeat(_latitude_bounds(self.lats, loc='s')),
            ]
        )

    def _cell_corner_longitudes(self):
        return np.array(
            [
                self._tile(_longitudes, loc='e'),
                self._tile(_longitudes, loc='w'),
                self._tile(_longitudes, loc='w'),
                self._tile(_longitudes, loc='e'),
            ]
        )

    def cell_corners(self):
        return np.array(
            [self._cell_corner_latitudes(), self._cell_corner_longitudes()]
        )

    def cell_areas(self):
        areas = 2*np.pi*EARTH_RADIUS**2*np.abs(
            np.sin(np.radians(_latitude_bounds(self.lats, loc='n')))
            - np.sin(np.radians(_latitude_bounds(self.lats, loc='s')))
        )/self.nlons
        return self._repeat(areas)
