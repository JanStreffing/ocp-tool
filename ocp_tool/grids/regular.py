import numpy as np

from .earth import RADIUS as EARTH_RADIUS


def _equidistant(start, end, N, *, first_at_start=False):
    """Divides the interval [start, end]  into N subintervals of equal size and
    returns the midpoints of the subintervals. If 'first_at_start' is true, the
    first interval is centered around the start point."""
    if first_at_start:
        return np.linspace(start, end, N+1)[:-1]
    else:
        return np.linspace(start, end, 2*N+1)[1::2]


def _interval_bounds(start, centers, end, *, loc, wrap=False):
    """Returns the upper or lower bounds of [start, end] subintervals defined by
    their centers. It is checked if the first center coincides with 'start' and
    it is in that case assumed that the first interval is centered around the
    starting point. If 'wrap' is True, the first lower bound (last upper) is
    returned the same as the last upper (first lower)."""
    center_midpoints = 0.5*(centers[:-1]+centers[1:])
    last_dx = 0.5*(end-centers[-1]) if np.isclose(start, centers[0]) else 0
    if loc in ('l', 'left', 'lower'):
        return np.array(
            (end-last_dx if wrap else start-last_dx, *center_midpoints)
        )
    elif loc in ('r', 'right', 'u', 'upper'):
        return np.array(
            (*center_midpoints, start-last_dx if wrap else end-last_dx)
        )
    raise ValueError(f"Invalid value for 'loc' argument: {loc}")


def _row_distribute(row, nrows):
    return np.tile(row, (nrows, 1))


def _col_distribute(col, ncols):
    return np.tile(col, (ncols, 1)).T


def _is_monotonic(array):
    """Checks for *strict* monotonicity"""
    darray = np.diff(array)
    return all(darray < 0) or all(darray > 0)


class LatLonGrid:

    def __init__(self, lats, lons, first_lat=-90):
        if not _is_monotonic((first_lat, *lats, -first_lat)):
            raise ValueError('Non-monotonic latitude values')
        if not _is_monotonic(lons):
            raise ValueError('Non-monotonic longitude values')

        self._OP = first_lat
        self.lats = np.array(lats)
        self.lons = np.array(lons)

    @property
    def nlats(self):
        return len(self.lats)

    @property
    def nlons(self):
        return len(self.lons)

    def cell_latitudes(self):
        return _col_distribute(self.lats, len(self.lons))

    def cell_longitudes(self):
        return _row_distribute(self.lons, len(self.lats))

    def _cell_corner_latitudes(self):
        upper_lats = _interval_bounds(self._OP, self.lats, -self._OP, loc='u')
        lower_lats = _interval_bounds(self._OP, self.lats, -self._OP, loc='l')
        return np.array(
            [
                _col_distribute(upper_lats, len(self.lons)),  # 1 ---- 0
                _col_distribute(upper_lats, len(self.lons)),  # |      |
                _col_distribute(lower_lats, len(self.lons)),  # |      |
                _col_distribute(lower_lats, len(self.lons)),  # 2 ---- 3
            ]
        )

    def _cell_corner_longitudes(self):
        left_lons = _interval_bounds(0, self.lons, 360, loc='l')
        right_lons = _interval_bounds(0, self.lons, 360, loc='r', wrap=True)
        return np.array(
            [
                _row_distribute(right_lons, len(self.lats)),  # 1 ---- 0
                _row_distribute(left_lons, len(self.lats)),   # |      |
                _row_distribute(left_lons, len(self.lats)),   # |      |
                _row_distribute(right_lons, len(self.lats)),  # 2 ---- 3
            ]
        )

    def cell_corners(self):
        return np.array(
            [self._cell_corner_latitudes(), self._cell_corner_longitudes()]
        )

    def cell_areas(self):
        upper_lats = _interval_bounds(self._OP, self.lats, -self._OP, loc='u')
        lower_lats = _interval_bounds(self._OP, self.lats, -self._OP, loc='l')
        return _col_distribute(
            2*np.pi*EARTH_RADIUS**2
            * np.abs(
                np.sin(np.radians(upper_lats))
                - np.sin(np.radians(lower_lats))
            )/len(self.lons),
            len(self.lons)
        )


class RegularLatLonGrid(LatLonGrid):
    def __init__(self, nlats, nlons, first_lat=-90):
        super().__init__(
            lats=_equidistant(first_lat, -first_lat, nlats),
            lons=_equidistant(0, 360, nlons),
            first_lat=first_lat
        )


class FullGaussianGrid(LatLonGrid):
    def __init__(self, lats, first_lat=-90):
        super().__init__(
            lats=lats,
            lons=_equidistant(0, 360, 2*len(lats), first_at_start=True),
            first_lat=first_lat
        )
