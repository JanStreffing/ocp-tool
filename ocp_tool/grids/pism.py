import numpy as np
import pyproj
from netCDF4 import Dataset


class PISM:
    """A PISM ice-sheet domain as an OASIS grid.

    Centres come from the file's own lat/lon; corners are inverse-projected
    from the cell corner (x, y). The grid is fixed for the whole run, only the
    ice inside it moves.
    """

    def __init__(self, grid_file, name="ismp"):
        if len(name) != 4:
            raise ValueError(f"OASIS grid name must be 4 characters: {name!r}")
        self.grid_file = grid_file
        self.name = name

        with Dataset(grid_file) as nc:
            if not {"x", "y"}.issubset(nc.dimensions):
                raise RuntimeError(f"Missing x/y dimensions in PISM file {grid_file}")
            if not {"lat", "lon", "mapping"}.issubset(nc.variables):
                raise RuntimeError(
                    f"Missing lat/lon/mapping variable(s) in PISM file {grid_file}"
                )
            mapping = {k: nc["mapping"].getncattr(k) for k in nc["mapping"].ncattrs()}

        gm = mapping.get("grid_mapping_name", "").lower()
        if gm != "polar_stereographic":
            raise RuntimeError(f"Unsupported PISM projection: {gm!r}")

        self.proj4 = (
            f"+proj=stere "
            f"+lat_0={mapping['latitude_of_projection_origin']} "
            f"+lat_ts={mapping['standard_parallel']} "
            f"+lon_0={mapping.get('straight_vertical_longitude_from_pole', 0.0)} "
            f"+x_0={mapping.get('false_easting', 0.0)} "
            f"+y_0={mapping.get('false_northing', 0.0)} "
            f"+ellps={mapping.get('ellipsoid', 'WGS84')} +units=m +no_defs"
        )

    def _xy(self):
        with Dataset(self.grid_file) as nc:
            return (np.array(nc["x"][:], dtype=float),
                    np.array(nc["y"][:], dtype=float))

    def cell_latitudes(self):
        with Dataset(self.grid_file) as nc:
            return np.array(nc["lat"][:], dtype=float)

    def cell_longitudes(self):
        with Dataset(self.grid_file) as nc:
            return np.array(nc["lon"][:], dtype=float)

    def cell_corners(self):
        """Corner lat/lon as (2, 4, ny, nx), ORCA convention (0=lat, 1=lon).

        Counterclockwise seen from outside the sphere, as SCRIP wants.
        """
        x, y = self._xy()
        dx = float(np.diff(x).mean())
        dy = float(np.diff(y).mean())
        xc = np.concatenate([x - dx / 2.0, [x[-1] + dx / 2.0]])
        yc = np.concatenate([y - dy / 2.0, [y[-1] + dy / 2.0]])
        XC, YC = np.meshgrid(xc, yc)

        transformer = pyproj.Transformer.from_crs(
            pyproj.CRS.from_proj4(self.proj4), "EPSG:4326", always_xy=True
        )
        lonc, latc = transformer.transform(XC, YC)

        corners = np.empty((2, 4, *self.cell_latitudes().shape))
        corners[0] = np.stack([latc[:-1, :-1], latc[:-1, 1:],
                               latc[1:, 1:], latc[1:, :-1]])
        corners[1] = np.stack([lonc[:-1, :-1], lonc[:-1, 1:],
                               lonc[1:, 1:], lonc[1:, :-1]])
        return corners

    def cell_areas(self, earth_radius=6371229.0):
        """Cell area [m2], spherical excess over two triangles.

        Not a plane formula: the domain wraps the pole, so cells on the +/-180
        meridian would get a lon span of ~2*pi instead of ~0.
        """
        corners = self.cell_corners()
        lat, lon = corners[0], corners[1]

        def unit(k):
            la, lo = np.radians(lat[k]), np.radians(lon[k])
            return np.stack([np.cos(la) * np.cos(lo),
                             np.cos(la) * np.sin(lo),
                             np.sin(la)], axis=0)

        v = [unit(k) for k in range(4)]

        def triangle(a, b, c):
            num = np.abs(np.einsum("i...,i...->...", a, np.cross(b, c, axis=0)))
            den = (1.0
                   + np.einsum("i...,i...->...", a, b)
                   + np.einsum("i...,i...->...", b, c)
                   + np.einsum("i...,i...->...", c, a))
            return 2.0 * np.arctan2(num, den) * earth_radius * earth_radius

        return triangle(v[0], v[1], v[2]) + triangle(v[0], v[2], v[3])

    def cell_masks(self):
        """OASIS mask, 1 = excluded, 0 = active. All active.

        Masking by ice cover would make the weights follow the ice geometry,
        which defeats having a fixed grid.
        """
        return np.zeros(self.cell_latitudes().shape, dtype=np.int32)
