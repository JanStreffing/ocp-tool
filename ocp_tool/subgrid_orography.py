"""
Subgrid-scale orography for OpenIFS, computed the way OpenIFS expects it.

OpenIFS carries four subgrid orography fields, all gridpoint on the reduced
Gaussian grid, all consumed by the Lott and Miller (1997) gravity-wave drag:

    sdor  160  standard deviation of the orography inside the grid box   [m]
    isor  161  anisotropy of the subgrid orography                       [0-1]
    anor  162  angle of the principal axis, from east                    [rad]
    slor  163  slope of the subgrid orography                            [-]

They are derived from a fine orography by accumulating, over the fine points
inside each model grid box, the gradient correlation tensor

    K = 1/2 (<hx^2> + <hy^2>)     L = 1/2 (<hx^2> - <hy^2>)     M = <hx hy>

and then

    slor = sqrt(K + sqrt(L^2 + M^2))
    isor = sqrt((K - sqrt(L^2+M^2)) / (K + sqrt(L^2+M^2)))
    anor = 1/2 atan2(M, L)

This module exists because ``paleo_subgrid_oro.py`` cannot be used for
OpenIFS. That path shells out to ``calnoro``, an ECHAM tool implementing a
different definition: fed ECMWF's own 9 km orography it returns a global mean
sdor of 43.3 against ECMWF's 60.4, and a point-by-point RMS difference of 74.8
against a field mean of 60.4. The gap is definitional, not resolution: moving
its input from ECHAM5 T511 (26 km, 2007) to ECMWF TCO1279 (9 km) closed only
about 9% of it. ``paleo_subgrid_oro.py`` stays as it is for ECHAM.

Validation target is what ECMWF ships beside the model, for example
``climate.v020/95_4/{stdgwd,isogwd,anggwd,slogwd}``, which is where the fields
in a stock ICMGG come from.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Optional, Tuple

import numpy as np

# Earth radius used by the IFS (YOMCST RA).
RA = 6371229.0

# GRIB paramIds of the four fields this module produces.
PARAM_SDOR = 160
PARAM_ISOR = 161
PARAM_ANOR = 162
PARAM_SLOR = 163

# Name each field is known by in the ECMWF climate file set, for validation.
CLIMATE_FILES = {
    PARAM_SDOR: "stdgwd",
    PARAM_ISOR: "isogwd",
    PARAM_ANOR: "anggwd",
    PARAM_SLOR: "slogwd",
}


# ---------------------------------------------------------------------------
# Target grid
# ---------------------------------------------------------------------------
@dataclass
class ReducedGaussian:
    """Row structure of a reduced Gaussian grid, read from a GRIB message.

    Reading it from the file rather than deriving it from N keeps this correct
    for octahedral grids, where the row lengths are not the classic table.
    """

    pl: np.ndarray        # points in each latitude row
    row_lat: np.ndarray   # latitude of each row, north to south
    start: np.ndarray     # index of the first point of each row
    size: int             # total number of points

    @property
    def nrows(self) -> int:
        return int(self.pl.size)


def read_target_grid(grib_file: Path) -> ReducedGaussian:
    """Read the reduced Gaussian row structure from the first GRIB message."""
    import eccodes as ec

    with open(grib_file, "rb") as f:
        gid = ec.codes_grib_new_from_file(f)
        if gid is None:
            raise ValueError(f"No GRIB message in {grib_file}")
        try:
            grid_type = ec.codes_get(gid, "gridType")
            if grid_type != "reduced_gg":
                raise ValueError(
                    f"{grib_file} is a {grid_type} grid; the subgrid orography "
                    "fields are gridpoint and this needs a reduced_gg message"
                )
            pl = np.asarray(ec.codes_get_array(gid, "pl"), dtype=np.int64)
            lats = np.asarray(ec.codes_get_array(gid, "latitudes"), dtype=float)
        finally:
            ec.codes_release(gid)

    start = np.concatenate([[0], np.cumsum(pl)[:-1]])
    return ReducedGaussian(pl=pl, row_lat=lats[start], start=start,
                           size=int(pl.sum()))


# ---------------------------------------------------------------------------
# Fine orography
# ---------------------------------------------------------------------------
@dataclass
class FineOrography:
    """A global orography on a regular latitude/longitude grid.

    ETOPO1 is the intended source (``ETOPO1_*_gmt4.grd``, 21600x10800, cell
    registered, 1 arc-minute). Use the **Ice** surface version, not Bed: the Bed
    version carries bedrock under Greenland and Antarctica, which is wrong for
    an atmospheric orography anywhere PISM does not supply the surface itself.
    """

    path: Path
    variable: str = "z"
    lat_name: str = "y"
    lon_name: str = "x"
    band_rows: int = 600

    def open(self):
        import netCDF4 as nc

        ds = nc.Dataset(self.path)
        lat = np.asarray(ds[self.lat_name][:], dtype=float)
        lon = np.asarray(ds[self.lon_name][:], dtype=float)
        if lat[1] < lat[0]:
            raise ValueError(
                f"{self.path} stores latitude descending; this reader assumes "
                "ascending (ETOPO1 gmt4 does)"
            )
        return ds, lat, lon


# A delta is called with (lat_band, lon) and returns metres to add to the fine
# orography on that band. This is how an ice sheet enters: the fine texture of
# the DEM is kept and only the ice surface change is added, so the 1-2 km
# variance that the DEM resolves and the ice model does not is not destroyed.
DeltaFn = Callable[[np.ndarray, np.ndarray], np.ndarray]


# ---------------------------------------------------------------------------
# The computation
# ---------------------------------------------------------------------------
def compute_sso(
    fine: FineOrography,
    grid: ReducedGaussian,
    delta: Optional[DeltaFn] = None,
    sea_level_clip: bool = True,
    verbose: bool = False,
) -> Dict[str, np.ndarray]:
    """Accumulate the Lott and Miller statistics onto ``grid``.

    Returns a dict with ``mean``, ``sdor``, ``isor``, ``anor``, ``slor`` and the
    per-box fine-point count ``n``.

    The fine grid is streamed in latitude bands with a one-row halo, so the
    memory cost is the band, not the DEM: ETOPO1 at the default band size is
    about 100 MB rather than 1.8 GB.
    """
    ds, flat, flon = fine.open()
    try:
        ny, nx = len(flat), len(flon)
        dlon = np.radians(360.0 / nx)
        dy = RA * np.radians(180.0 / ny)

        acc = {k: np.zeros(grid.size) for k in
               ("n", "h", "h2", "gx2", "gy2", "gxy")}

        for j0 in range(0, ny, fine.band_rows):
            j1 = min(j0 + fine.band_rows, ny)
            # one-row halo so the meridional gradient is centred everywhere
            a, b = max(j0 - 1, 0), min(j1 + 1, ny)
            h = np.asarray(ds[fine.variable][a:b, :], dtype=np.float64)
            if sea_level_clip:
                # The IFS orography is zero over sea, not bathymetry.
                h = np.maximum(h, 0.0)
            if delta is not None:
                h = h + delta(flat[a:b], flon)

            # Zonal spacing shrinks as cos(lat) and reaches zero at the pole;
            # clip so the polar rows cannot produce an infinite gradient.
            dx = (RA * np.cos(np.radians(flat[a:b]))[:, None] * dlon).clip(1.0)
            gx = (np.roll(h, -1, axis=1) - np.roll(h, 1, axis=1)) / (2.0 * dx)
            gy = np.empty_like(h)
            gy[1:-1] = (h[2:] - h[:-2]) / (2.0 * dy)
            gy[0], gy[-1] = gy[1], gy[-2]

            lo = j0 - a
            hi = lo + (j1 - j0)
            hb, gxb, gyb = h[lo:hi], gx[lo:hi], gy[lo:hi]

            rows = np.abs(
                grid.row_lat[None, :] - flat[j0:j1][:, None]
            ).argmin(axis=1)

            for k, r in enumerate(rows):
                col = np.floor(
                    (flon % 360.0) / (360.0 / grid.pl[r])
                ).astype(np.int64) % grid.pl[r]
                idx = grid.start[r] + col
                acc["n"] += np.bincount(idx, minlength=grid.size)
                acc["h"] += np.bincount(idx, weights=hb[k], minlength=grid.size)
                acc["h2"] += np.bincount(idx, weights=hb[k] ** 2,
                                         minlength=grid.size)
                acc["gx2"] += np.bincount(idx, weights=gxb[k] ** 2,
                                          minlength=grid.size)
                acc["gy2"] += np.bincount(idx, weights=gyb[k] ** 2,
                                          minlength=grid.size)
                acc["gxy"] += np.bincount(idx, weights=gxb[k] * gyb[k],
                                          minlength=grid.size)
            if verbose:
                print(f"    rows {j0}-{j1} of {ny}")
    finally:
        ds.close()

    empty = int((acc["n"] == 0).sum())
    if empty:
        raise RuntimeError(
            f"{empty} of {grid.size} target boxes received no fine point. The "
            "fine orography is too coarse for this target grid."
        )

    n = acc["n"]
    mean = acc["h"] / n
    sdor = np.sqrt(np.maximum(acc["h2"] / n - mean ** 2, 0.0))

    gx2, gy2, gxy = acc["gx2"] / n, acc["gy2"] / n, acc["gxy"] / n
    K = 0.5 * (gx2 + gy2)
    L = 0.5 * (gx2 - gy2)
    M = gxy
    R = np.sqrt(L ** 2 + M ** 2)

    slor = np.sqrt(np.maximum(K + R, 0.0))
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(K + R > 0.0, (K - R) / (K + R), 0.0)
    isor = np.sqrt(np.clip(ratio, 0.0, 1.0))
    anor = 0.5 * np.arctan2(M, L)

    return {"mean": mean, "sdor": sdor, "isor": isor, "anor": anor,
            "slor": slor, "n": n}


# ---------------------------------------------------------------------------
# Validation against what ECMWF ships
# ---------------------------------------------------------------------------
def read_climate_field(climate_dir: Path, param: int) -> np.ndarray:
    """Read one ECMWF climate field, e.g. ``climate.v020/95_4/stdgwd``."""
    import eccodes as ec

    path = Path(climate_dir) / CLIMATE_FILES[param]
    with open(path, "rb") as f:
        gid = ec.codes_grib_new_from_file(f)
        if gid is None:
            raise ValueError(f"No GRIB message in {path}")
        try:
            return np.asarray(ec.codes_get_values(gid), dtype=float)
        finally:
            ec.codes_release(gid)


def validate(
    computed: Dict[str, np.ndarray],
    climate_dir: Path,
    land_mask: Optional[np.ndarray] = None,
) -> Dict[str, Dict[str, float]]:
    """Compare computed fields against the ECMWF climate set.

    Reports the mean ratio, the correlation and the least-squares slope over
    the cells selected by ``land_mask``. Correlation is the number that matters:
    a ratio away from one is a calibration question, a low correlation means the
    formulation is wrong.
    """
    names = {PARAM_SDOR: "sdor", PARAM_ISOR: "isor",
             PARAM_ANOR: "anor", PARAM_SLOR: "slor"}
    if land_mask is None:
        land_mask = np.ones(computed["sdor"].size, dtype=bool)

    report = {}
    for param, key in names.items():
        mine = computed[key][land_mask]
        ref = read_climate_field(climate_dir, param)[land_mask]
        ok = np.isfinite(mine) & np.isfinite(ref)
        if param != PARAM_ANOR:
            # anor is a circular angle on (-pi/2, pi/2]; a ratio or a linear
            # correlation of it is not meaningful, so only its spread is shown.
            ok &= ref > 0
        mine, ref = mine[ok], ref[ok]
        entry = {"n": int(ok.sum()),
                 "mean_mine": float(mine.mean()),
                 "mean_ref": float(ref.mean())}
        if param != PARAM_ANOR:
            entry["ratio"] = float(mine.mean() / ref.mean())
            entry["corr"] = float(np.corrcoef(mine, ref)[0, 1])
            entry["slope"] = float(np.sum(mine * ref) / np.sum(ref * ref))
        report[key] = entry
    return report


def format_report(report: Dict[str, Dict[str, float]]) -> str:
    lines = [f"{'field':6s} {'mine':>10s} {'ECMWF':>10s} {'ratio':>7s} "
             f"{'corr':>7s} {'slope':>7s}  cells"]
    for key, e in report.items():
        if "corr" in e:
            lines.append(f"{key:6s} {e['mean_mine']:10.5g} {e['mean_ref']:10.5g} "
                         f"{e['ratio']:7.3f} {e['corr']:7.3f} {e['slope']:7.3f}  "
                         f"{e['n']}")
        else:
            lines.append(f"{key:6s} {e['mean_mine']:10.5g} {e['mean_ref']:10.5g} "
                         f"{'-':>7s} {'-':>7s} {'-':>7s}  {e['n']}  (circular)")
    return "\n".join(lines)
