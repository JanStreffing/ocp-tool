"""
Subgrid-scale orography for OpenIFS, computed the way OpenIFS expects it.

OpenIFS carries four subgrid orography fields, all gridpoint on the reduced
Gaussian grid, all consumed by the Lott and Miller (1997) gravity-wave drag:

    sdor  160  standard deviation of the orography inside the grid box   [m]
    isor  161  anisotropy of the subgrid orography                       [0-1]
    anor  162  angle of the principal axis, from east                    [rad]
    slor  163  slope of the subgrid orography                            [-]

ECMWF build them to represent scales between 5 km and the grid length, in three
steps (IFS orographic drag documentation, and Elvidge et al., ECMWF Tech Memo
893):

1. the 30 arc-second source is averaged to 2'30", about 5 km, because scales
   below that belong to TOFD and the effective roughness scheme, not here;
2. the grid-scale orography is subtracted from it, leaving a band-passed
   residual that carries only 5 km to grid length;
3. the statistics are accumulated from that residual, not from the orography
   relative to its box mean.

Skipping step 2 leaves the near-grid-scale ramp across each box inside the
gradient tensor. That barely moves the standard deviation or the slope but it
wrecks the anisotropy, which is a ratio of the tensor's eigenvalues and so is
ruined by a term common to both.

The statistics are the gradient correlation tensor of the residual

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
def interpolate_to_fine(values: np.ndarray, grid: ReducedGaussian,
                        flat: np.ndarray, flon: np.ndarray) -> np.ndarray:
    """Bilinearly interpolate a reduced Gaussian field onto a regular lat/lon grid.

    Bilinear and not nearest: a piecewise-constant grid-scale orography would put
    a cliff at every box edge, and those cliffs would then dominate the very
    gradient statistics this is subtracted from.
    """
    rl = grid.row_lat                      # north to south
    jj = np.clip(np.searchsorted(-rl, -flat), 1, grid.nrows - 1)
    j0, j1 = jj - 1, jj
    wlat = (rl[j0] - flat) / (rl[j0] - rl[j1])
    wlat = np.clip(wlat, 0.0, 1.0)

    def row_interp(j):
        step = (360.0 / grid.pl[j])[:, None]
        x = (flon[None, :] % 360.0) / step
        i0 = np.floor(x).astype(np.int64) % grid.pl[j][:, None]
        i1 = (i0 + 1) % grid.pl[j][:, None]
        w = x - np.floor(x)
        base = grid.start[j][:, None]
        return values[base + i0] * (1 - w) + values[base + i1] * w

    return row_interp(j0) * (1 - wlat)[:, None] + row_interp(j1) * wlat[:, None]


def compute_sso(
    fine: FineOrography,
    grid: ReducedGaussian,
    delta: Optional[DeltaFn] = None,
    sea_level_clip: bool = True,
    remove_mean_slope: bool = False,
    gradient_stencil: int = 1,
    prefilter_arcmin: float = 3.0,
    highpass_km: Optional[float] = None,   # None means 1.5 * grid spacing
    grid_scale_orography: Optional[np.ndarray] = None,
    verbose: bool = False,
) -> Dict[str, np.ndarray]:
    """Accumulate the Lott and Miller statistics onto ``grid``.

    Returns a dict with ``mean``, ``sdor``, ``isor``, ``anor``, ``slor`` and the
    per-box fine-point count ``n``.

    The fine grid is streamed in latitude bands, so the memory cost is the band
    and not the DEM: ETOPO1 at the default band size is about 100 MB rather than
    1.8 GB.

    ``sdor`` is the spread of the prefiltered orography about its box mean.
    ``slor``, ``isor`` and ``anor`` come from the band-passed residual: taking
    ``sdor`` from the residual too costs accuracy, ratio 0.68 against 0.78,
    for a correlation gain of 0.02.

    ``prefilter_arcmin`` is a physical choice and not a free parameter: at
    TCO319 refining it from 3' to 1' takes the slope ratio from 0.93 to 1.73,
    because scales below about 5 km belong to TOFD rather than here.

    Behaviour across resolution, TCO79 to TCO319, against climate.v020:
    sdor holds at 0.76-0.80 and slor at 0.91-0.93, both with correlations near
    0.85 and 0.90. isor is flat at 0.52 for TCO79 and TCO95 and then falls to
    0.39 by TCO319, for reasons not understood; it is not sample starvation,
    since nine times the fine points per box does not move it.
    """
    ds, flat, flon = fine.open()
    try:
        ny, nx = len(flat), len(flon)
        src_arcmin = 360.0 * 60.0 / nx
        blk = max(int(round(prefilter_arcmin / src_arcmin)), 1)
        if verbose:
            print(f"    source {src_arcmin:.2f}', prefilter block {blk} "
                  f"({blk * src_arcmin:.2f}')")
        # Step 1: average the source to about 5 km. Trim to a whole number of
        # blocks rather than padding, so no block is a partial average.
        ny2, nx2 = (ny // blk) * blk, (nx // blk) * blk
        h = np.asarray(ds[fine.variable][:ny2, :nx2], dtype=np.float32)
    finally:
        ds.close()

    if sea_level_clip:
        # The IFS orography is zero over sea, not bathymetry.
        h = np.maximum(h, 0.0)
    flat = flat[:ny2].reshape(-1, blk).mean(axis=1)
    flon = flon[:nx2].reshape(-1, blk).mean(axis=1)
    h = h.reshape(ny2 // blk, blk, nx2 // blk, blk).mean(axis=(1, 3))
    ny, nx = h.shape

    # The block average above makes the cell about 5 km north to south, but a
    # regular lat/lon cell shrinks as cos(lat) east to west: 5.5 km at the
    # equator is 1 km at 80S. Left alone that samples the gradient tensor far
    # more finely in longitude than in latitude, biasing the anisotropy toward
    # zonal structure in a pattern organised radially about the pole. It shows
    # as radial spokes in isor over Antarctica, which is precisely where an ice
    # sheet coupling needs the field. So bring the zonal resolution back to the
    # meridional one, per row.
    from scipy.ndimage import uniform_filter1d

    dy_km = 180.0 / ny * 111.0
    dx_km = 360.0 / nx * 111.0 * np.cos(np.radians(flat))
    zonal_width = np.clip(
        np.round(dy_km / np.maximum(dx_km, 1e-9)).astype(int), 1, max(nx // 4, 1))
    for w in np.unique(zonal_width):
        if w > 1:
            rows_w = zonal_width == w
            h[rows_w] = uniform_filter1d(h[rows_w], int(w), axis=1, mode="wrap")
    if verbose:
        print(f"    zonal prefilter width 1 at the equator, "
              f"{zonal_width.max()} at the polar row")

    if delta is not None:
        h = h + delta(flat, flon).astype(np.float32)

    # Box mean and standard deviation come from the prefiltered orography
    # itself: they describe the box, not the band-passed residual.
    rows = np.abs(grid.row_lat[None, :] - flat[:, None]).argmin(axis=1)
    cols = [np.floor((flon % 360.0) / (360.0 / grid.pl[r])).astype(np.int64)
            % grid.pl[r] for r in range(grid.nrows)]
    idx_of = lambda r: grid.start[r] + cols[r]

    # Weight by cos(lat): a regular lat/lon cell carries less area the further
    # poleward it is, so an unweighted box mean over-represents its poleward
    # edge. Small inside one box, but free to get right.
    wrow = np.cos(np.radians(flat)).clip(1e-9)

    acc = {k: np.zeros(grid.size) for k in
           ("n", "h", "h2", "gx", "gy", "gx2", "gy2", "gxy")}
    for k, r in enumerate(rows):
        idx = idx_of(r)
        acc["n"] += np.bincount(idx, weights=np.full(nx, wrow[k]),
                                minlength=grid.size)
        acc["h"] += np.bincount(idx, weights=h[k] * wrow[k], minlength=grid.size)
        acc["h2"] += np.bincount(idx, weights=h[k].astype(np.float64) ** 2 * wrow[k],
                                 minlength=grid.size)

    if (acc["n"] <= 0).any():
        raise RuntimeError(
            f"{int((acc['n'] == 0).sum())} of {grid.size} target boxes received "
            "no fine point. The prefiltered orography is too coarse for this "
            "target grid; lower prefilter_arcmin."
        )
    n = acc["n"]
    mean = acc["h"] / n
    sdor = np.sqrt(np.maximum(acc["h2"] / n - mean ** 2, 0.0))

    # Step 2: band-pass, so the tensor sees only 5 km to grid length.
    #
    # Two ways to remove the large scales. The IFS uses a smoothed version of
    # the prefiltered source, which is what highpass_km does; subtracting the
    # grid-scale orography interpolated back onto the fine grid is the other
    # form given in the same description. The smoothed version wins on the
    # anisotropy, 0.51 against 0.45, because it has no box-edge structure of
    # its own for the gradient to pick up.
    if highpass_km is None:
        # The documentation calls this a "1-dx filter", so the cutoff is the
        # grid length and is derived, not fitted. Confirmed across a fourfold
        # resolution range: held at a fixed 150 km the slor ratio drifts from
        # 0.910 at TCO79 to 0.955 at TCO319, while scaled with the grid it
        # stays inside 0.911 to 0.927.
        highpass_km = 1.5 * (360.0 / int(grid.pl.max())) * 111.0
    if highpass_km > 0:
        from scipy.ndimage import gaussian_filter1d

        row_km = 180.0 / ny * 111.0
        sm = gaussian_filter1d(h, highpass_km / row_km, axis=0, mode="nearest")
        # Zonal spacing shrinks as cos(lat), so a single sigma in grid points
        # would smooth over far too short a distance near the poles, which is
        # exactly where the ice sheets are. Apply it in latitude bands of a few
        # degrees, each with its own sigma.
        band = max(int(round(3.0 / (180.0 / ny))), 1)
        for b0 in range(0, ny, band):
            b1 = min(b0 + band, ny)
            coslat = max(np.cos(np.radians(flat[b0:b1])).mean(), 1e-3)
            sig = highpass_km / (360.0 / nx * 111.0 * coslat)
            sm[b0:b1] = gaussian_filter1d(sm[b0:b1], sig, axis=1, mode="wrap")
        resid = h - sm
    else:
        gso = mean if grid_scale_orography is None else np.asarray(
            grid_scale_orography, dtype=float)
        resid = h - interpolate_to_fine(gso, grid, flat, flon).astype(np.float32)

    # Step 3: the gradient correlation tensor of the residual.
    st = max(int(gradient_stencil), 1)
    dlon = np.radians(360.0 / nx)
    dy = RA * np.radians(180.0 / ny)
    dx = (RA * np.cos(np.radians(flat))[:, None] * dlon).clip(1.0)

    # Same reason as the zonal prefilter: the zonal difference has to span the
    # same physical distance as the meridional one, or the tensor is measuring
    # two different scales in the two directions and the anisotropy is an
    # artefact of the grid.
    gx = np.empty_like(resid)
    for w in np.unique(zonal_width):
        rows_w = zonal_width == w
        sx = max(int(w) * st, 1)
        block = resid[rows_w]
        gx[rows_w] = ((np.roll(block, -sx, axis=1) - np.roll(block, sx, axis=1))
                      / (2.0 * sx * dx[rows_w]))
    gy = np.empty_like(resid)
    gy[st:-st] = (resid[2 * st:] - resid[:-2 * st]) / (2.0 * st * dy)
    gy[:st], gy[-st:] = gy[st], gy[-st - 1]

    for k, r in enumerate(rows):
        idx = idx_of(r)
        acc["gx"] += np.bincount(idx, weights=gx[k] * wrow[k], minlength=grid.size)
        acc["gy"] += np.bincount(idx, weights=gy[k] * wrow[k], minlength=grid.size)
        acc["gx2"] += np.bincount(idx, weights=gx[k].astype(np.float64) ** 2 * wrow[k],
                                  minlength=grid.size)
        acc["gy2"] += np.bincount(idx, weights=gy[k].astype(np.float64) ** 2 * wrow[k],
                                  minlength=grid.size)
        acc["gxy"] += np.bincount(idx,
                                  weights=(gx[k] * gy[k]).astype(np.float64) * wrow[k],
                                  minlength=grid.size)

    gx2, gy2, gxy = acc["gx2"] / n, acc["gy2"] / n, acc["gxy"] / n
    if remove_mean_slope:
        gxm, gym = acc["gx"] / n, acc["gy"] / n
        gx2 = gx2 - gxm ** 2
        gy2 = gy2 - gym ** 2
        gxy = gxy - gxm * gym
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
