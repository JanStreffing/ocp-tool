#!/usr/bin/env python3
"""Remove the dateline seam of false land from a published ICMGG, in place.

THE DEFECT.  Before the polar-cap guard in read_fesom_grid_polygon, a single
triangle spanning more than 180 degrees after the dateline shift made the
trapezoid map unbuildable, the whole dateline pass was lost, and every point it
would have classified kept the "land" default.  On TCO95/CORE2 that put all 192
cells within half a degree of lon 180 on land against 9 % in the ECMWF input --
a seam of false land from pole to pole.  A bare except turned it into a warning,
so the files were written and published.

WHAT THIS FIXES, AND WHAT IT DELIBERATELY DOES NOT.  Only cells in the dateline
band whose land-sea value the corrected mask disagrees with.  Recomputing the
mask today also disagrees in a handful of places away from the dateline -- 1
cell on TCO95/CORE2, 36 on TL255/CORE2 -- because the published files do not all
date from the current mesh.nc.  Those are coastline drift, not this defect, and
folding them in would change a coastline silently under cover of a bug fix.  The
band is the fix; everything else is left exactly as published.

A cell that becomes ocean cannot keep a land surface column: it would carry soil
moisture, vegetation and a land skin temperature into a sea point, which is the
inconsistency that produces NaNs in the surface scheme.  So each repaired cell
is rebuilt from the nearest stable OCEAN neighbour, the same great-circle
nearest-neighbour rule ocp-tool uses for the flip itself.  Donors exclude the
cells being repaired, so a repaired cell is never seeded from another.

Usage:
    repair_dateline_mask.py PUBLISHED CORRECTED_MASK.npy OUTPUT [--dry-run]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import eccodes

# Half-width of the seam, in degrees of longitude. The defect is confined to the
# single column of cells centred exactly on 180; 0.5 covers it with no room to
# catch a genuine coastal cell.
DATELINE_HALFWIDTH_DEG = 0.5


def read_all(path):
    """[(shortName, values)], plus lat/lon off the lsm message."""
    out, lat, lon = [], None, None
    with open(path, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                sn = eccodes.codes_get(gid, 'shortName')
                out.append((sn, eccodes.codes_get_values(gid)))
                if sn == 'lsm' and lat is None:
                    lat = eccodes.codes_get_array(gid, 'latitudes')
                    lon = eccodes.codes_get_array(gid, 'longitudes')
            finally:
                eccodes.codes_release(gid)
    return out, lat, lon


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('published')
    ap.add_argument('corrected_mask')
    ap.add_argument('output')
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    fields, lat, lon = read_all(a.published)
    names = [sn for sn, _ in fields]
    if 'lsm' not in names:
        sys.exit('  no lsm in the published file')
    k_lsm = names.index('lsm')
    pub = np.array(fields[k_lsm][1], copy=True)

    fixed = np.load(a.corrected_mask)
    if len(fixed) != len(pub):
        sys.exit(f'  mask has {len(fixed)} cells, file has {len(pub)}')

    lon180 = np.where(lon > 180, lon - 360, lon)
    band = np.abs(np.abs(lon180) - 180) < DATELINE_HALFWIDTH_DEG
    wrong = band & (np.round(pub) != np.round(fixed))
    to_ocean = np.where(wrong & (fixed < 0.5))[0]
    to_land = np.where(wrong & (fixed >= 0.5))[0]

    print(f'  {a.published}')
    print(f'    cells in the dateline band          : {int(band.sum())}')
    print(f'    of those, land in the published file: {int((pub[band] >= 0.5).sum())}')
    print(f'    of those, land in the corrected mask: {int((fixed[band] >= 0.5).sum())}')
    print(f'    to repair: {len(to_ocean)} -> ocean, {len(to_land)} -> land')
    outside = int(((np.round(pub) != np.round(fixed)) & ~band).sum())
    print(f'    left alone, outside the band        : {outside} '
          f'(coastline drift, not this defect)')
    if len(to_ocean) == 0 and len(to_land) == 0:
        print('    nothing to do')
        return

    new = [np.array(v, copy=True) for _, v in fields]
    new[k_lsm][to_ocean] = 0.0
    new[k_lsm][to_land] = 1.0

    # Rebuild each repaired cell from the nearest stable neighbour of its NEW
    # type. Unit vectors on the sphere: nearest = largest dot product.
    rlat, rlon = np.deg2rad(np.asarray(lat, float)), np.deg2rad(np.asarray(lon, float))
    vx, vy, vz = (np.cos(rlat) * np.cos(rlon), np.cos(rlat) * np.sin(rlon),
                  np.sin(rlat))
    npts = len(pub)
    repaired = np.zeros(npts, bool)
    repaired[to_ocean] = True
    repaired[to_land] = True
    is_land = new[k_lsm] >= 0.5
    donor_ocean = np.where(~is_land & ~repaired)[0]
    donor_land = np.where(is_land & ~repaired)[0]

    def seed(cells, donors, label):
        if len(cells) and donors.size == 0:
            sys.exit(f'  no stable {label} donors available')
        for i in cells:
            dots = vx[donors] * vx[i] + vy[donors] * vy[i] + vz[donors] * vz[i]
            j = int(donors[int(np.argmax(dots))])
            for f in range(len(new)):
                if f == k_lsm or len(new[f]) != npts:
                    continue
                new[f][i] = new[f][j]

    seed(to_ocean, donor_ocean, 'ocean')
    seed(to_land, donor_land, 'land')

    # ---------------------------------------------------------------- verify
    touched = set(to_ocean.tolist()) | set(to_land.tolist())
    bad = []
    for k, name in enumerate(names):
        ch = np.where(fields[k][1] != new[k])[0]
        stray = [int(c) for c in ch if int(c) not in touched]
        if stray:
            bad.append(f'{name}: {len(stray)} cells changed outside the repair set')
    nb = int((new[k_lsm][band] >= 0.5).sum())
    print(f'    after: dateline band land = {nb}/{int(band.sum())} '
          f'({100 * nb / band.sum():.0f} %)')
    print(f'    land total {int((pub >= 0.5).sum())} -> {int((new[k_lsm] >= 0.5).sum())}')
    print(f'    changes confined to the repair set: {not bad}')
    if bad:
        for b in bad:
            print(f'      FAIL: {b}')
        sys.exit('  verification failed; nothing written')

    if a.dry_run:
        print('    --dry-run, nothing written')
        return

    out = Path(a.output)
    tmp = out.with_suffix(out.suffix + '.tmp')
    with open(a.published, 'rb') as fin, open(tmp, 'wb') as fout:
        k = 0
        while True:
            gid = eccodes.codes_grib_new_from_file(fin)
            if gid is None:
                break
            try:
                eccodes.codes_set_values(gid, new[k])
                eccodes.codes_write(gid, fout)
            finally:
                eccodes.codes_release(gid)
            k += 1
    tmp.replace(out)
    check, _, _ = read_all(out)
    assert [sn for sn, _ in check] == names, 'message order changed'
    cl_lsm = check[k_lsm][1]
    print(f'    wrote {out}  (read back: land {int((cl_lsm >= 0.5).sum())}, '
          f'dateline {int((cl_lsm[band] >= 0.5).sum())}/{int(band.sum())})')


if __name__ == '__main__':
    main()
