#!/usr/bin/env python3
"""Put the lake fields back into a published ICMGG, in place.

WHY NOT JUST RE-RUN OCP-TOOL.  Regenerating an ICMGG from mesh.nc recomputes
the land-sea mask, and the published files do not all date from the current
meshes -- for DARS2 a regeneration moved the coastline by 1989 cells.  Moving a
coastline is a far larger change than the defect being fixed, and it would also
sweep in every unrelated ocp-tool change since the file was published.  This is
the same argument that produced the _v2 files by repair rather than by
regeneration, and it applies unchanged here.

WHAT THIS DOES INSTEAD.  It never recomputes the mask.  The flipped cells are
recovered by comparing the published file's ``lsm`` against the pristine
input's -- a cell that was ocean and is now land is a cell the flip touched --
so no mesh, no polygon test, and zero coastline risk by construction.  At those
cells and nowhere else it restores the pristine ``cl``, ``dl`` and FLake
prognostics, which ``fill_flipped_from_nearest_neighbour`` had overwritten from
the nearest dry land neighbour.

WHAT IS DELIBERATELY NOT TOUCHED.  ``lsm`` and ``slt``.  The mask must stay
byte-identical for the file to be a drop-in, and the soil repair already in the
_v2 files is a separate, independently validated change that this must not
disturb.  Because ``slt`` is unchanged, the derived LPJ-GUESS soil file
(``grib_copy -w shortName=slt``) is unchanged too and does NOT need
republishing -- only the ICMGG does.

The result is verified before it is written: every field other than the lake
set must come out bit-identical, and the lake set must differ only at flipped
cells.  A failure of either check aborts rather than producing a file.

Usage:
    repair_lake_cover.py PRISTINE PUBLISHED OUTPUT [--dry-run]

    PRISTINE   the untouched ICMGG the published file was built from
    PUBLISHED  the _v2 file to repair
    OUTPUT     the file to write (conventionally _v3)
"""
import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import eccodes

# FLake integrates MIN(50, MAX(2, dl)); flakeene_mod.F90:211-212,238 and
# surftstp_ctl_mod.F90:882. Writing a depth outside that band claims a
# behaviour the model will not produce.
FLAKE_DEPTH_MIN_M, FLAKE_DEPTH_MAX_M = 2.0, 50.0

# The fields the nearest-neighbour fill destroys and this restores.
LAKE_FIELDS = ('cl', 'dl', 'lmlt', 'lblt', 'lmld', 'lict', 'licd', 'ltlt',
               'lshf')
# Must come through untouched: the mask so the file stays a drop-in, the soil
# type so the _v2 repair is preserved.
PROTECTED = ('lsm', 'slt')


def read_all(path):
    """[(shortName, values)] in file order, plus the coordinates."""
    out, lon = [], None
    with open(path, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                break
            try:
                sn = eccodes.codes_get(gid, 'shortName')
                out.append((sn, eccodes.codes_get_values(gid)))
                if sn == 'lsm' and lon is None:
                    lon = eccodes.codes_get_array(gid, 'longitudes')
            finally:
                eccodes.codes_release(gid)
    return out, lon


def first(fields, name):
    for sn, v in fields:
        if sn == name:
            return v
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('pristine')
    ap.add_argument('published')
    ap.add_argument('output')
    ap.add_argument('--dry-run', action='store_true',
                    help='report what would change and write nothing')
    a = ap.parse_args()

    pristine, _ = read_all(a.pristine)
    published, lon = read_all(a.published)
    print(f'  pristine  {a.pristine}  ({len(pristine)} messages)')
    print(f'  published {a.published} ({len(published)} messages)')

    p_lsm, q_lsm = first(pristine, 'lsm'), first(published, 'lsm')
    if p_lsm is None or q_lsm is None:
        sys.exit('  no lsm in one of the files')
    if len(p_lsm) != len(q_lsm):
        sys.exit(f'  grid mismatch: {len(p_lsm)} vs {len(q_lsm)} cells')

    # The flip, recovered from the two masks. No mesh is read.
    flipped = np.where((p_lsm < 0.5) & (q_lsm >= 0.5))[0]
    print(f'\n  cells flipped ocean -> land: {len(flipped)}')
    if len(flipped) == 0:
        sys.exit('  nothing was flipped; this file needs no lake repair')

    pmap = {}
    for i, (sn, v) in enumerate(pristine):
        pmap.setdefault(sn, v)

    # Build the repaired values.
    new = [np.array(v, copy=True) for _, v in published]
    names = [sn for sn, _ in published]
    touched = {}
    for name in LAKE_FIELDS:
        if name not in names or name not in pmap:
            print(f'    {name:6s} absent, skipped')
            continue
        k = names.index(name)
        src = pmap[name]
        if len(src) != len(new[k]):
            sys.exit(f'  {name}: length mismatch between the two files')
        before = np.array(new[k], copy=True)
        new[k][flipped] = src[flipped]
        if name == 'dl':
            new[k][flipped] = np.clip(new[k][flipped],
                                      FLAKE_DEPTH_MIN_M, FLAKE_DEPTH_MAX_M)
        n = int((before != new[k]).sum())
        touched[name] = n
        print(f'    {name:6s} restored at {n:5d} cells')

    cl = first(published, 'cl')
    print(f'\n  cl at the flipped cells: '
          f'{cl[flipped].mean():.4f} -> {new[names.index("cl")][flipped].mean():.4f}'
          f'   (majority-lake cells: '
          f'{int((new[names.index("cl")][flipped] >= 0.5).sum())})')

    # ---------------------------------------------------------------- verify
    print('\n  VERIFICATION')
    bad = []
    for k, name in enumerate(names):
        changed = np.where(published[k][1] != new[k])[0]
        if name in PROTECTED and changed.size:
            bad.append(f'{name} changed at {changed.size} cells (must not)')
        elif name not in LAKE_FIELDS and changed.size:
            bad.append(f'{name} changed at {changed.size} cells (not a lake field)')
        elif changed.size and not np.isin(changed, flipped).all():
            outside = int((~np.isin(changed, flipped)).sum())
            bad.append(f'{name} changed at {outside} cells outside the flip')
    for name in PROTECTED:
        if name in names:
            k = names.index(name)
            ok = np.array_equal(published[k][1], new[k])
            print(f'    {name:6s} bit-identical to the published file: {ok}')
    print(f'    changes confined to the {len(flipped)} flipped cells: '
          f'{not bad}')
    if bad:
        for b in bad:
            print(f'      FAIL: {b}')
        sys.exit('\n  verification failed; nothing written')

    if a.dry_run:
        print('\n  --dry-run, nothing written')
        return

    # ------------------------------------------------------------ write out
    # Copy the published file and overwrite values message by message, so every
    # GRIB header, packing choice and message order is inherited exactly.
    out = Path(a.output)
    shutil.copy2(a.published, out)
    tmp = out.with_suffix(out.suffix + '.tmp')
    with open(a.published, 'rb') as fin, open(tmp, 'wb') as fout:
        k = 0
        while True:
            gid = eccodes.codes_grib_new_from_file(fin)
            if gid is None:
                break
            try:
                if names[k] in touched:
                    eccodes.codes_set_values(gid, new[k])
                eccodes.codes_write(gid, fout)
            finally:
                eccodes.codes_release(gid)
            k += 1
    tmp.replace(out)
    print(f'\n  wrote {out}')

    # Read it back: the file on disk is the thing that matters, not the arrays.
    check, _ = read_all(out)
    cnames = [sn for sn, _ in check]
    assert cnames == names, 'message order changed'
    for k, name in enumerate(names):
        if name in PROTECTED:
            assert np.array_equal(check[k][1], published[k][1]), \
                f'{name} not preserved on disk'
    ncl = first(check, 'cl')
    print(f'  read back: cl >= 0.5 at {int((ncl >= 0.5).sum())} cells '
          f'(published had {int((cl >= 0.5).sum())})')


if __name__ == '__main__':
    main()
