"""restore_lake_cover: put back what the nearest-neighbour fill overwrote.

Driven by the real TCO95 ICMGG rather than a synthetic array, because the
property under test is a statement about the ECMWF field convention -- that
``cl`` is the lake share of a box and is capped by its water share -- and a
hand-built fixture would simply assert whatever the author already believed.
"""
from pathlib import Path

import numpy as np
import pytest

from ocp_tool.lsm import (
    FLAKE_DEPTH_MAX_M,
    FLAKE_DEPTH_MIN_M,
    _clip_flake_depth,
    read_field_ids,
    read_grib_fields,
    restore_lake_cover,
)

ICMGG = Path('/work/ab0246/a270092/software/ocp-tool/input/'
             'openifs_input_default/ICMGGab45INIT')

pytestmark = pytest.mark.skipif(
    not ICMGG.exists(), reason=f'{ICMGG} not available')


@pytest.fixture(scope='module')
def fields():
    gribfield, lsm_id, slt_id, cl_id, gids, _ = read_grib_fields(
        ICMGG, num_fields=None, verbose=False)
    import gribapi
    for g in gids:
        if g is not None:
            try:
                gribapi.grib_release(g)
            except Exception:                                   # noqa: BLE001
                pass
    return gribfield, read_field_ids(ICMGG), lsm_id


def test_depth_clip_matches_flake():
    """The band is FLake's, not ours: flakeene_mod.F90:211-212,238."""
    assert _clip_flake_depth(0.0) == FLAKE_DEPTH_MIN_M
    assert _clip_flake_depth(1.0) == FLAKE_DEPTH_MIN_M
    assert _clip_flake_depth(25.0) == 25.0
    assert _clip_flake_depth(4064.0) == FLAKE_DEPTH_MAX_M
    # A Caspian at 1000 m and one at 50 m are the same model.
    assert _clip_flake_depth(1000.0) == _clip_flake_depth(50.0)


def test_lake_cover_never_exceeds_water_share(fields):
    """
    The convention the restore relies on, asserted against the shipped file.

    If this ever fails, ``cl`` is not what this module assumes and the restore
    is unsafe -- it would be putting back a lake fraction larger than the box's
    water fraction.
    """
    gribfield, ids, lsm_id = fields
    cl = np.asarray(gribfield[ids['cl']], float)
    water = 1.0 - np.asarray(gribfield[lsm_id], float)
    assert np.all(cl <= water + 1e-4), float((cl - water).max())


def test_restore_recovers_the_pristine_values(fields):
    """A flipped cell gets its own cl back, not 1.0 and not a neighbour's."""
    gribfield, ids, _ = fields
    cl_id, dl_id = ids['cl'], ids['dl']
    pristine = {f: np.array(gribfield[f], copy=True)
                for f in (cl_id, dl_id)}

    # Pick real lake cells and real non-lake cells from the shipped file.
    cl0 = pristine[cl_id]
    lake_cells = np.where(cl0 >= 0.5)[0][:10].tolist()
    dry_cells = np.where(cl0 == 0.0)[0][:10].tolist()
    assert lake_cells and dry_cells

    # Simulate what the NN fill does: overwrite cl and dl from a dry neighbour.
    mod = [np.array(a, copy=True) for a in gribfield]
    for i in lake_cells + dry_cells:
        mod[cl_id][i] = 0.0
        mod[dl_id][i] = 7.0

    n = restore_lake_cover(
        gribfield_pristine=pristine,
        gribfield_mod=mod,
        lake_ids={'cl': cl_id, 'dl': dl_id},
        changed_to_land=lake_cells + dry_cells,
        verbose=False,
    )
    assert n == len(lake_cells) + len(dry_cells)

    for i in lake_cells:
        # the cell's OWN value, not a forced 1.0
        assert mod[cl_id][i] == pytest.approx(cl0[i])
        assert mod[dl_id][i] == pytest.approx(
            _clip_flake_depth(pristine[dl_id][i]))
    for i in dry_cells:
        # a non-lake cell is restored to its pristine zero, not invented into one
        assert mod[cl_id][i] == 0.0


def test_restore_writes_only_the_flipped_cells(fields):
    """Untouched cells must come out bit-identical."""
    gribfield, ids, _ = fields
    cl_id, dl_id = ids['cl'], ids['dl']
    pristine = {f: np.array(gribfield[f], copy=True) for f in (cl_id, dl_id)}
    mod = [np.array(a, copy=True) for a in gribfield]
    before = np.array(mod[cl_id], copy=True)

    targets = [11, 2222, 33333]
    for i in targets:
        mod[cl_id][i] = 0.123
    restore_lake_cover(pristine, mod, {'cl': cl_id, 'dl': dl_id}, targets)

    changed = np.where(before != mod[cl_id])[0]
    assert changed.size == 0, f'cl changed at unexpected cells: {changed[:20]}'


def test_no_flipped_cells_is_a_noop(fields):
    gribfield, ids, _ = fields
    cl_id, dl_id = ids['cl'], ids['dl']
    pristine = {f: np.array(gribfield[f], copy=True) for f in (cl_id, dl_id)}
    mod = [np.array(a, copy=True) for a in gribfield]
    assert restore_lake_cover(pristine, mod, {'cl': cl_id, 'dl': dl_id}, []) == 0
    assert np.array_equal(mod[cl_id], gribfield[cl_id])
