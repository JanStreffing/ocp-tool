"""
Shared fixtures for the paleo / model-cycle tests.

The tests run against real OpenIFS ICMGG initial files, because the whole point
of the code under test is the difference between how the cycles lay their
surface fields out — a synthetic GRIB file would encode the assumptions being
tested. Default paths point at the AWI Levante input trees and can be
overridden with OCP_TEST_ICMGG_48R1 / OCP_TEST_ICMGG_43R3. Tests skip when the
files are not reachable rather than fail, so the suite is safe to run anywhere.

The paleo functions take an OCPConfig, but only touch a handful of its
attributes. They are driven here with a small stand-in instead of a full
config, so the tests do not need the complete ocp-tool input tree (FESOM mesh,
Gaussian grid tables, albedo and CO2 GRIB files).
"""

import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from ocp_tool.cycles import get_cycle

ICMGG_48R1_DEFAULT = '/work/ab0246/a270092/input/oifs-48r1/TCO95L91/ICMGGab45INIT_AMIP'
ICMGG_43R3_DEFAULT = '/work/ab0246/a270179/runtime/input/oifs-43r3/TCO95L91/ICMGGaackINIT_MIDPLI1'


def _icmgg(env_var: str, default: str, label: str) -> Path:
    path = Path(os.environ.get(env_var, default))
    if not path.exists():
        pytest.skip(f"{label} ICMGG file not available: {path} (set {env_var})")
    return path


@pytest.fixture(scope='session')
def icmgg_48r1() -> Path:
    """A CY48R1 ICMGG file: multi-layer snow, 5 layers, kg m-2."""
    return _icmgg('OCP_TEST_ICMGG_48R1', ICMGG_48R1_DEFAULT, 'CY48R1')


@pytest.fixture(scope='session')
def icmgg_43r3() -> Path:
    """A CY43R3 ICMGG file: single-layer snow, m of water equivalent."""
    return _icmgg('OCP_TEST_ICMGG_43R3', ICMGG_43R3_DEFAULT, 'CY43R3')


@pytest.fixture
def work_48r1(icmgg_48r1, tmp_path) -> Path:
    """Writable copy of the CY48R1 file (tests modify it in place)."""
    work = tmp_path / 'ICMGGtestINIT_48r1'
    shutil.copy(icmgg_48r1, work)
    return work


@pytest.fixture
def work_43r3(icmgg_43r3, tmp_path) -> Path:
    """Writable copy of the CY43R3 file."""
    work = tmp_path / 'ICMGGtestINIT_43r3'
    shutil.copy(icmgg_43r3, work)
    return work


@pytest.fixture
def ice_mask_nc(tmp_path) -> Path:
    """
    Synthetic PlioMIP3-style ice mask: 1 = ice free, 2 = ice sheet.

    Ice is placed south of ~70S so it overlaps the Antarctic ice of the real
    input files, giving both already-glaciated and newly glaciated points.
    """
    nc = pytest.importorskip('netCDF4')
    path = tmp_path / 'TEST_icemask_v1.0.nc'
    ds = nc.Dataset(path, 'w')
    ds.createDimension('lat', 181)
    ds.createDimension('lon', 360)
    lat = ds.createVariable('lat', 'f4', ('lat',))
    lon = ds.createVariable('lon', 'f4', ('lon',))
    var = ds.createVariable('icemask', 'f4', ('lat', 'lon'))
    lat[:] = np.linspace(-90, 90, 181)
    lon[:] = np.linspace(0, 359, 360)
    values = np.ones((181, 360))
    values[:20, :] = 2
    var[:] = values
    ds.close()
    return path


def make_config(cycle_name, icmgg_file, ice_mask=None, verbose=False):
    """
    Minimal stand-in for OCPConfig covering what the paleo functions read.

    Deliberately not an OCPConfig: building one requires the full input tree,
    and these tests are about field handling, not configuration loading (see
    test_config.py for that).
    """
    paleo = SimpleNamespace(
        get_reconstruction_file=lambda key: ice_mask,
        get_modern_file=lambda key: None,
    )
    return SimpleNamespace(
        options=SimpleNamespace(verbose=verbose),
        cycle=get_cycle(cycle_name),
        paleo=paleo,
        get_icmgg_input_file=lambda: icmgg_file,
        get_icmgg_output_file=lambda: icmgg_file,
    )


def make_grid(icmgg_file):
    """Gaussian grid stand-in built from the coordinates of a GRIB message."""
    from ocp_tool.paleo_input import _get_grib_field_by_code

    _, lats, lons = _get_grib_field_by_code(icmgg_file, 172)
    return SimpleNamespace(lons_list=list(lons), center_lats=np.array(lats))
