"""
Tests for loading the model-cycle setting through the YAML config.

Only the atmosphere/cycle part of the config is exercised: a full config load
needs nothing more than the paths section to be syntactically present, but
resolving the cycle needs a real ICMGG file under openifs_default.
"""

import textwrap

import pytest

from ocp_tool.config import load_config

CONFIG_TEMPLATE = textwrap.dedent("""
    atmosphere:
      resolution_list: [95]
      truncation_type: "cubic-octahedral"
      experiment_name: "{expid}"
    {cycle_line}
    ocean:
      grid_name: "AMIP"
      has_ice_cavities: false
      mesh_file: null
    runoff:
      manual_basin_removal: []
    paths:
      input:
        fesom_mesh: "input/fesom_mesh/"
        gaussian_grids_full: "input/gaussian_grids_full/"
        gaussian_grids_octahedral_reduced: "input/gaussian_grids_octahedral_reduced/"
        gaussian_grids_linear_reduced: "input/gaussian_grids_linear_reduced/"
        openifs_default: "input/openifs_input_default/"
        runoff_default: "input/runoff_map_default/"
        lpj_guess: "input/lpj-guess/"
    options:
      verbose: false
      parallel_workers: 1
      use_dask: false
""")

EXPID = 'test'


@pytest.fixture
def write_config(tmp_path):
    """Build a config tree whose openifs_default holds the given ICMGG file."""
    def _write(icmgg_file, model_cycle=None):
        root = tmp_path / 'root'
        (root / 'configs').mkdir(parents=True, exist_ok=True)
        openifs = root / 'input' / 'openifs_input_default'
        openifs.mkdir(parents=True, exist_ok=True)

        target = openifs / f'ICMGG{EXPID}INIT'
        if target.exists() or target.is_symlink():
            target.unlink()
        target.symlink_to(icmgg_file)

        cycle_line = f'  model_cycle: "{model_cycle}"' if model_cycle else ''
        path = root / 'configs' / 'test.yaml'
        path.write_text(CONFIG_TEMPLATE.format(expid=EXPID, cycle_line=cycle_line))
        return path
    return _write


def test_default_is_auto_detection(write_config, icmgg_48r1):
    """model_cycle is optional; omitting it detects from the input file."""
    config = load_config(write_config(icmgg_48r1))
    assert config.atmosphere.model_cycle == 'auto'
    assert config.cycle.name == '48r1'


def test_auto_detects_43r3(write_config, icmgg_43r3):
    assert load_config(write_config(icmgg_43r3, 'auto')).cycle.name == '43r3'


def test_explicit_cycle_is_used(write_config, icmgg_43r3):
    assert load_config(write_config(icmgg_43r3, '43r3')).cycle.name == '43r3'


def test_alias_is_accepted(write_config, icmgg_48r1):
    assert load_config(write_config(icmgg_48r1, 'CY48R1')).cycle.name == '48r1'


def test_cycle_contradicting_the_input_file_raises(write_config, icmgg_43r3):
    config = load_config(write_config(icmgg_43r3, '48r1'))
    with pytest.raises(ValueError, match='but .* looks like'):
        config.cycle


def test_cycle_is_resolved_once(write_config, icmgg_48r1):
    config = load_config(write_config(icmgg_48r1))
    assert config.cycle is config.cycle


def test_icmgg_input_path_follows_experiment_name(write_config, icmgg_48r1):
    """
    The experiment id is stored in the GRIB header, so the input path is
    derived from experiment_name rather than overridable per file.
    """
    config = load_config(write_config(icmgg_48r1))
    assert config.get_icmgg_input_file().name == f'ICMGG{EXPID}INIT'
