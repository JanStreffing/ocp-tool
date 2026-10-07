"""
Tests for the truncation helpers: which Gaussian grid and which name prefix go
with a spectral truncation, and the computed full-grid description used for
quadratic truncation.
"""

import numpy as np
import pytest

from ocp_tool.config import gaussian_number, truncation_prefix
from ocp_tool.gaussian_grids import extract_grid_data, read_grid_file


@pytest.mark.parametrize("resolution, truncation_type, expected", [
    (159, "linear", 80),
    (255, "linear", 128),
    (21, "quadratic", 16),
    (42, "quadratic", 32),
    (63, "quadratic", 48),
    (95, "cubic-octahedral", 96),
    (319, "cubic-octahedral", 320),
])
def test_gaussian_number(resolution, truncation_type, expected):
    assert gaussian_number(resolution, truncation_type) == expected


def test_truncation_prefix():
    assert truncation_prefix("linear") == "TL"
    assert truncation_prefix("quadratic") == "TQ"
    assert truncation_prefix("cubic-octahedral") == "TCO"


def test_unknown_truncation_type_raises():
    with pytest.raises(ValueError):
        truncation_prefix("octahedral")
    with pytest.raises(ValueError):
        gaussian_number(21, "octahedral")


def test_quadratic_grid_is_computed_full_gaussian(tmp_path):
    lines, nn = read_grid_file(
        resolution=21,
        reduced_grid_path=tmp_path,
        full_grid_path=tmp_path,
        truncation_type="quadratic",
    )
    assert nn == 16
    assert (tmp_path / "f16_full.txt").exists()

    lons, lats, numlons, dlons, lat_list = extract_grid_data(lines)
    assert len(lat_list) == 32
    assert set(numlons) == {64}
    assert len(lons) == 2048

    # North to south, symmetric about the equator, first latitude of N16
    assert lat_list == sorted(lat_list, reverse=True)
    assert np.allclose(lat_list, -np.array(lat_list[::-1]))
    assert lat_list[0] == pytest.approx(85.7605871204438, abs=1e-6)
