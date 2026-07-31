"""Tests for the OpenIFS model-cycle switch."""

import pytest

from ocp_tool.cycles import (
    AUTO,
    CY43R3,
    CY48R1,
    ICE_SHEET_SNOW_MWE,
    detect_cycle,
    get_cycle,
    has_snow_layers,
    normalise_cycle_name,
    resolve_cycle,
)


class TestCycleSpec:
    def test_snow_field_layout_differs(self):
        assert not CY43R3.snow_multi_layer
        assert CY43R3.snow_mass_code == 141
        assert CY43R3.snow_layer_codes == ()

        assert CY48R1.snow_multi_layer
        assert CY48R1.snow_mass_code == 228141
        # sd, asn, rsn, tsn, lwcs
        assert set(CY48R1.snow_layer_codes) == {228141, 32, 33, 238, 228038}

    def test_snow_mass_units_differ_by_factor_1000(self):
        """
        CY43R3 stores m of water equivalent, CY48R1 kg m-2. Mixing the two is
        silent, which is why the scale lives in the cycle definition.
        """
        assert CY43R3.ice_sheet_snow_mass == pytest.approx(ICE_SHEET_SNOW_MWE)
        assert CY48R1.ice_sheet_snow_mass == pytest.approx(1000 * ICE_SHEET_SNOW_MWE)
        assert CY43R3.scale_snow_mass(2.5) == pytest.approx(2.5)
        assert CY48R1.scale_snow_mass(2.5) == pytest.approx(2500.0)


class TestNameHandling:
    @pytest.mark.parametrize('name,expected', [
        ('43r3', '43r3'), ('CY43R3', '43r3'), ('43', '43r3'),
        ('48r1', '48r1'), ('cy48r1', '48r1'), ('48', '48r1'),
        (' 48R1 ', '48r1'),
    ])
    def test_aliases_accepted(self, name, expected):
        assert normalise_cycle_name(name) == expected
        assert get_cycle(name).name == expected

    def test_auto_passes_through(self):
        assert normalise_cycle_name('auto') == AUTO
        assert normalise_cycle_name('') == AUTO

    def test_unknown_cycle_rejected(self):
        with pytest.raises(ValueError, match='Unknown OpenIFS cycle'):
            normalise_cycle_name('47r1')

    def test_get_cycle_refuses_auto(self):
        with pytest.raises(ValueError, match='explicit cycle'):
            get_cycle('auto')


class TestDetection:
    def test_detects_48r1_from_snow_layers(self, icmgg_48r1):
        assert has_snow_layers(icmgg_48r1)
        assert detect_cycle(icmgg_48r1).name == '48r1'

    def test_detects_43r3_without_snow_layers(self, icmgg_43r3):
        assert not has_snow_layers(icmgg_43r3)
        assert detect_cycle(icmgg_43r3).name == '43r3'


class TestResolve:
    def test_auto_matches_detection(self, icmgg_48r1, icmgg_43r3):
        assert resolve_cycle(AUTO, icmgg_48r1, verbose=False).name == '48r1'
        assert resolve_cycle(AUTO, icmgg_43r3, verbose=False).name == '43r3'

    def test_auto_without_file_is_an_error(self, tmp_path):
        with pytest.raises(ValueError, match='does not exist'):
            resolve_cycle(AUTO, tmp_path / 'absent', verbose=False)

    def test_explicit_agreeing_with_file(self, icmgg_48r1):
        assert resolve_cycle('48r1', icmgg_48r1, verbose=False).name == '48r1'

    @pytest.mark.parametrize('requested,fixture', [
        ('48r1', 'icmgg_43r3'),
        ('43r3', 'icmgg_48r1'),
    ])
    def test_explicit_contradicting_file_raises(self, requested, fixture, request):
        """
        The mistake this guard exists for otherwise surfaces much later, as
        IOSTREAM_MIX:GRID_IN - MISSING FIELD at model start-up.
        """
        icmgg = request.getfixturevalue(fixture)
        with pytest.raises(ValueError, match='but .* looks like'):
            resolve_cycle(requested, icmgg, verbose=False)

    def test_explicit_without_file_is_accepted(self, tmp_path):
        """No file to cross-check against; the config is taken at its word."""
        assert resolve_cycle('48r1', tmp_path / 'absent', verbose=False).name == '48r1'
