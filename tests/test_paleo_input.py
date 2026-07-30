"""
Tests for the cycle-dependent parts of the paleo ICMGG modifications.

Covers the three failure modes that the model-cycle work fixes, all of which
were silent: fields identified by paramId alone collapsing the snow layers
onto one entry, the snow mass being written in the wrong units or to no
message at all, and the lake variables never matching.
"""

import numpy as np
import pytest

from conftest import make_config, make_grid
from ocp_tool.cycles import get_cycle
from ocp_tool.paleo_input import (
    REMAINING_FIELD_CODES,
    _get_grib_field_by_code,
    _read_all_grib_fields,
    _read_total_snow_mass,
    _replace_grib_fields,
    create_paleo_masks,
    interpolate_remaining_fields,
    modify_ice_snow,
)

SNOW_CODES_48R1 = (228141, 32, 33, 238, 228038)


def count_messages(grib_file):
    """Number of GRIB messages in a file, duplicates included."""
    import eccodes
    n = 0
    with open(grib_file, 'rb') as f:
        while True:
            gid = eccodes.codes_grib_new_from_file(f)
            if gid is None:
                return n
            eccodes.codes_release(gid)
            n += 1


class TestFieldKeying:
    def test_48r1_snow_block_keyed_per_layer(self, icmgg_48r1):
        """
        Keying on the code alone collapsed five layers onto one entry, so four
        fifths of the snowpack was invisible.
        """
        fields = _read_all_grib_fields(icmgg_48r1)
        for code in SNOW_CODES_48R1:
            levels = sorted(lev for c, lev in fields if c == code)
            assert len(levels) == 5, f"code {code} has levels {levels}"

    def test_43r3_snow_is_single_level(self, icmgg_43r3):
        fields = _read_all_grib_fields(icmgg_43r3)
        for code in (141, 32, 33, 238):
            levels = sorted(lev for c, lev in fields if c == code)
            assert len(levels) == 1, f"code {code} has levels {levels}"

    def test_snow_albedo_resolves_despite_paramid_zero(self, icmgg_48r1):
        """
        IFS writes multi-layer snow albedo as GRIB2 (0, 19, 192), which some
        ecCodes builds report as paramId 0. It must still come out as code 32.
        """
        fields = _read_all_grib_fields(icmgg_48r1)
        assert sorted(lev for c, lev in fields if c == 32)
        assert not [key for key in fields if key[0] == 0]

    def test_get_field_by_code_can_select_a_level(self, icmgg_48r1):
        fields = _read_all_grib_fields(icmgg_48r1)
        levels = sorted(lev for c, lev in fields if c == 238)
        first, last = levels[0], levels[-1]
        values_first, _, _ = _get_grib_field_by_code(icmgg_48r1, 238, level=first)
        values_last, _, _ = _get_grib_field_by_code(icmgg_48r1, 238, level=last)
        np.testing.assert_allclose(values_first, fields[(238, first)])
        np.testing.assert_allclose(values_last, fields[(238, last)])
        assert not np.allclose(values_first, values_last)


class TestLakeCodes:
    def test_lake_variables_are_table_228_paramids(self, icmgg_48r1, icmgg_43r3):
        """
        mod_input_oifs.sh listed the lake fields as 8-14, correct for
        `cdo -selcode` which matches indicatorOfParameter within table 228.
        This port matches on paramId, where they are 228008-228014; as 8-14
        they matched nothing and were skipped without an error.
        """
        expected = set(range(228008, 228015))
        assert expected <= set(REMAINING_FIELD_CODES)
        assert not (set(range(8, 15)) & set(REMAINING_FIELD_CODES))

        for icmgg in (icmgg_48r1, icmgg_43r3):
            codes = {code for code, _ in _read_all_grib_fields(icmgg)}
            assert expected <= codes
            assert not (set(range(8, 15)) & codes)


class TestTotalSnowMass:
    def test_48r1_sums_over_layers(self, icmgg_48r1):
        """No single CY48R1 layer holds the column total."""
        cycle = get_cycle('48r1')
        fields = _read_all_grib_fields(icmgg_48r1)
        layers = {lev: fields[(cycle.snow_mass_code, lev)]
                  for c, lev in fields if c == cycle.snow_mass_code}
        total = _read_total_snow_mass(icmgg_48r1, cycle)
        np.testing.assert_allclose(total, sum(layers.values()))
        # kg m-2, so a 10 m w.e. ice sheet reads as 1e4
        assert total.max() == pytest.approx(cycle.ice_sheet_snow_mass)

    def test_43r3_uses_the_single_field(self, icmgg_43r3):
        cycle = get_cycle('43r3')
        total = _read_total_snow_mass(icmgg_43r3, cycle)
        values, _, _ = _get_grib_field_by_code(icmgg_43r3, 141)
        np.testing.assert_allclose(total, values)
        # m of water equivalent
        assert total.max() == pytest.approx(cycle.ice_sheet_snow_mass)

    def test_absent_snow_code_returns_none(self, icmgg_43r3):
        assert _read_total_snow_mass(icmgg_43r3, get_cycle('48r1')) is None


class TestReplaceGribFields:
    def test_bare_code_on_multi_level_field_is_rejected(self, work_48r1):
        """A bare 238 would write one profile into all five snow layers."""
        n = len(_read_all_grib_fields(work_48r1)[(172, 0)])
        with pytest.raises(ValueError, match='exists on levels'):
            _replace_grib_fields(work_48r1, work_48r1, {238: np.zeros(n)})

    def test_bare_code_on_single_level_field_is_accepted(self, work_43r3):
        """The CY43R3 contract: modify_ice_snow returns a bare 141 key."""
        n = len(_read_all_grib_fields(work_43r3)[(172, 0)])
        _replace_grib_fields(work_43r3, work_43r3, {141: np.full(n, 7.0)})
        assert _read_all_grib_fields(work_43r3)[(141, 0)].max() == pytest.approx(7.0)

    def test_per_level_write_leaves_other_levels_alone(self, work_48r1):
        before = _read_all_grib_fields(work_48r1)
        n = len(before[(172, 0)])
        _replace_grib_fields(work_48r1, work_48r1, {(238, 3): np.full(n, 260.0)})
        after = _read_all_grib_fields(work_48r1)
        assert after[(238, 3)].max() == pytest.approx(260.0)
        for level in (1, 2, 4, 5):
            np.testing.assert_allclose(after[(238, level)], before[(238, level)])

    def test_mixed_bare_and_per_level_keys(self, work_48r1):
        """
        modify_paleo_input merges bare keys from the surface steps with
        per-level keys from the snow steps into one replacement dict.
        """
        n = len(_read_all_grib_fields(work_48r1)[(172, 0)])
        _replace_grib_fields(work_48r1, work_48r1, {
            43: np.full(n, 2.0),
            (228141, 1): np.full(n, 5.0),
        })
        after = _read_all_grib_fields(work_48r1)
        assert after[(43, 0)].max() == pytest.approx(2.0)
        assert after[(228141, 1)].max() == pytest.approx(5.0)
        assert after[(228141, 2)].max() == pytest.approx(0.0)

    def test_message_count_is_preserved(self, work_48r1, icmgg_48r1):
        """
        Duplicated fields (lsm, cl, 200026 appear twice) must survive: the
        reader de-duplicates, the writer must not.
        """
        n = len(_read_all_grib_fields(work_48r1)[(172, 0)])
        expected = count_messages(icmgg_48r1)
        _replace_grib_fields(work_48r1, work_48r1, {(228141, 1): np.zeros(n)})
        assert count_messages(work_48r1) == expected


class TestModifyIceSnow:
    def test_43r3_writes_single_field_in_metres(self, work_43r3, ice_mask_nc):
        config = make_config('43r3', work_43r3, ice_mask=ice_mask_nc)
        masks = {'modern_has_ice': _read_total_snow_mass(work_43r3, config.cycle) > 0}
        result = modify_ice_snow(config, make_grid(work_43r3), masks, work_43r3)

        assert set(result) == {141}, "CY43R3 must keep the bare-code contract"
        assert set(np.unique(result[141])) == {0.0, 10.0}

    def test_48r1_puts_mass_in_top_layer_in_kg_per_m2(self, work_48r1, ice_mask_nc):
        config = make_config('48r1', work_48r1, ice_mask=ice_mask_nc)
        cycle = config.cycle
        masks = {'modern_has_ice': _read_total_snow_mass(work_48r1, cycle) > 0}
        result = modify_ice_snow(config, make_grid(work_48r1), masks, work_48r1)

        assert (cycle.snow_mass_code, 1) in result
        mass = result[(cycle.snow_mass_code, 1)]
        assert set(np.unique(mass)) == {0.0, cycle.ice_sheet_snow_mass}
        assert mass.max() == pytest.approx(10000.0)

    def test_48r1_total_column_matches_the_ice_mask(self, work_48r1, ice_mask_nc):
        config = make_config('48r1', work_48r1, ice_mask=ice_mask_nc)
        cycle = config.cycle
        masks = {'modern_has_ice': _read_total_snow_mass(work_48r1, cycle) > 0}
        result = modify_ice_snow(config, make_grid(work_48r1), masks, work_48r1)
        _replace_grib_fields(work_48r1, work_48r1, result)

        total = _read_total_snow_mass(work_48r1, cycle)
        expected_points = int((result[(cycle.snow_mass_code, 1)] > 0).sum())
        assert int((total > 0).sum()) == expected_points
        assert total.max() == pytest.approx(cycle.ice_sheet_snow_mass)

    def test_48r1_does_not_prescribe_snow_properties(self, work_48r1, ice_mask_nc):
        """
        The original workflow left albedo/density/temperature to the
        distance-weighted fill; only the mass field is set here.
        """
        config = make_config('48r1', work_48r1, ice_mask=ice_mask_nc)
        masks = {'modern_has_ice': _read_total_snow_mass(work_48r1, config.cycle) > 0}
        result = modify_ice_snow(config, make_grid(work_48r1), masks, work_48r1)
        assert {code for code, _ in result} == {config.cycle.snow_mass_code}

    def test_missing_ice_mask_is_skipped(self, work_48r1, tmp_path):
        config = make_config('48r1', work_48r1, ice_mask=tmp_path / 'absent.nc')
        masks = {'modern_has_ice': _read_total_snow_mass(work_48r1, config.cycle) > 0}
        assert modify_ice_snow(config, make_grid(work_48r1), masks, work_48r1) == {}


class TestCreatePaleoMasks:
    def test_modern_ice_comes_from_total_column_snow(self, work_48r1, ice_mask_nc):
        """
        Reading a single layer, or code 141 which does not exist on CY48R1,
        left this mask empty and the reduced-ice handling dead.
        """
        from types import SimpleNamespace

        config = make_config('48r1', work_48r1, ice_mask=ice_mask_nc)
        lsm, _, _ = _get_grib_field_by_code(work_48r1, 172)
        lsm_data = SimpleNamespace(lsm_binary_atm=lsm.copy())

        masks = create_paleo_masks(config, make_grid(work_48r1), lsm_data)

        expected = _read_total_snow_mass(work_48r1, config.cycle) > 0
        np.testing.assert_array_equal(masks['modern_has_ice'], expected)
        assert masks['modern_has_ice'].sum() > 0
        assert masks['reduced_ice'].sum() > 0


class TestInterpolateRemainingFields:
    @pytest.fixture
    def masks(self, work_48r1):
        """Synthetic added / drowned / reduced-ice points on the real LSM."""
        lsm, _, _ = _get_grib_field_by_code(work_48r1, 172)
        n = len(lsm)
        land = lsm >= 0.5
        rng = np.random.default_rng(0)
        added = np.zeros(n, bool)
        drowned = np.zeros(n, bool)
        reduced_ice = np.zeros(n, bool)
        added[rng.choice(np.where(~land)[0], 200, replace=False)] = True
        drowned[rng.choice(np.where(land)[0], 200, replace=False)] = True
        reduced_ice[rng.choice(np.where(land)[0], 50, replace=False)] = True
        paleo_land = (land | added) & ~drowned
        return dict(land=paleo_land, ocean=~paleo_land, added=added,
                    drowned=drowned, reduced_ice=reduced_ice,
                    modern_has_ice=np.zeros(n, bool))

    def test_each_snow_layer_is_filled_separately(self, work_48r1, ice_mask_nc, masks):
        config = make_config('48r1', work_48r1, ice_mask=ice_mask_nc)
        result = interpolate_remaining_fields(config, make_grid(work_48r1), masks)

        for code in (32, 33, 238):
            levels = sorted(lev for c, lev in result if c == code)
            assert levels == [1, 2, 3, 4, 5], f"code {code} filled on {levels}"

    def test_layer_profiles_stay_distinct(self, work_48r1, ice_mask_nc, masks):
        """A per-code fill would have flattened all layers to one profile."""
        config = make_config('48r1', work_48r1, ice_mask=ice_mask_nc)
        result = interpolate_remaining_fields(config, make_grid(work_48r1), masks)
        means = [result[(238, lev)].mean() for lev in (1, 2, 3, 4, 5)]
        assert len(set(np.round(means, 6))) > 1

    def test_only_fill_points_change(self, work_48r1, ice_mask_nc, masks):
        config = make_config('48r1', work_48r1, ice_mask=ice_mask_nc)
        before = _read_all_grib_fields(work_48r1)
        result = interpolate_remaining_fields(config, make_grid(work_48r1), masks)

        allowed = masks['added'] | masks['drowned'] | masks['reduced_ice']
        for key, values in result.items():
            changed = np.where(values != before[key])[0]
            assert allowed[changed].all(), f"{key} changed outside the fill masks"

    def test_43r3_still_fills_its_single_snow_level(self, work_43r3, ice_mask_nc):
        lsm, _, _ = _get_grib_field_by_code(work_43r3, 172)
        n = len(lsm)
        land = lsm >= 0.5
        rng = np.random.default_rng(1)
        added = np.zeros(n, bool)
        added[rng.choice(np.where(~land)[0], 100, replace=False)] = True
        paleo_land = land | added
        masks = dict(land=paleo_land, ocean=~paleo_land, added=added,
                     drowned=np.zeros(n, bool), reduced_ice=np.zeros(n, bool),
                     modern_has_ice=np.zeros(n, bool))

        config = make_config('43r3', work_43r3, ice_mask=ice_mask_nc)
        result = interpolate_remaining_fields(config, make_grid(work_43r3), masks)
        for code in (32, 33, 238):
            assert [lev for c, lev in result if c == code] == [0]
