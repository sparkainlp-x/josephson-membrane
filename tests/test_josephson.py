"""Tests for the Josephson laws (josephson.py): exact evaluations of
established physics. The laws are inputs, so these tests check that the
implementation reproduces the tabulated values exactly."""

import math

import pytest

from josephson_membrane.josephson import (
    E_CHARGE_C,
    FLUX_QUANTUM_WB,
    H_PLANCK_J_S,
    JOSEPHSON_CONSTANT_GHZ_PER_MV,
    ac_frequency_ghz,
    ac_frequency_hz,
    dc_cardinal_table,
    dc_current,
    squid_cardinal_table,
    squid_critical_current,
)

CODATA_KJ_GHZ_PER_V = 483597.848416984
CODATA_FLUX_QUANTUM_WB = 2.067833848e-15


def test_exact_si_constants():
    assert E_CHARGE_C == 1.602176634e-19
    assert H_PLANCK_J_S == 6.62607015e-34


def test_dc_cardinal_table_exact():
    table = dc_cardinal_table(1.0)
    assert table["0"] == pytest.approx(0.0, abs=1e-15)
    assert table["pi/2"] == pytest.approx(1.0, abs=1e-15)
    assert table["pi"] == pytest.approx(0.0, abs=1e-15)
    assert table["3pi/2"] == pytest.approx(-1.0, abs=1e-15)
    assert table["2pi"] == pytest.approx(0.0, abs=1e-15)


def test_dc_scales_with_critical_current():
    for i_c in (0.5, 2.5, 100.0):
        assert dc_current(math.pi / 2, i_c) == pytest.approx(i_c, abs=1e-12)
        assert dc_current(3 * math.pi / 2, i_c) == pytest.approx(-i_c, abs=1e-12)
        assert dc_current(0.0, i_c) == pytest.approx(0.0, abs=1e-12)


def test_dc_odd_and_periodic():
    for phi in (0.3, 1.0, 2.0, 3.0, 5.5):
        assert dc_current(-phi) == pytest.approx(-dc_current(phi), abs=1e-12)
        assert dc_current(phi + 2 * math.pi) == pytest.approx(
            dc_current(phi), abs=1e-12
        )


def test_dc_bounded_by_critical_current():
    for phi in [i * 0.37 for i in range(20)]:
        assert abs(dc_current(phi, 2.0)) <= 2.0 + 1e-12


def test_ac_constant_matches_codata():
    computed_ghz_per_v = JOSEPHSON_CONSTANT_GHZ_PER_MV * 1000.0
    rel = abs(computed_ghz_per_v - CODATA_KJ_GHZ_PER_V) / CODATA_KJ_GHZ_PER_V
    assert rel < 1e-12
    # The rounded figure quoted in the task (483.5979 GHz/mV) agrees to
    # 6 significant figures; the exact computed value is 483.597848...
    assert JOSEPHSON_CONSTANT_GHZ_PER_MV == pytest.approx(483.5978484, rel=1e-7)


def test_ac_sample_values():
    # K_J = 483597848416984 Hz/V, so 1 uV -> 483597848.417 Hz (483.6 MHz),
    # 1 mV -> 483597848416.984 Hz (483.598 GHz), 1 V -> 483.598 THz.
    assert ac_frequency_hz(1e-6) == pytest.approx(483597848.4169836, rel=1e-12)
    assert ac_frequency_hz(1e-3) == pytest.approx(483597848416.9836, rel=1e-12)
    assert ac_frequency_ghz(1e-3) == pytest.approx(483.5978484169836, rel=1e-12)
    assert ac_frequency_hz(1.0) == pytest.approx(483597848416983.6, rel=1e-12)
    assert ac_frequency_hz(0.0) == pytest.approx(0.0, abs=1e-12)


def test_ac_linear_in_voltage():
    assert ac_frequency_hz(2e-3) == pytest.approx(2 * ac_frequency_hz(1e-3), rel=1e-12)
    assert ac_frequency_hz(-1e-3) == pytest.approx(-ac_frequency_hz(1e-3), rel=1e-12)


def test_flux_quantum_matches_codata():
    rel = abs(FLUX_QUANTUM_WB - CODATA_FLUX_QUANTUM_WB) / CODATA_FLUX_QUANTUM_WB
    assert rel < 1e-9
    assert FLUX_QUANTUM_WB == pytest.approx(2.0678e-15, rel=1e-4)


def test_squid_cardinal_table_exact():
    table = squid_cardinal_table(1.0)
    assert table["0"] == pytest.approx(2.0, abs=1e-12)
    assert table["Phi0/4"] == pytest.approx(math.sqrt(2.0), abs=1e-12)
    assert table["Phi0/2"] == pytest.approx(0.0, abs=1e-12)
    assert table["Phi0"] == pytest.approx(2.0, abs=1e-12)


def test_squid_scales_with_critical_current():
    for i_c in (0.5, 3.0):
        assert squid_critical_current(0.0, i_c) == pytest.approx(2 * i_c, abs=1e-12)
        assert squid_critical_current(FLUX_QUANTUM_WB / 2, i_c) == pytest.approx(
            0.0, abs=1e-12
        )


def test_squid_periodic_and_nonnegative():
    for k in range(9):
        flux = k * FLUX_QUANTUM_WB / 8
        assert squid_critical_current(flux) >= 0.0
        assert squid_critical_current(flux + FLUX_QUANTUM_WB) == pytest.approx(
            squid_critical_current(flux), abs=1e-12
        )
    # Maximum is 2*I_c, minimum is 0.
    assert max(squid_critical_current(f) for f in
               [i * FLUX_QUANTUM_WB / 64 for i in range(65)]) == pytest.approx(
        2.0, abs=1e-9)
