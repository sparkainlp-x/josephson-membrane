"""Tests for the odd-function variant (odd_phase.py) and the ramp.

Three kinds of tests: (1) the quantities -- the two-sided imaginary
part is identically zero, the one-sided overlap Q is odd for identical
cells and not odd for distinct cells; (2) the measured outcome --
claim (e) FAILS on the primary setup (identical cells, linear
generator) by the pre-registered rule, the binary drive is an exact
sine, and the binary/distinct flow PASSES the rule; (3) the ramp --
claim (f) shapes: the v0.1.0 flow is a notch train under the linear
generator, the binary variant drive is a pure sine at the ramp rate.
"""

import math

import numpy as np
import pytest

from josephson_membrane.membrane_phase import phase_unitary
from josephson_membrane.odd_phase import (
    overlap_imag,
    overlap_imag_two_sided,
    run_odd_experiment,
    run_ramp_experiment,
)
from membrane_arena_3x3.core import DIM
from membrane_arena_3x3.lattice import LatticeConfig, MembraneLattice


def _cells(seed, identical):
    rng = np.random.default_rng(seed)
    lat = MembraneLattice(LatticeConfig(), rng)
    rho_a = lat.state[0, 0]
    rho_b = rho_a if identical else lat.state[0, 1]
    return rho_a, rho_b


@pytest.fixture(scope="module")
def odd_exp():
    return run_odd_experiment()


@pytest.fixture(scope="module")
def ramp_exp():
    return run_ramp_experiment()


def test_two_sided_imag_is_identically_zero():
    # Im Tr(rho_A rho_B(phi)) = 0: the trace of a product of two
    # Hermitian matrices is real. Asserted on the real lattice, both
    # generators, identical and distinct cells.
    for generator in ("linear", "binary"):
        for identical in (True, False):
            rho_a, rho_b = _cells(0, identical)
            for phi in (0.3, 1.0, 2.5):
                U = phase_unitary(phi, DIM, generator)
                rho_b_phi = U @ rho_b @ U.conj().T
                assert abs(overlap_imag_two_sided(rho_a, rho_b_phi)) < 1e-12


def test_overlap_imag_zero_at_zero_phase():
    for generator in ("linear", "binary"):
        for identical in (True, False):
            rho_a, rho_b = _cells(0, identical)
            assert abs(overlap_imag(rho_a, rho_b, 0.0, generator)) < 1e-12


def test_overlap_imag_odd_for_identical_cells():
    # Q(-phi) = -Q(phi) for identical cells, both generators.
    for generator in ("linear", "binary"):
        rho_a, rho_b = _cells(0, True)
        for phi in (0.3, 1.0, 2.5):
            assert abs(
                overlap_imag(rho_a, rho_b, phi, generator)
                + overlap_imag(rho_a, rho_b, -phi, generator)
            ) < 1e-9


def test_overlap_imag_changes_sign_identical():
    # Q(pi/2) and Q(3pi/2) have opposite signs (identical, linear):
    # the variant drive changes sign, unlike the v0.1.0 drive.
    rho_a, rho_b = _cells(0, True)
    q_half = overlap_imag(rho_a, rho_b, math.pi / 2, "linear")
    q_three = overlap_imag(rho_a, rho_b, 3 * math.pi / 2, "linear")
    assert q_half > 0 and q_three < 0
    assert abs(q_half + q_three) < 1e-9


def test_overlap_imag_not_odd_for_distinct_cells():
    # REGRESSION: Q is odd for identical cells only. For distinct
    # cells the deviation is ~1e-3, well above machine precision.
    worst = 0.0
    for generator in ("linear", "binary"):
        rho_a, rho_b = _cells(1, False)
        for phi in (0.3, 1.0, 2.5):
            worst = max(
                worst,
                abs(
                    overlap_imag(rho_a, rho_b, phi, generator)
                    + overlap_imag(rho_a, rho_b, -phi, generator)
                ),
            )
    assert worst > 1e-4


def test_binary_identical_drive_is_exact_sine(odd_exp):
    # For the binary generator the odd drive is a single-frequency
    # sine: Q(phi) = W sin(phi) with W = Im Tr(rho_A P_even rho_B).
    # Measured on the full experiment: R^2 = 1.0000, delta = 0.
    fit = odd_exp["setups"]["binary_identical"]["drive_fit"]
    assert fit["r2"] > 0.999
    assert abs(fit["delta"]) < 0.01
    assert fit["max_residual_fraction"] < 0.01


def test_linear_identical_drive_is_multi_frequency(odd_exp):
    # For the linear generator the odd drive is a Fourier series over
    # 31 frequencies, not a single sine: R^2 = 0.0266.
    fit = odd_exp["setups"]["linear_identical"]["drive_fit"]
    assert fit["r2"] < 0.5


def test_claim_e_primary_fails_and_binary_distinct_passes(odd_exp):
    # The pre-registered verdict, full run (5 seeds, 64 points).
    # Primary (identical, linear): FAIL by a wide margin (R^2 = 0.037,
    # max|residual|/|A| = 17.9). Binary/distinct: PASSES the rule
    # (R^2 = 0.9999, max|residual|/|A| = 0.009, delta = -0.008).
    # Binary/identical flow is a near-miss (R^2 = 0.92), distorted
    # by the distance notch; linear/distinct FAILS.
    setups = odd_exp["setups"]
    primary = setups["linear_identical"]
    assert primary["verdict"]["passed"] is False
    assert primary["flow_fit"]["r2"] < 0.5
    assert primary["flow_fit"]["max_residual_fraction"] > 0.5
    assert setups["linear_distinct"]["verdict"]["passed"] is False
    assert setups["binary_identical"]["verdict"]["passed"] is False
    assert 0.85 < setups["binary_identical"]["flow_fit"]["r2"] < 0.95
    passing = setups["binary_distinct"]
    assert passing["verdict"]["passed"] is True
    assert passing["flow_fit"]["r2"] > 0.95
    assert passing["flow_fit"]["max_residual_fraction"] < 0.05
    assert abs(passing["flow_fit"]["delta"]) < 0.05


def test_ramp_v010_flow_is_notch_train_linear_identical(ramp_exp):
    # Claim (f), primary setup: the v0.1.0 flow under a linear phase
    # ramp is a notch train -- nearly constant at ~-0.0016, dipping
    # to 0 at each cycle boundary -- not an oscillation at the ramp
    # rate (13.5% of power at the ramp frequency; sine fit at the
    # ramp rate R^2 = -32.75, the same number as claim (d)).
    s = ramp_exp["setups"]["linear_identical"]
    series = np.asarray(s["mean_flow_v010"])
    assert abs(series.max()) < 1e-6  # the notch: S = 0 at phi = 0
    assert series.mean() < -0.001
    assert s["fft"]["flow_v010"]["power_fraction_at_ramp"] < 0.2
    assert s["sine_at_ramp"]["flow_v010"]["r2"] < 0


def test_ramp_linear_odd_flow_not_at_ramp_rate(ramp_exp):
    # The linear variant flow and drive do not oscillate at the ramp
    # rate: their power sits at high harmonics (peak bins 120 and 96
    # of 128), with under 4% of power at the ramp frequency.
    s = ramp_exp["setups"]["linear_identical"]
    assert s["fft"]["flow_odd"]["power_fraction_at_ramp"] < 0.1
    assert s["fft"]["drive_q"]["power_fraction_at_ramp"] < 0.1
    assert s["fft"]["flow_odd"]["peak_bin"] > 64
    assert s["fft"]["drive_q"]["peak_bin"] > 64


def test_ramp_binary_drive_is_pure_sine_at_ramp_rate(ramp_exp):
    # Claim (f), binary generator: the variant drive Q(t) is a pure
    # sine at the ramp rate (100% of power, R^2 = 1.0000, delta = 0),
    # and the binary/distinct variant flow is a clean sine too
    # (R^2 = 0.9999). The binary/identical flow oscillates at the
    # ramp rate but is distorted by the distance notch (R^2 = 0.92).
    s = ramp_exp["setups"]["binary_identical"]
    assert s["fft"]["drive_q"]["power_fraction_at_ramp"] > 0.99
    assert s["sine_at_ramp"]["drive_q"]["r2"] > 0.999
    assert abs(s["sine_at_ramp"]["drive_q"]["delta"]) < 0.01
    assert 0.85 < s["sine_at_ramp"]["flow_odd"]["r2"] < 0.95
    s2 = ramp_exp["setups"]["binary_distinct"]
    assert s2["fft"]["flow_odd"]["power_fraction_at_ramp"] > 0.99
    assert s2["sine_at_ramp"]["flow_odd"]["r2"] > 0.95
