# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP)

"""Tests for the membrane phase response (membrane_phase.py).

Two kinds of tests: (1) mechanics -- the phase unitary, the pairwise
flow, the sine fit -- verified on synthetic inputs and on the real
lattice; (2) the measured outcome -- the sine-law verdict FAILs on the
real lattice, asserted with wide margins (measured R^2 ~ -30 and
max|residual|/|A| ~ 12 on the primary setup, versus the bar R^2 >= 0.95
and max|residual|/|A| <= 0.05).
"""

import math

import numpy as np
import pytest

from josephson_membrane.membrane_phase import (
    SEEDS,
    SINE_MAX_RESIDUAL_FRACTION,
    SINE_R2_MIN,
    lattice_provenance,
    pairwise_flow,
    phase_unitary,
    run_membrane_experiment,
    sine_fit,
    sine_law_verdict,
    sweep,
)
from membrane_arena_3x3.core import DIM, fidelity, is_valid_density_matrix
from membrane_arena_3x3.lattice import LatticeConfig, MembraneLattice


def test_phase_unitary_is_unitary():
    for gen in ("linear", "binary"):
        for phi in (0.0, 0.7, math.pi, 2.0):
            U = phase_unitary(phi, DIM, gen)
            assert np.allclose(U.conj().T @ U, np.eye(DIM), atol=1e-12)
            assert np.allclose(U @ U.conj().T, np.eye(DIM), atol=1e-12)


def test_phase_unitary_zero_is_identity():
    for gen in ("linear", "binary"):
        assert np.allclose(phase_unitary(0.0, DIM, gen), np.eye(DIM), atol=1e-12)


def test_phase_unitary_preserves_populations_and_validity():
    rng = np.random.default_rng(0)
    lat = MembraneLattice(LatticeConfig(), rng)
    rho = lat.state[0, 0]
    for gen in ("linear", "binary"):
        for phi in (0.5, 1.3, 3.0):
            U = phase_unitary(phi, DIM, gen)
            rotated = U @ rho @ U.conj().T
            # Populations (the diagonal) are exactly unchanged.
            assert np.allclose(np.diag(rotated), np.diag(rho), atol=1e-12)
            assert is_valid_density_matrix(rotated)
            # Trace and purity are preserved by unitary conjugation.
            assert np.trace(rotated).real == pytest.approx(1.0, abs=1e-12)
            assert np.trace(rotated @ rotated).real == pytest.approx(
                np.trace(rho @ rho).real, abs=1e-12
            )


def test_phase_unitary_binary_touches_only_even_levels():
    U = phase_unitary(0.9, DIM, "binary")
    diag = np.diag(U)
    assert np.allclose(diag[1::2], 1.0, atol=1e-12)
    assert np.allclose(diag[0::2], np.exp(1j * 0.9), atol=1e-12)


def test_phase_unitary_unknown_generator():
    with pytest.raises(ValueError):
        phase_unitary(1.0, DIM, "cubic")


def test_pairwise_flow_identical_at_zero_phase():
    rng = np.random.default_rng(0)
    lat = MembraneLattice(LatticeConfig(), rng)
    rho = lat.state[0, 0]
    flow = pairwise_flow(rho, rho, lat.config.coupling, lat.config.resonance_threshold)
    assert flow["fidelity"] == pytest.approx(1.0, abs=1e-12)
    assert flow["distance"] == pytest.approx(0.0, abs=1e-12)
    assert flow["signed_flow"] == pytest.approx(0.0, abs=1e-12)
    assert flow["magnitude"] == pytest.approx(0.0, abs=1e-12)
    # Drive is positive at fidelity 1: coupling * (1 - 0.90).
    assert flow["drive"] == pytest.approx(0.02 * 0.10, abs=1e-12)


def test_pairwise_flow_negative_below_resonance():
    rng = np.random.default_rng(0)
    lat = MembraneLattice(LatticeConfig(), rng)
    rho_a, rho_b = lat.state[0, 0], lat.state[0, 1]
    flow = pairwise_flow(rho_a, rho_b, lat.config.coupling,
                         lat.config.resonance_threshold)
    # Two independent random states sit far below the 0.90 threshold.
    assert flow["fidelity"] < 0.90
    assert flow["drive"] < 0.0
    assert flow["signed_flow"] < 0.0


def test_sine_fit_recovers_known_sine():
    phis = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    y = 2.5 * np.sin(phis + 0.7)
    fit = sine_fit(phis, y)
    assert fit["amplitude"] == pytest.approx(2.5, abs=1e-9)
    assert fit["delta"] == pytest.approx(0.7, abs=1e-9)
    assert fit["r2"] == pytest.approx(1.0, abs=1e-9)
    assert fit["max_residual"] < 1e-9
    assert sine_law_verdict(fit)["passed"] is True


def test_sine_fit_on_constant_fails():
    phis = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    y = np.full_like(phis, -0.00164)
    fit = sine_fit(phis, y)
    assert fit["amplitude"] < 1e-12
    assert fit["r2"] < 0.5
    assert sine_law_verdict(fit)["passed"] is False


def test_sine_fit_with_dc_recovers_offset_sine():
    phis = np.linspace(0, 2 * np.pi, 64, endpoint=False)
    y = -0.001 + 0.0008 * np.sin(phis + 1.571)
    fit = sine_fit(phis, y, with_dc=True)
    assert fit["amplitude"] == pytest.approx(0.0008, abs=1e-9)
    assert fit["delta"] == pytest.approx(1.571, abs=1e-3)
    assert fit["dc"] == pytest.approx(-0.001, abs=1e-9)
    assert fit["r2"] == pytest.approx(1.0, abs=1e-9)


def test_sweep_shapes_and_ranges():
    sw = sweep(0, "linear", n_phi=32, identical=True, step_change_grid=8)
    assert sw["phis"].shape == (32,)
    assert sw["signed_flow"].shape == (32,)
    assert sw["magnitude"].shape == (32,)
    assert sw["lattice_delta"].shape == (32,)
    assert np.all((sw["fidelity"] >= 0.0) & (sw["fidelity"] <= 1.0))
    assert np.all(sw["magnitude"] >= 0.0)
    # Identical cells: fidelity 1 and zero flow at delta_phi = 0.
    assert sw["fidelity"][0] == pytest.approx(1.0, abs=1e-12)
    assert sw["signed_flow"][0] == pytest.approx(0.0, abs=1e-12)
    assert sw["step_change"].shape == (8,)
    assert np.all(sw["step_change"] > 0.0)


def test_signed_flow_even_in_phase():
    # Structural property: the flow depends on the fidelity (an even
    # function of delta_phi), so S(phi) == S(2*pi - phi). The Josephson
    # current is odd -- this evenness alone forbids a sine law.
    sw = sweep(1, "linear", n_phi=64, identical=True)
    s = sw["signed_flow"]
    assert np.allclose(s[1:], s[:0:-1], atol=1e-12)
    sw2 = sweep(1, "binary", n_phi=64, identical=True)
    s2 = sw2["signed_flow"]
    assert np.allclose(s2[1:], s2[:0:-1], atol=1e-12)


def test_lattice_provenance_points_at_real_tree():
    prov = lattice_provenance()
    assert "membrane-arena-3x3-upgrade" in prov["path"]
    assert prov["branch"] == "rounds-on-main"
    assert prov["head"] is not None and len(prov["head"]) == 40


def test_measured_sine_law_fails_with_margin():
    # Small run (2 seeds, 32 points): the primary verdict must FAIL by a
    # wide margin, as measured in the exploratory probes (R^2 ~ -30,
    # max|residual|/|A| ~ 12, versus the bar 0.95 / 0.05).
    exp = run_membrane_experiment(seeds=(0, 1), n_phi=32, step_change_grid=0)
    fit = exp["linear_identical"]["fit"]
    assert fit["r2"] < 0.5
    assert fit["max_residual_fraction"] > 0.5
    assert exp["linear_identical"]["verdict"]["passed"] is False
    # The best-fit phase is ~pi/2 (a cosine), not the Josephson delta=0.
    assert abs(fit["delta"] - math.pi / 2) < 0.1
    # The binary generator fails too.
    assert exp["binary_identical"]["verdict"]["passed"] is False


def test_step_change_nearly_phase_independent():
    # The one-step state change is dominated by the unitary /
    # dissipative / mirror stages, not by the interaction: it barely
    # moves with delta_phi.
    sw = sweep(0, "linear", n_phi=16, identical=True, step_change_grid=16)
    vals = sw["step_change"]
    assert (vals.max() - vals.min()) / vals.mean() < 0.01


def test_distinct_flow_not_even():
    # REGRESSION (review v0.2.0): the evenness argument is scoped to
    # identical cells. For distinct cells the fidelity is not even in
    # general -- S(pi/2) != S(3pi/2). Measured deviation 3.5e-04
    # (linear, seed 1); assert it is well above machine precision so
    # a future claim of "even for every setup" fails loudly.
    # Analytic 2-level counterexample (arena fidelity = squared
    # fidelity F^2): rho_A = |+><+|, rho_B = |psi><psi| with
    # |psi> = (|0> + i|1>)/sqrt(2), U = diag(1, e^{i phi}). Then
    # F^2(phi) = (1 - sin phi)/2, so F^2(pi/2) = 0 but F^2(3pi/2) = 1.
    # The evenness for identical cells rests on unitary invariance
    # plus the symmetry of the fidelity, which distinct cells lack.
    sw = sweep(1, "linear", n_phi=64, identical=False)
    s = sw["signed_flow"]
    dev = abs(s[16] - s[48])  # pi/2 vs 3pi/2 on a 64-point grid
    assert dev > 1e-4


def test_binary_identical_positive_lobe():
    # REGRESSION (review v0.2.0): the "drive negative for every
    # delta_phi != 0" claim holds for the linear generator only.
    # The binary generator has a positive lobe at small phase:
    # mean S(2pi/64) = +2.45e-05 across seeds (20 of 63 grid points
    # positive), because the binary fidelity stays above the 0.90
    # threshold at small phase.
    vals = []
    for seed in SEEDS:
        sw = sweep(seed, "binary", n_phi=64, identical=True)
        vals.append(sw["signed_flow"][1])  # phi = 2pi/64
    assert sum(vals) / len(vals) > 0


def test_binary_identical_fidelity_above_threshold_small_phi():
    # REGRESSION (review v0.2.0): the binary fidelity stays above the
    # 0.90 resonance threshold at small phase (mean F = 0.9990 at
    # 2pi/64, 0.96 at pi/4) and only falls below it at larger phase
    # (mean F = 0.79 at pi/2) -- so the threshold CAN be crossed by a
    # phase rotation under the binary generator.
    for idx, floor in ((1, 0.95), (8, 0.90)):
        vals = []
        for seed in SEEDS:
            sw = sweep(seed, "binary", n_phi=64, identical=True)
            vals.append(sw["fidelity"][idx])
        assert sum(vals) / len(vals) > floor
    vals = []
    for seed in SEEDS:
        sw = sweep(seed, "binary", n_phi=64, identical=True)
        vals.append(sw["fidelity"][16])  # pi/2
    assert sum(vals) / len(vals) < 0.90
