# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP)

"""Tests for the tunneling measurement (claims (h)-(k)).

The full experiment runs once per module (module-scoped fixture,
~90 s); the verdict tests read the measured numbers from it, and the
structural tests use small direct calls. Measured values are pinned
with tolerance, following the package pattern (test_harmonics.py).
"""

import math

import numpy as np
import pytest

from josephson_membrane.tunneling import (
    HEIGHTS,
    WIDTHS,
    barrier_links,
    link_coupling,
    log_linear_fit,
    measure_transmission,
    run_tunneling_experiment,
)


@pytest.fixture(scope="module")
def tun_exp():
    return run_tunneling_experiment()


# --- Wiring: the barrier arguments reach their named slots -----------

def test_barrier_links_centred():
    assert barrier_links(4, 0) == frozenset()
    assert barrier_links(4, 1) == frozenset({1})
    assert barrier_links(4, 2) == frozenset({0, 1})
    assert barrier_links(4, 3) == frozenset({0, 1, 2})
    with pytest.raises(ValueError):
        barrier_links(4, 4)


def test_link_coupling_scales_only_barrier_links():
    assert link_coupling(True, 0.75, 0.02) == pytest.approx(0.005)
    assert link_coupling(False, 0.75, 0.02) == pytest.approx(0.02)
    assert link_coupling(True, 0.0, 0.02) == pytest.approx(0.02)
    assert link_coupling(True, 1.0, 0.02) == pytest.approx(0.0)


def test_log_linear_fit_exact_line():
    fit = log_linear_fit(np.array([0.0, 1.0, 2.0, 3.0]),
                         np.log(np.array([1.0, 0.25, 0.0625, 0.015625])))
    assert fit["slope"] == pytest.approx(math.log(0.25))
    assert fit["r2"] == pytest.approx(1.0)


def test_log_linear_fit_flat_line_has_no_r2():
    fit = log_linear_fit(np.array([0.0, 1.0, 2.0]),
                         np.log(np.array([1.0, 1.0, 1.0])))
    assert fit["slope"] == pytest.approx(0.0)
    assert math.isnan(fit["r2"])


# --- Structural checks (not claims): the apparatus behaves ----------

def test_fully_cut_chain_transmits_nothing():
    # Height 1.0 cuts every barrier link: the detector's evolution is
    # bit-identical with and without injection, so T is exactly 0.
    r = measure_transmission(0, 1.0, 2)
    assert r["transmission"] == 0.0
    assert r["injected_deviation"] > 0.0


def test_zero_height_barrier_is_a_noop(tun_exp):
    # At height 0 the barrier links are full links, so the height
    # sweep's h=0 row and the width sweep's w=0 row are the same
    # computation and must agree exactly.
    h0 = tun_exp["height_sweep"][0]["mean_transmission"]
    w0 = tun_exp["width_sweep"][0]["mean_transmission"]
    assert h0 == w0
    assert h0 > 0.0


def test_injection_is_nonzero_and_detector_deviation_finite(tun_exp):
    for row in tun_exp["height_sweep"] + tun_exp["width_sweep"]:
        for t in row["per_seed_transmission"]:
            assert math.isfinite(t) and t >= 0.0
    assert tun_exp["height_sweep"][0]["mean_transmission"] > 1e-6


# --- Directional assertions (skill: freeze direction before measuring)

def test_transmission_falls_with_height_on_every_seed(tun_exp):
    for s in range(5):
        curve = [row["per_seed_transmission"][s]
                 for row in tun_exp["height_sweep"]]
        assert all(curve[i + 1] < curve[i] for i in range(len(curve) - 1)), (
            f"seed {s} height curve not monotone: {curve}"
        )


def test_transmission_falls_with_width_on_every_seed(tun_exp):
    for s in range(5):
        curve = [row["per_seed_transmission"][s]
                 for row in tun_exp["width_sweep"]]
        assert all(curve[i + 1] < curve[i] for i in range(len(curve) - 1)), (
            f"seed {s} width curve not monotone: {curve}"
        )


def test_a_rise_in_height_reaches_the_coupling_slot():
    # Same seed, same width, two heights: the higher barrier must
    # transmit less. This is the wiring check that the height argument
    # reaches the per-link coupling and not another slot.
    lo = measure_transmission(1, 0.0, 2)
    hi = measure_transmission(1, 0.99, 2)
    assert hi["transmission"] < lo["transmission"]
    assert hi["transmission"] / lo["transmission"] < 0.01


# --- Claims (h)-(k): frozen thresholds, measured numbers -------------

def test_claim_h_height_suppression(tun_exp):
    v = tun_exp["verdicts"]["h"]
    assert v["passed"] is True
    assert v["monotone_decreasing"] is True
    assert v["suppression_ratio"] == pytest.approx(1.12e-04, rel=0.05)
    assert v["suppression_ratio"] <= 0.20


def test_claim_i_width_exponential(tun_exp):
    v = tun_exp["verdicts"]["i"]
    assert v["passed"] is True
    assert v["log_linear_fit"]["slope"] == pytest.approx(-1.353, abs=0.01)
    assert v["log_linear_fit"]["r2"] == pytest.approx(0.999996, abs=0.001)
    assert v["log_linear_fit"]["slope"] < -0.10
    assert v["log_linear_fit"]["r2"] >= 0.90


def test_claim_j_graded_not_a_step(tun_exp):
    v = tun_exp["verdicts"]["j"]
    assert v["passed"] is True
    assert v["graded_ratio"] == pytest.approx(0.260, abs=0.01)
    assert 0.20 < v["graded_ratio"] < 0.80


def test_claim_k_wide_barrier_suppression(tun_exp):
    v = tun_exp["verdicts"]["k"]
    assert v["passed"] is True
    assert v["suppression_ratio"] == pytest.approx(1.21e-06, rel=0.05)
    assert v["suppression_ratio"] <= 0.20
    # The wide-barrier point sits at the numerical floor (detector
    # deviations ~1e-11), so its value is an upper bound on the true
    # transmission -- which only strengthens the <= 0.20 verdict.
    assert tun_exp["wide_barrier"]["mean_transmission"] < 1e-9


def test_transmission_factorises_as_one_minus_height_to_the_width(tun_exp):
    # Diagnostic, not a claim: T(h, w) ~= T(0) * (1 - h)^w to within
    # 15% on every mean point. The exponential-in-width shape is a
    # direct consequence of the drive being linear in the per-link
    # coupling -- a property of the designed rule, not a discovery.
    h0 = tun_exp["height_sweep"][0]["mean_transmission"]
    for row in tun_exp["height_sweep"]:
        pred = h0 * (1.0 - row["height"]) ** 2
        assert row["mean_transmission"] / pred == pytest.approx(1.0, rel=0.15)
    w0 = tun_exp["width_sweep"][0]["mean_transmission"]
    for row in tun_exp["width_sweep"]:
        pred = w0 * 0.25 ** row["width"]
        assert row["mean_transmission"] / pred == pytest.approx(1.0, rel=0.15)


def test_peak_time_never_advances_with_the_barrier(tun_exp):
    # Diagnostic: within each seed the arrival step t_peak never gets
    # earlier as the barrier strengthens or widens -- for seeds 0, 1,
    # 3, 4 it is barrier-independent (97, 100, 95, 63), and for seed 2
    # it is delayed, 48 -> 66, at the stronger/wider barrier. In this
    # window the barrier mostly scales the amplitude, and where it
    # moves the peak it delays it.
    for s in range(5):
        peaks_h = [row["per_seed_t_peak"][s]
                   for row in tun_exp["height_sweep"]]
        peaks_w = [row["per_seed_t_peak"][s]
                   for row in tun_exp["width_sweep"]]
        assert all(peaks_h[i + 1] >= peaks_h[i]
                   for i in range(len(peaks_h) - 1))
        assert all(peaks_w[i + 1] >= peaks_w[i]
                   for i in range(len(peaks_w) - 1))
