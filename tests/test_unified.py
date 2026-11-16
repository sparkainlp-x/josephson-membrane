# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP)

"""Tests for the unified claim (unified.py): six statements (five
verdicts), no headline; the membrane verdict is FAIL with the measured
numbers; the odd variant FAILS on its primary setup but PASSES for
binary/distinct; the ramp carries no verdict."""

import math

import pytest

from josephson_membrane.unified import HEADLINE_ALLOWED, unified_claim


@pytest.fixture(scope="module")
def result():
    return unified_claim()


def test_unified_structure(result):
    assert result["n_claims"] == 5
    assert result["n_passed"] == 3
    assert result["overall_pass"] is False
    assert result["headline_allowed"] is False
    assert HEADLINE_ALLOWED is False
    assert [c["id"] for c in result["claims"]] == ["a", "b", "c", "d", "e", "f"]


def test_law_claims_pass(result):
    by_id = {c["id"]: c for c in result["claims"]}
    assert by_id["a"]["passed"] is True
    assert by_id["a"]["max_error"] < 1e-12
    assert by_id["b"]["passed"] is True
    assert by_id["b"]["constant_relative_error"] < 1e-12
    assert by_id["c"]["passed"] is True
    assert by_id["c"]["max_error"] < 1e-12


def test_membrane_claim_fails_with_measured_numbers(result):
    claim = {c["id"]: c for c in result["claims"]}["d"]
    assert claim["passed"] is False
    primary = claim["primary"]
    fit = primary["fit"]
    # The bar is R^2 >= 0.95 and max|residual|/|A| <= 0.05. Measured:
    # R^2 is hugely negative and the residual is ~12x the amplitude.
    assert fit["r2"] < 0.5
    assert fit["max_residual_fraction"] > 0.5
    assert primary["verdict"]["passed"] is False
    # The best-fit phase is ~pi/2: the best sinusoid is a cosine, not
    # the Josephson sine (delta = 0).
    assert abs(fit["delta"] - math.pi / 2) < 0.1
    # Every per-seed fit fails too.
    assert all(r2 < 0.5 for r2 in primary["per_seed_r2"])
    assert all(
        frac > 0.5 for frac in primary["per_seed_max_residual_fraction"]
    )


def test_membrane_cardinal_flows_never_positive(result):
    # The Josephson law swings +I_c at pi/2 and -I_c at 3*pi/2. The
    # membrane flow is zero at 0 and negative everywhere else -- no
    # sign reversal, no positive lobe. (Primary setup: identical
    # cells, linear generator.)
    claim = {c["id"]: c for c in result["claims"]}["d"]
    card = claim["primary"]["cardinal_signed_flow"]
    assert card["0"] == pytest.approx(0.0, abs=1e-12)
    assert card["pi/2"] < 0.0
    assert card["pi"] < 0.0
    assert card["3pi/2"] < 0.0
    assert card["pi/2"] == pytest.approx(card["3pi/2"], abs=1e-9)
    # Even in delta_phi: pi/2 and 3*pi/2 are identical, as the
    # structural evenness requires. Fidelity is 1 at delta_phi = 0
    # (identical cells) and stays above 0.5 everywhere else.
    assert claim["primary"]["mean_fidelity_range"][0] > 0.5
    assert claim["primary"]["mean_fidelity_range"][1] == pytest.approx(1.0, abs=1e-9)


def test_supporting_setups_also_fail(result):
    claim = {c["id"]: c for c in result["claims"]}["d"]
    for key, entry in claim["supporting"].items():
        assert entry["verdict"]["passed"] is False, key


def test_odd_variant_primary_fails_with_measured_numbers(result):
    # Claim (e), primary setup: FAIL by a wide margin, but with the
    # Josephson phase (delta = 0), unlike (d)'s cosine.
    claim = {c["id"]: c for c in result["claims"]}["e"]
    assert claim["passed"] is False
    fit = claim["primary"]["flow_fit"]
    assert fit["r2"] < 0.5
    assert fit["max_residual_fraction"] > 0.5
    assert abs(fit["delta"]) < 0.1
    # The literal two-sided imaginary part is zero (machine precision).
    assert claim["primary"]["max_two_sided_imag"] < 1e-12


def test_odd_variant_binary_distinct_passes(result):
    # Claim (e), supporting setup: the binary/distinct flow PASSES
    # the pre-registered rule, and the binary drive is an exact sine.
    claim = {c["id"]: c for c in result["claims"]}["e"]
    passing = claim["supporting"]["binary_distinct"]
    assert passing["verdict"]["passed"] is True
    assert passing["flow_fit"]["r2"] > 0.95
    assert passing["flow_fit"]["max_residual_fraction"] < 0.05
    assert abs(passing["flow_fit"]["delta"]) < 0.05
    binary_drive = claim["supporting"]["binary_identical"]["drive_fit"]
    assert binary_drive["r2"] > 0.999
    assert abs(binary_drive["delta"]) < 0.01


def test_ramp_carries_no_verdict_and_reports_shape(result):
    # Claim (f): no verdict is forced; the shape is reported. The
    # v0.1.0 flow under the linear ramp is a notch train (mean
    # ~-0.0016, max ~0, little power at the ramp frequency), while
    # the binary variant drive is a pure sine at the ramp rate.
    claim = {c["id"]: c for c in result["claims"]}["f"]
    assert claim["passed"] is None
    lin = claim["setups"]["linear_identical"]["shape"]
    assert abs(lin["flow_v010"]["max"]) < 1e-6
    assert lin["flow_v010"]["mean"] < -0.001
    assert lin["flow_v010"]["power_fraction_at_ramp"] < 0.2
    assert lin["flow_v010"]["sine_r2_at_ramp"] < 0
    binary = claim["setups"]["binary_identical"]["shape"]
    assert binary["drive_q"]["power_fraction_at_ramp"] > 0.99
    assert binary["drive_q"]["sine_r2_at_ramp"] > 0.999
