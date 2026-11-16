# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP)

"""Tests for the harmonics test (harmonics.py) and claim (g).

Two kinds of tests: (1) the spectrum/classifier on synthetic signals
-- a pure sine is a single tone, a sine plus a 2nd or 3rd harmonic is
"harmonics", white noise is "incoherent", a DC offset changes nothing,
and the sine-part power fraction separates even from odd functions;
(2) the measured outcome on the real lattice -- the failing corner
(identical cells, linear generator) is INCOHERENT, not hiding a
dominant harmonic; the binary/distinct control is a single sine; the
v0.1.0 control is exactly even.
"""

import numpy as np
import pytest

from josephson_membrane.harmonics import (
    classify,
    run_harmonics_experiment,
    shuffled_spectrum,
    spectrum,
)

PHIS = np.linspace(0.0, 2 * np.pi, 256, endpoint=False)


@pytest.fixture(scope="module")
def harm_exp():
    return run_harmonics_experiment()


def test_spectrum_pure_sine_is_single_tone():
    spec = spectrum(np.sin(PHIS))
    assert spec["p1"] == pytest.approx(1.0, abs=1e-9)
    assert spec["p2"] == pytest.approx(0.0, abs=1e-9)
    assert spec["sine_power_fraction"] == pytest.approx(1.0, abs=1e-9)
    assert classify(spec, shuffled_spectrum(np.sin(PHIS)))["verdict"] == "neither"


def test_spectrum_pure_cosine_is_even():
    # An even function has (numerically) zero sine-part power.
    spec = spectrum(np.cos(PHIS))
    assert spec["p1"] == pytest.approx(1.0, abs=1e-9)
    assert spec["sine_power_fraction"] == pytest.approx(0.0, abs=1e-9)


def test_spectrum_ignores_dc_offset():
    base = spectrum(np.sin(PHIS))
    shifted = spectrum(np.sin(PHIS) + 5.0)
    for key in ("p1", "p2", "p3", "p_top3", "peak_fraction"):
        assert shifted[key] == pytest.approx(base[key], abs=1e-12)


def test_second_harmonic_is_harmonics():
    # sin(phi) + 0.7 sin(2phi): 2nd-harmonic power 0.49/1.49 = 0.329.
    y = np.sin(PHIS) + 0.7 * np.sin(2 * PHIS)
    spec = spectrum(y)
    assert spec["p2"] == pytest.approx(0.49 / 1.49, abs=0.01)
    assert classify(spec, shuffled_spectrum(y))["verdict"] == "harmonics"


def test_third_harmonic_is_harmonics():
    y = np.sin(PHIS) + 0.7 * np.sin(3 * PHIS)
    spec = spectrum(y)
    assert spec["p3"] == pytest.approx(0.49 / 1.49, abs=0.01)
    assert classify(spec, shuffled_spectrum(y))["verdict"] == "harmonics"


def test_white_noise_is_incoherent():
    # Broadband by construction: no bin above 10%, top-3 below 30%.
    for seed in (0, 1):
        y = np.random.default_rng(seed).normal(size=256)
        spec = spectrum(y)
        assert spec["peak_fraction"] < 0.10
        assert spec["p_top3"] < 0.30
        assert classify(spec, shuffled_spectrum(y))["verdict"] == "incoherent"


def test_shuffled_keeps_total_power():
    # Shuffling redistributes power across bins but cannot change the
    # total demeaned power (Parseval: the sum of squared samples is
    # invariant under permutation). The harmonic-bin total alone is
    # NOT conserved -- the Nyquist bin (excluded from the harmonic
    # fractions) absorbs the difference -- so the comparison uses
    # total_power_all_bins (2x the interior bins + Nyquist, the
    # Parseval total in rfft terms).
    y = np.sin(PHIS) + 0.7 * np.sin(2 * PHIS)
    assert shuffled_spectrum(y)["total_power_all_bins"] == pytest.approx(
        spectrum(y)["total_power_all_bins"], rel=1e-12
    )


def test_spectrum_rejects_too_few_samples_and_constant():
    with pytest.raises(ValueError):
        spectrum(np.ones(4))
    with pytest.raises(ValueError):
        spectrum(np.full(256, 3.0))


def test_failing_corner_flow_is_incoherent(harm_exp):
    # The question: higher harmonics, or incoherent phases? Measured:
    # no harmonic dominates -- peak bin 2 carries 7.3%, top-3 16.8%,
    # 84% of the power sits above the 3rd harmonic, spread out.
    flow = harm_exp["setups"]["linear_identical"]["flow"]
    assert flow["classification"]["verdict"] == "incoherent"
    spec = flow["spectrum"]
    assert spec["p1"] < 0.10
    assert spec["p2_plus_p3"] < 0.20
    assert spec["peak_fraction"] < 0.10
    assert spec["p_top3"] < 0.30
    assert spec["p_high"] > 0.80
    # It is still a sine series (odd function): sine-part power 1.0.
    assert spec["sine_power_fraction"] > 0.99


def test_failing_corner_drive_is_incoherent_high_harmonics(harm_exp):
    # The drive's strongest bins are high (24th, 30th, 29th -- the
    # same high harmonics the ramp experiment already hinted at), but
    # each carries only ~5% of the power: no hidden dominant harmonic.
    drive = harm_exp["setups"]["linear_identical"]["drive"]
    assert drive["classification"]["verdict"] == "incoherent"
    spec = drive["spectrum"]
    assert min(spec["top_bins"]) >= 9
    assert spec["p_high"] > 0.85
    assert spec["p1"] < 0.05


def test_linear_distinct_also_incoherent(harm_exp):
    flow = harm_exp["setups"]["linear_distinct"]["flow"]
    assert flow["classification"]["verdict"] == "incoherent"
    assert flow["spectrum"]["p_high"] > 0.85


def test_binary_distinct_flow_is_single_sine(harm_exp):
    # The control: the passing regime is a single sine, p1 = 1.0000.
    flow = harm_exp["setups"]["binary_distinct"]["flow"]
    assert flow["classification"]["verdict"] == "neither"
    assert flow["spectrum"]["p1"] > 0.99


def test_binary_drives_are_exact_sines(harm_exp):
    for key in ("binary_identical", "binary_distinct"):
        drive = harm_exp["setups"][key]["drive"]
        assert drive["spectrum"]["p1"] > 0.999
        assert drive["classification"]["verdict"] == "neither"


def test_binary_identical_flow_near_miss_spectrum(harm_exp):
    # The R^2 = 0.92 near-miss is a single sine plus a small 2nd
    # harmonic from the distance notch: p1 = 0.9222, p2 = 0.0753.
    flow = harm_exp["setups"]["binary_identical"]["flow"]
    assert flow["classification"]["verdict"] == "neither"
    spec = flow["spectrum"]
    assert spec["p1"] > 0.90
    assert 0.05 < spec["p2"] < 0.10


def test_v010_control_is_even(harm_exp):
    # The v0.1.0 flow is EVEN in phi for identical cells, so its
    # sine-part power is exactly zero; its spectrum sits on the odd
    # harmonics 1, 3, 5 (a notch train, not a cosine).
    for part in ("flow", "fidelity"):
        spec = harm_exp["v010_control"][part]["spectrum"]
        assert spec["sine_power_fraction"] < 0.01
        assert spec["top_bins"] == [1, 3, 5]


def test_shuffled_control_is_incoherent(harm_exp):
    shuf = harm_exp["setups"]["linear_identical"]["flow"]["shuffled_spectrum"]
    assert classify(shuf, shuf)["verdict"] == "incoherent"
    assert shuf["peak_fraction"] < 0.10
    assert shuf["p_top3"] < 0.30
