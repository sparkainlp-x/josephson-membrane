# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP)

"""Unified claim: six statements, measured, no headline.

(a) DC law reproduces the cardinal table exactly: I_s = 0 at
    delta_phi = 0, pi, 2*pi; +I_c at pi/2; -I_c at 3*pi/2.
(b) AC frequency f = 2eV/h: the constant 2e/h reproduces the CODATA
    Josephson constant K_J = 483597.848416984 GHz/V to machine
    precision, and sample voltages give the tabulated frequencies.
(c) SQUID flux modulation: I_c(Phi) = 2*I_c*|cos(pi*Phi/Phi_0)| gives
    2*I_c at Phi = 0 and Phi_0, 0 at Phi_0/2, sqrt(2)*I_c at Phi_0/4.
(d) Membrane phase response vs the sine law (pre-registered rule in
    membrane_phase.py): the pairwise signed flow S(delta_phi) on
    identical cells, mean over seeds 0..4, fitted by A*sin(delta_phi
    + delta) with no DC term, passes iff R^2 >= 0.95 AND
    max|residual| <= 0.05*|A|.
(e) Odd-function variant vs the sine law (same pre-registered rule,
    odd_phase.py): the variant flow S_odd(delta_phi), whose drive is
    the imaginary part of the one-sided overlap Im Tr(rho_A U rho_B),
    fitted the same way on identical cells with the linear generator.
(f) Phase-ramp analogue (odd_phase.py): the phase driven linearly in
    time, phi_t = 2 pi t / 64; the flow's time series is measured and
    its shape reported -- FFT power at the ramp frequency and a sine
    fit at the ramp rate. No verdict is forced for (f); "passed" is
    None.

These are measurements on synthetic classical code plus exact
evaluations of established laws. They are not a proof about physical
reality; see REPORT.md, "Honest limits". No headline claim is emitted.
"""

from __future__ import annotations

import json
import math

from .josephson import (
    FLUX_QUANTUM_WB,
    JOSEPHSON_CONSTANT_GHZ_PER_MV,
    ac_frequency_ghz,
    ac_frequency_hz,
    dc_cardinal_table,
    squid_cardinal_table,
)
from .membrane_phase import run_membrane_experiment
from .odd_phase import run_odd_experiment, run_ramp_experiment

HEADLINE_ALLOWED = False


def claim_dc_law() -> dict:
    """(a) DC law cardinal table, exact."""
    table = dc_cardinal_table(1.0)
    expected = {"0": 0.0, "pi/2": 1.0, "pi": 0.0, "3pi/2": -1.0, "2pi": 0.0}
    max_err = max(abs(table[k] - v) for k, v in expected.items())
    # Oddness and 2*pi-periodicity, exact to machine precision.
    from .josephson import dc_current
    odd_err = max(abs(dc_current(p) + dc_current(-p)) for p in (0.3, 1.0, 2.0, 3.0))
    per_err = max(
        abs(dc_current(p) - dc_current(p + 2 * math.pi))
        for p in (0.3, 1.0, 2.0, 3.0)
    )
    passed = max_err < 1e-12 and odd_err < 1e-12 and per_err < 1e-12
    return {
        "id": "a",
        "name": "DC law cardinal table",
        "table": table,
        "max_error": max_err,
        "oddness_error": odd_err,
        "periodicity_error": per_err,
        "passed": bool(passed),
    }


def claim_ac_frequency() -> dict:
    """(b) AC frequency f = 2eV/h: constant and sample values."""
    # CODATA 2018 Josephson constant, GHz/V.
    codata_kj_ghz_per_v = 483597.848416984
    computed = JOSEPHSON_CONSTANT_GHZ_PER_MV * 1000.0
    rel_err = abs(computed - codata_kj_ghz_per_v) / codata_kj_ghz_per_v
    # Sample expectations are computed from the CODATA constant
    # (f = K_J * V), not hand-typed: K_J in Hz/V is 483597848416984,
    # so 1 uV -> 483597848.417 Hz, 1 mV -> 483.597848 GHz.
    codata_kj_hz_per_v = codata_kj_ghz_per_v * 1e9
    samples = {
        "1 uV": (1e-6, ac_frequency_hz(1e-6), codata_kj_hz_per_v * 1e-6),
        "10 uV": (1e-5, ac_frequency_hz(1e-5), codata_kj_hz_per_v * 1e-5),
        "1 mV": (1e-3, ac_frequency_hz(1e-3), codata_kj_hz_per_v * 1e-3),
        "1 V": (1.0, ac_frequency_hz(1.0), codata_kj_hz_per_v * 1.0),
    }
    sample_rows = {
        label: {"voltage_V": v, "frequency_Hz": f, "expected_Hz": e,
                "rel_error": abs(f - e) / e}
        for label, (v, f, e) in samples.items()
    }
    max_sample_err = max(r["rel_error"] for r in sample_rows.values())
    passed = rel_err < 1e-12 and max_sample_err < 1e-12
    return {
        "id": "b",
        "name": "AC frequency f = 2eV/h",
        "josephson_constant_GHz_per_mV": JOSEPHSON_CONSTANT_GHZ_PER_MV,
        "codata_KJ_GHz_per_V": codata_kj_ghz_per_v,
        "constant_relative_error": rel_err,
        "samples": sample_rows,
        "max_sample_relative_error": max_sample_err,
        "passed": bool(passed),
    }


def claim_squid() -> dict:
    """(c) SQUID flux modulation cardinal points."""
    table = squid_cardinal_table(1.0)
    expected = {"0": 2.0, "Phi0/4": math.sqrt(2.0), "Phi0/2": 0.0, "Phi0": 2.0}
    max_err = max(abs(table[k] - v) for k, v in expected.items())
    passed = max_err < 1e-12
    return {
        "id": "c",
        "name": "SQUID flux modulation",
        "flux_quantum_Wb": FLUX_QUANTUM_WB,
        "table": table,
        "max_error": max_err,
        "passed": bool(passed),
    }


def claim_membrane_phase() -> dict:
    """(d) Membrane phase response vs the sine law (pre-registered)."""
    exp = run_membrane_experiment()
    primary = exp["linear_identical"]
    fit = primary["fit"]
    verdict = primary["verdict"]
    supporting = {}
    for key in ("binary_identical", "linear_distinct", "binary_distinct"):
        entry = exp[key]
        supporting[key] = {
            "fit": entry["fit"],
            "fit_with_dc": entry["fit_with_dc"],
            "verdict": entry["verdict"],
            "cardinal_signed_flow": entry["cardinal_signed_flow"],
        }
    return {
        "id": "d",
        "name": "Membrane phase response vs sine law",
        "primary": {
            "setup": "identical cells, linear generator, mean over seeds 0..4",
            "fit": fit,
            "fit_with_dc": primary["fit_with_dc"],
            "verdict": verdict,
            "cardinal_signed_flow": primary["cardinal_signed_flow"],
            "per_seed_r2": primary["per_seed_r2"],
            "per_seed_delta": primary["per_seed_delta"],
            "per_seed_max_residual_fraction":
                primary["per_seed_max_residual_fraction"],
            "mean_fidelity_range": [
                float(primary["mean_fidelity"].min()),
                float(primary["mean_fidelity"].max()),
            ],
            "mean_lattice_delta_range": [
                float(primary["mean_lattice_delta"].min()),
                float(primary["mean_lattice_delta"].max()),
            ],
            "step_change_mean": primary["step_change_mean"],
            "step_change_range": primary["step_change_range"],
        },
        "supporting": supporting,
        "provenance": exp["provenance"],
        "passed": bool(verdict["passed"]),
    }


def claim_odd_phase() -> dict:
    """(e) Odd-function variant vs the sine law (pre-registered rule)."""
    exp = run_odd_experiment()
    primary = exp["setups"]["linear_identical"]
    verdict = primary["verdict"]
    supporting = {}
    for key in ("binary_identical", "linear_distinct", "binary_distinct"):
        entry = exp["setups"][key]
        supporting[key] = {
            "flow_fit": entry["flow_fit"],
            "drive_fit": entry["drive_fit"],
            "verdict": entry["verdict"],
            "max_two_sided_imag": entry["max_two_sided_imag"],
        }
    return {
        "id": "e",
        "name": "Odd-function variant vs sine law",
        "primary": {
            "setup": "identical cells, linear generator, mean over seeds 0..4",
            "flow_fit": primary["flow_fit"],
            "drive_fit": primary["drive_fit"],
            "verdict": verdict,
            "max_two_sided_imag": primary["max_two_sided_imag"],
        },
        "supporting": supporting,
        "verdict_rule": exp["verdict_rule"],
        "provenance": exp["lattice"],
        "passed": bool(verdict["passed"]),
    }


def claim_phase_ramp() -> dict:
    """(f) Phase-ramp analogue: measured shape, no forced verdict."""
    exp = run_ramp_experiment()
    setups = {}
    for key, entry in exp["setups"].items():
        setups[key] = {
            "shape": entry["shape"],
            "fft": entry["fft"],
            "sine_at_ramp": entry["sine_at_ramp"],
        }
    return {
        "id": "f",
        "name": "Phase-ramp analogue (AC)",
        "note": (
            "No verdict is forced: the shape is reported. The v0.1.0 "
            "flow under a linear phase ramp is a notch train (the same "
            "shape as claim (d), repeated); the binary variant drive "
            "is a pure sine at the ramp rate."
        ),
        "n_steps": exp["n_steps"],
        "points_per_cycle": exp["points_per_cycle"],
        "setups": setups,
        "provenance": exp["lattice"],
        "passed": None,
    }


def unified_claim() -> dict:
    """Run all six statements; return measured numbers, no headline."""
    claims = [
        claim_dc_law(),
        claim_ac_frequency(),
        claim_squid(),
        claim_membrane_phase(),
        claim_odd_phase(),
        claim_phase_ramp(),
    ]
    verdicts = [c for c in claims if c["passed"] is not None]
    return {
        "claims": claims,
        "n_passed": sum(1 for c in verdicts if c["passed"]),
        "n_claims": len(verdicts),
        "overall_pass": all(c["passed"] for c in verdicts),
        "headline_allowed": HEADLINE_ALLOWED,
    }


def _cli() -> None:
    """Print the unified claim as JSON (python -m josephson_membrane.unified)."""
    result = unified_claim()
    print(json.dumps(result, indent=2, default=float))


if __name__ == "__main__":
    _cli()
