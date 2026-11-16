"""Odd-function variant (claim (e)) and phase-ramp analogue (claim (f)).

The v0.1.0 result was a structural FAIL: the membrane flow depends on
the fidelity, an overlap magnitude, so it is EVEN in delta_phi, while
the Josephson current is ODD. A Josephson-like response needs the
interaction to depend on a signed, odd function of phase. This module
implements that variant.

The drive. The literal quantity in the task brief, Im Tr(rho_A
rho_B(phi)) with rho_B(phi) = U rho_B U^dagger, is identically zero:
the trace of a product of two Hermitian matrices is real. (The test
suite asserts this.) The quantity that is actually odd is the
ONE-SIDED overlap

    Q(phi) = Im Tr(rho_A U(phi) rho_B),

with rho_B unrotated. Q is odd in phi for identical cells
(rho_A = rho_B): then Q(phi) = sum_ij |rho_ij|^2 sin(j phi) for the
linear generator, manifestly odd, zero at phi = 0, and changing sign.
For distinct cells Q is NOT odd in general (measured deviation
~1e-3; the evenness proof fails because Tr(rho_A U rho_B) is not
symmetric under phi -> -phi when rho_A != rho_B). The variant flow is

    S_odd(phi) = drive(phi) * ||rho_B(phi) - rho_A||_F,
    drive(phi) = coupling * Q(phi),

with the same coupling (0.02), lattice, seeds 0-4, and 64-point sweep
grid as v0.1.0, so the variant is directly comparable to claim (d).

Pre-registered verdict rule (identical to v0.1.0, recorded here before
measuring): fit A sin(delta_phi + delta) with no DC term; the sine law
holds iff R^2 >= 0.95 AND max|residual| <= 0.05 |A|; report delta
(Josephson is delta = 0). The primary setup is identical cells with
the linear generator, mirroring claim (d). The other three setups
(binary generator, distinct cells) are supporting measurements, and
the drive fit is reported as a diagnostic alongside the flow fit.

The ramp (claim (f)). The phase is driven linearly in time,
phi_t = 2 pi t / 64 (one cycle per 64 steps, 256 steps = 4 cycles),
and the flow's time series is measured: the v0.1.0 flow S(t), the
variant flow S_odd(t), and the variant drive Q(t). No verdict is
forced for (f): the measured shape is reported -- FFT power at the
ramp frequency and a sine fit at the ramp rate, as diagnostics.
"""

from __future__ import annotations

import numpy as np

from membrane_arena_3x3.core import DIM, fidelity
from membrane_arena_3x3.lattice import LatticeConfig, MembraneLattice

from .membrane_phase import (
    COUPLING,
    N_PHI,
    RESONANCE_THRESHOLD,
    SEEDS,
    lattice_provenance,
    pairwise_flow,
    phase_unitary,
    sine_fit,
    sine_law_verdict,
)

# Phase ramp: one cycle per 64 steps, 256 steps = 4 cycles.
RAMP_POINTS_PER_CYCLE = 64
RAMP_N_STEPS = 256


def overlap_imag_two_sided(rho_a: np.ndarray, rho_b_phi: np.ndarray) -> float:
    """Im Tr(rho_A rho_B(phi)) with rho_B(phi) = U rho_B U^dagger.

    Identically zero: the trace of a product of two Hermitian matrices
    is real. This is the literal quantity in the task brief; it cannot
    drive anything, which is why the variant uses the one-sided
    overlap below.
    """
    return float(np.imag(np.trace(rho_a @ rho_b_phi)))


def overlap_imag(
    rho_a: np.ndarray, rho_b: np.ndarray, phi: float, generator: str = "linear"
) -> float:
    """Q(phi) = Im Tr(rho_A U(phi) rho_B), the one-sided overlap.

    Odd in phi for identical cells (rho_A = rho_B), zero at phi = 0,
    changes sign. Not odd in general for distinct cells.
    """
    U = phase_unitary(phi, DIM, generator)
    return float(np.imag(np.trace(rho_a @ U @ rho_b)))


def odd_pairwise_flow(
    rho_a: np.ndarray,
    rho_b: np.ndarray,
    phi: float,
    coupling: float = COUPLING,
    generator: str = "linear",
) -> dict:
    """Variant interaction: drive from the imaginary part of the overlap.

    Returns Q(phi), the drive, the distance, the signed variant flow,
    the (identically zero) two-sided imaginary part, and the fidelity
    for reference. The signed flow keeps the same structure as the
    v0.1.0 flow -- drive times distance -- with only the drive
    changed from coupling * (fidelity - 0.90) to coupling * Q(phi).
    """
    U = phase_unitary(phi, DIM, generator)
    rho_b_phi = U @ rho_b @ U.conj().T
    q = overlap_imag(rho_a, rho_b, phi, generator)
    drive = coupling * q
    dist = float(np.linalg.norm(rho_b_phi - rho_a))
    return {
        "q": q,
        "drive": drive,
        "distance": dist,
        "signed_flow": drive * dist,
        "two_sided_imag": overlap_imag_two_sided(rho_a, rho_b_phi),
        "fidelity": fidelity(rho_a, rho_b_phi),
    }


def odd_sweep(
    seed: int,
    generator: str = "linear",
    identical: bool = True,
    n_phi: int = N_PHI,
) -> dict:
    """Sweep the variant flow over delta_phi on the real lattice."""
    rng = np.random.default_rng(seed)
    lat = MembraneLattice(LatticeConfig(), rng)
    rho_a = lat.state[0, 0]
    rho_b = rho_a if identical else lat.state[0, 1]
    phis = np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False)
    out = {
        "phis": phis,
        "q": np.zeros(n_phi),
        "drive": np.zeros(n_phi),
        "signed_flow": np.zeros(n_phi),
        "distance": np.zeros(n_phi),
        "two_sided_imag": np.zeros(n_phi),
        "fidelity": np.zeros(n_phi),
    }
    for i, phi in enumerate(phis):
        res = odd_pairwise_flow(rho_a, rho_b, float(phi), generator=generator)
        for key in ("q", "drive", "signed_flow", "distance", "two_sided_imag", "fidelity"):
            out[key][i] = res[key]
    return out


def run_odd_experiment(seeds=SEEDS, n_phi: int = N_PHI) -> dict:
    """Run the variant sweep on all four setups.

    Returns the mean curves, the sine fits (flow and drive), and the
    pre-registered verdicts. Primary: linear generator, identical
    cells. The verdict rule is the v0.1.0 rule, imported, not
    redefined.
    """
    setups = {}
    for generator in ("linear", "binary"):
        for identical in (True, False):
            key = f"{generator}_{'identical' if identical else 'distinct'}"
            sweeps = [odd_sweep(s, generator, identical, n_phi) for s in seeds]
            mean_flow = np.mean([sw["signed_flow"] for sw in sweeps], axis=0)
            mean_drive = np.mean([sw["drive"] for sw in sweeps], axis=0)
            mean_q = np.mean([sw["q"] for sw in sweeps], axis=0)
            mean_dist = np.mean([sw["distance"] for sw in sweeps], axis=0)
            mean_fid = np.mean([sw["fidelity"] for sw in sweeps], axis=0)
            max_two_sided = max(
                float(np.max(np.abs(sw["two_sided_imag"]))) for sw in sweeps
            )
            phis = sweeps[0]["phis"]
            flow_fit = sine_fit(phis, mean_flow)
            drive_fit = sine_fit(phis, mean_drive)
            setups[key] = {
                "generator": generator,
                "identical": identical,
                "phis": phis.tolist(),
                "mean_flow": mean_flow.tolist(),
                "mean_drive": mean_drive.tolist(),
                "mean_q": mean_q.tolist(),
                "mean_distance": mean_dist.tolist(),
                "mean_fidelity": mean_fid.tolist(),
                "max_two_sided_imag": max_two_sided,
                "flow_fit": flow_fit,
                "drive_fit": drive_fit,
                "verdict": sine_law_verdict(flow_fit),
            }
    return {
        "lattice": lattice_provenance(),
        "seeds": list(seeds),
        "n_phi": n_phi,
        "coupling": COUPLING,
        "verdict_rule": (
            "fit A sin(delta_phi + delta), no DC term; sine law holds iff "
            "R^2 >= 0.95 and max|residual| <= 0.05 |A|; Josephson is delta = 0"
        ),
        "setups": setups,
    }


def phase_ramp(
    seed: int,
    generator: str = "linear",
    identical: bool = True,
    n_steps: int = RAMP_N_STEPS,
    points_per_cycle: int = RAMP_POINTS_PER_CYCLE,
) -> dict:
    """Drive the phase linearly in time and measure the flow series.

    phi_t = 2 pi t / points_per_cycle. Measures the v0.1.0 flow S(t),
    the variant flow S_odd(t), and the variant drive Q(t) on the real
    lattice, plus FFT power at the ramp frequency and a sine fit at
    the ramp rate (diagnostics, not a verdict).
    """
    rng = np.random.default_rng(seed)
    lat = MembraneLattice(LatticeConfig(), rng)
    rho_a = lat.state[0, 0]
    rho_b = rho_a if identical else lat.state[0, 1]
    phis = 2 * np.pi * np.arange(n_steps) / points_per_cycle
    s = np.zeros(n_steps)
    s_odd = np.zeros(n_steps)
    q = np.zeros(n_steps)
    for t, phi in enumerate(phis):
        U = phase_unitary(float(phi), DIM, generator)
        rho_b_phi = U @ rho_b @ U.conj().T
        s[t] = pairwise_flow(rho_a, rho_b_phi, COUPLING, RESONANCE_THRESHOLD)["signed_flow"]
        res = odd_pairwise_flow(rho_a, rho_b, float(phi), generator=generator)
        s_odd[t] = res["signed_flow"]
        q[t] = res["q"]
    return {
        "phis": phis,
        "flow_v010": s,
        "flow_odd": s_odd,
        "drive_q": q,
        "fft": {
            name: _ramp_spectrum(series, n_steps, points_per_cycle)
            for name, series in (
                ("flow_v010", s),
                ("flow_odd", s_odd),
                ("drive_q", q),
            )
        },
        "sine_at_ramp": {
            name: sine_fit(phis, series)
            for name, series in (
                ("flow_v010", s),
                ("flow_odd", s_odd),
                ("drive_q", q),
            )
        },
    }


def _ramp_spectrum(series: np.ndarray, n_steps: int, points_per_cycle: int) -> dict:
    """FFT power at the ramp frequency vs total power (demeaned)."""
    centered = series - series.mean()
    spectrum = np.abs(np.fft.rfft(centered)) ** 2
    ramp_bin = n_steps // points_per_cycle
    total = float(spectrum[1:].sum())
    peak_bin = int(np.argmax(spectrum[1:]) + 1)
    return {
        "ramp_bin": ramp_bin,
        "power_at_ramp": float(spectrum[ramp_bin]),
        "total_power": total,
        "power_fraction_at_ramp": float(spectrum[ramp_bin] / total) if total > 0 else 0.0,
        "peak_bin": peak_bin,
        "power_fraction_at_peak": float(spectrum[peak_bin] / total) if total > 0 else 0.0,
    }


def run_ramp_experiment(
    seeds=SEEDS,
    n_steps: int = RAMP_N_STEPS,
    points_per_cycle: int = RAMP_POINTS_PER_CYCLE,
) -> dict:
    """Run the phase ramp on all four setups. No verdict is forced.

    Returns the mean time series, the FFT power at the ramp frequency,
    and the sine fit at the ramp rate for each series. The "shape"
    field summarises what the series actually does.
    """
    setups = {}
    for generator in ("linear", "binary"):
        for identical in (True, False):
            key = f"{generator}_{'identical' if identical else 'distinct'}"
            ramps = [
                phase_ramp(s, generator, identical, n_steps, points_per_cycle)
                for s in seeds
            ]
            mean = {
                name: np.mean([r[name] for r in ramps], axis=0)
                for name in ("flow_v010", "flow_odd", "drive_q")
            }
            phis = ramps[0]["phis"]
            fft = {
                name: _ramp_spectrum(series, n_steps, points_per_cycle)
                for name, series in mean.items()
            }
            sine = {name: sine_fit(phis, series) for name, series in mean.items()}
            setups[key] = {
                "generator": generator,
                "identical": identical,
                "mean_flow_v010": mean["flow_v010"].tolist(),
                "mean_flow_odd": mean["flow_odd"].tolist(),
                "mean_drive_q": mean["drive_q"].tolist(),
                "fft": fft,
                "sine_at_ramp": sine,
                "shape": _describe_ramp_shape(mean, fft, sine),
            }
    return {
        "lattice": lattice_provenance(),
        "seeds": list(seeds),
        "n_steps": n_steps,
        "points_per_cycle": points_per_cycle,
        "setups": setups,
    }


def _describe_ramp_shape(mean: dict, fft: dict, sine: dict) -> dict:
    """Plain-language summary of what each ramp series actually does."""
    out = {}
    for name in ("flow_v010", "flow_odd", "drive_q"):
        series = mean[name]
        out[name] = {
            "mean": float(series.mean()),
            "min": float(series.min()),
            "max": float(series.max()),
            "power_fraction_at_ramp": fft[name]["power_fraction_at_ramp"],
            "peak_bin": fft[name]["peak_bin"],
            "sine_r2_at_ramp": sine[name]["r2"],
            "sine_delta_at_ramp": sine[name]["delta"],
        }
    return out
