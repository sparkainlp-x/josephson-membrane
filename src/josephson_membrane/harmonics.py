# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP)

"""Harmonics vs incoherence in the failing corner (claim (g)).

The question (Jean, to Dr. Sinem Salva Ersoy): in the failing corner
-- the identical/linear regime, where the variant sine fit gives
R^2 = 0.037 -- is the response hiding higher harmonics, or are the
phases there just incoherent? A single sine fit cannot tell those
apart. This module builds the test that can: a Fourier decomposition
of the drive/response curve, densely sampled over phi in [0, 2*pi).

Why the rule looks the way it does. The linear generator's drive is
analytically a Fourier series over 31 frequencies (CHANGELOG v0.2.0:
"the odd drive is a 31-frequency Fourier series, not a single
sine"), so harmonic power can sit at ANY harmonic k <= 31, not only
at k = 2 or 3. The ramp experiment (claim (f)) already hints at this:
the linear variant drive's ramp-spectrum peaks at bin 96 of 128 and
the flow's at bin 120, i.e. the 24th and 30th harmonics of phi. A
rule that only looked at 2phi and 3phi would miss that. So the
verdict keys on CONCENTRATION: harmonics means the power is
concentrated in a few discrete bins (wherever they sit); incoherent
means it is spread broadband with no bin dominant.

PRE-REGISTERED VERDICT RULE (fixed before measuring; the thresholds
below are not moved after seeing the numbers):
  Let p1, p2, p3 be the power fractions at harmonics 1, 2, 3 (of the
  demeaned curve, bins 1..127 of a 256-point rfft), p_nonfund = 1 - p1
  the non-fundamental fraction, p_top3 the sum of the three largest
  single-bin fractions, and peak the largest single-bin fraction.
  The shuffled control is the same mean curve with its samples
  shuffled (fixed seed 12345), decomposed the same way.
  - "harmonics" iff p_top3 >= 0.50 AND p_nonfund >= 0.20 AND
    p_top3 >= 3.0 * (shuffled p_top3).
  - "incoherent" iff peak <= 0.10 AND p_top3 <= 0.30.
  - "neither" otherwise (this is where a single-sine-dominated curve
    lands: p1 ~= 1, so p_nonfund < 0.20 and peak = p1 > 0.10).
  The p2+p3 fraction is reported as a diagnostic alongside, but the
  verdict does not rest on it alone, for the reason above.

Claim (g) is a claim about the TEST's discriminating power, not
about the sine law. It passes iff (i) the failing corner
(linear/identical flow) gets a definite verdict, "harmonics" or
"incoherent"; (ii) the binary/distinct control flow is "neither"
with p1 >= 0.90 (a single sine, as its R^2 = 0.9999 requires); and
(iii) the shuffled control on the failing corner is "incoherent".
If the failing corner is "neither", (g) FAILS: the test cannot tell
them apart there.

Controls. The v0.1.0 flow (claim (d)) is EVEN in phi for identical
cells, so its spectrum should be cosine-like: the sine-part power
fraction (imaginary part of the rfft coefficients) should be small,
and the fundamental should dominate. The binary/distinct variant
flow is the single-sine control (R^2 = 0.9999).
"""

from __future__ import annotations

import numpy as np

from .membrane_phase import SEEDS, sweep
from .odd_phase import odd_sweep

# --- Pre-registered verdict thresholds (see module docstring) ---
N_HARM = 256  # dense grid points over [0, 2*pi)
SHUFFLE_SEED = 12345
HARMONICS_TOP3_MIN = 0.50
HARMONICS_NONFUND_MIN = 0.20
HARMONICS_FLOOR_FACTOR = 3.0
INCOHERENT_MAX_BIN = 0.10
INCOHERENT_TOP3_MAX = 0.30


def spectrum(y: np.ndarray) -> dict:
    """Fourier power fractions of a curve sampled on [0, 2*pi).

    The curve is demeaned (the DC term carries no harmonic
    information), then rfft'd. Bins 1..127 are harmonics 1..127 of
    phi; bin 128 is the Nyquist frequency and is reported separately,
    not counted in the harmonic fractions. Returns the per-harmonic
    fractions p1, p2, p3, the non-fundamental fraction, the top-3 bin
    fractions, the peak bin, and the sine-part power fraction (the
    fraction of power in the imaginary part of the coefficients --
    0 for an even function, 1 for a pure sine).
    """
    y = np.asarray(y, dtype=float)
    n = y.size
    if n < 8:
        raise ValueError("need at least 8 samples for a spectrum")
    centered = y - y.mean()
    coeffs = np.fft.rfft(centered)
    power = np.abs(coeffs) ** 2
    # Harmonic bins 1..n//2 - 1; bin n//2 is Nyquist.
    harm = power[1 : n // 2]
    total = float(harm.sum())
    nyquist = float(power[n // 2]) if n % 2 == 0 else 0.0
    if total <= 0:
        raise ValueError("curve has no AC power")
    frac = harm / total
    order = np.argsort(frac)[::-1]
    top3_idx = order[:3]
    sine_part = float((coeffs[1 : n // 2].imag ** 2).sum() / total)
    return {
        "n": int(n),
        "total_power": total,
        # Parseval total (up to the 1/N factor): interior rfft bins
        # each appear twice in the full FFT, the Nyquist bin once.
        # This is the quantity a permutation cannot change; the
        # harmonic-bin total alone is not conserved, because the
        # Nyquist power moves when samples are shuffled.
        "total_power_all_bins": 2.0 * total + nyquist,
        "nyquist_fraction": (
            nyquist / (2.0 * total + nyquist) if total + nyquist > 0 else 0.0
        ),
        "p1": float(frac[0]),
        "p2": float(frac[1]),
        "p3": float(frac[2]),
        "p2_plus_p3": float(frac[1] + frac[2]),
        "p_high": float(frac[3:].sum()),
        "p_nonfund": float(1.0 - frac[0]),
        "p_top3": float(frac[top3_idx].sum()),
        "top_bins": [int(i + 1) for i in top3_idx],
        "top_fractions": [float(frac[i]) for i in top3_idx],
        "peak_bin": int(np.argmax(frac) + 1),
        "peak_fraction": float(frac.max()),
        "sine_power_fraction": sine_part,
    }


def shuffled_spectrum(y: np.ndarray, seed: int = SHUFFLE_SEED) -> dict:
    """The spectrum of the same samples with their order shuffled.

    Shuffling destroys all phase structure while keeping the value
    distribution (hence total power) identical, so this is the
    broadband noise floor: power spread roughly uniformly over bins.
    """
    rng = np.random.default_rng(seed)
    shuffled = rng.permutation(np.asarray(y, dtype=float))
    return spectrum(shuffled)


def classify(spec: dict, shuffled_spec: dict) -> dict:
    """Apply the pre-registered rule to a spectrum and its control.

    Returns the verdict ("harmonics" / "incoherent" / "neither") and
    which conditions held, so a verdict can be audited, not trusted.
    """
    harmonics = bool(
        spec["p_top3"] >= HARMONICS_TOP3_MIN
        and spec["p_nonfund"] >= HARMONICS_NONFUND_MIN
        and spec["p_top3"] >= HARMONICS_FLOOR_FACTOR * shuffled_spec["p_top3"]
    )
    incoherent = bool(
        spec["peak_fraction"] <= INCOHERENT_MAX_BIN
        and spec["p_top3"] <= INCOHERENT_TOP3_MAX
    )
    if harmonics:
        verdict = "harmonics"
    elif incoherent:
        verdict = "incoherent"
    else:
        verdict = "neither"
    return {
        "verdict": verdict,
        "harmonics_conditions": {
            "p_top3 >= 0.50": bool(spec["p_top3"] >= HARMONICS_TOP3_MIN),
            "p_nonfund >= 0.20": bool(
                spec["p_nonfund"] >= HARMONICS_NONFUND_MIN
            ),
            "p_top3 >= 3x shuffled": bool(
                spec["p_top3"]
                >= HARMONICS_FLOOR_FACTOR * shuffled_spec["p_top3"]
            ),
        },
        "incoherent_conditions": {
            "peak <= 0.10": bool(
                spec["peak_fraction"] <= INCOHERENT_MAX_BIN
            ),
            "p_top3 <= 0.30": bool(spec["p_top3"] <= INCOHERENT_TOP3_MAX),
        },
    }


def _analyse_curve(y: np.ndarray) -> dict:
    """Spectrum + shuffled control + verdict for one mean curve."""
    spec = spectrum(y)
    shuf = shuffled_spectrum(y)
    return {
        "spectrum": spec,
        "shuffled_spectrum": shuf,
        "classification": classify(spec, shuf),
    }


def run_harmonics_experiment(
    seeds: tuple[int, ...] = SEEDS, n_phi: int = N_HARM
) -> dict:
    """Fourier-decompose the variant drive and flow on a dense grid.

    Mean curves over seeds 0..4 on n_phi points over [0, 2*pi), all
    four setups (linear/binary x identical/distinct), plus the v0.1.0
    control (mean signed flow and mean fidelity, linear/identical).
    The verdict is computed on the variant FLOW (the response); the
    drive spectrum is reported alongside as the explanation.
    """
    setups = {}
    for generator in ("linear", "binary"):
        for identical in (True, False):
            key = f"{generator}_{'identical' if identical else 'distinct'}"
            sweeps = [
                odd_sweep(s, generator, identical, n_phi) for s in seeds
            ]
            mean_q = np.mean([sw["q"] for sw in sweeps], axis=0)
            mean_flow = np.mean([sw["signed_flow"] for sw in sweeps], axis=0)
            setups[key] = {
                "generator": generator,
                "identical": identical,
                "drive": _analyse_curve(mean_q),
                "flow": _analyse_curve(mean_flow),
            }
    # v0.1.0 control: the (even) flow and the fidelity on the same
    # dense grid, linear generator, identical cells.
    v010 = [
        sweep(s, "linear", n_phi, True, step_change_grid=0) for s in seeds
    ]
    mean_flow_v010 = np.mean([sw["signed_flow"] for sw in v010], axis=0)
    mean_fid = np.mean([sw["fidelity"] for sw in v010], axis=0)
    control = {
        "flow": _analyse_curve(mean_flow_v010),
        "fidelity": _analyse_curve(mean_fid),
    }
    return {
        "n_phi": n_phi,
        "seeds": list(seeds),
        "shuffle_seed": SHUFFLE_SEED,
        "verdict_rule": {
            "harmonics": (
                "p_top3 >= 0.50 AND p_nonfund >= 0.20 AND "
                "p_top3 >= 3.0 * shuffled p_top3"
            ),
            "incoherent": "peak bin <= 0.10 AND p_top3 <= 0.30",
            "neither": "otherwise",
        },
        "setups": setups,
        "v010_control": control,
    }
