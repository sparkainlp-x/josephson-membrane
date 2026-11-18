# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP)

"""Tunneling: transmission across a barrier in the membrane lattice (claims (h)-(k)).

The question. No transmission-vs-barrier measurement exists anywhere
in this work yet. A Josephson junction IS a tunneling barrier: the
supercurrent crosses an insulating gap that a classical particle
cannot cross, and the transmission falls exponentially with the
barrier's width and with the square root of its height. This module
asks whether the membrane lattice shows the same SHAPE: does the
fraction of an injected membrane disturbance that reaches a detector
on the far side of a barrier fall exponentially with barrier width
and height (tunneling-like), or linearly, or as an on/off step
(classical blocking)?

Why a chain. The 3x3 lattice's neighbour graph is dense: a cell's
neighbours are its whole row and its whole column, so every pair of
cells is directly connected or shares a neighbour, and no single
region can choke the flow -- a barrier column would be bypassed by
the direct links around it. A barrier is only definable on a graph
with a chokepoint, so the experiment runs on a 1-D chain of cells
taken from the REAL lattice: the states, Hamiltonians, configuration
and every evolution stage come from membrane_arena_3x3, and the
interaction term is the lattice's own pairwise formula
(drive = coupling * (fidelity - threshold), delta += drive *
(other - rho), averaged over neighbours -- the inner loop of
MembraneLattice.interaction_delta). What is new is only the neighbour
graph (a chain instead of the 3x3 grid) and a per-link coupling, so a
barrier can be a set of consecutive suppressed links. Nothing about
the dynamics is reimplemented: unitary_stage, dissipative_stage and
mirror_stage are the lattice's own methods, composed exactly as
propose_step composes them in legacy mode.

Setup. A chain of CHAIN_LENGTH cells (4) from one real
MembraneLattice (seeded; states and Hamiltonians are the lattice's
own, flattened row-major), run at coupling 0.10 -- the top of the
lattice's own validated range (LatticeConfig.coupling_max). The
default 0.02 is unusable here: the pilot showed the per-link
attenuation there is ~200x, so the injected deviation is ~2e-3 at
one hop, ~8e-6 at two hops, and below the numerical floor (~1e-11
in detector deviation) by the detector, and no barrier shape is
measurable at all. The chain is short (4 cells, 3 hops) for the
same reason: each hop attenuates the signal ~50x at coupling 0.10,
so a longer chain puts even the open-chain transmission below the
floor. Cell 0 is the source, cell 3 the detector. The barrier is a
block of WIDTH consecutive centred links whose coupling is scaled by
(1 - height): height 0 is no barrier, height 1 cuts the links
entirely. At t = 0 the source cell is rotated by the linear phase
unitary U(pi/2) (the same phase twist membrane_phase.py uses to
impose a phase difference), which preserves populations and trace
and changes only coherences. The lattice then runs N_STEPS steps of
dt = 0.05 (the config default). A control run with the same seed and
settings but no injection isolates the injected signal:
D(t) = ||rho_detector(t) - rho_detector_control(t)||_F, and the
transmission is T = max_t D(t) / D_inj, where D_inj is the injected
deviation ||rho_source(0) - rho_source_control(0)||_F. The peak
(rather than a fixed time) is used because the arrival time grows
with barrier width; the arrival step t_peak is reported alongside.
The coupling 0.10, chain length 4 and N_STEPS = 100 were fixed in
the pilot phase, before the measurement run, because they are what
puts the open-chain transmission measurably above the floor; every
claim below is a ratio or a shape, so no threshold depends on them.

PRE-REGISTERED VERDICT RULES (fixed before measuring; the thresholds
below are not moved after seeing the numbers). All claims use the
mean transmission over seeds 0..4. Every claim is a ratio or a
shape, so no threshold depends on the absolute scale of T.
  (h) HEIGHT: at width 2, the mean T strictly decreases with height
      over HEIGHTS = (0, 0.25, 0.5, 0.75, 0.9, 0.99), AND the
      suppression ratio T(0.99) / T(0) <= 0.20 (a nearly-closed
      barrier transmits at most a fifth of the open chain).
  (i) WIDTH EXPONENTIAL: at height 0.75, over widths 0..3, the
      least-squares fit of log(mean T) against width has slope
      < -0.10 per link AND R^2 >= 0.90 -- exponential (tunneling-like)
      falloff, not linear and not a step.
  (j) GRADED, NOT A STEP: at width 2, 0.20 < T(0.5) / T(0) < 0.80 --
      an intermediate height gives intermediate transmission, not an
      on/off step.
  (k) WIDE-BARRIER SUPPRESSION: at width 3, T(0.99) / T(0) <= 0.20 --
      the suppression of (h) also holds for the widest barrier.
Structural checks (asserted in tests, not claims): height 1.0 gives
T = 0 exactly (the chain is cut, so the detector never sees the
injection); height 0 gives the same T for every width (the barrier
is a no-op when its height is zero).

Honest limits. This is an engineered analogue on a synthetic
classical emulation, not a physical junction and not a new law: the
chain is a restriction of the lattice, not the 3x3 lattice itself;
the barrier is a designed coupling suppression, not a discovered
potential; the exponential shape, if found, is a property of this
interaction rule (drive linear in coupling, fidelity-gated), not
evidence that the membrane tunnels. The measurement is what it is:
four directional claims with frozen thresholds, reported PASS/FAIL
with the measured numbers.
"""

from __future__ import annotations

import numpy as np

from membrane_arena_3x3.core import DIM, fidelity, project_valid
from membrane_arena_3x3.lattice import LatticeConfig, MembraneLattice

from .membrane_phase import lattice_provenance, phase_unitary

# --- Pre-registered design constants (see module docstring) ---
CHAIN_LENGTH = 4
SOURCE = 0
DETECTOR = CHAIN_LENGTH - 1
N_STEPS = 100
SEEDS = (0, 1, 2, 3, 4)
HEIGHTS = (0.0, 0.25, 0.5, 0.75, 0.9, 0.99)
WIDTHS = (0, 1, 2, 3)
HEIGHT_SWEEP_WIDTH = 2
WIDTH_SWEEP_HEIGHT = 0.75
WIDE_WIDTH = 3
INJECTION_PHI = float(np.pi / 2)
# Coupling for the chain. The lattice default is 0.02, but the pilot
# run showed the per-link attenuation there is ~200x: the injected
# deviation is ~2e-3 at one hop, ~8e-6 at two hops, and below the
# numerical floor (~1e-11 in detector deviation) by the detector, so
# no barrier shape is measurable at all at the default. 0.10 is the
# top of the lattice's own validated range
# (LatticeConfig.coupling_max), and it lifts the signal into the
# measurable range. Fixed in the pilot phase, before measuring; every
# claim is a ratio or a shape, so no threshold depends on this value.
TUNNEL_COUPLING = 0.10

# --- Pre-registered verdict thresholds (frozen before measuring) ---
SUPPRESSION_RATIO_MAX = 0.20
WIDTH_LOG_SLOPE_MAX = -0.10
WIDTH_LOG_R2_MIN = 0.90
GRADED_RATIO_MIN = 0.20
GRADED_RATIO_MAX = 0.80


def barrier_links(length: int = CHAIN_LENGTH, width: int = 0) -> frozenset:
    """Return the set of link indices (k = link between cell k and k+1)
    that form the barrier: WIDTH consecutive links, centred in the
    chain. The barrier is the region between source and detector, so
    its links may include the entry and exit links.
    """
    n_links = length - 1
    if width <= 0:
        return frozenset()
    if width > n_links:
        raise ValueError(f"width {width} exceeds the {n_links} links")
    start = (n_links - width) // 2
    return frozenset(range(start, start + width))


def link_coupling(is_barrier: bool, height: float, coupling: float) -> float:
    """Coupling on one link: barrier links are scaled by (1 - height)."""
    return coupling * (1.0 - height) if is_barrier else coupling


def _chain_step(lat: MembraneLattice, states: list,
                flat_hamiltonians: np.ndarray, barriers: frozenset,
                height: float) -> list:
    """One chain step, composing the lattice's own stages exactly as
    propose_step does in legacy mode (unitary Euler step + dephasing
    + interaction + projection + mirror), except the interaction sums
    over chain neighbours with per-link coupling instead of the 3x3
    neighbour graph. The interaction term is the lattice's own
    pairwise formula: drive = coupling_link * (fidelity - threshold),
    delta += drive * (other - rho), averaged over neighbours. Each
    link's fidelity and drive are computed once and shared by its two
    endpoint cells, exactly as the two directed terms in
    interaction_delta share one fidelity evaluation per pair.
    """
    cfg = lat.config
    n = len(states)
    # Per-link interaction terms, computed once per link.
    link_drive = {}
    for link in range(n - 1):
        rho, other = states[link], states[link + 1]
        score = fidelity(rho, other)
        c_link = link_coupling(link in barriers, height, cfg.coupling)
        link_drive[link] = c_link * (score - cfg.resonance_threshold)
    nxt = []
    for i, rho in enumerate(states):
        dt_cell = cfg.dt * cfg.dt_multipliers[i] + cfg.phase_offsets[i]
        h = flat_hamiltonians[i]
        unitary = lat.unitary_stage(rho, h, dt_cell)
        dissipative = lat.dissipative_stage(rho, cfg.dephasing * cfg.dt)
        delta = (unitary - rho) + (dissipative - rho)
        interaction = np.zeros_like(rho)
        n_neigh = 0
        if i > 0:
            interaction += link_drive[i - 1] * (states[i - 1] - rho)
            n_neigh += 1
        if i < n - 1:
            interaction += link_drive[i] * (states[i + 1] - rho)
            n_neigh += 1
        interaction /= max(n_neigh, 1)
        delta = delta + interaction * cfg.dt
        cand = project_valid(rho + delta)
        nxt.append(lat.mirror_stage(cand, cfg.mirror_strength))
    return nxt


def _run_chain(seed: int, height: float, width: int,
               n_steps: int = N_STEPS) -> dict:
    """Run one injected chain against its control (same seed and
    settings, no injection) and return the detector deviation series.
    """
    rng = np.random.default_rng(seed)
    lat = MembraneLattice(LatticeConfig(coupling=TUNNEL_COUPLING), rng)
    flat_state = lat.state.reshape(-1, DIM, DIM)
    flat_hamiltonians = lat.hamiltonians.reshape(-1, DIM, DIM)
    barriers = barrier_links(CHAIN_LENGTH, width)
    base = [flat_state[i].copy() for i in range(CHAIN_LENGTH)]
    runs = {}
    for tag, do_inject in (("injected", True), ("control", False)):
        states = [rho.copy() for rho in base]
        if do_inject:
            U = phase_unitary(INJECTION_PHI, DIM, "linear")
            states[SOURCE] = U @ states[SOURCE] @ U.conj().T
        series = [[s.copy() for s in states]]
        for _ in range(n_steps):
            states = _chain_step(lat, states, flat_hamiltonians,
                                 barriers, height)
            series.append([s.copy() for s in states])
        runs[tag] = series
    inj = runs["injected"]
    ctrl = runs["control"]
    d_inj = float(np.linalg.norm(inj[0][SOURCE] - ctrl[0][SOURCE]))
    deviations = np.array([
        float(np.linalg.norm(inj[t][DETECTOR] - ctrl[t][DETECTOR]))
        for t in range(n_steps + 1)
    ])
    return {
        "seed": seed,
        "height": height,
        "width": width,
        "n_steps": n_steps,
        "injected_deviation": d_inj,
        "detector_deviations": deviations,
        "transmission": float(np.max(deviations) / d_inj) if d_inj > 0 else 0.0,
        "t_peak": int(np.argmax(deviations)),
    }


def measure_transmission(seed: int, height: float, width: int,
                         n_steps: int = N_STEPS) -> dict:
    """Transmission of one injected run against its control (seed, height, width)."""
    return _run_chain(seed, height, width, n_steps)


def _mean_transmission(results: list[dict]) -> float:
    return float(np.mean([r["transmission"] for r in results]))


def log_linear_fit(x: np.ndarray, log_y: np.ndarray) -> dict:
    """Least-squares fit log_y ~ a + b*x; return slope, intercept, R^2."""
    x = np.asarray(x, dtype=float)
    log_y = np.asarray(log_y, dtype=float)
    design = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(design, log_y, rcond=None)
    pred = design @ coef
    ss_res = float(np.sum((log_y - pred) ** 2))
    ss_tot = float(np.sum((log_y - np.mean(log_y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"intercept": float(coef[0]), "slope": float(coef[1]), "r2": r2}


def run_tunneling_experiment(seeds: tuple[int, ...] = SEEDS,
                             heights: tuple[float, ...] = HEIGHTS,
                             widths: tuple[int, ...] = WIDTHS,
                             n_steps: int = N_STEPS) -> dict:
    """Run the height sweep, the width sweep and the wide-barrier check,
    and judge claims (h)-(k) by the pre-registered rules.
    """
    results: dict = {
        "provenance": lattice_provenance(),
        "seeds": list(seeds),
        "heights": list(heights),
        "widths": list(widths),
        "n_steps": n_steps,
        "chain_length": CHAIN_LENGTH,
        "coupling": TUNNEL_COUPLING,
        "verdict_rule": {
            "suppression_ratio_max": SUPPRESSION_RATIO_MAX,
            "width_log_slope_max": WIDTH_LOG_SLOPE_MAX,
            "width_log_r2_min": WIDTH_LOG_R2_MIN,
            "graded_ratio_min": GRADED_RATIO_MIN,
            "graded_ratio_max": GRADED_RATIO_MAX,
        },
    }
    # (h) height sweep at fixed width.
    height_rows = []
    for h in heights:
        runs = [measure_transmission(s, h, HEIGHT_SWEEP_WIDTH, n_steps)
                for s in seeds]
        height_rows.append({
            "height": h,
            "width": HEIGHT_SWEEP_WIDTH,
            "mean_transmission": _mean_transmission(runs),
            "per_seed_transmission": [r["transmission"] for r in runs],
            "per_seed_t_peak": [r["t_peak"] for r in runs],
        })
    results["height_sweep"] = height_rows
    # (i) width sweep at fixed height.
    width_rows = []
    for w in widths:
        runs = [measure_transmission(s, WIDTH_SWEEP_HEIGHT, w, n_steps)
                for s in seeds]
        width_rows.append({
            "width": w,
            "height": WIDTH_SWEEP_HEIGHT,
            "mean_transmission": _mean_transmission(runs),
            "per_seed_transmission": [r["transmission"] for r in runs],
            "per_seed_t_peak": [r["t_peak"] for r in runs],
        })
    results["width_sweep"] = width_rows
    # (k) wide-barrier suppression: height 0.99 at width 3. At height 0
    # the barrier is a no-op, so T(0, width 3) equals T(0, width 2)
    # from the height sweep (verified structurally in tests).
    wide_runs = [measure_transmission(s, 0.99, WIDE_WIDTH, n_steps)
                 for s in seeds]
    results["wide_barrier"] = {
        "width": WIDE_WIDTH,
        "height": 0.99,
        "mean_transmission": _mean_transmission(wide_runs),
        "per_seed_transmission": [r["transmission"] for r in wide_runs],
        "per_seed_t_peak": [r["t_peak"] for r in wide_runs],
    }

    # --- Verdicts (frozen thresholds, applied to the mean curves) ---
    h_mean = [row["mean_transmission"] for row in height_rows]
    h_ratio = h_mean[-1] / h_mean[0] if h_mean[0] > 0 else float("inf")
    h_monotone = all(
        h_mean[i + 1] < h_mean[i] for i in range(len(h_mean) - 1)
    )
    verdict_h = {
        "passed": bool(h_monotone and h_ratio <= SUPPRESSION_RATIO_MAX),
        "monotone_decreasing": bool(h_monotone),
        "suppression_ratio": float(h_ratio),
        "suppression_ratio_max": SUPPRESSION_RATIO_MAX,
        "reason": (
            "PASS" if (h_monotone and h_ratio <= SUPPRESSION_RATIO_MAX)
            else f"FAIL: monotone={h_monotone}, "
                 f"T(0.99)/T(0)={h_ratio:.4f} "
                 f"(max {SUPPRESSION_RATIO_MAX})"
        ),
    }
    w_mean = np.array([row["mean_transmission"] for row in width_rows])
    w_x = np.array(widths, dtype=float)
    fit = log_linear_fit(w_x, np.log(w_mean))
    verdict_i = {
        "passed": bool(
            fit["slope"] < WIDTH_LOG_SLOPE_MAX
            and fit["r2"] >= WIDTH_LOG_R2_MIN
        ),
        "log_linear_fit": fit,
        "slope_max": WIDTH_LOG_SLOPE_MAX,
        "r2_min": WIDTH_LOG_R2_MIN,
        "reason": (
            "PASS" if (
                fit["slope"] < WIDTH_LOG_SLOPE_MAX
                and fit["r2"] >= WIDTH_LOG_R2_MIN
            ) else f"FAIL: log-linear slope={fit['slope']:.4f} "
                    f"(max {WIDTH_LOG_SLOPE_MAX}), R^2={fit['r2']:.4f} "
                    f"(min {WIDTH_LOG_R2_MIN})"
        ),
    }
    graded_ratio = h_mean[2] / h_mean[0] if h_mean[0] > 0 else float("inf")
    verdict_j = {
        "passed": bool(GRADED_RATIO_MIN < graded_ratio < GRADED_RATIO_MAX),
        "graded_ratio": float(graded_ratio),
        "graded_ratio_range": (GRADED_RATIO_MIN, GRADED_RATIO_MAX),
        "reason": (
            "PASS" if GRADED_RATIO_MIN < graded_ratio < GRADED_RATIO_MAX
            else f"FAIL: T(0.5)/T(0)={graded_ratio:.4f} outside "
                 f"({GRADED_RATIO_MIN}, {GRADED_RATIO_MAX})"
        ),
    }
    k_ratio = (
        results["wide_barrier"]["mean_transmission"] / h_mean[0]
        if h_mean[0] > 0 else float("inf")
    )
    verdict_k = {
        "passed": bool(k_ratio <= SUPPRESSION_RATIO_MAX),
        "suppression_ratio": float(k_ratio),
        "suppression_ratio_max": SUPPRESSION_RATIO_MAX,
        "reason": (
            "PASS" if k_ratio <= SUPPRESSION_RATIO_MAX
            else f"FAIL: T(0.99, w=3)/T(0)={k_ratio:.4f} "
                 f"(max {SUPPRESSION_RATIO_MAX})"
        ),
    }
    results["verdicts"] = {"h": verdict_h, "i": verdict_i,
                           "j": verdict_j, "k": verdict_k}
    return results
