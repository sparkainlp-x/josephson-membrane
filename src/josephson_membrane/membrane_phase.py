"""Membrane phase response: does the 3x3 lattice follow the Josephson sine law?

PRE-REGISTERED VERDICT RULE (fixed before the measurement run below):
the membrane's phase response follows the Josephson sine law iff, on
the least-squares fit y = A * sin(delta_phi + delta) with NO DC term,
    R^2 >= 0.95  AND  max |residual| <= 0.05 * |A|.
The rule is applied to the pairwise signed flow S(delta_phi) on
identical cells, averaged over seeds 0..4, for each phase generator.
A fit with a DC term (y = C + A*sin(delta_phi+delta)) is reported as a
diagnostic only: the Josephson law has no DC term, so a good DC fit does
not rescue a failed Josephson fit. The fitted phase delta is also
reported: the Josephson law is specifically sin(delta_phi), i.e.
delta = 0; a fit with delta ~ pi/2 is a cosine, not the Josephson law.

How a phase difference is imposed (the honest way). The Josephson
delta_phi is the relative phase of two wavefunctions. The lattice's own
phase_offsets are clock advances in the unitary stage (a time
difference, already measured as an additive rate shift in round 2), not
a quantum phase. So delta_phi is imposed directly on the state: cell B's
density matrix is conjugated by the diagonal phase unitary
U(delta_phi), which rotates the relative phases of B's coherences and
leaves its populations exactly unchanged. Two generators are measured:
  "linear": U = diag(exp(i*k*delta_phi)), k = 0..31 (phase grows with
            level index -- the standard multi-level phase twist);
  "binary": U = diag(exp(i*delta_phi)) on even levels, 1 on odd levels
            (a clean two-group phase difference).
Both use the REAL lattice: fidelity() and the interaction formula come
from membrane_arena_3x3.core / .lattice, and the full-lattice check
calls MembraneLattice.interaction_delta() and propose_step() directly.
The lattice is never reimplemented here.

Observables, for delta_phi sweeping [0, 2*pi):
  S(delta_phi): pairwise signed flow between neighbour cells A=(0,0)
    and B=(0,1): coupling * (F(rho_A, rho_B(phi)) - 0.90)
    * ||rho_B(phi) - rho_A||_F. Positive = A flows toward B.
  M(delta_phi): |S(delta_phi)|, the interaction-term magnitude.
  D(delta_phi): ||interaction_delta(A)||_F with B rotated in place --
    the full-lattice interaction magnitude for A (real method).
  C(delta_phi): ||state_A after one propose_step - state_A before||_F
    -- the one-step state change (measured on a coarser grid; the step
    is dominated by the unitary/dissipative/mirror stages).
  F(delta_phi): fidelity(rho_A, rho_B(phi)) -- the driver. Note the
    arena's fidelity() returns the SQUARED fidelity F^2 (its
    docstring: "squared root fidelity"), so the 0.90 resonance
    threshold below is a threshold on F^2 (F >= 0.949), and all
    fidelity numbers in this package are F^2.
Two setups: "identical" (rho_B = rho_A, so the cells differ ONLY by
delta_phi -- the clean Josephson analogue) and "distinct" (rho_B is a
different random cell, the generic case).

Structural note, recorded before measuring: for identical cells
(rho_A = rho_B) the membrane flow depends on the fidelity, which is
a function of the overlap magnitude and hence EVEN in delta_phi
(S(phi) = S(-phi)); the Josephson current is ODD
(I(-phi) = -I(phi)). An even function can equal an odd one only at
zero. For distinct cells the fidelity is NOT even in general (the
overlap Tr(rho_A U rho_B) is not symmetric under phi -> -phi when
rho_A != rho_B), so the evenness argument is scoped to identical
cells. The measurement below tests whether that structural mismatch
shows up in the numbers, and by how much.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np

# --- Pre-registered verdict thresholds (see module docstring) ---
SINE_R2_MIN = 0.95
SINE_MAX_RESIDUAL_FRACTION = 0.05

# --- Lattice provenance: use the REAL code, never a reimplementation ---
_LATTICE_CANDIDATES = (
    Path.home() / "workspace" / "membrane-arena-3x3-upgrade",
    Path.home() / "workspace" / "membrane-arena-3x3",
)

_lattice_root: Path | None = None
for _candidate in _LATTICE_CANDIDATES:
    if (_candidate / "membrane_arena_3x3" / "lattice.py").is_file():
        _lattice_root = _candidate
        break
if _lattice_root is None:
    raise ImportError(
        "membrane-arena-3x3 lattice tree not found in "
        + ", ".join(str(c) for c in _LATTICE_CANDIDATES)
    )
if str(_lattice_root) not in sys.path:
    sys.path.insert(0, str(_lattice_root))

from membrane_arena_3x3.core import DIM, fidelity  # noqa: E402
from membrane_arena_3x3.lattice import LatticeConfig, MembraneLattice  # noqa: E402

SEEDS = (0, 1, 2, 3, 4)
N_PHI = 64

# The interaction constants, taken from the lattice defaults
# (LatticeConfig) so the variant module shares the same values.
COUPLING = LatticeConfig().coupling  # 0.02
RESONANCE_THRESHOLD = LatticeConfig().resonance_threshold  # 0.90


def lattice_provenance() -> dict:
    """Return the lattice tree path, branch and HEAD commit used."""
    info: dict = {"path": str(_lattice_root), "branch": None, "head": None}
    try:
        out = subprocess.run(
            ["git", "-C", str(_lattice_root), "rev-parse", "--abbrev-ref", "HEAD"],
            capture_output=True, text=True, timeout=10, check=True,
        )
        info["branch"] = out.stdout.strip()
        out = subprocess.run(
            ["git", "-C", str(_lattice_root), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=10, check=True,
        )
        info["head"] = out.stdout.strip()
    except Exception:
        pass
    return info


def phase_unitary(phi: float, dim: int = DIM, generator: str = "linear") -> np.ndarray:
    """Diagonal phase unitary U(phi) for imposing a phase difference.

    "linear": U = diag(exp(i*k*phi)), k = 0..dim-1.
    "binary": U = diag(exp(i*phi)) on even levels, 1 on odd levels.
    Both are unitary, preserve populations (the diagonal of rho), and
    rotate only the relative phases of the coherences.
    """
    if generator == "linear":
        idx = np.arange(dim)
        return np.diag(np.exp(1j * idx * phi))
    if generator == "binary":
        diag = np.ones(dim, dtype=complex)
        diag[0::2] = np.exp(1j * phi)
        return np.diag(diag)
    raise ValueError(f"unknown generator: {generator!r}")


def pairwise_flow(rho_a: np.ndarray, rho_b_phi: np.ndarray,
                  coupling: float, threshold: float) -> dict:
    """Pairwise interaction between A and phase-rotated B, real formula.

    drive = coupling * (fidelity(rho_A, rho_B(phi)) - threshold), exactly
    the lattice's interaction rule for one neighbour pair; the signed
    flow S = drive * ||rho_B(phi) - rho_A||_F is the coefficient of
    (rho_B - rho_A) times its length -- the current between the cells.
    """
    f = fidelity(rho_a, rho_b_phi)
    drive = coupling * (f - threshold)
    dist = float(np.linalg.norm(rho_b_phi - rho_a))
    signed = drive * dist
    return {
        "fidelity": float(f),
        "drive": float(drive),
        "distance": dist,
        "signed_flow": float(signed),
        "magnitude": float(abs(signed)),
    }


def sine_fit(phis: np.ndarray, y: np.ndarray, with_dc: bool = False) -> dict:
    """Least-squares fit y ~ A*sin(phi + delta) [+ C if with_dc].

    Returns amplitude A (>= 0), phase delta, R^2 (against the data mean),
    max |residual|, and max|residual|/|A|. R^2 can be negative: the fit
    has no intercept unless with_dc, so a constant signal scores far
    below a constant mean -- which is the point of the Josephson test.
    """
    y = np.asarray(y, dtype=float)
    cols = [np.sin(phis), np.cos(phis)]
    if with_dc:
        cols.append(np.ones_like(phis))
    design = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(design, y, rcond=None)
    a, b = float(coef[0]), float(coef[1])
    amp = float(np.hypot(a, b))
    delta = float(np.arctan2(b, a))
    pred = design @ coef
    resid = y - pred
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    max_res = float(np.max(np.abs(resid)))
    return {
        "amplitude": amp,
        "delta": delta,
        "r2": r2,
        "max_residual": max_res,
        "max_residual_fraction": max_res / amp if amp > 0 else float("inf"),
        "dc": float(coef[2]) if with_dc else 0.0,
    }


def sine_law_verdict(fit: dict) -> dict:
    """Apply the pre-registered rule to a no-DC sine fit."""
    passed = bool(
        fit["r2"] >= SINE_R2_MIN
        and fit["max_residual_fraction"] <= SINE_MAX_RESIDUAL_FRACTION
    )
    return {
        "passed": passed,
        "r2_min": SINE_R2_MIN,
        "max_residual_fraction_max": SINE_MAX_RESIDUAL_FRACTION,
        "reason": (
            "PASS" if passed else
            f"FAIL: R^2={fit['r2']:.4f} (< {SINE_R2_MIN}) and/or "
            f"max|residual|/|A|={fit['max_residual_fraction']:.3f} "
            f"(> {SINE_MAX_RESIDUAL_FRACTION})"
        ),
    }


def sweep(seed: int, generator: str = "linear", n_phi: int = N_PHI,
          identical: bool = True, step_change_grid: int = 0) -> dict:
    """Sweep delta_phi over [0, 2*pi) on the real lattice for one seed.

    step_change_grid: if > 0, also measure the one-step state change C
    on that many evenly spaced grid points (propose_step is the slow
    stage, so C runs on a coarser grid than S/M/D/F).
    """
    rng = np.random.default_rng(seed)
    lat = MembraneLattice(LatticeConfig(), rng)
    rho_a = lat.state[0, 0]
    rho_b = rho_a if identical else lat.state[0, 1]
    phis = np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False)
    out = {
        "seed": seed, "generator": generator, "identical": identical,
        "phis": phis,
        "fidelity": np.empty(n_phi),
        "signed_flow": np.empty(n_phi),
        "magnitude": np.empty(n_phi),
        "lattice_delta": np.empty(n_phi),
        "step_change_phis": np.empty(0),
        "step_change": np.empty(0),
    }
    cfg = lat.config
    step_idx = (
        set(np.linspace(0, n_phi - 1, step_change_grid, dtype=int).tolist())
        if step_change_grid > 0 else set()
    )
    step_phis: list[float] = []
    step_vals: list[float] = []
    for n, phi in enumerate(phis):
        U = phase_unitary(float(phi), DIM, generator)
        rho_b_phi = U @ rho_b @ U.conj().T
        flow = pairwise_flow(rho_a, rho_b_phi, cfg.coupling, cfg.resonance_threshold)
        out["fidelity"][n] = flow["fidelity"]
        out["signed_flow"][n] = flow["signed_flow"]
        out["magnitude"][n] = flow["magnitude"]
        saved = lat.state[0, 1].copy()
        lat.state[0, 1] = rho_b_phi
        out["lattice_delta"][n] = float(np.linalg.norm(lat.interaction_delta(0, 0)))
        if n in step_idx:
            before = lat.state[0, 0].copy()
            cand = lat.propose_step()
            step_phis.append(float(phi))
            step_vals.append(float(np.linalg.norm(cand[0, 0] - before)))
        lat.state[0, 1] = saved
    out["step_change_phis"] = np.asarray(step_phis)
    out["step_change"] = np.asarray(step_vals)
    return out


def run_membrane_experiment(seeds: tuple[int, ...] = SEEDS,
                            n_phi: int = N_PHI,
                            step_change_grid: int = 16) -> dict:
    """Run the full phase-response experiment and judge the sine law.

    Primary: identical cells, linear generator, mean over seeds.
    Supporting: binary generator, distinct cells, full-lattice delta,
    one-step state change. The verdict uses the pre-registered rule on
    the primary mean curve (and is reported per seed too).
    """
    phis = np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False)
    results: dict = {"phis": phis, "provenance": lattice_provenance(),
                     "seeds": list(seeds), "n_phi": n_phi,
                     "verdict_rule": {
                         "r2_min": SINE_R2_MIN,
                         "max_residual_fraction_max": SINE_MAX_RESIDUAL_FRACTION,
                     }}
    for generator in ("linear", "binary"):
        for identical in (True, False):
            key = f"{generator}_{'identical' if identical else 'distinct'}"
            curves = {"signed_flow": [], "magnitude": [],
                      "lattice_delta": [], "fidelity": []}
            per_seed_fit = []
            step_all: list[float] = []
            for seed in seeds:
                sw = sweep(seed, generator, n_phi, identical,
                           step_change_grid if identical else 0)
                for name in curves:
                    curves[name].append(sw[name])
                per_seed_fit.append(sine_fit(phis, sw["signed_flow"]))
                step_all.extend(sw["step_change"].tolist())
            mean_flow = np.mean(curves["signed_flow"], axis=0)
            fit = sine_fit(phis, mean_flow)
            fit_dc = sine_fit(phis, mean_flow, with_dc=True)
            verdict = sine_law_verdict(fit)
            results[key] = {
                "generator": generator,
                "identical": identical,
                "mean_signed_flow": mean_flow,
                "mean_magnitude": np.mean(curves["magnitude"], axis=0),
                "mean_lattice_delta": np.mean(curves["lattice_delta"], axis=0),
                "mean_fidelity": np.mean(curves["fidelity"], axis=0),
                "fit": fit,
                "fit_with_dc": fit_dc,
                "verdict": verdict,
                "per_seed_r2": [f["r2"] for f in per_seed_fit],
                "per_seed_delta": [f["delta"] for f in per_seed_fit],
                "per_seed_max_residual_fraction": [
                    f["max_residual_fraction"] for f in per_seed_fit
                ],
                "cardinal_signed_flow": {
                    "0": float(mean_flow[0]),
                    "pi/2": float(mean_flow[n_phi // 4]),
                    "pi": float(mean_flow[n_phi // 2]),
                    "3pi/2": float(mean_flow[3 * n_phi // 4]),
                },
                "step_change_mean": float(np.mean(step_all)) if step_all else None,
                "step_change_range": (
                    [float(np.min(step_all)), float(np.max(step_all))]
                    if step_all else None
                ),
            }
    return results
