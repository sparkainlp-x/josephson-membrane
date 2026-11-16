# Changelog

## 0.2.0 — 2026-11-15

- Independent review of v0.1.0 (see REPORT.md §3): all four claims
  reproduce exactly; claim (d)'s FAIL stands. Two overgeneralisations
  corrected and regression-tested: the evenness argument is scoped to
  identical cells only (for distinct cells S(pi/2) - S(3pi/2) reaches
  3.5e-04; analytic counterexample F^2(phi) = (1 - sin phi)/2), and
  the "drive negative for every delta_phi != 0" claim is scoped to
  the linear generator only (the binary generator has a positive lobe:
  mean S(2pi/64) = +2.45e-05, 20 of 63 grid points positive; binary
  fidelity stays above the 0.90 threshold at small phase). Clarified
  that the arena's fidelity() returns the squared fidelity F^2.
- `odd_phase.py` (new): odd-function variant, claim (e). The literal
  Im Tr(rho_A rho_B(phi)) is identically zero (trace of two Hermitian
  matrices is real); the variant uses the one-sided overlap
  Q(phi) = Im Tr(rho_A U(phi) rho_B), odd for identical cells. Same
  pre-registered verdict rule as v0.1.0. Primary (identical, linear):
  FAIL (flow R^2 = 0.0366, max|residual|/|A| = 17.86, delta = 0 --
  the odd drive is a 31-frequency Fourier series, not a single sine).
  Binary drive is an exact sine (R^2 = 1.0000, delta = 0).
  Binary/distinct flow PASSES the rule (R^2 = 0.9999,
  max|residual|/|A| = 0.009, delta = -0.0079). Binary/identical flow
  is a near-miss (R^2 = 0.9222), distorted by the distance notch.
- `odd_phase.py`: phase-ramp analogue, claim (f), no forced verdict.
  Under a linear phase ramp the v0.1.0 flow is a notch train
  (mean -0.0016, dips to 0 at each cycle boundary; 13.5% of power at
  the ramp frequency; sine fit R^2 = -32.75, the same number as (d)).
  The binary variant drive is a pure sine at the ramp rate
  (100% of power, R^2 = 1.0000, delta = 0); the binary/distinct
  variant flow is a clean sine at the ramp rate (R^2 = 0.9999).
- `unified.py`: six statements; (f) carries no verdict ("passed" is
  null). Overall 3 of 5 verdicts pass. headline_allowed stays False.
- The honest-limits section gains the variant caveat: the variant
  interaction is a designed rule, not a discovered law.

## 0.1.0 — 2026-11-15

- Initial release. SYNTHETIC classical emulation only.
- `josephson.py`: DC law I_s = I_c sin(delta_phi) with the cardinal
  table; AC frequency f = 2eV/h from the exact SI constants e and h
  (2e/h = 483.597848 GHz/mV, CODATA K_J to machine precision); SQUID
  flux modulation I_c(Phi) = 2 I_c |cos(pi Phi / Phi_0)| with
  Phi_0 = h/2e = 2.067833848e-15 Wb.
- `membrane_phase.py`: phase difference imposed on the REAL
  membrane-arena-3x3 lattice (branch rounds-on-main) by conjugating a
  cell's state with a diagonal phase unitary; pairwise signed flow,
  interaction magnitude, full-lattice delta and one-step state change
  swept over delta_phi in [0, 2*pi), seeds 0..4, two generators
  (linear, binary), identical and distinct cells.
- `unified.py`: four claims. (a) DC cardinal table exact. (b) AC
  constant and sample values exact. (c) SQUID cardinals exact.
  (d) Membrane phase response vs the sine law: FAIL by the
  pre-registered rule (R^2 = -32.75, max|residual|/|A| = 12.51,
  best-fit phase pi/2 -- a cosine, not the Josephson sine).
- No headline claim is emitted. headline_allowed is always False.
