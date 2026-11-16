# Josephson–Membrane

**SYNTHETIC classical emulation only. No claim about physical reality.
No headline claim is emitted.**

The Josephson laws — DC current, AC frequency, SQUID flux modulation —
measured against the membrane-arena-3x3 lattice's response to a phase
difference. The laws are inputs, not discoveries. The central question
is whether the lattice follows the Josephson sine law. It does not.

## Results (v0.2.1, measured)

| Claim | Result |
|---|---|
| (a) DC law I_s = I_c sin(Δφ), cardinal table | **PASS** — exact (max error 2.4e-16) |
| (b) AC frequency f = 2eV/h | **PASS** — 2e/h = 483.597848 GHz/mV, CODATA K_J to 7e-16 relative |
| (c) SQUID I_c(Φ) = 2I_c\|cos(πΦ/Φ₀)\| | **PASS** — cardinals exact (max error 1.2e-16) |
| (d) Membrane phase response vs sine law | **FAIL** — R² = −32.75, max\|residual\|/\|A\| = 12.51, best-fit phase π/2 (a cosine) |
| (e) Odd-function variant vs sine law | **FAIL** on the primary (identical, linear): R² = 0.037, δ = 0 — but **PASS** for binary/distinct (R² = 0.9999, δ = −0.008); binary drive is an exact sine (R² = 1.0000) |
| (f) Phase-ramp analogue (AC) | measured, no verdict — v0.1.0 flow is a notch train; binary variant drive is a pure sine at the ramp rate (R² = 1.0000) |
| (g) Harmonics vs incoherence in the failing corner | **PASS** — the failing corner is INCOHERENT, not hiding a dominant harmonic: peak bin 7.3% of power, top-3 bins 16.8%, 84% of power above the 3rd harmonic; binary/distinct control is a single sine (p1 = 1.0000) |

4 of 6 verdicts pass (a–c, g). `headline_allowed` is always False.

## The FAIL in one paragraph

On identical cells differing only by the phase difference Δφ, the
membrane's pairwise signed flow (linear generator) is 0 at Δφ = 0
and negative everywhere else (−0.00167 at π/2, −0.00177 at π,
−0.00167 at 3π/2): no sign reversal, no positive lobe. The Josephson
law wants +I_c at π/2 and −I_c at 3π/2. The structural reason: for
identical cells the membrane flow depends on the fidelity, an overlap
magnitude, so it is **even** in Δφ (S(φ) = S(−φ)); the Josephson
current is **odd** (I(−φ) = −I(φ)). An even function can equal an odd
one only at zero. The best-fitting sinusoid is a cosine (δ = π/2),
and even it explains only 13.6% of the variance. (Scope: the
"negative everywhere" and "even" claims hold for identical cells
with the linear generator. The binary generator has a small positive
lobe at small phase, and for distinct cells the flow is not even —
see REPORT.md §3. Neither changes the FAIL.)

## The v0.2.0 follow-up in one paragraph

The recorded direction was to make the interaction depend on a
signed, odd function of phase. The literal quantity Im Tr(ρ_A
ρ_B(φ)) is identically zero (the trace of two Hermitian matrices is
real), so the variant uses the one-sided overlap Q(φ) = Im Tr(ρ_A
U(φ) ρ_B), odd for identical cells. On the primary setup (identical,
linear) it still FAILS (R² = 0.037) — the odd drive is a Fourier
series over 31 frequencies, not a single sine — but with the right
phase (δ = 0, the Josephson phase). The binary generator makes the
drive an exact sine (R² = 1.0000), and for distinct cells the variant
flow PASSES the pre-registered rule (R² = 0.9999, δ = −0.008): a
Josephson-like sine law in the membrane, in that regime, by design.
Under a linear phase ramp the v0.1.0 flow is a notch train, not an
oscillation; the binary variant drive is a pure sine at the ramp
rate. See REPORT.md §4–§5. The variant is a designed rule, not a
discovered law — no physical junction is implied.

## The v0.2.1 follow-up in one paragraph

Dr. Sinem Salva Ersoy asked the next question: in the failing
corner, is the response hiding higher harmonics, or are the phases
just incoherent? A single sine fit cannot tell those apart, so the
test is a Fourier decomposition of the drive and flow on a dense
256-point grid, with a shuffled-phase control setting the noise
floor. The answer is incoherent: the failing corner's flow has no
dominant harmonic — peak bin 7.3% of power, top-3 bins 16.8%, 84%
of the power spread above the 3rd harmonic. The drive's strongest
bins are high (24th, 30th, 29th), but each carries only ~5% of the
power. The binary/distinct control confirms the test works: there
the flow is a single sine (p1 = 1.0000). See REPORT.md §7. A Fourier
peak is a property of the sampled curve, not of a device.

## Run it

```bash
python3 -m pytest -q            # 66 tests
PYTHONPATH=src python3 -m josephson_membrane.unified   # full claim, JSON
```

Requires Python ≥ 3.10 and NumPy. The membrane experiment uses the
REAL membrane-arena-3x3 lattice at
`~/workspace/membrane-arena-3x3-upgrade` (branch `rounds-on-main`,
head `e7d8fea70cbc0429ad9acb1cea8d55e989cfca70`); it is imported, never
reimplemented. If that tree is absent the module raises ImportError.

## Files

- `src/josephson_membrane/josephson.py` — DC law, AC frequency, SQUID.
- `src/josephson_membrane/membrane_phase.py` — Δφ sweep on the real
  lattice, sine fit, pre-registered verdict rule.
- `src/josephson_membrane/odd_phase.py` — odd-function variant
  (claim e) and phase ramp (claim f).
- `src/josephson_membrane/harmonics.py` — Fourier harmonics test
  (claim g) with shuffled-phase control.
- `src/josephson_membrane/unified.py` — the seven statements.
- `REPORT.md` — full measured numbers and honest limits.

## Credit

The Josephson laws are Brian Josephson's (Nobel Prize in Physics,
1973). Brisson is the package author. No Lee connection.

## Licence

Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP). AGPL-3.0-only. See LICENSE.
