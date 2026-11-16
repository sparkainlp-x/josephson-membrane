# REPORT — Josephson–Membrane v0.2.0

**SYNTHETIC classical emulation only. No claim about physical reality.
No headline claim is emitted (`headline_allowed` is always False).**

## 1. Verdicts

| Claim | Verdict | Key numbers |
|---|---|---|
| (a) DC law cardinal table | PASS | max error 2.4e-16 |
| (b) AC frequency f = 2eV/h | PASS | 2e/h = 483.597848 GHz/mV; CODATA K_J reproduced to 7.2e-16 relative |
| (c) SQUID flux modulation | PASS | max error 1.2e-16 |
| (d) Membrane phase response vs sine law | **FAIL** | R² = −32.75; max\|residual\|/\|A\| = 12.51; best-fit phase δ = π/2 |
| (e) Odd-function variant vs sine law | **FAIL** (primary) | R² = 0.037; max\|residual\|/\|A\| = 17.86; δ = 0 — but binary/distinct PASSES (R² = 0.9999, δ = −0.008) |
| (f) Phase-ramp analogue (AC) | measured, no verdict | v0.1.0 flow is a notch train at the ramp rate; binary variant drive is a pure sine at the ramp rate (R² = 1.0000) |

Overall: 3 of 5 verdicts pass (a–c). The failures are the point of
the package: the lattice does not follow the Josephson sine law under
its own interaction rule, and the measurement says exactly how it
differs — and exactly which designed change does produce a sine.

## 2. The laws (claims a–c): exact reproductions

These are established physics, taken as inputs. The package checks
that the implementation reproduces them; it does not discover them.

### (a) DC law: I_s = I_c sin(Δφ)

Cardinal table (I_c = 1), the infographic's table:

| Δφ | I_s (measured) | I_s (law) |
|---|---|---|
| 0 | 0.0 | 0 |
| π/2 | 1.0 | +I_c |
| π | 1.22e-16 | 0 |
| 3π/2 | −1.0 | −I_c |
| 2π | −2.45e-16 | 0 |

Max error 2.4e-16 (floating-point sin near π and 2π). Oddness
error 0.0, 2π-periodicity error 3.9e-16. The law is odd and periodic,
exactly as stated.

### (b) AC frequency: f = 2eV/h

Constants are the exact SI values fixed in 2019 (e and h are exact by
definition), so 2e/h is computed, not copied:

- 2e/h = **483.59784841698365 GHz/mV** (483597.848416984 GHz/V).
- CODATA 2018 K_J = 483597.848416984 GHz/V. Relative error 7.2e-16
  (machine precision).
- The rounded figure 483.5979 GHz/mV agrees to 6 significant figures;
  the exact computed value is 483.597848.

Sample voltages (f = K_J × V):

| V | f (measured) |
|---|---|
| 1 µV | 483,597,848.417 Hz (483.6 MHz) |
| 10 µV | 4,835,978,484.170 Hz (4.836 GHz) |
| 1 mV | 483,597,848,416.984 Hz (483.598 GHz) |
| 1 V | 483,597,848,416,983.6 Hz (483.6 THz) |

Max sample relative error < 1e-12. Linearity in V exact.

### (c) SQUID: I_c(Φ) = 2 I_c |cos(πΦ/Φ₀)|

Flux quantum Φ₀ = h/2e = **2.0678338484619295e-15 Wb** (CODATA
2.067833848e-15 Wb; relative error 2.2e-10, limited by the rounded
CODATA figure quoted).

| Φ | I_c(Φ) (measured) | I_c(Φ) (law) |
|---|---|---|
| 0 | 2.0 | 2 I_c |
| Φ₀/4 | 1.4142135623730951 | √2 I_c |
| Φ₀/2 | 1.22e-16 | 0 |
| Φ₀ | 2.0 | 2 I_c |

Max error 1.2e-16. Non-negative everywhere, Φ₀-periodic, maximum
2 I_c, minimum 0.

## 3. The membrane experiment (claim d): the sine law FAILS

### Method

The Josephson Δφ is the relative phase of two wavefunctions. The
lattice's own `phase_offsets` are clock advances in the unitary stage
(a time difference, already measured as an additive rate shift in
round 2), not a quantum phase — so Δφ is imposed directly on the
state: cell B's density matrix is conjugated by the diagonal phase
unitary U(Δφ), which rotates the relative phases of B's coherences
and leaves its populations exactly unchanged. Two generators:

- **linear**: U = diag(exp(i·k·Δφ)), k = 0..31.
- **binary**: U = diag(exp(i·Δφ)) on even levels, 1 on odd levels.

Both use the REAL lattice: `fidelity()` and the interaction formula
from `membrane_arena_3x3.core` / `.lattice`; the full-lattice check
calls `MembraneLattice.interaction_delta()` and `propose_step()`
directly. Lattice tree: `~/workspace/membrane-arena-3x3-upgrade`,
branch `rounds-on-main`, head `e7d8fea70cbc0429ad9acb1cea8d55e989cfca70`.

Observables, Δφ sweeping [0, 2π) in 64 points, seeds 0..4:

- **S(Δφ)** — pairwise signed flow between neighbours A=(0,0) and
  B=(0,1): coupling × (F(ρ_A, ρ_B(φ)) − 0.90) × ‖ρ_B(φ) − ρ_A‖_F.
  Positive = A flows toward B.
- **M(Δφ)** = |S(Δφ)|, the interaction-term magnitude.
- **D(Δφ)** — ‖interaction_delta(A)‖_F with B rotated in place.
- **C(Δφ)** — one-step state change ‖state_A after propose_step −
  state_A before‖_F (16-point grid; the step is slow).
- **F(Δφ)** — the fidelity driving the interaction.

Two setups: **identical** (ρ_B = ρ_A, so the cells differ ONLY by
Δφ — the clean Josephson analogue) and **distinct** (ρ_B a different
random cell, the generic case).

### Pre-registered verdict rule

Fixed in `membrane_phase.py` before the measurement run: the sine law
holds iff, on the least-squares fit y = A·sin(Δφ + δ) with NO DC term,
**R² ≥ 0.95 AND max|residual| ≤ 0.05·|A|**. The fit is applied to the
mean S(Δφ) over seeds 0..4. A fit with a DC term is reported as a
diagnostic only (the Josephson law has no DC term). The fitted phase δ
is reported: the Josephson law is specifically sin(Δφ), i.e. δ = 0.

### Measured result — primary (identical cells, linear generator)

| Quantity | Measured |
|---|---|
| Fit A | 1.442150e-04 |
| Fit δ | 1.5708 (= π/2) |
| Fit R² | **−32.7501** (bar: ≥ 0.95) |
| max\|residual\| | 1.803523e-03 |
| max\|residual\|/\|A\| | **12.506** (bar: ≤ 0.05) |
| Fit with DC: R² / max\|residual\|/\|A\| / C | 0.1355 / 10.017 / −1.588872e-03 |
| Per-seed R² | −31.88, −30.50, −28.52, −31.56, −31.35 |
| Per-seed δ | 1.571 in all 5 seeds |
| Per-seed max\|residual\|/\|A\| | 12.56, 15.67, 14.07, 12.17, 10.86 |

Cardinal flows (mean over seeds):

| Δφ | S(Δφ) | Josephson I_s/I_c |
|---|---|---|
| 0 | 0.0 | 0 |
| π/2 | −0.00167416 | +1 |
| π | −0.00176917 | 0 |
| 3π/2 | −0.00167416 | −1 |

**Verdict: FAIL.** The flow is 0 at Δφ = 0 and negative everywhere
else. There is no positive lobe and no sign reversal. The best-fitting
sinusoid is a cosine (δ = π/2), and even it explains only 13.6% of the
variance with a DC term (R² = 0.1355) — the actual shape is a sharp
notch at Δφ = 0 (S = 0 there because the distance ‖ρ_B − ρ_A‖ is 0)
dropping immediately to ≈ −0.0017 and staying there, not a sine.

### Supporting setups — all FAIL

| Setup | Fit R² | max\|residual\|/\|A\| | δ | Fit+DC R² | Fit+DC max\|residual\|/\|A\| |
|---|---|---|---|---|---|
| binary, identical | 0.1158 | 1.000 | 1.571 | 0.9392 | 0.338 |
| linear, distinct | −2290.99 | 218.14 | −0.338 | 0.0263 | 9.539 |
| binary, distinct | −5025.61 | 50.26 | 0.842 | 0.9971 | 0.061 |

Notes:

- **binary, identical**: the best no-DC fit has R² = 0.12 and a
  residual exactly equal to the amplitude (the residual at Δφ = 0 is
  A, because S(0) = 0 while the fitted cosine peaks there). The DC
  fit reaches R² = 0.94 but its residual is still 34% of the
  amplitude — and δ = π/2, a cosine, not the Josephson sine.
- **linear, distinct**: the flow is nearly constant (−0.00164 ± 4%)
  with irregular wiggles; the sine fit is worse than a constant mean
  by three orders of magnitude (R² = −2291).
- **binary, distinct**: the DC fit explains 99.7% of the variance of
  the tiny wiggles — but those wiggles are ~4% of a constant −0.0016
  offset, the residual is still 6.1% of the (tiny) amplitude, and
  δ = 0.842 ≠ 0. A 4% wiggle on a constant is not a sine law by any
  reading; the no-DC fit fails catastrophically (R² = −5026).

### The structural reason

For **identical cells** (ρ_A = ρ_B, the primary setup) the membrane
flow is **even** in Δφ: S(φ) = S(2π − φ) exactly (measured:
S(π/2) = S(3π/2) to machine precision in every seed and both
generators; the test suite asserts this). It has to be: the flow
depends on the fidelity, an overlap magnitude, and on the distance
‖ρ_B(φ) − ρ_A‖ — both even functions of the phase when ρ_A = ρ_B.
The Josephson current is **odd**: I(−φ) = −I(φ). An even function
can equal an odd one only at zero. The mismatch is structural, not
parametric: no coupling, threshold, or generator choice can fix it,
because the interaction rule itself — drive = coupling ×
(fidelity − 0.90) — takes the fidelity as its only phase-dependent
input, and the fidelity is an even function of Δφ for identical cells.

Two scoping corrections (independent review, v0.2.0). First, the
evenness argument holds for identical cells **only**. For **distinct**
cells the fidelity is not even in general: S(π/2) − S(3π/2) reaches
3.5×10⁻⁴ (linear generator, seed 1) — a genuine counterexample, not
noise (an analytic 2-level counterexample is in the test suite). The
distinct setups fail the sine law on their own numbers (R² = −2291
and −5026), but not by the evenness argument. Second, the "drive
negative for every Δφ ≠ 0" claim holds for the **linear** generator
only. The binary generator has a **positive lobe**: for identical
cells the mean flow is positive at small phase (mean S(2π/64) =
+2.45×10⁻⁵ across seeds; 20 of 63 grid points positive), because
the binary fidelity stays above the 0.90 threshold at small phase
(mean F = 0.9990 at 2π/64, 0.96 at π/4) and only falls below it at
larger phase (mean F = 0.79 at π/2). The binary flow still fails the
sine law (R² = 0.12) — a small positive lobe is not a sine — but the
threshold mechanism is not "never crosses zero under phase
rotation": under the binary generator the fidelity does cross the
threshold.

A second, mechanistic reason (linear generator): the fidelity stays
in [0.55, 1.0] and is below the 0.90 resonance threshold for every
Δφ ≠ 0 (mean fidelity at π/2 ≈ 0.57, at π ≈ 0.55 for the linear
generator), so the drive is negative for every Δφ ≠ 0. Under the
linear generator the flow can never change sign under a phase
rotation, because a phase rotation can never push the fidelity
across the threshold. The Josephson current changes sign because
sin changes sign; the membrane flow cannot change sign because its
driver cannot cross zero. (All fidelity numbers here are the arena's
squared fidelity F² — the arena's `fidelity()` returns (Σ√λ)² — so
the 0.90 threshold is a threshold on F², i.e. F ≥ 0.949.)

### The other observables

- **Full-lattice delta D(Δφ)**: mean range [0.001039, 0.001360] —
  varies by ~±15% around its mean, mostly a dip at Δφ = 0 where the
  rotated neighbour coincides with A and contributes nothing. No sine.
- **One-step state change C(Δφ)**: mean 0.126144 across seeds
  (0.121–0.133); within any one sweep it varies by < 1% with Δφ.
  The step is dominated by the unitary / dissipative / mirror stages,
  not by the interaction — the phase difference is invisible at the
  one-step level to ~1e-6.

## 4. The odd-function variant (claim e): FAIL on the primary, PASS in one supporting setup

### The correction that had to come first

The brief said "a variant interaction whose drive depends on the
imaginary part of the overlap, Im Tr(ρ_A ρ_B(φ))". The literal
quantity — with ρ_B(φ) = U ρ_B U†, both states rotated into the same
frame — is **identically zero**: the trace of a product of two
Hermitian matrices is real. Measured max |Im| ≈ 1e-18 (machine zero)
on every setup. A drive built from it would be exactly zero
everywhere. The variant therefore uses the **one-sided** overlap

    Q(φ) = Im Tr(ρ_A U(φ) ρ_B),

with ρ_B unrotated. Q is odd in φ for identical cells (ρ_A = ρ_B):
then Q(φ) = Σ_ij |ρ_ij|² sin(jφ) for the linear generator —
manifestly odd, zero at φ = 0, changing sign (measured Q(π/2) =
+4.34×10⁻³, Q(3π/2) = −4.34×10⁻³). For distinct cells Q is **not**
odd in general (measured deviation ~1.3×10⁻³; the oddness proof
fails because Tr(ρ_A U ρ_B) is not symmetric under φ → −φ when
ρ_A ≠ ρ_B). The variant flow keeps the v0.1.0 structure with only
the drive changed:

    S_odd(φ) = drive(φ) · ‖ρ_B(φ) − ρ_A‖_F,
    drive(φ) = 0.02 · Q(φ).

Same lattice, same seeds 0–4, same 64-point grid, same coupling.

### Pre-registered verdict rule

Identical to v0.1.0, recorded in `odd_phase.py` before measuring:
fit A·sin(Δφ + δ) with no DC term; the sine law holds iff R² ≥ 0.95
AND max|residual| ≤ 0.05·|A|; report δ (Josephson is δ = 0).
Primary setup: identical cells, linear generator — mirroring (d).

### Measured result — primary (identical cells, linear generator)

| Quantity | Value |
|---|---|
| Flow fit R² | 0.0366 |
| Flow max\|residual\|/\|A\| | 17.86 |
| Flow δ | 0.0000 |
| Drive fit R² | 0.0266 |
| Drive max\|residual\|/\|A\| | 22.10 |
| Drive δ | 0.0000 |

**Verdict: FAIL.** The odd drive exists and has the right phase
(δ = 0, the Josephson phase — unlike (d)'s cosine), but under the
linear generator it is a Fourier series over 31 frequencies
(Q(φ) = Σ_ij |ρ_ij|² sin(jφ)), not a single sine: the drive fit
explains only 2.7% of the variance. Making the drive odd was
necessary but not sufficient — a single-frequency structure is
needed too.

### Supporting setups

| Setup | Flow R² | Flow max\|residual\|/\|A\| | Flow δ | Drive R² | Drive δ | Verdict |
|---|---|---|---|---|---|---|
| binary, identical | 0.9222 | 0.320 | 0.0000 | 1.0000 | −0.0000 | FAIL (near-miss) |
| linear, distinct | 0.0237 | 23.29 | 0.0294 | 0.0238 | 0.0322 | FAIL |
| binary, distinct | 0.9999 | 0.009 | −0.0079 | 0.9999 | −0.0079 | **PASS** |

Three findings:

- **The binary drive is an exact sine.** For the binary generator
  only the even levels carry the phase, so Q(φ) = W·sin(φ) with a
  single weight W = Im Tr(ρ_A P_even ρ_B): drive fit R² = 1.0000,
  δ = 0, for identical cells; R² = 0.9999, δ = −0.008, for distinct
  cells (the small deviation is the constant and cosine terms from
  the odd levels, which do not vanish exactly for distinct cells).
- **The binary/distinct flow PASSES the pre-registered rule**
  (R² = 0.9999, max|residual|/|A| = 0.009, δ = −0.008). That is a
  real result: in the distinct-cells, binary-generator regime the
  variant flow follows the Josephson sine law to within 1%. The
  reason is mechanistic — for distinct cells the distance
  ‖ρ_B(φ) − ρ_A‖ is nearly constant in φ (the two states are
  already far apart), so S_odd ≈ (pure sine) × (constant) stays a
  sine.
- **The binary/identical flow is a near-miss** (R² = 0.92,
  residual 32% of amplitude): the drive is an exact sine, but the
  distance has the notch at φ = 0 (distance 0 there, rising to a
  plateau), which distorts the product. The distance factor, not
  the drive, is what breaks the sine law for identical cells.

The design lesson, now measured rather than guessed: a
Josephson-like response needs **three** things — an odd drive, a
single-frequency phase structure, and a non-notched distance
factor. The linear generator fails the second; identical cells fail
the third; the binary/distinct combination satisfies all three.

## 5. The phase ramp (claim f): the notch train, and a pure sine where the binary generator rules

### Method

The phase is driven linearly in time, φ_t = 2πt/64 (one cycle per
64 steps, 256 steps = 4 cycles), and three time series are measured
on the real lattice: the v0.1.0 flow S(t), the variant flow
S_odd(t), and the variant drive Q(t). Two diagnostics per series:
FFT power at the ramp frequency (on the demeaned series) and a sine
fit at the ramp rate. No verdict is forced — the shape is reported.
Note the two diagnostics measure different things: the FFT is on
the demeaned series (it sees the wiggle), while the sine fit has no
DC term (it sees the offset too). A series can carry 99% of its
demeaned power at the ramp frequency and still fail the sine fit,
if it is a small wiggle on a large constant.

### Measured shapes

| Setup | Series | Mean | Range | Power at ramp freq. | Sine R² at ramp rate | δ |
|---|---|---|---|---|---|---|
| linear, identical | S(t) v0.1.0 | −1.59e-03 | [−1.77e-03, 0] | 0.135 | −32.75 | π/2 |
| linear, identical | S_odd(t) | ~0 | ±1.50e-04 | 0.037 | 0.0366 | 0.0000 |
| linear, identical | Q(t) | ~0 | ±3.91e-02 | 0.027 | 0.0266 | 0.0000 |
| linear, distinct | S(t) v0.1.0 | −1.64e-03 | [−1.72e-03, −1.57e-03] | 0.026 | −2290.99 | −0.34 |
| binary, identical | S(t) v0.1.0 | −6.15e-04 | [−1.77e-03, +9.6e-05] | 0.939 | 0.1158 | π/2 |
| binary, identical | S_odd(t) | ~0 | ±1.23e-04 | 0.922 | 0.9222 | 0.0000 |
| binary, identical | Q(t) | ~0 | ±3.16e-02 | 1.000 | 1.0000 | 0.0000 |
| binary, distinct | S(t) v0.1.0 | −1.65e-03 | [−1.68e-03, −1.62e-03] | 0.997 | −5025.61 | 0.84 |
| binary, distinct | S_odd(t) | ~0 | ±8.15e-05 | 1.000 | 0.9999 | −0.0079 |
| binary, distinct | Q(t) | 1.3e-04 | ±1.65e-02 | 1.000 | 0.9999 | −0.0079 |

What the time series actually do:

- **Primary (linear, identical): the v0.1.0 flow is a notch train.**
  S(t) sits at ≈ −0.0016 and dips to 0 at each cycle boundary
  (φ ≡ 0 mod 2π), four times in 256 steps. It does not oscillate at
  the ramp rate: only 13.5% of its demeaned power is at the ramp
  frequency, and the sine fit at the ramp rate is R² = −32.75 —
  the same number as claim (d), as it must be: the ramp samples the
  same φ-curve four times. The ramp reproduces (d) as a time series.
  The variant flow and drive under the linear generator also do not
  oscillate at the ramp rate (3.7% and 2.7% of power there; the
  spectra peak at high harmonics, bins 120 and 96 of 128) — the
  linear odd drive is a multi-frequency signal, not a ramp-rate
  oscillation.
- **Binary generator: the variant drive is a pure sine at the ramp
  rate.** Q(t) carries 100% of its power at the ramp frequency with
  sine-fit R² = 1.0000, δ = 0 (identical cells). The binary/identical
  variant flow oscillates at the ramp rate (92% of power, R² = 0.92,
  the distance-notch distortion again), and the binary/distinct
  variant flow is a clean sine at the ramp rate (100% of power,
  R² = 0.9999, δ = −0.008).
- **The v0.1.0 flow under the binary generator is periodic at the
  ramp rate but not a sine.** Binary/identical: 94% of demeaned
  power at the ramp frequency, yet sine R² = 0.12 (an offset,
  distorted wave, δ = π/2). Binary/distinct: 99.7% of demeaned
  power at the ramp frequency, yet sine R² = −5026 — a 4% wiggle on
  a constant −0.00165 offset. The FFT and the sine fit disagree
  exactly as the offset explains.

## 6. What the FAIL means (and does not mean)

- It does **not** mean the Josephson laws are wrong. They are
  established physics (Nobel 1973) and reproduce exactly here
  (claims a–c). The failure is on the membrane side: the lattice's
  phase response is not Josephson-like.
- It does **not** mean phase is irrelevant to the membrane. Phase
  offsets are a real, first-class parameter (round 2: an additive
  rate shift). What fails is specifically the sine law: the response
  is even (for identical cells), threshold-gated, and essentially
  constant, not odd and sinusoidal.
- It **does** mean the "phase difference drives a flow" bridge
  between the Josephson effect and the membrane is an analogy of
  shape, not of law. The shared shape is "a difference drives a
  flow"; the differing law is "sine (odd, sign-reversing)" versus
  "fidelity-gated repulsion (even, sign-fixed)".
- A falsification this clean is useful: any future membrane design
  that wants a Josephson-like response must make the interaction
  depend on a signed, odd function of the phase — e.g. on the
  imaginary part of the overlap ⟨ρ_A, ρ_B⟩ rather than on the
  fidelity |overlap|. That direction was built and measured in v0.2.0
  (§4–§5): it works exactly where the design lesson predicts
  (binary generator, distinct cells) and fails where it predicts
  (linear generator, identical cells).

## 7. Honest limits

1. **Classical emulation.** Everything here is NumPy on classical
   hardware. The lattice is NOT a superconductor; no Cooper pairs,
   no tunnelling, no 2e charge, no real Josephson junction was
   measured. The Josephson laws are inputs evaluated in code.
2. **The Josephson law is an input, not a discovery.** Claims (a)–(c)
   verify the implementation against tabulated values. They cannot
   validate the physics.
3. **The phase unitary is a convention.** U(Δφ) = diag(exp(i·k·Δφ))
   and the binary variant are two natural choices for "a phase
   difference between two states"; other generators (different level
   weightings, a single-coherence rotation) would give different
   wiggle shapes. Both generators measured here fail the sine law by
   wide margins, and the evenness argument covers every diagonal
   phase unitary for identical cells (for distinct cells the
   fidelity is not even; see §3) — but the exact wiggle amplitudes are
   generator-dependent.
4. **Random states.** The lattice cells are Ginibre random density
   matrices (32 dimensions). The fidelity between two independent
   random states sits near 0.57, far below the 0.90 resonance
   threshold. A lattice of nearly-identical states (fidelity near
   1) would have a different regime; the identical-cells setup
   covers the "differ only by phase" case, which is the Josephson
   analogue, and still fails.
5. **One-step and pairwise scope.** The sweep measures the
   instantaneous interaction, not a long-time steady state. A
   multi-step evolution could in principle develop structure the
   one-step response lacks — but the one-step response is the
   mechanism, and it shows no sine.
6. **No hardware data.** No SQUID, no junction, no cryostat. The
   SQUID numbers are evaluations of the formula, not measurements
   of a device.
7. **The verdict rule is a choice.** R² ≥ 0.95 and max|residual| ≤
   5% of amplitude is a strict bar, chosen before measuring. The
   measured numbers miss it by orders of magnitude (R² = −32.75,
   residual 12.5× amplitude), so the verdict does not hinge on the
   exact bar — any reasonable bar fails. The δ = π/2 finding (cosine,
   not sine) is independent of the bar.
8. **The variant interaction is a designed rule, not a discovered
   law.** The odd-function variant (claim e) was constructed to be
   odd; its binary/distinct pass is a property of that designed rule
   acting on the lattice, not evidence of a physical Josephson
   effect in the membrane. The one-sided overlap Q is a
   bookkeeping quantity (the imaginary part of an overlap between
   two states in different frames), not a current. No physical
   junction is implied.

## 8. Reproduction

```bash
cd ~/workspace/josephson-membrane
python3 -m pytest -q            # 49 tests, ~45 s
PYTHONPATH=src python3 -m josephson_membrane.unified   # full claim as JSON
```

Lattice provenance is recorded by the measurement itself:
`membrane_phase.lattice_provenance()` returns the tree path, branch
(`rounds-on-main`) and HEAD commit (`e7d8fea70cbc0429ad9acb1cea8d55e989cfca70`).
Seeds 0..4, 64 phase points, identical and distinct cells, linear and
binary generators. All numbers in this report come from that run.
