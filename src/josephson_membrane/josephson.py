# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2026 Jean-Francois Brisson (Spark AI NLP)

"""Josephson laws: DC current, AC frequency, SQUID flux modulation.

Established physics. Brian Josephson predicted the DC and AC effects
in 1962; the Nobel Prize in Physics followed in 1973. These are laws
taken as INPUT to this package, not discoveries of it -- see REPORT.md,
"Honest limits".

Constants are the exact SI values fixed in 2019 (elementary charge e
and Planck constant h are exact by definition), so 2e/h and h/2e are
computed here, not copied from a table. They reproduce the CODATA 2018
values (K_J = 483597.848416984 GHz/V, Phi_0 = 2.067833848e-15 Wb) to
machine precision.
"""

from __future__ import annotations

import math

# Exact SI constants (2019 redefinition).
E_CHARGE_C = 1.602176634e-19  # elementary charge, coulomb
H_PLANCK_J_S = 6.62607015e-34  # Planck constant, joule-second

# Josephson constant K_J = 2e/h, in Hz per volt. Computed, not tabulated.
JOSEPHSON_CONSTANT_HZ_PER_V = 2.0 * E_CHARGE_C / H_PLANCK_J_S
# Same constant expressed as GHz per mV (the unit the infographic uses).
JOSEPHSON_CONSTANT_GHZ_PER_MV = JOSEPHSON_CONSTANT_HZ_PER_V / 1e12

# Magnetic flux quantum Phi_0 = h / 2e, in weber. Computed, not tabulated.
FLUX_QUANTUM_WB = H_PLANCK_J_S / (2.0 * E_CHARGE_C)


def dc_current(delta_phi: float, i_c: float = 1.0) -> float:
    """DC Josephson current I_s = I_c * sin(delta_phi).

    delta_phi is the phase difference across the junction in radians;
    i_c is the junction critical current (same units as the result).
    The law is odd in delta_phi and 2*pi-periodic.
    """
    return float(i_c * math.sin(delta_phi))


def dc_cardinal_table(i_c: float = 1.0) -> dict[str, float]:
    """The five cardinal points of the DC law (the infographic's table).

    delta_phi = 0   -> I_s = 0;  pi/2 -> +I_c;  pi -> 0;
    3*pi/2 -> -I_c;  2*pi -> 0.
    """
    return {
        "0": dc_current(0.0, i_c),
        "pi/2": dc_current(math.pi / 2.0, i_c),
        "pi": dc_current(math.pi, i_c),
        "3pi/2": dc_current(3.0 * math.pi / 2.0, i_c),
        "2pi": dc_current(2.0 * math.pi, i_c),
    }


def ac_frequency_hz(voltage_v: float) -> float:
    """AC Josephson frequency f = 2eV / h, in hertz.

    A DC voltage V across the junction makes the supercurrent oscillate
    at f = K_J * V with K_J = 2e/h. Voltage in volts, frequency in Hz.
    """
    return float(JOSEPHSON_CONSTANT_HZ_PER_V * voltage_v)


def ac_frequency_ghz(voltage_v: float) -> float:
    """AC Josephson frequency f = 2eV / h, in gigahertz."""
    return ac_frequency_hz(voltage_v) / 1e9


def squid_critical_current(flux_wb: float, i_c: float = 1.0) -> float:
    """SQUID critical current I_c(Phi) = 2 * I_c * |cos(pi * Phi / Phi_0)|.

    A superconducting loop interrupted by two identical junctions has a
    flux-modulated critical current. flux_wb is the magnetic flux through
    the loop in weber; i_c is one junction's critical current.
    """
    return float(2.0 * i_c * abs(math.cos(math.pi * flux_wb / FLUX_QUANTUM_WB)))


def squid_cardinal_table(i_c: float = 1.0) -> dict[str, float]:
    """Cardinal points of the SQUID modulation.

    Phi = 0      -> 2*I_c;  Phi_0/4 -> sqrt(2)*I_c;  Phi_0/2 -> 0;
    Phi_0 -> 2*I_c.
    """
    phi0 = FLUX_QUANTUM_WB
    return {
        "0": squid_critical_current(0.0, i_c),
        "Phi0/4": squid_critical_current(phi0 / 4.0, i_c),
        "Phi0/2": squid_critical_current(phi0 / 2.0, i_c),
        "Phi0": squid_critical_current(phi0, i_c),
    }
