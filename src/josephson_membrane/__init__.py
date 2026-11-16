"""Josephson-membrane: the Josephson laws measured against the 3x3 membrane lattice.

SYNTHETIC classical emulation only. No claim about physical reality.
No headline claim is emitted (see unified.unified_claim()).

Import the functions from their modules, e.g.::

    from josephson_membrane.josephson import dc_current, ac_frequency_hz, squid_critical_current
    from josephson_membrane.membrane_phase import run_membrane_experiment, sine_fit
    from josephson_membrane.unified import unified_claim
"""

from . import josephson, membrane_phase, unified

__version__ = "0.1.0"

__all__ = ["josephson", "membrane_phase", "unified"]
