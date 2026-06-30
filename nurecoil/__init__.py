"""
NuRecoil — CEvNS event-rate calculator for reactor-neutrino experiments.

Public API
----------
Nucleus                 : target nucleus dataclass
compute_rate            : backward-compatible differential-rate integrator
compute_binned_spectrum : binned event-rate integrator
constants               : physical constants and unit conversions
"""

from nurecoil.nucleus import Nucleus
from nurecoil.rate import (
    compute_binned_spectrum,
    compute_cevns_spectrum,
    compute_differential_rate,
    compute_eves_spectrum,
    compute_rate,
    target_count_from_mass,
)
from nurecoil import constants

__all__ = [
    "Nucleus",
    "compute_binned_spectrum",
    "compute_cevns_spectrum",
    "compute_differential_rate",
    "compute_eves_spectrum",
    "compute_rate",
    "target_count_from_mass",
    "constants",
]
