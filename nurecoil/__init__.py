"""
NuRecoil — CEvNS event-rate calculator for reactor-neutrino experiments.

Public API
----------
Nucleus          : target nucleus dataclass
compute_rate     : top-level event-rate integrator
constants        : physical constants and unit conversions
"""

from nurecoil.nucleus import Nucleus
from nurecoil.rate import compute_rate
from nurecoil import constants

__all__ = ["Nucleus", "compute_rate", "constants"]
