"""nurecoil.detector — Detector-response module."""

from nurecoil.detector.quenching import QuenchingBase, LindhardQuenching, ConstantQuenching
from nurecoil.detector.resolution import ResolutionBase, CDEXResolution, ConstantResolution

__all__ = [
    "QuenchingBase",
    "LindhardQuenching",
    "ConstantQuenching",
    "ResolutionBase",
    "CDEXResolution",
    "ConstantResolution",
]
