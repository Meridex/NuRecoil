"""nurecoil.flux — Neutrino flux module."""

from nurecoil.flux.base import FluxBase
from nurecoil.flux.phenomenological import PhenomenologicalFlux
from nurecoil.flux.interpolated import InterpolatedFlux

#: Short alias for :class:`PhenomenologicalFlux`.
PhenoFlux = PhenomenologicalFlux

__all__ = ["FluxBase", "PhenomenologicalFlux", "PhenoFlux", "InterpolatedFlux"]
