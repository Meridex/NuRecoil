"""nurecoil.flux — Neutrino flux module."""

from nurecoil.flux.base import FluxBase
from nurecoil.flux.reactor_mix import ReactorMix
from nurecoil.flux.phenomenological import PhenomenologicalFlux
from nurecoil.flux.interpolated import InterpolatedFlux, IsotopeInterpolatedFlux, make_flux
from nurecoil.flux.data_loader import load_spectrum, load_all_isotopes, list_sources, list_isotopes

#: Short alias for :class:`PhenomenologicalFlux`.
PhenoFlux = PhenomenologicalFlux

__all__ = [
    "FluxBase",
    "ReactorMix",
    "PhenomenologicalFlux",
    "PhenoFlux",
    "InterpolatedFlux",
    "IsotopeInterpolatedFlux",
    "make_flux",
    "load_spectrum",
    "load_all_isotopes",
    "list_sources",
    "list_isotopes",
]
