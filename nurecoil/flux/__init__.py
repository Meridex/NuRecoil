"""nurecoil.flux — Neutrino flux module."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from nurecoil.flux.base import FluxBase
from nurecoil.flux.reactor_mix import ReactorMix
from nurecoil.flux.phenomenological import PhenomenologicalFlux
from nurecoil.flux.interpolated import InterpolatedFlux, IsotopeInterpolatedFlux, make_flux
from nurecoil.flux.tabulated import TabulatedFlux
from nurecoil.flux.data_loader import load_spectrum, load_all_isotopes, list_sources, list_isotopes

#: Short alias for :class:`PhenomenologicalFlux`.
PhenoFlux = PhenomenologicalFlux


def load_flux(spec: str, **kwargs: Any) -> FluxBase:
    """
    Unified flux factory — the recommended public entry point.

    *spec* is a string with a prefix that selects the underlying implementation:

    ========================  =============================================
    spec format               Result
    ========================  =============================================
    ``"pheno:<model>"``       :class:`PhenomenologicalFlux`
    ``"interp:<dataset>"``    :func:`make_flux` (auto-detects composite vs
                              isotope mode)
    ``"file:<path>"``         :class:`TabulatedFlux` loaded from *path*
    ``"file:built-in"``       :class:`TabulatedFlux` using the bundled
                              example file
    ``"conflux"``             :class:`~nurecoil.flux.conflux.ConfluxFlux`
    ========================  =============================================

    Extra keyword arguments are forwarded to the underlying constructor.
    Unsupported kwargs for a given prefix raise :class:`ValueError`.

    Parameters
    ----------
    spec : str
        Flux specification string (see table above).
    **kwargs
        Constructor keyword arguments relevant to the chosen backend:

        * ``pheno:`` accepts ``fission_fractions``, ``energy_per_fission``
        * ``interp:`` accepts ``fission_fractions``, ``energy_per_fission``,
          ``extrapolate``
        * ``file:`` accepts no extra kwargs
        * ``conflux`` accepts no extra kwargs

    Returns
    -------
    FluxBase
        A ready-to-call flux object.

    Examples
    --------
    >>> flux = load_flux("pheno:huber_mueller")
    >>> flux = load_flux("pheno:huber_mueller",
    ...                  fission_fractions="fractions:typical",
    ...                  energy_per_fission="energy:ma_2013")
    >>> flux = load_flux("interp:estienne2019")
    >>> flux = load_flux("interp:kopeikin2012",
    ...                  fission_fractions="fractions:daya_bay")
    >>> flux = load_flux("file:built-in")
    >>> flux = load_flux("file:/abs/path/to/my_flux.csv")
    """
    if ":" not in spec and spec != "conflux":
        raise ValueError(
            f"load_flux: spec must contain a prefix, e.g. 'pheno:huber_mueller'. "
            f"Got: {spec!r}"
        )

    prefix, _, body = spec.partition(":")

    if prefix == "pheno":
        _check_kwargs(kwargs, {"fission_fractions", "energy_per_fission"}, prefix)
        return PhenomenologicalFlux(
            spectrum_model=f"spectrum:{body}" if body else "huber_mueller",
            fission_fractions=kwargs.get("fission_fractions"),
            energy_per_fission=kwargs.get("energy_per_fission", "ma_2013"),
        )

    if prefix == "interp":
        _check_kwargs(kwargs, {"fission_fractions", "energy_per_fission", "extrapolate"}, prefix)
        return make_flux(
            source=body,
            fission_fractions=kwargs.get("fission_fractions"),
            energy_per_fission=kwargs.get("energy_per_fission", "ma_2013"),
            extrapolate=kwargs.get("extrapolate", False),
        )

    if prefix == "file":
        _check_kwargs(kwargs, set(), prefix)
        source = None if body == "built-in" else body
        return TabulatedFlux(source=source)

    if spec == "conflux":
        _check_kwargs(kwargs, set(), "conflux")
        from nurecoil.flux.conflux import ConfluxFlux
        return ConfluxFlux()

    raise ValueError(
        f"load_flux: unknown prefix '{prefix}'. "
        f"Valid prefixes: 'pheno', 'interp', 'file', or the bare spec 'conflux'."
    )


def _check_kwargs(
    kwargs: dict[str, Any],
    allowed: set[str],
    prefix: str,
) -> None:
    """Raise ValueError if any key in *kwargs* is not in *allowed*."""
    extra = set(kwargs) - allowed
    if extra:
        raise ValueError(
            f"load_flux '{prefix}:': unexpected keyword argument(s) {sorted(extra)}. "
            f"Allowed: {sorted(allowed) or '(none)'}."
        )


__all__ = [
    "FluxBase",
    "ReactorMix",
    "PhenomenologicalFlux",
    "PhenoFlux",
    "InterpolatedFlux",
    "IsotopeInterpolatedFlux",
    "TabulatedFlux",
    "make_flux",
    "load_flux",
    "load_spectrum",
    "load_all_isotopes",
    "list_sources",
    "list_isotopes",
]
