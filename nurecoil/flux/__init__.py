"""
nurecoil.flux — Neutrino flux module.

Quick Reference for :func:`load_flux`
======================================

``load_flux(spec, **kwargs)`` selects a flux backend via the *spec* string.
All valid *spec* values and their accepted keyword arguments are listed below.

------------------------------------------------------------------------------
1. ``"pheno:<model>"``  →  PhenomenologicalFlux
------------------------------------------------------------------------------
Analytical exp-polynomial spectrum; valid energy range: 2–8 MeV for all models.

    ``"huber_mueller"``  *(default)*
        Huber 2011 for U235/Pu239/Pu241, Mueller 2011 for U238.
        The standard combination used in the reactor-neutrino community.
    ``"mueller_2011"``
        Mueller et al. 2011 — all four isotopes.
    ``"huber_2011"``
        Huber 2011 — U235/Pu239/Pu241 only (no U238 coefficients;
        raises an error if fission_fractions includes U238).
    ``"vogel_1985"``
        Vogel & Engel 1989 — simplified quadratic fit for all four isotopes.

Accepted kwargs:

    fission_fractions : str | dict | None
        Per-isotope fission fractions.  Built-in preset names:
            ``"typical"``    Typical commercial PWR (equilibrium)  *(default)*
            ``"ksnps"``      Kuo-Sheng Nuclear Power Station
            ``"conus"``      CONUS experiment
            ``"daya_bay"``   Daya Bay experiment
        A custom dict is also accepted, e.g. ``{"U235": 0.56, "U238": 0.08, ...}``.

    energy_per_fission : str | dict | None
        Effective thermal energy per fission [MeV].  Built-in preset names:
            ``"ma_2013"``        Ma et al. 2013          *(default)*
            ``"kopeikin_2004"``  Kopeikin et al. 2004
            ``"james_1969"``     James 1969
        A custom dict is also accepted, e.g. ``{"U235": 202.36, ...}``.

------------------------------------------------------------------------------
2. ``"interp:<dataset>"``  →  InterpolatedFlux / IsotopeInterpolatedFlux
------------------------------------------------------------------------------
Numerically interpolated tabulated spectra.  Mode (composite vs. isotope) is
selected automatically based on the bundled data files.

Isotope mode — per-actinide spectra (U235, U238, Pu239, Pu241):

    ``"estienne2019"``   Estienne et al. 2019   0.05 – 10.05 MeV
    ``"mueller2011"``    Mueller et al. 2011     2.0  –  8.0  MeV
    ``"vogel1989"``      Vogel & Engel 1989      ~0.008 – 2.0 MeV
    ``"CEA2023"``        CEA 2023                ~0.01 – 12.5 MeV

Composite mode — pre-weighted total spectrum:

    ``"kopeikin2012"``   Kopeikin et al. 2012   0.01 –  9.0 MeV
    ``"kopeikin1999"``   Kopeikin et al. 1999   ~0.005 – 1.5 MeV

Accepted kwargs:

    fission_fractions : str | dict | None
        Same preset names as above.
        (Only meaningful in isotope mode; ignored for composite datasets.)

    energy_per_fission : str | dict | None
        Same preset names as above.

    extrapolate : bool
        Allow extrapolation outside the tabulated energy range (default: False).

------------------------------------------------------------------------------
3. ``"file:<path>"``  →  TabulatedFlux
------------------------------------------------------------------------------
Load an arbitrary two-column CSV file (E_nu [MeV], dN_dE [# / MeV / fission]).
Energy range is determined by the file content.

    ``"file:built-in"``           Use the bundled example spectrum file. 0.0005 - 9.9995 MeV
    ``"file:/absolute/path.csv"`` Load a user-supplied CSV file.

Accepted kwargs: none.

------------------------------------------------------------------------------
4. ``"conflux"``  →  ConfluxFlux
------------------------------------------------------------------------------
Interface to the external ``conflux`` package (must be installed separately).

Accepted kwargs: none.

------------------------------------------------------------------------------
"""

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
