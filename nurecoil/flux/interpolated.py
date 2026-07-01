r"""
Interpolated antineutrino flux from user-supplied tabulated data.

Two internal modes
------------------
**Composite mode** (``InterpolatedFlux``)
    A single pre-weighted total spectrum ``S_total(E)`` [# / MeV / fission].
    Suitable for datasets that already combine all actinides (e.g.
    ``kopeikin2012``).  The effective energy per fission ``e_bar`` is supplied
    by the user (or via a :class:`~nurecoil.flux.reactor_mix.ReactorMix`
    preset) because the spectrum itself does not encode that information.

**Isotope mode** (``IsotopeInterpolatedFlux``)
    Per-actinide spectra ``S_i(E)`` [# / MeV / fission], combined on the fly
    with user-specified fission fractions.  Suitable for datasets that provide
    separate tables for U235, U238, Pu239, Pu241 (e.g. ``estienne2019``,
    ``mueller2011``).

Unified constructor
-------------------
:func:`make_flux` inspects the data directory and automatically picks the
right class::

    flux = make_flux("estienne2019")   # → IsotopeInterpolatedFlux
    flux = make_flux("kopeikin2012")   # → InterpolatedFlux (composite)

Units (API boundary)
--------------------
    E_nu : MeV
    P    : GW    converted to fissions/s on entry
    L    : m     converted to cm on entry

Units (internal / output)
--------------------------
    d\Phi/dE_\nu : # / MeV / cm²
"""

from __future__ import annotations

from typing import Union

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import PchipInterpolator

from nurecoil import constants
from nurecoil.flux.base import FluxBase
from nurecoil.flux.reactor_mix import ReactorMix
from nurecoil.flux.data_loader import load_spectrum, load_all_isotopes, list_isotopes


# ---------------------------------------------------------------------------
# InterpolatedFlux — composite (pre-weighted) spectrum
# ---------------------------------------------------------------------------

class InterpolatedFlux(FluxBase):
    r"""
    Antineutrino flux from a single pre-weighted total spectrum table.

    Use this class when the data source already provides a combined
    ``S_total(E)`` [# / MeV / fission] (e.g. Kopeikin 2012 Table 3).
    For per-actinide datasets use :class:`IsotopeInterpolatedFlux` or the
    convenience constructor :func:`make_flux`.

    The flux is:

    .. math::

        \frac{d\Phi}{dE_\nu} = \frac{1}{4\pi L^2}
            \frac{P_{\rm th}}{\bar{e}} \, S_{\rm total}(E_\nu)

    where :math:`\bar{e}` is the fission-fraction-weighted average thermal
    energy per fission supplied via *fission_fractions* and
    *energy_per_fission*.

    Parameters
    ----------
    E_nu_table : array-like, shape (n,)
        Neutrino energy grid [MeV].  Must be strictly increasing.
    spectrum_table : array-like, shape (n,)
        Total dN/dE_ν per fission [# / MeV / fission].
    fission_fractions : str or dict[str, float], optional
        Used only to compute ``e_bar``.  Preset name or explicit dict.
        Defaults to ``"typical"``.
    energy_per_fission : str or dict[str, float], optional
        Model name or explicit dict [MeV / fission].  Defaults to
        ``"ma_2013"``.
    extrapolate : bool
        Allow extrapolation outside the table range.  Defaults to ``False``.
    """

    E_min: float = 0.0
    E_max: float = 0.0

    def __init__(
        self,
        E_nu_table: ArrayLike,
        spectrum_table: ArrayLike,
        fission_fractions: Union[dict[str, float], str, None] = None,
        energy_per_fission: Union[dict[str, float], str] = "ma_2013",
        extrapolate: bool = False,
    ) -> None:
        E_arr = np.asarray(E_nu_table, dtype=float)
        S_arr = np.asarray(spectrum_table, dtype=float)
        if E_arr.shape != S_arr.shape or E_arr.ndim != 1:
            raise ValueError(
                "E_nu_table and spectrum_table must be 1-D arrays of equal length."
            )
        if not np.all(np.diff(E_arr) > 0):
            raise ValueError("E_nu_table must be strictly increasing.")

        mix = ReactorMix(
            fission_fractions=fission_fractions,
            energy_per_fission=energy_per_fission,
        )
        self._e_bar = mix.e_bar

        self._interp      = PchipInterpolator(E_arr, S_arr, extrapolate=extrapolate)
        self.E_min        = float(E_arr[0])
        self.E_max        = float(E_arr[-1])
        self._extrapolate = extrapolate

    def _flux(self, E_nu: NDArray, P: float, L: float) -> NDArray:
        r"""Evaluate dΦ/dE_ν [# / MeV / cm²] for in-range energies."""
        L_cm         = L * constants.cm_per_m
        fission_rate = P * constants.GW_to_MeV_per_s / self._e_bar

        spectrum = self._interp(E_nu)
        spectrum = np.maximum(spectrum, 0.0)  # PCHIP can dip slightly below 0
        return fission_rate * spectrum / (4.0 * np.pi * L_cm ** 2)

    @classmethod
    def from_composite(
        cls,
        source: str,
        fission_fractions: Union[dict[str, float], str, None] = None,
        energy_per_fission: Union[dict[str, float], str] = "ma_2013",
        extrapolate: bool = False,
    ) -> "InterpolatedFlux":
        """
        Construct from a bundled composite spectrum file.

        Parameters
        ----------
        source : str
            Citation key, e.g. ``"kopeikin2012"``.
        fission_fractions : str or dict, optional
            Used to compute ``e_bar``.  Defaults to ``"typical"``.
        energy_per_fission : str or dict, optional
            Defaults to ``"ma_2013"``.
        extrapolate : bool
            Defaults to ``False``.

        Examples
        --------
        >>> flux = InterpolatedFlux.from_composite("kopeikin2012")
        """
        E, S = load_spectrum(source)
        return cls(E, S,
                   fission_fractions=fission_fractions,
                   energy_per_fission=energy_per_fission,
                   extrapolate=extrapolate)


# ---------------------------------------------------------------------------
# Type alias
# ---------------------------------------------------------------------------
_SpectrumEntry = tuple[ArrayLike, ArrayLike]  # (E_nu_table, spectrum_table)


# ---------------------------------------------------------------------------
# IsotopeInterpolatedFlux — per-actinide spectra combined on the fly
# ---------------------------------------------------------------------------

class IsotopeInterpolatedFlux(FluxBase):
    r"""
    Reactor antineutrino flux from per-actinide tabulated spectra.

    Each actinide's spectrum is independently interpolated, then combined
    using fission fractions:

    .. math::

        \frac{d\Phi}{dE_\nu} = \frac{1}{4\pi L^2}
            \frac{P_{\rm th}}{\bar{e}}
            \sum_i \frac{f_i}{F} \, S_i(E_\nu)

    where :math:`\bar{e} = \sum_i (f_i/F)\,e_i`.

    Parameters
    ----------
    spectra : dict[str, tuple[array-like, array-like]]
        Mapping ``isotope → (E_nu_table, spectrum_table)``.
        Keys must be a subset of ``{"U235", "U238", "Pu239", "Pu241"}``.
    fission_fractions : str or dict[str, float], optional
        Preset name or explicit dict.  Must sum to 1.  Defaults to
        ``"typical"``.
    energy_per_fission : str or dict[str, float], optional
        Model name or explicit dict [MeV / fission].  Defaults to
        ``"ma_2013"``.
    extrapolate : bool, optional
        Allow spline extrapolation.  Defaults to ``False``.

    Notes
    -----
    ``E_min`` / ``E_max`` span the **union** of all isotope table ranges.
    Energies outside an individual isotope's range contribute zero from
    that isotope but may still receive contributions from others.
    """

    E_min: float = 0.0
    E_max: float = 0.0

    def __init__(
        self,
        spectra: dict[str, _SpectrumEntry],
        fission_fractions: Union[dict[str, float], str, None] = None,
        energy_per_fission: Union[dict[str, float], str] = "ma_2013",
        extrapolate: bool = False,
    ) -> None:
        known = {"U235", "U238", "Pu239", "Pu241"}
        unknown = set(spectra) - known
        if unknown:
            raise ValueError(
                f"Unknown actinide(s) in spectra: {unknown}. "
                f"Allowed keys: {known}"
            )
        if not spectra:
            raise ValueError("spectra must contain at least one actinide.")

        mix = ReactorMix(
            fission_fractions=fission_fractions,
            energy_per_fission=energy_per_fission,
        )

        self._interps: dict[str, PchipInterpolator] = {}
        self._e_ranges: dict[str, tuple[float, float]] = {}
        e_min_global = float("inf")
        e_max_global = float("-inf")

        for iso, (E_table, S_table) in spectra.items():
            E_arr = np.asarray(E_table, dtype=float)
            S_arr = np.asarray(S_table, dtype=float)
            if E_arr.shape != S_arr.shape or E_arr.ndim != 1:
                raise ValueError(
                    f"spectra['{iso}']: E_nu_table and spectrum_table must be "
                    "1-D arrays of equal length."
                )
            if not np.all(np.diff(E_arr) > 0):
                raise ValueError(
                    f"spectra['{iso}']: E_nu_table must be strictly increasing."
                )
            self._interps[iso] = PchipInterpolator(E_arr, S_arr, extrapolate=extrapolate)
            self._e_ranges[iso] = (float(E_arr[0]), float(E_arr[-1]))
            e_min_global = min(e_min_global, float(E_arr[0]))
            e_max_global = max(e_max_global, float(E_arr[-1]))

        self.E_min = e_min_global
        self.E_max = e_max_global
        self._extrapolate = extrapolate
        self._fission_fractions = mix.fractions
        self._e_bar = mix.e_bar

    # --- total flux ---------------------------------------------------------

    def _flux(self, E_nu: NDArray, P: float, L: float) -> NDArray:
        r"""Evaluate dΦ/dE_ν [# / MeV / cm²] for in-range energies."""
        L_cm         = L * constants.cm_per_m
        fission_rate = P * constants.GW_to_MeV_per_s / self._e_bar

        combined = np.zeros_like(E_nu)
        for iso, interp in self._interps.items():
            frac = self._fission_fractions.get(iso, 0.0)
            if frac == 0.0:
                continue
            e_lo, e_hi = self._e_ranges[iso]
            mask = (E_nu >= e_lo) & (E_nu <= e_hi)
            if not np.any(mask):
                continue
            combined[mask] += frac * np.maximum(interp(E_nu[mask]), 0.0)

        return fission_rate * combined / (4.0 * np.pi * L_cm ** 2)

    # --- per-isotope interface ----------------------------------------------

    @property
    def isotopes(self) -> list[str]:
        """List of actinide names available in this instance."""
        return list(self._interps)

    def isotope_spectrum(
        self,
        isotope: str,
        E_nu: ArrayLike,
    ) -> NDArray:
        """
        Return the raw per-fission spectrum ``S_i(E_ν)`` for one actinide.

        This is the un-weighted, un-normalised spectrum table interpolated at
        *E_nu*, in units of **# / MeV / fission**.  Values outside the
        isotope's table range are returned as zero.

        Parameters
        ----------
        isotope : str
            One of :attr:`isotopes`.
        E_nu : array-like
            Neutrino energies [MeV].

        Returns
        -------
        ndarray
            ``S_i(E_nu)`` [# / MeV / fission].

        Examples
        --------
        >>> flux = IsotopeInterpolatedFlux.from_source("estienne2019")
        >>> S_u235 = flux.isotope_spectrum("U235", np.linspace(2, 8, 200))
        """
        if isotope not in self._interps:
            raise ValueError(
                f"Isotope '{isotope}' not in this instance. "
                f"Available: {self.isotopes}"
            )
        E_arr = np.asarray(E_nu, dtype=float)
        e_lo, e_hi = self._e_ranges[isotope]
        result = np.zeros_like(E_arr)
        mask = (E_arr >= e_lo) & (E_arr <= e_hi)
        if np.any(mask):
            result[mask] = np.maximum(self._interps[isotope](E_arr[mask]), 0.0)
        return result

    def isotope_flux(
        self,
        isotope: str,
        E_nu: ArrayLike,
        P: float,
        L: float,
    ) -> NDArray:
        """
        Return the flux contribution from a single actinide.

        Applies the fission fraction weight, the power-to-fission-rate
        conversion, and the geometric factor — giving the partial flux that
        this isotope contributes to the total:

        .. math::

            \\frac{d\\Phi_i}{dE_\\nu} = \\frac{f_i}{4\\pi L^2}
                \\frac{P_{\\rm th}}{\\bar{e}} \\, S_i(E_\\nu)

        Parameters
        ----------
        isotope : str
            One of :attr:`isotopes`.
        E_nu : array-like
            Neutrino energies [MeV].
        P : float
            Reactor thermal power [GW].
        L : float
            Baseline distance [m].

        Returns
        -------
        ndarray
            Partial flux [# / MeV / cm²].

        Examples
        --------
        >>> flux = IsotopeInterpolatedFlux.from_source("estienne2019")
        >>> phi_u235 = flux.isotope_flux("U235", E_nu, P=3.0, L=500.0)
        """
        if isotope not in self._interps:
            raise ValueError(
                f"Isotope '{isotope}' not in this instance. "
                f"Available: {self.isotopes}"
            )
        frac = self._fission_fractions.get(isotope, 0.0)
        S = self.isotope_spectrum(isotope, E_nu)
        L_cm = L * constants.cm_per_m
        fission_rate = P * constants.GW_to_MeV_per_s / self._e_bar
        return frac * fission_rate * S / (4.0 * np.pi * L_cm ** 2)

    # --- constructors -------------------------------------------------------

    @classmethod
    def from_source(
        cls,
        source: str,
        fission_fractions: Union[dict[str, float], str, None] = None,
        energy_per_fission: Union[dict[str, float], str] = "ma_2013",
        extrapolate: bool = False,
    ) -> "IsotopeInterpolatedFlux":
        """
        Construct by loading all per-actinide spectra for *source*.

        Parameters
        ----------
        source : str
            Citation key, e.g. ``"estienne2019"``, ``"mueller2011"``.
        fission_fractions : str or dict, optional
            Defaults to ``"typical"``.
        energy_per_fission : str or dict, optional
            Defaults to ``"ma_2013"``.
        extrapolate : bool
            Defaults to ``False``.

        Examples
        --------
        >>> flux = IsotopeInterpolatedFlux.from_source("estienne2019")
        >>> flux = IsotopeInterpolatedFlux.from_source(
        ...     "mueller2011", fission_fractions="daya_bay"
        ... )
        """
        spectra = load_all_isotopes(source)
        return cls(
            spectra,
            fission_fractions=fission_fractions,
            energy_per_fission=energy_per_fission,
            extrapolate=extrapolate,
        )


# ---------------------------------------------------------------------------
# make_flux — unified constructor that auto-detects mode
# ---------------------------------------------------------------------------

def make_flux(
    source: str,
    fission_fractions: Union[dict[str, float], str, None] = None,
    energy_per_fission: Union[dict[str, float], str] = "ma_2013",
    extrapolate: bool = False,
) -> Union[IsotopeInterpolatedFlux, InterpolatedFlux]:
    """
    Build a flux object from a bundled data source, auto-detecting mode.

    If the source provides per-actinide (``isospec``) files, returns an
    :class:`IsotopeInterpolatedFlux`.  If only a composite total spectrum is
    available, returns an :class:`InterpolatedFlux`.

    Parameters
    ----------
    source : str
        Citation key, e.g. ``"estienne2019"`` or ``"kopeikin2012"``, or
        with the ``"dataset:<name>"`` prefix form.
    fission_fractions : str or dict, optional
        Passed to the underlying class.  Defaults to ``"typical"``.
        Supports the ``"fractions:<name>"`` prefix form.
    energy_per_fission : str or dict, optional
        Passed to the underlying class.  Defaults to ``"ma_2013"``.
        Supports the ``"energy:<name>"`` prefix form.
    extrapolate : bool
        Passed to the underlying class.  Defaults to ``False``.

    Returns
    -------
    :class:`IsotopeInterpolatedFlux` or :class:`InterpolatedFlux`

    Examples
    --------
    >>> flux = make_flux("estienne2019")               # bare name
    >>> flux = make_flux("dataset:estienne2019")       # prefix form
    >>> flux = make_flux("kopeikin2012")               # InterpolatedFlux
    >>> phi  = flux(E_nu, P=3.0, L=500.0)
    """
    # Support "dataset:<name>" prefix as well as bare source name.
    source = source.removeprefix("dataset:")
    if list_isotopes(source):
        return IsotopeInterpolatedFlux.from_source(
            source,
            fission_fractions=fission_fractions,
            energy_per_fission=energy_per_fission,
            extrapolate=extrapolate,
        )
    else:
        return InterpolatedFlux.from_composite(
            source,
            fission_fractions=fission_fractions,
            energy_per_fission=energy_per_fission,
            extrapolate=extrapolate,
        )
