r"""
Interpolated antineutrino flux from user-supplied tabulated data.

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
from scipy.interpolate import CubicSpline

from nurecoil import constants
from nurecoil.flux.base import FluxBase
from nurecoil.flux.phenomenological import (
    FISSION_FRACTION_PRESETS,
    ENERGY_PER_FISSION_MODELS,
)


class InterpolatedFlux(FluxBase):
    r"""
    Antineutrino flux built from a user-supplied discrete spectrum table.

    The table gives the *normalised* spectrum per fission S(E_\nu) in
    units of  **# / MeV / fission**.  The class then applies the
    power-to-fission-rate and geometric 1/L² factors at call time.

    Parameters
    ----------
    E_nu_table : array-like, shape (n,)
        Neutrino energy grid [MeV].  Must be strictly increasing.
    spectrum_table : array-like, shape (n,)
        dN/dE_\nu per fission  [# / MeV / fission].
    extrapolate : bool
        If True, allow extrapolation outside the table range
        (the spline will extrapolate naturally).  If False (default),
        the flux is zero outside the table range.
    """

    # Placeholders; overwritten in __init__ with actual table bounds.
    E_min: float = 0.0
    E_max: float = 0.0

    def __init__(
        self,
        E_nu_table: ArrayLike,
        spectrum_table: ArrayLike,
        extrapolate: bool = False,
    ) -> None:
        E_arr = np.asarray(E_nu_table, dtype=float)
        S_arr = np.asarray(spectrum_table, dtype=float)
        if E_arr.shape != S_arr.shape or E_arr.ndim != 1:
            raise ValueError("E_nu_table and spectrum_table must be 1-D arrays of equal length.")
        if not np.all(np.diff(E_arr) > 0):
            raise ValueError("E_nu_table must be strictly increasing.")

        self._spline      = CubicSpline(E_arr, S_arr, extrapolate=extrapolate)
        self.E_min        = float(E_arr[0])
        self.E_max        = float(E_arr[-1])
        self._extrapolate = extrapolate

    def _flux(self, E_nu: NDArray, P: float, L: float) -> NDArray:
        r"""
        Evaluate d\Phi/dE_\nu [# / MeV / cm^2] for in-range energies.

        Parameters
        ----------
        E_nu : NDArray
            Neutrino energies [MeV], guaranteed within [E_min, E_max].
        P : float
            Reactor thermal power [GW].
        L : float
            Baseline distance [m].
        """
        # --- API boundary: convert P [GW] → fissions/s and L [m] → cm ---
        L_cm         = L * constants.cm_per_m
        fission_rate = P * constants.GW_to_MeV_per_s / constants.MeV_per_fission

        spectrum = self._spline(E_nu)  # # / MeV / fission

        # Clip to non-negative (spline can overshoot near table boundaries)
        spectrum = np.maximum(spectrum, 0.0)

        return fission_rate * spectrum / (4.0 * np.pi * L_cm ** 2)


# ---------------------------------------------------------------------------
# Type alias used in IsotopeInterpolatedFlux
# ---------------------------------------------------------------------------
_SpectrumEntry = tuple[ArrayLike, ArrayLike]  # (E_nu_table, spectrum_table)


class IsotopeInterpolatedFlux(FluxBase):
    r"""
    Reactor antineutrino flux from per-actinide tabulated spectra.

    Each actinide's spectrum is independently interpolated via a cubic spline,
    then combined using user-supplied (or preset) fission fractions:

    .. math::

        \frac{d\Phi}{dE_\nu} = \frac{1}{4\pi L^2}
            \frac{P_{\rm th}}{\bar{e}}
            \sum_i \frac{f_i}{F} S_i(E_\nu)

    where :math:`\bar{e} = \sum_i (f_i/F)\,e_i`.

    Parameters
    ----------
    spectra : dict[str, tuple[array-like, array-like]]
        Mapping of actinide name → ``(E_nu_table, spectrum_table)``.
        Keys must be a subset of ``{"U235", "U238", "Pu239", "Pu241"}``.
        Each ``spectrum_table`` gives :math:`dN/dE_\nu` per fission
        [# / MeV / fission].
    fission_fractions : str or dict[str, float], optional
        Preset name (``"typical"``, ``"ksnps"``, ``"conus"``, ``"daya_bay"``)
        or explicit dict mapping actinide name → fraction.  Must sum to 1.
        Defaults to ``"typical"``.
    energy_per_fission : str or dict[str, float], optional
        Preset name (``"ma_2013"``, ``"bemporad_2002"``, ``"meulenberg_1969"``)
        or explicit dict [MeV / fission].  Defaults to ``"ma_2013"``.
    extrapolate : bool, optional
        Allow spline extrapolation outside each table's range.  Defaults to
        ``False`` (values outside range are zeroed by the base class).

    Notes
    -----
    ``E_min`` and ``E_max`` are set to the **union** of all actinide table
    ranges (i.e. the widest interval covered by any isotope).  Energies
    outside an individual isotope's range contribute zero from that isotope
    but may still receive contributions from others.
    """

    # Placeholders; overwritten in __init__.
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

        # --- resolve fission fractions ---
        if fission_fractions is None:
            fission_fractions = "typical"
        if isinstance(fission_fractions, str):
            if fission_fractions not in FISSION_FRACTION_PRESETS:
                raise ValueError(
                    f"Unknown fission fraction preset '{fission_fractions}'. "
                    f"Choose from: {list(FISSION_FRACTION_PRESETS)}"
                )
            fission_fractions = dict(FISSION_FRACTION_PRESETS[fission_fractions])
        total = sum(fission_fractions.values())
        if not np.isclose(total, 1.0, atol=1e-3):
            raise ValueError(
                f"Fission fractions must sum to 1, got {total:.4f}."
            )

        # --- resolve energy per fission ---
        if isinstance(energy_per_fission, str):
            if energy_per_fission not in ENERGY_PER_FISSION_MODELS:
                raise ValueError(
                    f"Unknown energy-per-fission model '{energy_per_fission}'. "
                    f"Choose from: {list(ENERGY_PER_FISSION_MODELS)}"
                )
            energy_per_fission = dict(ENERGY_PER_FISSION_MODELS[energy_per_fission])

        # --- build per-actinide splines and determine energy range ---
        self._splines: dict[str, CubicSpline] = {}
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
            self._splines[iso] = CubicSpline(E_arr, S_arr, extrapolate=extrapolate)
            self._e_ranges[iso] = (float(E_arr[0]), float(E_arr[-1]))
            e_min_global = min(e_min_global, float(E_arr[0]))
            e_max_global = max(e_max_global, float(E_arr[-1]))

        self.E_min = e_min_global
        self.E_max = e_max_global
        self._extrapolate = extrapolate

        # --- pre-compute weighted average energy per fission ---
        self._fission_fractions = fission_fractions
        self._e_bar = sum(
            fission_fractions.get(iso, 0.0) * energy_per_fission[iso]
            for iso in known
            if iso in energy_per_fission
        )
        if self._e_bar <= 0:
            raise ValueError("Weighted average energy per fission is zero or negative.")

    def _flux(self, E_nu: NDArray, P: float, L: float) -> NDArray:
        r"""Evaluate d\Phi/dE_\nu [# / MeV / cm^2] for in-range energies."""
        L_cm = L * constants.cm_per_m
        fission_rate = P * constants.GW_to_MeV_per_s / self._e_bar

        # Weighted sum of per-actinide spectra
        combined = np.zeros_like(E_nu)
        for iso, spline in self._splines.items():
            frac = self._fission_fractions.get(iso, 0.0)
            if frac == 0.0:
                continue
            e_lo, e_hi = self._e_ranges[iso]
            in_range = (E_nu >= e_lo) & (E_nu <= e_hi)
            if not np.any(in_range):
                continue
            vals = spline(E_nu[in_range])
            combined[in_range] += frac * np.maximum(vals, 0.0)

        return fission_rate * combined / (4.0 * np.pi * L_cm ** 2)
