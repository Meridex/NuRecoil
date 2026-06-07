"""
Interpolated antineutrino flux from user-supplied tabulated data.

Units (API boundary)
--------------------
    E_nu : MeV
    P    : GW   → converted to fissions/s on entry
    L    : m    → converted to cm on entry

Units (internal / output)
--------------------------
    dΦ/dE_ν : # / MeV / cm²
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import CubicSpline

from nurecoil import constants
from nurecoil.flux.base import FluxBase


class InterpolatedFlux(FluxBase):
    """
    Antineutrino flux built from a user-supplied discrete spectrum table.

    The table gives the *normalised* spectrum per fission S(E_ν) in
    units of  **# / MeV / fission**.  The class then applies the
    power-to-fission-rate and geometric 1/L² factors at call time.

    Parameters
    ----------
    E_nu_table : array-like, shape (n,)
        Neutrino energy grid [MeV].  Must be strictly increasing.
    spectrum_table : array-like, shape (n,)
        dN/dE_ν per fission  [# / MeV / fission].
    extrapolate : bool
        If True, allow extrapolation outside the table range
        (the spline will extrapolate naturally).  If False (default),
        the flux is zero outside the table range.
    """

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
        self._E_min       = float(E_arr[0])
        self._E_max       = float(E_arr[-1])
        self._extrapolate = extrapolate

    def __call__(self, E_nu: ArrayLike, P: float, L: float) -> NDArray:
        """
        Evaluate dΦ/dE_ν [# / MeV / cm²].

        Parameters
        ----------
        E_nu : array-like
            Neutrino energies [MeV].
        P : float
            Reactor thermal power [GW].
        L : float
            Baseline distance [m].
        """
        E_nu_arr = np.asarray(E_nu, dtype=float)

        # --- API boundary: convert P [GW] → fissions/s and L [m] → cm ---
        L_cm         = L * constants.cm_per_m
        fission_rate = P * constants.GW_to_MeV_per_s / constants.MeV_per_fission

        spectrum = self._spline(E_nu_arr)  # # / MeV / fission

        if not self._extrapolate:
            outside = (E_nu_arr < self._E_min) | (E_nu_arr > self._E_max)
            spectrum = np.where(outside, 0.0, spectrum)

        # Clip to non-negative (spline can overshoot)
        spectrum = np.maximum(spectrum, 0.0)

        return fission_rate * spectrum / (4.0 * np.pi * L_cm ** 2)
