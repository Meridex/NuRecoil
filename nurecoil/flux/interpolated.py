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

import numpy as np
from numpy.typing import ArrayLike, NDArray  # ArrayLike kept for CubicSpline input annotation
from scipy.interpolate import CubicSpline

from nurecoil import constants
from nurecoil.flux.base import FluxBase


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
