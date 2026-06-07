"""
Huber–Mueller phenomenological antineutrino flux.

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

from nurecoil import constants
from nurecoil.flux.base import FluxBase


# Huber (2011) + Mueller et al. (2011) polynomial coefficients
# S_k(E_nu) = exp(Σ a_i^(k) E_nu^i)  with E_nu in MeV
# Keys: fissile isotope labels
_HUBER_MUELLER_COEFFS: dict[str, list[float]] = {
    # 235U — Huber 2011 (table III)
    "U235": [
         4.367,  -4.577, 2.100, -5.294e-1,  6.186e-2,  -2.777e-3,
    ],
    # 239Pu — Huber 2011 (table III)
    "Pu239": [
         4.757,  -5.392, 2.563, -6.596e-1,  7.820e-2,  -3.536e-3,
    ],
    # 241Pu — Huber 2011 (table III)
    "Pu241": [
         2.990,  -2.882, 1.278, -3.343e-1,  3.905e-2,  -1.754e-3,
    ],
    # 238U — Mueller et al. 2011 (table VI)
    "U238": [
         4.833e-1, 1.927e-1, -1.283e-1, -6.762e-3,  2.233e-3, -1.536e-4,
    ],
}

# Default fission fractions (PWR core, approximate equilibrium)
_DEFAULT_FISSION_FRACTIONS: dict[str, float] = {
    "U235":  0.570,
    "Pu239": 0.295,
    "Pu241": 0.070,
    "U238":  0.065,
}


class PhenomenologicalFlux(FluxBase):
    """
    Phenomenological antineutrino flux using Huber–Mueller parametrisation.

    Parameters
    ----------
    fission_fractions : dict[str, float] | None
        Per-isotope fission fractions {isotope: fraction}.
        Must sum to 1.  Defaults to ``_DEFAULT_FISSION_FRACTIONS``.
    """

    def __init__(
        self,
        fission_fractions: dict[str, float] | None = None,
    ) -> None:
        if fission_fractions is None:
            fission_fractions = dict(_DEFAULT_FISSION_FRACTIONS)
        total = sum(fission_fractions.values())
        if not np.isclose(total, 1.0, atol=1e-3):
            raise ValueError(
                f"Fission fractions must sum to 1, got {total:.4f}"
            )
        self._fractions = fission_fractions

    def _spectrum_per_fission(self, E_nu: NDArray) -> NDArray:
        """Combined antineutrino spectrum [# / MeV / fission]."""
        spectrum = np.zeros_like(E_nu)
        for isotope, fraction in self._fractions.items():
            coeffs = _HUBER_MUELLER_COEFFS[isotope]
            exponent = sum(
                a * E_nu ** i for i, a in enumerate(coeffs)
            )
            spectrum += fraction * np.exp(exponent)
        return spectrum

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

        spectrum = self._spectrum_per_fission(E_nu_arr)  # # / MeV / fission

        # dΦ/dE_ν = fission_rate × spectrum × 1 / (4π L²)
        return fission_rate * spectrum / (4.0 * np.pi * L_cm ** 2)
