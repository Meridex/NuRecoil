"""
Energy-resolution smearing classes.

Units
-----
    E_ee  : MeV  (electron-equivalent energy = f_Q · E_R, internal canonical unit)
    E_det : MeV  (detected energy axis)
    ΔE    : MeV  (energy resolution, stored internally in MeV)
    f_res : MeV⁻¹  (normalised Gaussian kernel)

Literature formulas that use eV or keV must convert:
    Input  E_ee [MeV] → keV  via  constants.keV_per_MeV
    Output ΔE   [eV]  → MeV  via  1 / constants.eV_per_MeV
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike, NDArray

from nurecoil import constants


class ResolutionBase(ABC):
    """Abstract base class for energy-resolution models."""

    @abstractmethod
    def sigma(self, E_ee: ArrayLike) -> NDArray:
        """
        Energy resolution ΔE at electron-equivalent energy E_ee.

        Parameters
        ----------
        E_ee : array-like
            Electron-equivalent energy [MeV].

        Returns
        -------
        NDArray
            ΔE(E_ee)  [MeV].
        """
        ...

    def smear(self, E_ee: ArrayLike, E_det: ArrayLike) -> NDArray:
        """
        Gaussian resolution kernel f_res(E_ee, E_det).

        f_res = (1 / (√2π ΔE)) · exp[−(E_ee − E_det)² / (2 ΔE²)]

        All quantities in MeV; output integrates to 1 over E_det [MeV].

        Parameters
        ----------
        E_ee : array-like, shape (n,)
            Electron-equivalent energies [MeV].
        E_det : array-like, shape (m,)
            Detected energy grid [MeV].

        Returns
        -------
        NDArray, shape (n, m)
            f_res [MeV⁻¹].
        """
        E_ee_arr  = np.asarray(E_ee,  dtype=float)[:, np.newaxis]  # (n, 1)
        E_det_arr = np.asarray(E_det, dtype=float)[np.newaxis, :]  # (1, m)
        dE = self.sigma(E_ee_arr)                                   # (n, 1) [MeV]

        return (
            np.exp(-0.5 * ((E_ee_arr - E_det_arr) / dE) ** 2)
            / (np.sqrt(2.0 * np.pi) * dE)
        )


class CDEXResolution(ResolutionBase):
    """
    CDEX energy-resolution model.

    Literature formula (in experimental units):
        ΔE [eV] = 35.8 + 16.6 × (E [keV])^{1/2}

    where E is the electron-equivalent energy E^{ee} = f_Q · E_R.

    Unit conversion inside ``sigma``:
        E_keV  = E_ee_MeV × constants.keV_per_MeV
        ΔE_MeV = (35.8 + 16.6 × E_keV^{0.5}) / constants.eV_per_MeV
    """

    def sigma(self, E_ee: ArrayLike) -> NDArray:
        E_ee_arr = np.asarray(E_ee, dtype=float)
        E_keV    = E_ee_arr * constants.keV_per_MeV
        dE_eV    = 35.8 + 16.6 * np.sqrt(np.maximum(E_keV, 0.0))
        return dE_eV / constants.eV_per_MeV  # [MeV]


class ConstantResolution(ResolutionBase):
    """
    Fixed energy resolution.

    Parameters
    ----------
    sigma_MeV : float
        Constant resolution ΔE [MeV].
    """

    def __init__(self, sigma_MeV: float) -> None:
        self._sigma = float(sigma_MeV)

    def sigma(self, E_ee: ArrayLike) -> NDArray:
        return np.full_like(np.asarray(E_ee, dtype=float), self._sigma)
