"""
Abstract base class and implementations for quenching factors.

Units
-----
    E_R  : MeV  (nuclear recoil energy, internal canonical unit)
    f_Q  : dimensionless
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike, NDArray

from nurecoil import constants
from nurecoil.nucleus import Nucleus


class QuenchingBase(ABC):
    """Abstract base class for quenching-factor models."""

    @abstractmethod
    def __call__(self, E_R: ArrayLike, nucleus: Nucleus) -> NDArray:
        """
        Evaluate the quenching factor f_Q(E_R).

        Parameters
        ----------
        E_R : array-like
            Nuclear recoil energy [MeV].
        nucleus : Nucleus
            Target nucleus (provides Z and A).

        Returns
        -------
        NDArray
            f_Q, dimensionless, same shape as *E_R*.  Values are in [0, 1].
        """
        ...


class LindhardQuenching(QuenchingBase):
    """
    Lindhard model quenching factor.

    The Lindhard formula is written in keV; E_R [MeV] is converted inside
    this class using ``constants.keV_per_MeV``.

    Formula
    -------
    ε = 11.5 · (E_R / keV) · Z^{-7/3}
    g = 3ε^{0.15} + 0.7ε^{0.6} + ε
    κ = 0.133 · Z^{2/3} · A^{-1/2}
    f_Q = κg / (1 + κg)

    Valid for E_R ≳ a few eV; returns f_Q = 0 below ``E_R_threshold_MeV``.
    """

    #: Minimum recoil energy below which f_Q is set to zero [MeV]
    E_R_threshold_MeV: float = 1.0e-6  # 1 eV

    def __call__(self, E_R: ArrayLike, nucleus: Nucleus) -> NDArray:
        E_R_arr = np.asarray(E_R, dtype=float)
        Z, A = nucleus.Z, nucleus.A

        # Convert MeV → keV for the Lindhard formula
        E_R_keV = E_R_arr * constants.keV_per_MeV

        epsilon = 11.5 * E_R_keV * Z ** (-7.0 / 3.0)
        g       = 3.0 * epsilon ** 0.15 + 0.7 * epsilon ** 0.6 + epsilon
        kappa   = 0.133 * Z ** (2.0 / 3.0) * A ** (-0.5)

        kg  = kappa * g
        f_Q = kg / (1.0 + kg)

        # Zero below threshold
        f_Q = np.where(E_R_arr >= self.E_R_threshold_MeV, f_Q, 0.0)
        return f_Q


class ConstantQuenching(QuenchingBase):
    """
    Trivial quenching factor: f_Q = *value* everywhere.

    Useful for electron-recoil targets (value=1) or unit tests.

    Parameters
    ----------
    value : float
        Constant quenching factor.  Default is 1.0.
    """

    def __init__(self, value: float = 1.0) -> None:
        self.value = float(value)

    def __call__(self, E_R: ArrayLike, nucleus: Nucleus) -> NDArray:
        return np.full_like(np.asarray(E_R, dtype=float), self.value)
