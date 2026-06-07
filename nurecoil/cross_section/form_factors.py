"""
Nuclear form-factor models.

Units
-----
    q (momentum transfer) : MeV  (ℏ = c = 1)
    R, s (nuclear radii)  : fm
    qR argument           : dimensionless  (q [MeV] × R [fm] / ℏc [MeV·fm])
    F(q)                  : dimensionless,  F(0) = 1
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from scipy.special import spherical_jn
from numpy.typing import ArrayLike, NDArray

from nurecoil import constants
from nurecoil.nucleus import Nucleus


class FormFactorBase(ABC):
    """Abstract base class for nuclear form factors."""

    @abstractmethod
    def __call__(self, q: ArrayLike, nucleus: Nucleus) -> NDArray:
        """
        Evaluate F(q) for a given nucleus.

        Parameters
        ----------
        q : array-like
            3-momentum transfer |q⃗| [MeV].
        nucleus : Nucleus
            Provides A (and Z if proton/neutron distinction is needed).

        Returns
        -------
        NDArray
            F(q), dimensionless, same shape as *q*.
        """
        ...


class HelmFormFactor(FormFactorBase):
    """
    Helm form factor — standard parametrisation used in most CEvNS analyses.

    F(q) = [3 j₁(q R₀) / (q R₀)] · exp(−(q s)² / 2)

    with
        j₁  : spherical Bessel function of order 1
        R₀  = √(R² − 5s²),  R ≈ 1.2 A^{1/3} fm,  s ≈ 1 fm  (skin thickness)

    The dimensionless argument is  q_fm = q [MeV] × R [fm] / ℏc [MeV·fm].

    Parameters
    ----------
    s_fm : float
        Nuclear skin thickness [fm].  Default: 1.0 fm.
    R_coeff : float
        Coefficient in  R = R_coeff × A^{1/3} [fm].  Default: 1.2 fm.
    """

    def __init__(self, s_fm: float = 1.0, R_coeff: float = 1.2) -> None:
        self.s_fm    = s_fm
        self.R_coeff = R_coeff

    def __call__(self, q: ArrayLike, nucleus: Nucleus) -> NDArray:
        q_arr = np.asarray(q, dtype=float)
        A = nucleus.A

        R_fm  = self.R_coeff * A ** (1.0 / 3.0)        # [fm]
        s_fm  = self.s_fm                               # [fm]
        R0_fm = np.sqrt(max(R_fm ** 2 - 5.0 * s_fm ** 2, 0.0))  # [fm]

        # Convert q [MeV] to dimensionless argument: q_fm = q * R0 / ℏc
        qR0 = q_arr * R0_fm / constants.hbar_c  # dimensionless
        qs  = q_arr * s_fm  / constants.hbar_c  # dimensionless

        # Avoid division by zero at q = 0
        with np.errstate(invalid="ignore", divide="ignore"):
            j1_over_x = np.where(
                qR0 == 0.0,
                1.0 / 3.0,
                spherical_jn(1, qR0) / qR0,
            )

        F = 3.0 * j1_over_x * np.exp(-0.5 * qs ** 2)
        return F


class GaussianFormFactor(FormFactorBase):
    """
    Simple Gaussian form factor — fast, useful for sanity checks.

    F(q) = exp(−q² R²  / (6 ℏc²))

    with  R ≈ 1.2 A^{1/3} fm  (RMS radius of a uniform sphere).

    Parameters
    ----------
    R_coeff : float
        Coefficient in  R = R_coeff × A^{1/3} [fm].  Default: 1.2 fm.
    """

    def __init__(self, R_coeff: float = 1.2) -> None:
        self.R_coeff = R_coeff

    def __call__(self, q: ArrayLike, nucleus: Nucleus) -> NDArray:
        q_arr = np.asarray(q, dtype=float)
        R_fm  = self.R_coeff * nucleus.A ** (1.0 / 3.0)  # [fm]
        # q² R² / (ℏc)²  (dimensionless)
        qR_over_hbarc_sq = (q_arr * R_fm / constants.hbar_c) ** 2
        return np.exp(-qR_over_hbarc_sq / 6.0)
