"""
Nuclear form-factor models.

Units
-----
    q (momentum transfer) : MeV  (hbar = c = 1)
    R, s (nuclear radii)  : fm
    qR argument           : dimensionless  (q [MeV] * R [fm] / hbar_c [MeV*fm])
    F(q)                  : dimensionless,  F(0) = 1

References
----------
arXiv:2203.07361  -- general weak form factor definition
arXiv:1902.07398  -- form-factor uncertainties for CEvNS
Phys. Rev. 104 (1956) 1466  -- Helm
Lewin & Smith, Astropart. Phys. 6 (1996) 87  -- Helm parametrisation
arXiv:hep-ph/0608035  -- Helm R1 formula
Phys. Rev. C 60 (1999) 014903 (arXiv:hep-ph/9902259)  -- Klein-Nystrand
arXiv:2104.01811  -- Klein-Nystrand parameters
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from scipy.special import spherical_jn
from numpy.typing import ArrayLike, NDArray

from nurecoil import constants
from nurecoil.nucleus import Nucleus


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _j1_over_x(x: NDArray) -> NDArray:
    """Return j1(x)/x, handling x=0 via the limit 1/3."""
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(x == 0.0, 1.0 / 3.0, spherical_jn(1, x) / x)


# ---------------------------------------------------------------------------
# Abstract base
# ---------------------------------------------------------------------------

class FormFactorBase(ABC):
    """Abstract base class for nuclear form factors."""

    @abstractmethod
    def __call__(self, q: ArrayLike, nucleus: Nucleus) -> NDArray:
        """
        Evaluate F(q) for a given nucleus.

        Parameters
        ----------
        q : array-like
            3-momentum transfer $|\\vec{q}|$ [MeV].
        nucleus : Nucleus
            Provides A, Z, N as needed.

        Returns
        -------
        NDArray
            $F(q)$, dimensionless, same shape as *q*.
        """
        ...


# ---------------------------------------------------------------------------
# Helm form factor
# ---------------------------------------------------------------------------

class HelmFormFactor(FormFactorBase):
    r"""
    Helm form factor (Lewin-Smith parametrisation).

    The nucleon distribution is the convolution of a uniform sphere of radius
    $R_1$ with a Gaussian surface of width $s$:

    .. math::

        F_{\mathrm{Helm}}(q) = \frac{3\,j_1(q R_1)}{q R_1}
                               \exp\!\left(-\frac{q^2 s^2}{2}\right)

    where $j_1$ is the spherical Bessel function of order 1 and the effective
    nuclear radius $R_1$ is obtained from the half-density radius $c$, the
    surface diffuseness $a$, and the skin thickness $s$ via

    .. math::

        R_1 = \sqrt{c^2 + \tfrac{7}{3}\pi^2 a^2 - 5s^2}, \quad
        c \approx 1.23\,A^{1/3} - 0.60\ \mathrm{fm}

    Default parameter values follow arXiv:hep-ph/0608035 and
    Lewin & Smith (1996):

    Parameters
    ----------
    s_fm : float
        Skin (surface) thickness [fm].  Default: 0.9 fm.
    a_fm : float
        Surface diffuseness [fm].  Default: 0.52 fm.
    c_coeff : float
        Coefficient in $c = c_{\mathrm{coeff}} \times A^{1/3} - 0.60$ [fm].
        Default: 1.23 fm.
    """

    def __init__(
        self,
        s_fm: float = 0.9,
        a_fm: float = 0.52,
        c_coeff: float = 1.23,
    ) -> None:
        self.s_fm    = s_fm
        self.a_fm    = a_fm
        self.c_coeff = c_coeff

    def __call__(self, q: ArrayLike, nucleus: Nucleus) -> NDArray:
        q_arr = np.asarray(q, dtype=float)
        A = nucleus.A

        c_fm  = self.c_coeff * A ** (1.0 / 3.0) - 0.60          # [fm]
        R1_fm = np.sqrt(
            max(c_fm ** 2 + (7.0 / 3.0) * np.pi ** 2 * self.a_fm ** 2
                - 5.0 * self.s_fm ** 2, 0.0)
        )                                                         # [fm]

        # Dimensionless arguments
        qR1 = q_arr * R1_fm  / constants.hbar_c  # q R1 / (hbar c)
        qs  = q_arr * self.s_fm / constants.hbar_c

        F = 3.0 * _j1_over_x(qR1) * np.exp(-0.5 * qs ** 2)
        return F


# ---------------------------------------------------------------------------
# Klein-Nystrand form factor
# ---------------------------------------------------------------------------

class KleinNystrandFormFactor(FormFactorBase):
    r"""
    Klein-Nystrand form factor.

    Obtained by folding a hard-sphere distribution of radius $R_A$ with a
    short-range Yukawa potential of range $a_k$:

    .. math::

        F_{\mathrm{KN}}(q) = \frac{3\,j_1(q R_A)}{q R_A}
                             \frac{1}{1 + q^2 a_k^2 / (\hbar c)^2}

    Default parameters follow arXiv:2104.01811:
    $r_0 = 1.3\ \mathrm{fm}$, $a_k = 0.7\ \mathrm{fm}$.

    Parameters
    ----------
    r0_fm : float
        Radius coefficient in $R_A = r_0\,A^{1/3}$ [fm].  Default: 1.3 fm.
    a_fm : float
        Yukawa range [fm].  Default: 0.7 fm.

    References
    ----------
    Phys. Rev. C 60 (1999) 014903 (arXiv:hep-ph/9902259)
    """

    def __init__(self, r0_fm: float = 1.3, a_fm: float = 0.7) -> None:
        self.r0_fm = r0_fm
        self.a_fm  = a_fm

    def __call__(self, q: ArrayLike, nucleus: Nucleus) -> NDArray:
        q_arr = np.asarray(q, dtype=float)

        R_A_fm = self.r0_fm * nucleus.A ** (1.0 / 3.0)  # [fm]

        qRA = q_arr * R_A_fm / constants.hbar_c           # dimensionless
        # Yukawa suppression: 1 / (1 + (q a_k / hbar_c)^2)
        yukawa = 1.0 / (1.0 + (q_arr * self.a_fm / constants.hbar_c) ** 2)

        F = 3.0 * _j1_over_x(qRA) * yukawa
        return F


# ---------------------------------------------------------------------------
# Gaussian form factor (sanity-check / fast approximation)
# ---------------------------------------------------------------------------

class GaussianFormFactor(FormFactorBase):
    """
    Simple Gaussian form factor -- fast, useful for sanity checks.

    .. math::

        F(q) = \\exp\\!\\left(-\\frac{q^2 R^2}{6\\,(\\hbar c)^2}\\right)

    with $R \\approx 1.2\\,A^{1/3}$ fm (RMS radius of a uniform sphere).

    Parameters
    ----------
    R_coeff : float
        Coefficient in $R = R_{\\mathrm{coeff}} \\times A^{1/3}$ [fm].  Default: 1.2 fm.
    """

    def __init__(self, R_coeff: float = 1.2) -> None:
        self.R_coeff = R_coeff

    def __call__(self, q: ArrayLike, nucleus: Nucleus) -> NDArray:
        q_arr = np.asarray(q, dtype=float)
        R_fm  = self.R_coeff * nucleus.A ** (1.0 / 3.0)  # [fm]
        qR_over_hbarc_sq = (q_arr * R_fm / constants.hbar_c) ** 2
        return np.exp(-qR_over_hbarc_sq / 6.0)
