r"""
Abstract base class for differential cross-section models.

Units
-----
    E_nu : MeV  (neutrino energy)
    E_R  : MeV  (nuclear recoil energy)
    $d\sigma/dE_R$ : cm^{2} / MeV
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike, NDArray

from nurecoil.nucleus import Nucleus


class CrossSectionBase(ABC):
    """
    Abstract base class for CEvNS (and BSM) differential cross sections.

    Parameters
    ----------
    nucleus : Nucleus
        Target nucleus.
    """

    def __init__(self, nucleus: Nucleus) -> None:
        self.nucleus = nucleus

    @abstractmethod
    def __call__(self, E_nu: ArrayLike, E_R: ArrayLike) -> NDArray:
        r"""
        Evaluate $d\sigma/dE_R$.

        Parameters
        ----------
        E_nu : array-like
            Neutrino energy [MeV].
        E_R : array-like
            Nuclear recoil energy [MeV].

        Returns
        -------
        NDArray
            $d\sigma/dE_R$  [cm^2 / MeV].  Must be >= 0 everywhere; return 0 (not
            NaN) outside the valid kinematic region.
        """
        ...

    def E_R_max(self, E_nu: float) -> float:
        """
        Kinematic upper bound on the recoil energy.

        Parameters
        ----------
        E_nu : float
            Neutrino energy [MeV].

        Returns
        -------
        float
            $E_R^{\\max}$ [MeV] $= 2 E_\\nu^2 / (M + 2 E_\\nu)$.
        """
        M = self.nucleus.M
        return 2.0 * E_nu ** 2 / (M + 2.0 * E_nu)
