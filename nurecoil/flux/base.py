"""
Abstract base class for neutrino flux implementations.

Units (API boundary)
--------------------
    E_nu : MeV   (neutrino energy; MeV is the standard experimental convention)
    P    : GW    (reactor thermal power)
    L    : m     (baseline distance from reactor core to detector)

Units (internal / output)
--------------------------
    dΦ/dE_ν : # / MeV / cm²

Implementations must convert P and L on entry:
    L_cm          = L_m  * constants.cm_per_m
    fission_rate  = P_GW * constants.GW_to_MeV_per_s / constants.MeV_per_fission
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike, NDArray


class FluxBase(ABC):
    """
    Abstract base class for antineutrino flux models.

    All subclasses must implement ``__call__`` and return the differential
    flux dΦ/dE_ν in units of  **# / MeV / cm²**.
    """

    @abstractmethod
    def __call__(self, E_nu: ArrayLike, P: float, L: float) -> NDArray:
        """
        Evaluate the differential antineutrino flux.

        Parameters
        ----------
        E_nu : array-like
            Neutrino energies [MeV].
        P : float
            Reactor thermal power [GW].  Converted to fissions/s internally.
        L : float
            Baseline distance [m].  Converted to cm internally.

        Returns
        -------
        NDArray
            dΦ/dE_ν  [# / MeV / cm²],  same shape as *E_nu*.
        """
        ...
