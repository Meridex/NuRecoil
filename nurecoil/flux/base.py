r"""
Abstract base class for neutrino flux implementations.

Units (API boundary)
--------------------
    E_nu : MeV   (neutrino energy; MeV is the standard experimental convention)
    P    : GW    (reactor thermal power)
    L    : m     (baseline distance from reactor core to detector)

Units (internal / output)
--------------------------
    $d\Phi/dE_\nu$ : # / MeV / cm^2

Implementations must convert P and L on entry:
    L_cm          = L_m  * constants.cm_per_m
    fission_rate  = P_GW * constants.GW_to_MeV_per_s / constants.MeV_per_fission

Subclass contract
-----------------
Every concrete subclass must define class or instance attributes:
    E_min : float   - lower bound of valid neutrino energy range [MeV]
    E_max : float   - upper bound of valid neutrino energy range [MeV]

Missing either attribute raises ``TypeError`` at class-definition time.
Subclasses implement ``_flux`` (not ``__call__``); the base ``__call__``
handles range-checking and zeroing automatically.
"""

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import ArrayLike, NDArray


class FluxBase(ABC):
    r"""
    Abstract base class for antineutrino flux models.

    Subclasses must:

    * Define ``E_min`` and ``E_max`` class (or instance) attributes [MeV].
    * Implement ``_flux``, which receives only energies already validated to
      be within ``[E_min, E_max]`` and returns $d\Phi/dE_\nu$
      [# / MeV / cm^2].

    ``__call__`` is provided by this base class. It masks out-of-range
    energies to zero. If ``warn_out_of_range`` is True, a ``UserWarning``
    is emitted when any out-of-range value is found.
    """

    E_min: float  # MeV - lower bound of valid neutrino energy range
    E_max: float  # MeV - upper bound of valid neutrino energy range
    warn_out_of_range: bool = False

    def __init_subclass__(cls, **kwargs: object) -> None:
        """Enforce that every concrete subclass declares E_min and E_max."""
        super().__init_subclass__(**kwargs)
        # Skip the check for abstract subclasses (they still have abstractmethods).
        if not getattr(cls, "__abstractmethods__", None):
            for attr in ("E_min", "E_max"):
                if not hasattr(cls, attr):
                    raise TypeError(
                        f"{cls.__name__} must define '{attr}' "
                        f"(valid neutrino energy range in MeV)."
                    )

    def __call__(self, E_nu: ArrayLike, P: float, L: float) -> NDArray:
        r"""
        Evaluate the differential antineutrino flux.

        Energies outside ``[E_min, E_max]`` are set to zero.
        If ``warn_out_of_range`` is True, a ``UserWarning`` is issued if
        any such values are present.

        Parameters
        ----------
        E_nu : array-like
            Neutrino energies [MeV].
        P : float
            Reactor thermal power [GW]. Converted to fissions/s internally.
        L : float
            Baseline distance [m]. Converted to cm internally.

        Returns
        -------
        NDArray
            $d\Phi/dE_\nu$  [# / MeV / cm^2], same shape as *E_nu*.
        """
        E_nu_arr = np.asarray(E_nu, dtype=float)
        in_range = (E_nu_arr >= self.E_min) & (E_nu_arr <= self.E_max)

        if self.warn_out_of_range and not np.all(in_range):
            warnings.warn(
                f"{type(self).__name__}: some E_nu values lie outside the "
                f"valid range [{self.E_min}, {self.E_max}] MeV. "
                "Flux set to 0 for those energies.",
                UserWarning,
                stacklevel=2,
            )

        result = np.zeros_like(E_nu_arr)
        if np.any(in_range):
            result[in_range] = self._flux(E_nu_arr[in_range], P, L)
        return result

    @abstractmethod
    def _flux(self, E_nu: NDArray, P: float, L: float) -> NDArray:
        r"""
        Compute $d\Phi/dE_\nu$ for energies guaranteed to be in ``[E_min, E_max]``.

        Parameters
        ----------
        E_nu : NDArray
            Neutrino energies [MeV], all within ``[E_min, E_max]``.
        P : float
            Reactor thermal power [GW]. Converted to fissions/s internally.
        L : float
            Baseline distance [m]. Converted to cm internally.

        Returns
        -------
        NDArray
            $d\Phi/dE_\nu$ [# / MeV / cm^2], same shape as *E_nu*.
        """
        ...
