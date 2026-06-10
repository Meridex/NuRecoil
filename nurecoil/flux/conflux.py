"""
CONFLUX library wrapper — optional flux backend.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from nurecoil.flux.base import FluxBase


class ConfluxFlux(FluxBase):
    """
    Thin wrapper around the CONFLUX library.

    Raises ``ImportError`` with a helpful message if CONFLUX is not installed.
    Install with:  pip install nurecoil[conflux]

    ``E_min`` and ``E_max`` are set to 0 and infinity as placeholders;
    update them once the CONFLUX integration is implemented.
    """

    E_min: float = 0.0          # MeV — placeholder; update when implemented
    E_max: float = float("inf") # MeV — placeholder; update when implemented

    def __init__(self, *args, **kwargs) -> None:
        try:
            import conflux  # noqa: F401
        except ImportError as exc:
            raise ImportError(
                "The 'conflux' package is required for ConfluxFlux. "
                "Install it with: pip install nurecoil[conflux]"
            ) from exc
        self._args   = args
        self._kwargs = kwargs

    def _flux(self, E_nu: NDArray, P: float, L: float) -> NDArray:
        raise NotImplementedError(
            "ConfluxFlux._flux is not yet implemented. "
            "Contribute at https://github.com/your-org/nurecoil"
        )
