"""
CONFLUX library wrapper — optional flux backend.
"""

from __future__ import annotations

from numpy.typing import ArrayLike, NDArray

from nurecoil.flux.base import FluxBase


class ConfluxFlux(FluxBase):
    """
    Thin wrapper around the CONFLUX library.

    Raises ``ImportError`` with a helpful message if CONFLUX is not installed.
    Install with:  pip install nurecoil[conflux]
    """

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

    def __call__(self, E_nu: ArrayLike, P: float, L: float) -> NDArray:
        raise NotImplementedError(
            "ConfluxFlux.__call__ is not yet implemented. "
            "Contribute at https://github.com/your-org/nurecoil"
        )
