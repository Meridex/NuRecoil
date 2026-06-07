"""
Shared numerical-integration helpers.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def make_gauss_legendre_grid(a: float, b: float, n: int) -> tuple[NDArray, NDArray]:
    """
    Return nodes and weights for Gauss–Legendre quadrature on [a, b].

    Parameters
    ----------
    a, b : float
        Integration limits.
    n : int
        Number of quadrature points.

    Returns
    -------
    nodes : NDArray, shape (n,)
    weights : NDArray, shape (n,)
        ∫_a^b f(x) dx ≈ Σ weights[i] * f(nodes[i])
    """
    xi, wi = np.polynomial.legendre.leggauss(n)
    # Map from [-1, 1] to [a, b]
    nodes   = 0.5 * (b - a) * xi + 0.5 * (b + a)
    weights = 0.5 * (b - a) * wi
    return nodes, weights
