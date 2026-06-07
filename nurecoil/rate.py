"""
Top-level event-rate integrator for CEvNS calculations.

The central formula is

    dR/dE_det = N_T ∫∫ dΦ/dE_ν · dσ/dE_R · f_res(f_Q(E_R)·E_R, E_det) dE_ν dE_R

Units
-----
    E_det      : MeV   (detected energy grid, internal canonical unit)
    E_nu_range : MeV
    P          : GW    (passed through to FluxBase; converted there)
    L          : m     (passed through to FluxBase; converted there)
    N_T        : dimensionless (number of target atoms)
    dR/dE_det  : counts / MeV  (per the exposure implicit in N_T)
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import dblquad

from nurecoil.flux.base import FluxBase
from nurecoil.cross_section.base import CrossSectionBase
from nurecoil.detector.quenching import QuenchingBase
from nurecoil.detector.resolution import ResolutionBase
from nurecoil.nucleus import Nucleus


def compute_rate(
    E_det: ArrayLike,
    flux: FluxBase,
    cross_section: CrossSectionBase,
    quenching: QuenchingBase,
    resolution: ResolutionBase,
    nucleus: Nucleus,
    N_T: float,
    P: float,
    L: float,
    E_nu_range: tuple[float, float] = (0.1, 10.0),
) -> NDArray:
    """
    Compute the differential event rate dR/dE_det.

    Parameters
    ----------
    E_det : array-like
        Detected energy bins [MeV].
    flux : FluxBase
        Neutrino flux model.  Accepts P [GW] and L [m] and returns
        dΦ/dE_ν [# / MeV / cm²].
    cross_section : CrossSectionBase
        Differential cross-section model.  Returns dσ/dE_R [cm² / MeV].
    quenching : QuenchingBase
        Quenching factor model.  Returns f_Q [dimensionless].
    resolution : ResolutionBase
        Energy-resolution model.  Returns f_res [MeV⁻¹].
    nucleus : Nucleus
        Target nucleus.
    N_T : float
        Number of target atoms [dimensionless].
    P : float
        Reactor thermal power [GW].  Passed to *flux*.
    L : float
        Baseline distance [m].  Passed to *flux*.
    E_nu_range : tuple[float, float]
        (E_nu_min, E_nu_max) integration range [MeV].  Default: (0.1, 10.0).

    Returns
    -------
    NDArray
        dR/dE_det  [counts / MeV],  same length as *E_det*.
    """
    E_det_arr = np.asarray(E_det, dtype=float)
    E_nu_min, E_nu_max = float(E_nu_range[0]), float(E_nu_range[1])

    result = np.zeros_like(E_det_arr)

    for idx, E_d in enumerate(E_det_arr):

        def integrand(E_R: float, E_nu: float) -> float:
            phi  = float(flux(E_nu, P, L))
            dsig = float(cross_section(E_nu, E_R))
            fq   = float(quenching(E_R, nucleus))
            E_ee = fq * E_R
            fres = float(resolution.smear([E_ee], [E_d])[0, 0])
            return phi * dsig * fres

        def E_R_upper(E_nu: float) -> float:
            return cross_section.E_R_max(E_nu)

        val, _ = dblquad(
            integrand,
            E_nu_min,
            E_nu_max,
            0.0,
            E_R_upper,
            limit=50,
        )
        result[idx] = N_T * val

    return result
