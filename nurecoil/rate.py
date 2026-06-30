"""
Top-level event-rate integrators for CEvNS-like calculations.

The central formula is

    dR/dE_det = N_T ∫∫ dΦ/dE_ν · dσ/dE_R · f_res(f_Q(E_R)·E_R, E_det) dE_ν dE_R

Units
-----
    E_det      : MeV   (detected energy grid, internal canonical unit)
    E_nu_range : MeV
    P          : GW    (passed through to FluxBase; converted there)
    L          : m     (passed through to FluxBase; converted there)
    N_T        : dimensionless (number of target atoms)
    dR/dE_det  : counts / s / MeV
    bin output : counts / s / bin, or counts / bin when exposure_s is supplied
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import dblquad

from nurecoil import constants
from nurecoil.flux.base import FluxBase
from nurecoil.cross_section.base import CrossSectionBase
from nurecoil.cross_section.sm_cevns import SMCEvNS
from nurecoil.cross_section.sm_eves import SMEvES
from nurecoil.detector.quenching import QuenchingBase
from nurecoil.detector.quenching import ConstantQuenching, LindhardQuenching
from nurecoil.detector.resolution import ResolutionBase
from nurecoil.nucleus import Nucleus


def target_count_from_mass(target_mass_g: float, nucleus: Nucleus) -> float:
    """
    Convert target mass to number of target atoms.

    Parameters
    ----------
    target_mass_g : float
        Target mass [g].
    nucleus : Nucleus
        Target nucleus. The molar mass is approximated as ``A`` [g/mol].

    Returns
    -------
    float
        Number of target atoms.
    """
    target_mass = float(target_mass_g)
    if target_mass <= 0.0:
        raise ValueError("target_mass_g must be positive.")
    return target_mass / float(nucleus.A) * constants.N_A


def _resolve_target_count(
    *,
    nucleus: Nucleus,
    N_T: float | None = None,
    target_mass_g: float | None = None,
) -> float:
    if N_T is None and target_mass_g is None:
        raise ValueError("Provide either N_T or target_mass_g.")
    if N_T is not None and target_mass_g is not None:
        raise ValueError("Provide only one of N_T or target_mass_g.")
    if N_T is not None:
        target_count = float(N_T)
        if target_count <= 0.0:
            raise ValueError("N_T must be positive.")
        return target_count
    assert target_mass_g is not None
    return target_count_from_mass(target_mass_g, nucleus)


def _validate_bin_edges(bin_edges: ArrayLike) -> NDArray:
    edges = np.asarray(bin_edges, dtype=float)
    if edges.ndim != 1:
        raise ValueError("bin_edges must be one-dimensional.")
    if edges.size < 2:
        raise ValueError("bin_edges must contain at least two edges.")
    if not np.all(np.isfinite(edges)):
        raise ValueError("bin_edges must be finite.")
    if not np.all(np.diff(edges) > 0.0):
        raise ValueError("bin_edges must be strictly increasing.")
    return edges


def compute_differential_rate(
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
        dR/dE_det  [counts / s / MeV],  same length as *E_det*.
    """
    E_det_arr = np.asarray(E_det, dtype=float)
    E_nu_min, E_nu_max = float(E_nu_range[0]), float(E_nu_range[1])
    if E_det_arr.ndim == 0:
        E_det_arr = E_det_arr.reshape(1)
    if E_nu_min >= E_nu_max:
        raise ValueError("E_nu_range must be strictly increasing.")

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
        )
        result[idx] = N_T * val

    return result


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
    """Backward-compatible alias for :func:`compute_differential_rate`."""
    return compute_differential_rate(
        E_det,
        flux,
        cross_section,
        quenching,
        resolution,
        nucleus,
        N_T,
        P,
        L,
        E_nu_range=E_nu_range,
    )


def compute_binned_spectrum(
    bin_edges: ArrayLike,
    flux: FluxBase,
    cross_section: CrossSectionBase,
    quenching: QuenchingBase,
    resolution: ResolutionBase,
    nucleus: Nucleus,
    N_T: float,
    P: float,
    L: float,
    E_nu_range: tuple[float, float] = (0.1, 10.0),
    exposure_s: float | None = None,
    samples_per_bin: int = 5,
) -> NDArray:
    """
    Compute detected-energy bin rates or counts.

    Parameters
    ----------
    bin_edges : array-like
        Detected-energy bin edges [MeV]. Bins are physical intervals
        ``[edge_i, edge_{i+1})``.
    exposure_s : float or None
        If provided, multiply rates by exposure time [s] and return counts/bin.
        If omitted, return counts/s/bin.
    samples_per_bin : int
        Number of detected-energy samples used for trapezoid integration in
        each bin. Must be at least 2.
    """
    edges = _validate_bin_edges(bin_edges)
    if samples_per_bin < 2:
        raise ValueError("samples_per_bin must be at least 2.")
    if exposure_s is not None and float(exposure_s) <= 0.0:
        raise ValueError("exposure_s must be positive.")

    result = np.zeros(edges.size - 1, dtype=float)
    for idx, (left, right) in enumerate(zip(edges[:-1], edges[1:])):
        E_det_grid = np.linspace(left, right, int(samples_per_bin))
        d_rate = compute_differential_rate(
            E_det_grid,
            flux,
            cross_section,
            quenching,
            resolution,
            nucleus,
            N_T,
            P,
            L,
            E_nu_range=E_nu_range,
        )
        result[idx] = np.trapezoid(d_rate, E_det_grid)

    if exposure_s is not None:
        result = result * float(exposure_s)
    return result


def _default_cevns_quenching() -> LindhardQuenching:
    return LindhardQuenching()


def _default_eves_quenching() -> ConstantQuenching:
    return ConstantQuenching(1.0)


def compute_cevns_spectrum(
    bin_edges: ArrayLike,
    flux: FluxBase,
    nucleus: Nucleus,
    P: float,
    L: float,
    *,
    N_T: float | None = None,
    target_mass_g: float | None = None,
    cross_section: CrossSectionBase | None = None,
    quenching: QuenchingBase | None = None,
    resolution: ResolutionBase,
    E_nu_range: tuple[float, float] = (0.1, 10.0),
    exposure_s: float | None = None,
    samples_per_bin: int = 5,
) -> NDArray:
    """Compute a binned SM CEvNS detected-energy spectrum."""
    target_count = _resolve_target_count(
        nucleus=nucleus,
        N_T=N_T,
        target_mass_g=target_mass_g,
    )
    xs = cross_section if cross_section is not None else SMCEvNS(nucleus)
    q = quenching if quenching is not None else _default_cevns_quenching()
    return compute_binned_spectrum(
        bin_edges,
        flux,
        xs,
        q,
        resolution,
        nucleus,
        target_count,
        P,
        L,
        E_nu_range=E_nu_range,
        exposure_s=exposure_s,
        samples_per_bin=samples_per_bin,
    )


def compute_eves_spectrum(
    bin_edges: ArrayLike,
    flux: FluxBase,
    nucleus: Nucleus,
    P: float,
    L: float,
    *,
    N_T: float | None = None,
    target_mass_g: float | None = None,
    cross_section: CrossSectionBase | None = None,
    quenching: QuenchingBase | None = None,
    resolution: ResolutionBase,
    E_nu_range: tuple[float, float] = (0.1, 10.0),
    exposure_s: float | None = None,
    samples_per_bin: int = 5,
) -> NDArray:
    """
    Compute a binned SM EvES detected-energy spectrum.

    ``N_T`` is the number of target atoms. Do not multiply it by ``Z``; the
    electron multiplicity is already included in ``SMEvES`` through ``Z_eff``.
    """
    target_count = _resolve_target_count(
        nucleus=nucleus,
        N_T=N_T,
        target_mass_g=target_mass_g,
    )
    xs = cross_section if cross_section is not None else SMEvES(nucleus)
    q = quenching if quenching is not None else _default_eves_quenching()
    return compute_binned_spectrum(
        bin_edges,
        flux,
        xs,
        q,
        resolution,
        nucleus,
        target_count,
        P,
        L,
        E_nu_range=E_nu_range,
        exposure_s=exposure_s,
        samples_per_bin=samples_per_bin,
    )


compute_cevns_spectrum.default_quenching = _default_cevns_quenching
compute_eves_spectrum.default_quenching = _default_eves_quenching
