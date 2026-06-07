"""
SM CEvNS differential cross section.

Units
-----
    E_nu, E_R, M  : MeV
    G_F           : MeV⁻²  (= constants.G_F)
    |q⃗|           : MeV  (natural units, ℏ = c = 1)
    dσ/dE_R       : cm² / MeV

The cross-section formula is

    dσ/dE_R = (G_F² M / π) · K(E_ν, E_R) · |Q_W(|q⃗|)|²

with

    K = 1 − E_R/E_ν + ½(E_R/E_ν)² − M E_R / (2 E_ν²)
    |q⃗| ≈ √(2 M E_R)   [MeV, natural units]
    Q_W = g_V^n N F_N(|q⃗|) + g_V^p Z F_Z(|q⃗|)

The output is in cm² / MeV.  The conversion from natural-unit MeV⁻² comes
from multiplying by (ℏc)² [MeV·fm]² / (ℏc / fm_per_cm)² — concretely:

    σ [cm²] = σ [MeV⁻²] × (ℏc [MeV·fm])² × (fm_per_cm [fm/cm])⁻²
             = σ [MeV⁻²] × hbar_c² × cm_per_fm²

References
----------
Freedman 1974; Drukier & Stodolsky 1984; Barranco et al. 2005.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from nurecoil import constants
from nurecoil.nucleus import Nucleus
from nurecoil.cross_section.base import CrossSectionBase
from nurecoil.cross_section.form_factors import FormFactorBase, HelmFormFactor


class SMCEvNS(CrossSectionBase):
    """
    Standard Model CEvNS differential cross section dσ/dE_R.

    Parameters
    ----------
    nucleus : Nucleus
        Target nucleus.
    form_factor : FormFactorBase, optional
        Nuclear form factor applied to both neutron and proton distributions.
        Defaults to ``HelmFormFactor()``.
    """

    def __init__(
        self,
        nucleus: Nucleus,
        form_factor: FormFactorBase | None = None,
    ) -> None:
        super().__init__(nucleus)
        self.form_factor: FormFactorBase = (
            form_factor if form_factor is not None else HelmFormFactor()
        )

        # Vector couplings (SM tree-level)
        # g_V^p = T₃ − 2 Q sin²θ_W = 1/2 − 2·(+1)·sin²θ_W
        # g_V^n = T₃ − 2 Q sin²θ_W = −1/2 − 2·(0)·sin²θ_W  = −1/2
        s2w = constants.sin2_theta_W
        self._g_V_p: float =  0.5 - 2.0 * s2w
        self._g_V_n: float = -0.5

        # Unit conversion factor: MeV⁻² → cm²
        # σ [cm²] = σ_nat [MeV⁻²] × (ℏc)² [MeV²·fm²] × (cm_per_fm)² [cm²/fm²]
        self._unit: float = constants.hbar_c ** 2 * constants.cm_per_fm ** 2

    def __call__(self, E_nu: ArrayLike, E_R: ArrayLike) -> NDArray:
        """
        Evaluate dσ/dE_R.

        Parameters
        ----------
        E_nu : array-like
            Neutrino energy [MeV].
        E_R : array-like
            Recoil energy [MeV].

        Returns
        -------
        NDArray
            dσ/dE_R  [cm² / MeV].  Returns 0 outside the valid kinematic region.
        """
        E_nu_arr = np.asarray(E_nu, dtype=float)
        E_R_arr  = np.asarray(E_R,  dtype=float)

        M = self.nucleus.M   # [MeV]
        N = self.nucleus.N
        Z = self.nucleus.Z

        # Kinematic factor K (dimensionless)
        r = E_R_arr / E_nu_arr
        K = 1.0 - r + 0.5 * r ** 2 - M * E_R_arr / (2.0 * E_nu_arr ** 2)

        # Momentum transfer |q⃗| [MeV]
        q = np.sqrt(np.maximum(2.0 * M * E_R_arr, 0.0))

        # Form factors (evaluated once; assumed same for n and p here)
        F = self.form_factor(q, self.nucleus)

        # Weak charge Q_W
        Q_W = self._g_V_n * N * F + self._g_V_p * Z * F

        # Prefactor G_F² M / π  [MeV⁻³]  (G_F in MeV⁻², M in MeV)
        prefactor = constants.G_F ** 2 * M / np.pi

        # dσ/dE_R in natural units [MeV⁻³], then convert to [cm²/MeV]
        dsigma = prefactor * K * Q_W ** 2 * self._unit  # [cm²/MeV]

        # Zero outside kinematic region
        E_R_max = self.E_R_max(E_nu_arr)
        kinematic_mask = (E_R_arr >= 0.0) & (E_R_arr <= E_R_max) & (K >= 0.0)
        return np.where(kinematic_mask, dsigma, 0.0)
