"""
SM CEvNS differential cross section.

Units
-----
    E_nu, E_R, M  : MeV
    G_F           : $\\mathrm{MeV}^{-2}$  (= constants.G_F)
    |q|           : MeV  (natural units, $\\hbar = c = 1$)
    $d\\sigma/dE_R$  : cm^2 / MeV

The cross-section formula is

.. math::

    \\frac{d\\sigma}{dE_R} = \\frac{G_F^2\\,M}{\\pi}
    \\left(1 - \\frac{E_R}{E_\\nu}
         + \\frac{1}{2}\\!\\left(\\frac{E_R}{E_\\nu}\\right)^{\\!2}
         - \\frac{M E_R}{2 E_\\nu^2}\\right)
    \\left|Q_W\\!\\left(|\\vec{q}|\\right)\\right|^2

where $|\\vec{q}| \\approx \\sqrt{2 M E_R}$ [MeV] and
$Q_W = g_V^n N F_N + g_V^p Z F_Z$.

Unit conversion from natural units ($\\mathrm{MeV}^{-2}$) to $\\mathrm{cm}^2$:

.. math::

    \\sigma\\,[\\mathrm{cm}^2]
    = \\sigma\\,[\\mathrm{MeV}^{-2}]
      \\times (\\hbar c)^2\\,[\\mathrm{MeV}^2 \\cdot \\mathrm{fm}^2]
      \\times (10^{-13}\\,\\mathrm{cm/fm})^2

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
    Standard Model CEvNS differential cross section $d\\sigma/dE_R$.

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
        # $g_V^p = T_3 - 2Q\sin^2\theta_W = +1/2 - 2(+1)\sin^2\theta_W$
        # $g_V^n = T_3 - 2Q\sin^2\theta_W = -1/2 - 2(0)\sin^2\theta_W = -1/2$
        s2w = constants.sin2_theta_W
        self._g_V_p: float =  0.5 - 2.0 * s2w
        self._g_V_n: float = -0.5

        # Unit conversion: MeV^{-2} -> cm^2
        # $\sigma\ [\mathrm{cm}^2] = \sigma\ [\mathrm{MeV}^{-2}] \times (\hbar c)^2 \times (\mathrm{cm/fm})^2$
        self._unit: float = constants.hbar_c ** 2 * constants.cm_per_fm ** 2

    def __call__(self, E_nu: ArrayLike, E_R: ArrayLike) -> NDArray:
        r"""
        Evaluate $d\sigma/dE_R$.

        Parameters
        ----------
        E_nu : array-like
            Neutrino energy [MeV].
        E_R : array-like
            Recoil energy [MeV].

        Returns
        -------
        NDArray
            $d\sigma/dE_R$  [cm^2 / MeV].  Returns 0 outside the valid kinematic region.
        """
        E_nu_arr = np.asarray(E_nu, dtype=float)
        E_R_arr  = np.asarray(E_R,  dtype=float)

        M = self.nucleus.M   # [MeV]
        N = self.nucleus.N
        Z = self.nucleus.Z

        # Kinematic factor K (dimensionless)
        r = E_R_arr / E_nu_arr
        K = 1.0 - r + 0.5 * r ** 2 - M * E_R_arr / (2.0 * E_nu_arr ** 2)

        # Momentum transfer $|\vec{q}| \approx \sqrt{2 M E_R}$ [MeV]
        q = np.sqrt(np.maximum(2.0 * M * E_R_arr, 0.0))

        # Form factors (evaluated once; assumed same for n and p here)
        F = self.form_factor(q, self.nucleus)

        # Weak charge: $Q_W = g_V^n N F + g_V^p Z F$
        Q_W = self._g_V_n * N * F + self._g_V_p * Z * F

        # Prefactor $G_F^2 M / \pi$ [MeV^{-3}]  (G_F in MeV^{-2}, M in MeV)
        prefactor = constants.G_F ** 2 * M / np.pi

        # d sigma/dE_R in natural units [MeV^{-3}], then convert to [cm^2/MeV]
        dsigma = prefactor * K * Q_W ** 2 * self._unit  # [cm^2/MeV]

        # Zero outside kinematic region
        E_R_max = self.E_R_max(E_nu_arr)
        kinematic_mask = (E_R_arr >= 0.0) & (E_R_arr <= E_R_max) & (K >= 0.0)
        return np.where(kinematic_mask, dsigma, 0.0)
