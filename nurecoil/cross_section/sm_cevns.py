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
$Q_W = g_V^n N F + g_V^p Z F$ (same form factor for n and p).

Convention note
---------------
This uses the $G_F^2 M / \\pi$ prefactor convention, where the factor of 4
relative to the $G_F^2 M / 4\\pi$ convention is absorbed into $Q_W$.
The two conventions are equivalent:
    $Q_w^{(4\\pi)} = Z(1 - 4s_W^2) - N  \\Leftrightarrow  2 Q_W^{(\\pi)}$.

Flavor / radiative corrections
-------------------------------
``flavor="tree"`` (default): tree-level SM couplings
    $g_V^n = -0.5$,  $g_V^p = 0.5 - 2 s_W^2$

``flavor="nu_e" | "nu_mu" | "nu_tau"``: include radiative corrections
    $g_V^n = -0.5117$
    $g_V^p$: 0.0382 (nu_e), 0.0300 (nu_mu), 0.0256 (nu_tau)
    Sources: arXiv:2411.03122; Eur. Phys. J. C 83(7):683 (2023).

Unit conversion from natural units ($\\mathrm{MeV}^{-2}$) to $\\mathrm{cm}^2$:

.. math::

    \\sigma\\,[\\mathrm{cm}^2]
    = \\sigma\\,[\\mathrm{MeV}^{-2}]
      \\times (\\hbar c)^2\\,[\\mathrm{MeV}^2 \\cdot \\mathrm{fm}^2]
      \\times (10^{-13}\\,\\mathrm{cm/fm})^2

References
----------
Freedman 1974; Barranco et al. 2005; arXiv:2203.07361; arXiv:2411.03122.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from nurecoil import constants
from nurecoil.nucleus import Nucleus
from nurecoil.cross_section.base import CrossSectionBase
from nurecoil.cross_section.form_factors import FormFactorBase, HelmFormFactor

# Radiative-corrected proton vector couplings per neutrino flavor
# (arXiv:2411.03122; Eur. Phys. J. C 83(7):683 (2023))
_G_V_P_RC: dict[str, float] = {
    "nu_e":   0.0382,
    "nu_mu":  0.0300,
    "nu_tau": 0.0256,
}
_G_V_N_RC: float = -0.5117  # same for all flavors

FlavorRC = Literal["tree", "nu_e", "nu_mu", "nu_tau"]


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
    flavor : {"tree", "nu_e", "nu_mu", "nu_tau"}, optional
        Neutrino flavor used to select vector couplings.  ``"tree"`` (default)
        uses SM tree-level values; the flavor strings include radiative
        corrections from arXiv:2411.03122.
    """

    def __init__(
        self,
        nucleus: Nucleus,
        form_factor: FormFactorBase | None = None,
        flavor: FlavorRC = "tree",
    ) -> None:
        super().__init__(nucleus)
        self.form_factor: FormFactorBase = (
            form_factor if form_factor is not None else HelmFormFactor()
        )
        self.flavor = flavor

        if flavor == "tree":
            # SM tree-level: g_V^p = 1/2 - 2 sin²θ_W,  g_V^n = -1/2
            s2w = constants.sin2_theta_W
            self._g_V_p: float =  0.5 - 2.0 * s2w
            self._g_V_n: float = -0.5
        elif flavor in _G_V_P_RC:
            self._g_V_p = _G_V_P_RC[flavor]
            self._g_V_n = _G_V_N_RC
        else:
            raise ValueError(
                f"Unknown flavor {flavor!r}. Choose from 'tree', 'nu_e', 'nu_mu', 'nu_tau'."
            )

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
