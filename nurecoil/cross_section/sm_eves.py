"""
SM EvES (Elastic neutrino-Electron Scattering) differential cross section.

Units
-----
    E_nu, E_R, m_e : MeV
    G_F            : MeV^{-2}  (= constants.G_F)
    d sigma/dE_R   : cm^2 / MeV

The cross-section formula is

.. math::

    \\frac{d\\sigma^{\\mathrm{EvES}}}{dT} =
        Z_{\\mathrm{eff}}^{\\mathcal{A}}(T)\\,
        \\frac{G_F^2 m_e}{2\\pi}
        \\Bigl[
            (g_V + g_A)^2
            + (g_V - g_A)^2 \\left(1 - \\frac{T}{E_\\nu}\\right)^{\\!2}
            - (g_V^2 - g_A^2)\\,\\frac{m_e T}{E_\\nu^2}
        \\Bigr]

where $T = E_R$ is the electron recoil energy and $Z_{\\mathrm{eff}}(T)$ is the
effective number of ionisable electrons at recoil energy $T$ (step function
from atomic binding energies).

Kinematic upper bound: $T_{\\max} = 2 E_\\nu^2 / (m_e + 2 E_\\nu)$.

Couplings (including radiative corrections, $s_W^2 = 0.23857$)
--------------------------------------------------------------
Source: arXiv:2207.05036; arXiv:2411.03122.

+----------+---------+---------+
| flavor   |  g_V    |  g_A    |
+==========+=========+=========+
| nu_e     | +0.9521 | +0.4938 |
+----------+---------+---------+
| nu_mu    | -0.0397 | -0.5062 |
+----------+---------+---------+
| nu_tau   | -0.0353 | -0.5062 |
+----------+---------+---------+

For anti-neutrinos the axial coupling changes sign: $g_A \\to -g_A$.

Z_eff step functions
--------------------
Implemented for Ge (Z=32), Cs (Z=55), I (Z=53).
Data from JHEP 09 (2022) 164; binding energies from https://xdb.lbl.gov.
Other elements raise ``NotImplementedError``.

References
----------
arXiv:1907.03379; arXiv:2207.05036; arXiv:2411.03122; JHEP 09 (2022) 164.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from nurecoil import constants
from nurecoil.nucleus import Nucleus
from nurecoil.cross_section.base import CrossSectionBase

# ---------------------------------------------------------------------------
# Coupling constants (with radiative corrections)
# ---------------------------------------------------------------------------
# Source: arXiv:2207.05036 (s_W^2 = 0.23857 at q^2 -> 0)
_COUPLINGS: dict[str, tuple[float, float]] = {
    #          g_V      g_A
    "nu_e":   ( 0.9521,  0.4938),
    "nu_mu":  (-0.0397, -0.5062),
    "nu_tau": (-0.0353, -0.5062),
}

FlavorEvES = Literal["nu_e", "nu_mu", "nu_tau"]

# ---------------------------------------------------------------------------
# Z_eff step-function data
# ---------------------------------------------------------------------------
# Each entry: list of (threshold_MeV, Z_eff_above_threshold), ordered from
# highest threshold to lowest.  Below the lowest threshold the last Z_eff
# applies.  Thresholds taken from electron binding energies at xdb.lbl.gov;
# compiled in JHEP 09 (2022) 164.

# keV -> MeV conversion applied inline
_keV = 1e-3  # 1 keV in MeV

# Ge (Z=32): thresholds in MeV, Z_eff when T > threshold
_Z_EFF_GE: list[tuple[float, int]] = [
    (11.103 * _keV, 32),
    ( 1.4146 * _keV, 30),
    ( 1.2481 * _keV, 28),
    ( 1.217  * _keV, 26),
    ( 0.1801 * _keV, 22),
    ( 0.1249 * _keV, 20),
    ( 0.1208 * _keV, 18),
    ( 0.0298 * _keV, 14),
    ( 0.0292 * _keV, 10),
]
_Z_EFF_GE_BELOW: int = 4  # T < 0.0292 keV

# Cs (Z=55)
_Z_EFF_CS: list[tuple[float, int]] = [
    (35.99  * _keV, 55),
    ( 5.71  * _keV, 53),
    ( 5.36  * _keV, 51),
    ( 5.01  * _keV, 49),
    ( 1.21  * _keV, 45),
    ( 1.07  * _keV, 43),
    ( 1.00  * _keV, 41),
    ( 0.74  * _keV, 37),
    ( 0.73  * _keV, 33),
    ( 0.23  * _keV, 27),
    ( 0.17  * _keV, 25),
    ( 0.16  * _keV, 23),
]
_Z_EFF_CS_BELOW: int = 19  # T < 0.16 keV

# I (Z=53)
_Z_EFF_I: list[tuple[float, int]] = [
    (33.17  * _keV, 53),
    ( 5.19  * _keV, 51),
    ( 4.86  * _keV, 49),
    ( 4.56  * _keV, 47),
    ( 1.07  * _keV, 43),
    ( 0.93  * _keV, 41),
    ( 0.88  * _keV, 39),
    ( 0.63  * _keV, 35),
    ( 0.62  * _keV, 31),
    ( 0.19  * _keV, 25),
    ( 0.124 * _keV, 23),
    ( 0.123 * _keV, 21),
]
_Z_EFF_I_BELOW: int = 17  # T < 0.123 keV

_Z_EFF_TABLE: dict[int, tuple[list[tuple[float, int]], int]] = {
    32: (_Z_EFF_GE, _Z_EFF_GE_BELOW),
    55: (_Z_EFF_CS, _Z_EFF_CS_BELOW),
    53: (_Z_EFF_I,  _Z_EFF_I_BELOW),
}


def _z_eff_scalar(T: float, steps: list[tuple[float, int]], below: int) -> int:
    """Return Z_eff for a single recoil energy T [MeV]."""
    for threshold, z in steps:
        if T > threshold:
            return z
    return below


def _z_eff_array(
    T_arr: NDArray,
    steps: list[tuple[float, int]],
    below: int,
) -> NDArray:
    """Vectorised Z_eff over an array of recoil energies [MeV]."""
    result = np.full(T_arr.shape, below, dtype=float)
    # Apply thresholds from lowest to highest so highest wins
    for threshold, z in reversed(steps):
        result = np.where(T_arr > threshold, z, result)
    return result


# ---------------------------------------------------------------------------
# SMEvES class
# ---------------------------------------------------------------------------

class SMEvES(CrossSectionBase):
    r"""
    Standard Model EvES (Elastic neutrino-Electron Scattering)
    differential cross section $d\sigma/dE_R$.

    Parameters
    ----------
    nucleus : Nucleus
        Target atom.  The atomic number ``nucleus.Z`` determines which
        $Z_{\mathrm{eff}}(T)$ step function is used.  Currently supported:
        Ge (Z=32), Cs (Z=55), I (Z=53).  Other elements raise
        ``NotImplementedError``.
    flavor : {"nu_e", "nu_mu", "nu_tau"}
        Neutrino flavor, selects radiatively-corrected couplings.
        Default ``"nu_e"``.
    antineutrino : bool, optional
        If ``True``, flip the sign of the axial coupling ($g_A \to -g_A$).
        Default ``False``.
    """

    def __init__(
        self,
        nucleus: Nucleus,
        flavor: FlavorEvES = "nu_e",
        antineutrino: bool = False,
    ) -> None:
        super().__init__(nucleus)

        if nucleus.Z not in _Z_EFF_TABLE:
            raise NotImplementedError(
                f"Z_eff is not implemented for Z={nucleus.Z}. "
                f"Supported: Ge (Z=32), Cs (Z=55), I (Z=53)."
            )
        if flavor not in _COUPLINGS:
            raise ValueError(
                f"Unknown flavor {flavor!r}. Choose from 'nu_e', 'nu_mu', 'nu_tau'."
            )

        self.flavor = flavor
        self.antineutrino = antineutrino

        g_V, g_A = _COUPLINGS[flavor]
        self._g_V: float = g_V
        self._g_A: float = -g_A if antineutrino else g_A

        self._z_steps, self._z_below = _Z_EFF_TABLE[nucleus.Z]

        # Unit conversion: MeV^{-2} -> cm^2
        self._unit: float = constants.hbar_c ** 2 * constants.cm_per_fm ** 2

    def E_R_max(self, E_nu: float) -> float:  # type: ignore[override]
        r"""
        Kinematic upper bound on the electron recoil energy.

        .. math::
            T_{\max} = \frac{2 E_\nu^2}{m_e + 2 E_\nu}

        Parameters
        ----------
        E_nu : float
            Neutrino energy [MeV].

        Returns
        -------
        float
            Maximum recoil energy [MeV].
        """
        m_e = constants.m_e
        return 2.0 * E_nu ** 2 / (m_e + 2.0 * E_nu)

    def __call__(self, E_nu: ArrayLike, E_R: ArrayLike) -> NDArray:
        r"""
        Evaluate $d\sigma/dE_R$.

        Parameters
        ----------
        E_nu : array-like
            Neutrino energy [MeV].
        E_R : array-like
            Electron recoil energy [MeV].

        Returns
        -------
        NDArray
            $d\sigma/dE_R$  [cm^2 / MeV].  Returns 0 outside the valid
            kinematic region.
        """
        E_nu_arr = np.asarray(E_nu, dtype=float)
        E_R_arr  = np.asarray(E_R,  dtype=float)

        m_e  = constants.m_e
        g_V  = self._g_V
        g_A  = self._g_A

        # Kinematic upper bound
        E_R_max = 2.0 * E_nu_arr ** 2 / (m_e + 2.0 * E_nu_arr)

        # Bracketed factor
        t = E_R_arr / E_nu_arr
        bracket = (
            (g_V + g_A) ** 2
            + (g_V - g_A) ** 2 * (1.0 - t) ** 2
            - (g_V ** 2 - g_A ** 2) * m_e * E_R_arr / E_nu_arr ** 2
        )

        # Prefactor G_F^2 m_e / (2 pi) [MeV^{-3}]
        prefactor = constants.G_F ** 2 * m_e / (2.0 * np.pi)

        # Z_eff step function
        Z_eff = _z_eff_array(E_R_arr, self._z_steps, self._z_below)

        # d sigma / dE_R  [cm^2 / MeV]
        dsigma = Z_eff * prefactor * bracket * self._unit

        # Zero outside kinematic region
        kinematic_mask = (E_R_arr >= 0.0) & (E_R_arr <= E_R_max)
        return np.where(kinematic_mask, dsigma, 0.0)
