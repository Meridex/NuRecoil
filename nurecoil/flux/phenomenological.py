r"""
Phenomenological antineutrino flux from reactor fission.

The differential flux is:

    d\Phi/dE_\nu = 1/(4\pi L^2) * P_th/e_bar * \sum_i (f_i/F) * dN_i/dE_\nu

where e_bar = \sum_i (f_i/F) * e_i is the fission-fraction-weighted average
thermal energy per fission.

Units (API boundary)
--------------------
    E_nu : MeV
    P    : GW   (converted to MeV/s on entry via GW_to_MeV_per_s)
    L    : m    (converted to cm on entry)

Units (internal / output)
--------------------------
    d\Phi/dE_\nu : # / MeV / cm^2 / s
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from nurecoil import constants
from nurecoil.flux.base import FluxBase

# ---------------------------------------------------------------------------
# Neutrino spectrum coefficients  dN_k/dE_nu = exp(sum_p alpha_p * E_nu^(p-1))
# ---------------------------------------------------------------------------

# Mueller et al. 2011 [arXiv:1101.2663, Table VI]  —  valid 2–8 MeV
# All four isotopes available.
SPECTRUM_MUELLER_2011: dict[str, list[float]] = {
    "U235":  [ 3.217, -3.111,  1.395, -3.690e-1,  4.445e-2, -2.053e-3],
    "U238":  [ 4.833e-1,  1.927e-1, -1.283e-1, -6.762e-3,  2.233e-3, -1.536e-4],
    "Pu239": [ 6.413, -7.432,  3.535, -8.820e-1,  1.025e-1, -4.550e-3],
    "Pu241": [ 3.251, -3.204,  1.428, -3.675e-1,  4.254e-2, -1.896e-3],
}

# Schreckenbach et al. 1985 [Phys.Rev.D 39 (1989) 3378, Table I]  —  valid 2–8 MeV
# Simple quadratic exponent: dN/dE = exp(a0 + a1*E + a2*E^2)
# Note: overestimates by factor 2–3 for E > 8 MeV.
SPECTRUM_SCHRECKENBACH_1985: dict[str, list[float]] = {
    "U235":  [0.870, -0.160,  -0.0910],
    "U238":  [0.976, -0.162,  -0.0790],
    "Pu239": [0.896, -0.239,  -0.0981],
    "Pu241": [0.793, -0.080,  -0.1085],
}

# Huber 2011 [arXiv:1106.0687, Table III]  —  valid 2–8 MeV
# U238 is NOT available in this reference; cannot be used alone for a mixed spectrum.
SPECTRUM_HUBER_2011: dict[str, list[float]] = {
    "U235":  [ 4.367, -4.577,  2.100, -5.294e-1,  6.186e-2, -2.777e-3],
    "Pu239": [ 4.757, -5.392,  2.563, -6.596e-1,  7.820e-2, -3.536e-3],
    "Pu241": [ 2.990, -2.882,  1.278, -3.343e-1,  3.905e-2, -1.754e-3],
    # U238 not available in this reference.
}

# Combined "Huber-Mueller": Huber 2011 for U235/Pu239/Pu241, Mueller 2011 for U238.
# This is the most commonly used combination in the reactor neutrino community.
SPECTRUM_HUBER_MUELLER: dict[str, list[float]] = {
    **{k: SPECTRUM_HUBER_2011[k] for k in ("U235", "Pu239", "Pu241")},
    "U238": SPECTRUM_MUELLER_2011["U238"],
}

# ---------------------------------------------------------------------------
# Effective thermal energy per fission e_i  [MeV / fission]
# ---------------------------------------------------------------------------

# Meulenberg 1969  [Journal of Nuclear Energy 23 (1969) 517, Table A.2]
ENERGY_PER_FISSION_MEULENBERG_1969: dict[str, float] = {
    "U235":  201.7,
    "U238":  205.0,
    "Pu239": 210.0,
    "Pu241": 212.4,
}

# Bemporad et al. 2002  [arXiv:hep-ph/0410100, Table 4]
ENERGY_PER_FISSION_BEMPORAD_2002: dict[str, float] = {
    "U235":  201.92,
    "U238":  205.52,
    "Pu239": 209.99,
    "Pu241": 213.60,
}

# Ma et al. 2013  [Phys. Rev. C 88 (2013) 014605, Table IX]
ENERGY_PER_FISSION_MA_2013: dict[str, float] = {
    "U235":  202.36,
    "U238":  205.99,
    "Pu239": 211.12,
    "Pu241": 214.26,
}

# ---------------------------------------------------------------------------
# Fission fractions  f_i / F
# ---------------------------------------------------------------------------

# Kuo-Sheng Nuclear Power Station  [arXiv:2402.06416]
FISSION_FRACTIONS_KSNPS: dict[str, float] = {
    "U235":  0.55,
    "U238":  0.07,
    "Pu239": 0.32,
    "Pu241": 0.06,
}

# CONUS  [arXiv:2401.07684]
FISSION_FRACTIONS_CONUS: dict[str, float] = {
    "U235":  0.491,
    "U238":  0.074,
    "Pu239": 0.361,
    "Pu241": 0.074,
}

# Typical commercial reactor (PWR equilibrium)  [arXiv:2310.113070]
FISSION_FRACTIONS_TYPICAL: dict[str, float] = {
    "U235":  0.58,
    "U238":  0.08,
    "Pu239": 0.29,
    "Pu241": 0.05,
}

# Daya Bay  [arXiv:2210.01068]
FISSION_FRACTIONS_DAYA_BAY: dict[str, float] = {
    "U235":  0.564,
    "U238":  0.076,
    "Pu239": 0.304,
    "Pu241": 0.056,
}

# ---------------------------------------------------------------------------
# Convenience name maps exposed for user selection
# ---------------------------------------------------------------------------
SPECTRUM_MODELS: dict[str, dict[str, list[float]]] = {
    "huber_mueller":      SPECTRUM_HUBER_MUELLER,
    "mueller_2011":       SPECTRUM_MUELLER_2011,
    "huber_2011":         SPECTRUM_HUBER_2011,
    "schreckenbach_1985": SPECTRUM_SCHRECKENBACH_1985,
}

ENERGY_PER_FISSION_MODELS: dict[str, dict[str, float]] = {
    "ma_2013":         ENERGY_PER_FISSION_MA_2013,
    "bemporad_2002":   ENERGY_PER_FISSION_BEMPORAD_2002,
    "meulenberg_1969": ENERGY_PER_FISSION_MEULENBERG_1969,
}

FISSION_FRACTION_PRESETS: dict[str, dict[str, float]] = {
    "typical":  FISSION_FRACTIONS_TYPICAL,
    "ksnps":    FISSION_FRACTIONS_KSNPS,
    "conus":    FISSION_FRACTIONS_CONUS,
    "daya_bay": FISSION_FRACTIONS_DAYA_BAY,
}

# Valid energy range (MeV) for each named spectrum model
_SPECTRUM_ENERGY_RANGE: dict[str, tuple[float, float]] = {
    "huber_mueller":      (2.0, 8.0),
    "mueller_2011":       (2.0, 8.0),
    "huber_2011":         (2.0, 8.0),
    "schreckenbach_1985": (2.0, 8.0),
}


def _eval_spectrum(coeffs: list[float], E_nu: NDArray) -> NDArray:
    """
    Evaluate exp-polynomial spectrum for a single isotope.

    S(E) = exp(a0 + a1*E + a2*E^2 + ...)   (works for both 3- and 6-term fits)
    """
    exponent = np.zeros_like(E_nu)
    for p, a in enumerate(coeffs):
        exponent = exponent + a * E_nu ** p
    return np.exp(exponent)


class PhenomenologicalFlux(FluxBase):
    r"""
    Phenomenological reactor antineutrino flux.

    Computes:

        d\Phi/dE_\nu = (P_th / e_bar) * \sum_i (f_i/F) * S_i(E_\nu)  /  (4\pi L^2)

    where e_bar = \sum_i (f_i/F) * e_i is the fission-fraction-weighted average thermal
    energy per fission.

    Parameters
    ----------
    fission_fractions : dict[str, float] | str | None
        Per-isotope fission fractions, or a preset name:
        ``"typical"`` (default), ``"ksnps"``, ``"conus"``, ``"daya_bay"``.
        Must sum to 1.
    spectrum_model : dict[str, list[float]] | str
        Polynomial coefficients per isotope, or a model name:
        ``"huber_mueller"`` (default), ``"mueller_2011"``,
        ``"huber_2011"``, ``"schreckenbach_1985"``.
    energy_per_fission : dict[str, float] | str
        Effective thermal energy per fission [MeV] per isotope, or a model name:
        ``"ma_2013"`` (default), ``"bemporad_2002"``, ``"meulenberg_1969"``.

    Notes
    -----
    ``"huber_2011"`` does not include U238.  Passing a non-zero U238 fission
    fraction with this model raises a ``ValueError``.
    """

    # Overwritten in __init__ based on the chosen spectrum model.
    E_min: float = 2.0
    E_max: float = 8.0

    def __init__(
        self,
        fission_fractions: dict[str, float] | str | None = None,
        spectrum_model: dict[str, list[float]] | str = "huber_mueller",
        energy_per_fission: dict[str, float] | str = "ma_2013",
    ) -> None:
        # --- resolve fission fractions ---
        if fission_fractions is None:
            fission_fractions = "typical"
        if isinstance(fission_fractions, str):
            if fission_fractions not in FISSION_FRACTION_PRESETS:
                raise ValueError(
                    f"Unknown fission fraction preset '{fission_fractions}'. "
                    f"Choose from: {list(FISSION_FRACTION_PRESETS)}"
                )
            fission_fractions = dict(FISSION_FRACTION_PRESETS[fission_fractions])
        total = sum(fission_fractions.values())
        if not np.isclose(total, 1.0, atol=1e-3):
            raise ValueError(
                f"Fission fractions must sum to 1, got {total:.4f}."
            )

        # --- resolve spectrum model ---
        if isinstance(spectrum_model, str):
            model_name = spectrum_model
            if model_name not in SPECTRUM_MODELS:
                raise ValueError(
                    f"Unknown spectrum model '{model_name}'. "
                    f"Choose from: {list(SPECTRUM_MODELS)}"
                )
            spectrum_coeffs = SPECTRUM_MODELS[model_name]
            self.E_min, self.E_max = _SPECTRUM_ENERGY_RANGE[model_name]
        else:
            spectrum_coeffs = spectrum_model

        # Guard: huber_2011 lacks U238
        if "U238" not in spectrum_coeffs and fission_fractions.get("U238", 0.0) > 0:
            raise ValueError(
                "The chosen spectrum model does not include U238, but the "
                f"U238 fission fraction is {fission_fractions['U238']:.3f}. "
                "Use 'huber_mueller' or 'mueller_2011', or set U238 fraction to 0."
            )

        # --- resolve energy per fission ---
        if isinstance(energy_per_fission, str):
            if energy_per_fission not in ENERGY_PER_FISSION_MODELS:
                raise ValueError(
                    f"Unknown energy-per-fission model '{energy_per_fission}'. "
                    f"Choose from: {list(ENERGY_PER_FISSION_MODELS)}"
                )
            energy_per_fission = ENERGY_PER_FISSION_MODELS[energy_per_fission]

        self._fractions       = fission_fractions
        self._spectrum_coeffs = spectrum_coeffs
        self._e_per_fission   = energy_per_fission

        # e_bar = \sum_i (f_i/F) * e_i  [MeV / fission]
        self._e_bar: float = sum(
            frac * self._e_per_fission[iso]
            for iso, frac in self._fractions.items()
            if iso in self._e_per_fission
        )

    def _spectrum_per_fission(self, E_nu: NDArray) -> NDArray:
        r"""Weighted spectrum \sum_i (f_i/F) S_i(E_\nu)  [# / MeV / fission]."""
        spectrum = np.zeros_like(E_nu)
        for iso, frac in self._fractions.items():
            if iso not in self._spectrum_coeffs:
                continue
            spectrum += frac * _eval_spectrum(self._spectrum_coeffs[iso], E_nu)
        return spectrum

    def _flux(self, E_nu: NDArray, P: float, L: float) -> NDArray:
        r"""
        Evaluate d\Phi/dE_\nu [# / MeV / cm^2 / s] for in-range energies.

        Parameters
        ----------
        E_nu : NDArray
            Neutrino energies [MeV], guaranteed within [E_min, E_max].
        P : float
            Reactor thermal power [GW].
        L : float
            Baseline distance [m].
        """
        L_cm = L * constants.cm_per_m

        # Fission rate [fissions/s] = P [MeV/s] / ē [MeV/fission]
        fission_rate = P * constants.GW_to_MeV_per_s / self._e_bar

        spectrum = self._spectrum_per_fission(E_nu)  # # / MeV / fission

        return fission_rate * spectrum / (4.0 * np.pi * L_cm ** 2)
