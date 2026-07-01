r"""
Reactor fission mixture: fission fractions and thermal energy per fission.

This module centralises all tabulated data and the bookkeeping logic that is
shared by every flux model that combines per-isotope spectra into a total
antineutrino spectrum:

    S_total(E) = \sum_i (f_i / F) * S_i(E)
    e_bar      = \sum_i (f_i / F) * e_i

Both :mod:`nurecoil.flux.phenomenological` and
:mod:`nurecoil.flux.interpolated` delegate to :class:`ReactorMix` for this
bookkeeping.

Supported isotopes
------------------
``"U235"``, ``"U238"``, ``"Pu239"``, ``"Pu241"``
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Effective thermal energy per fission  e_i  [MeV / fission]
# ---------------------------------------------------------------------------

# James 1969  [Journal of Nuclear Energy 23 (1969) 517, Table A.2]
ENERGY_PER_FISSION_JAMES_1969: dict[str, float] = {
    "U235":  201.7,
    "U238":  205.0,
    "Pu239": 210.0,
    "Pu241": 212.4,
}

# Kopeikin et al. 2004  [arXiv:hep-ph/0410100, Table 4]
ENERGY_PER_FISSION_KOPEIKIN_2004: dict[str, float] = {
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

ENERGY_PER_FISSION_MODELS: dict[str, dict[str, float]] = {
    "ma_2013":         ENERGY_PER_FISSION_MA_2013,
    "kopeikin_2004":   ENERGY_PER_FISSION_KOPEIKIN_2004,
    "james_1969":      ENERGY_PER_FISSION_JAMES_1969,
}

FISSION_FRACTION_PRESETS: dict[str, dict[str, float]] = {
    "typical":  FISSION_FRACTIONS_TYPICAL,
    "ksnps":    FISSION_FRACTIONS_KSNPS,
    "conus":    FISSION_FRACTIONS_CONUS,
    "daya_bay": FISSION_FRACTIONS_DAYA_BAY,
}

# ---------------------------------------------------------------------------
# ReactorMix
# ---------------------------------------------------------------------------

class ReactorMix:
    """
    Resolved fission-fraction mixture for a reactor core.

    Validates inputs, resolves preset names, and pre-computes the
    fission-fraction-weighted average thermal energy per fission:

        e_bar = sum_i (f_i / F) * e_i   [MeV / fission]

    Parameters
    ----------
    fission_fractions : dict[str, float] | str | None
        Per-isotope fission fractions summing to 1, or a preset name:
        ``"typical"`` (default), ``"ksnps"``, ``"conus"``, ``"daya_bay"``.
    energy_per_fission : dict[str, float] | str
        Effective thermal energy per fission [MeV] per isotope, or a model
        name: ``"ma_2013"`` (default), ``"kopeikin_2004"``,
        ``"james_1969"``.
    """

    def __init__(
        self,
        fission_fractions: dict[str, float] | str | None = None,
        energy_per_fission: dict[str, float] | str = "ma_2013",
    ) -> None:
        # --- resolve fission fractions ---
        if fission_fractions is None:
            fission_fractions = "typical"
        if isinstance(fission_fractions, str):
            # Support "fractions:<name>" prefix as well as bare preset name.
            key = fission_fractions.removeprefix("fractions:")
            if key not in FISSION_FRACTION_PRESETS:
                raise ValueError(
                    f"Unknown fission fraction preset '{key}'. "
                    f"Choose from: {list(FISSION_FRACTION_PRESETS)}. "
                    f"You can also use the 'fractions:<name>' prefix form."
                )
            fission_fractions = dict(FISSION_FRACTION_PRESETS[key])
        total = sum(fission_fractions.values())
        if not np.isclose(total, 1.0, atol=1e-3):
            raise ValueError(
                f"Fission fractions must sum to 1, got {total:.4f}."
            )
        self.fractions: dict[str, float] = fission_fractions

        # --- resolve energy per fission ---
        if isinstance(energy_per_fission, str):
            # Support "energy:<name>" prefix as well as bare model name.
            key = energy_per_fission.removeprefix("energy:")
            if key not in ENERGY_PER_FISSION_MODELS:
                raise ValueError(
                    f"Unknown energy-per-fission model '{key}'. "
                    f"Choose from: {list(ENERGY_PER_FISSION_MODELS)}. "
                    f"You can also use the 'energy:<name>' prefix form."
                )
            energy_per_fission = dict(ENERGY_PER_FISSION_MODELS[key])
        self.energy_per_fission: dict[str, float] = energy_per_fission

        # e_bar = sum_i (f_i/F) * e_i  [MeV / fission]
        self.e_bar: float = sum(
            frac * self.energy_per_fission[iso]
            for iso, frac in self.fractions.items()
            if iso in self.energy_per_fission
        )
        if self.e_bar <= 0:
            raise ValueError(
                "Weighted average energy per fission is zero or negative. "
                "Check that fission fractions and energy_per_fission share "
                "at least one common isotope."
            )
