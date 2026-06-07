"""
Nucleus dataclass and built-in isotope table.

Units
-----
All masses are in MeV (internal canonical unit).
"""

from __future__ import annotations

from dataclasses import dataclass

from nurecoil import constants


@dataclass(frozen=True)
class Nucleus:
    """
    A nuclear target characterised by atomic number Z and mass number A.

    Parameters
    ----------
    Z : int
        Atomic number (number of protons).
    A : int
        Mass number (protons + neutrons).

    Properties
    ----------
    N   : neutron number  (A - Z)
    M   : nuclear mass [MeV]  ≈ A × u_to_MeV
    """

    Z: int
    A: int

    def __post_init__(self) -> None:
        if self.Z < 1 or self.A < 1 or self.A < self.Z:
            raise ValueError(f"Invalid nucleus: Z={self.Z}, A={self.A}")

    @property
    def N(self) -> int:
        """Neutron number."""
        return self.A - self.Z

    @property
    def M(self) -> float:
        """Nuclear mass in MeV  (approximation $M \\approx A \\cdot u$, < 0.1 % for A > 10)."""
        return self.A * constants.u_to_MeV

    def __repr__(self) -> str:
        return f"Nucleus(Z={self.Z}, A={self.A}, N={self.N}, M={self.M:.1f} MeV)"


# ---------------------------------------------------------------------------
# Built-in isotope table
# Keys follow the convention  "<symbol><A>", e.g. "Ge76", "Si28".
# ---------------------------------------------------------------------------
ISOTOPE_TABLE: dict[str, Nucleus] = {
    # --- Germanium (CDEX, TEXONO targets) ---
    "Ge70": Nucleus(Z=32, A=70),
    "Ge72": Nucleus(Z=32, A=72),
    "Ge73": Nucleus(Z=32, A=73),
    "Ge74": Nucleus(Z=32, A=74),
    "Ge76": Nucleus(Z=32, A=76),
    # --- Silicon ---
    "Si28": Nucleus(Z=14, A=28),
    "Si29": Nucleus(Z=14, A=29),
    "Si30": Nucleus(Z=14, A=30),
    # --- Caesium and Iodine (CsI[Na] detectors) ---
    "Cs133": Nucleus(Z=55, A=133),
    "I127":  Nucleus(Z=53, A=127),
    # --- Argon (liquid-argon detectors) ---
    "Ar40": Nucleus(Z=18, A=40),
    # --- Xenon (LXe detectors) ---
    "Xe132": Nucleus(Z=54, A=132),
    "Xe134": Nucleus(Z=54, A=134),
    "Xe136": Nucleus(Z=54, A=136),
}


def get_nucleus(name: str) -> Nucleus:
    """
    Look up a nucleus by name from the built-in isotope table.

    Parameters
    ----------
    name : str
        Key in ``ISOTOPE_TABLE``, e.g. ``"Ge76"``.

    Returns
    -------
    Nucleus

    Raises
    ------
    KeyError
        If *name* is not found.  Available keys are listed in
        ``ISOTOPE_TABLE``.
    """
    try:
        return ISOTOPE_TABLE[name]
    except KeyError:
        available = ", ".join(sorted(ISOTOPE_TABLE))
        raise KeyError(
            f"Unknown nucleus {name!r}. Available: {available}"
        ) from None
