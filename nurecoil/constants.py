"""
Physical constants and unit-conversion factors for NuRecoil.

Design rules
------------
* All *internal* calculations use the canonical units defined in
  ``misc/plans/units_convention.md``:
    - Energy          : MeV
    - Nuclear length  : fm
    - Geometric length: cm
    - Cross section   : cm²
    - Momentum transfer: MeV  (ℏ = c = 1)
    - Nuclear mass    : MeV
    - Flux            : # / MeV / cm²

* Conversions happen **only at the API boundary** (user input → internal,
  internal → user output).  Every conversion must reference a named constant
  below — no bare numeric literals elsewhere in the package.

PDG values: 2024 edition (https://pdg.lbl.gov).
CODATA values: 2018 CODATA recommended values.
"""

# ---------------------------------------------------------------------------
# Fermi constant
# ---------------------------------------------------------------------------
# PDG: G_F / (ℏc)³ = 1.1663788(6) × 10⁻⁵ GeV⁻²
G_F_GeV2: float = 1.1663788e-5   # GeV⁻²  — PDG reference value; do NOT use in formulas
G_F: float      = 1.1663788e-11  # MeV⁻²  — use this in all cross-section formulas
                                  # (= G_F_GeV2 × 10⁻⁶, since 1 GeV⁻² = 10⁻⁶ MeV⁻²)

# ---------------------------------------------------------------------------
# Electroweak mixing
# ---------------------------------------------------------------------------
sin2_theta_W: float = 0.23122  # sin²θ_W, MS-bar scheme at M_Z (PDG 2024)

# ---------------------------------------------------------------------------
# Natural unit bridge
# ---------------------------------------------------------------------------
hbar_c: float = 197.3269804  # MeV · fm   (ℏc; CODATA 2018)
# Useful derived quantity: (ℏc)² in MeV²·fm²
hbar_c_sq: float = hbar_c ** 2  # MeV²·fm²

# ---------------------------------------------------------------------------
# Energy unit conversions  (internal energy unit: MeV)
# ---------------------------------------------------------------------------
keV_per_MeV: float = 1.0e3   # 1 MeV = 1 000 keV
eV_per_MeV: float  = 1.0e6   # 1 MeV = 1 000 000 eV

# ---------------------------------------------------------------------------
# Length unit conversions
# ---------------------------------------------------------------------------
fm_per_cm: float  = 1.0e13   # 1 cm = 10¹³ fm  (used for σ unit conversion via (ℏc)²)
cm_per_fm: float  = 1.0e-13  # 1 fm = 10⁻¹³ cm
cm_per_m: float   = 1.0e2    # 1 m  = 100 cm   (baseline L: m → cm at API boundary)

# ---------------------------------------------------------------------------
# Atomic mass unit
# ---------------------------------------------------------------------------
u_to_MeV: float = 931.494    # 1 u = 931.494 MeV/c²  (CODATA 2018)
# Nuclear mass approximation: M ≈ A × u_to_MeV  (< 0.1 % error for A > 10)

# ---------------------------------------------------------------------------
# Reactor / flux conversions
# ---------------------------------------------------------------------------
MeV_per_fission: float = 205.0       # average energy released per fission [MeV]
                                      # (Huber & Mueller consensus value)

# 1 GW = 1×10⁹ J/s; 1 MeV = 1.602176634×10⁻¹³ J
GW_to_MeV_per_s: float = 1.0e9 / 1.602176634e-13  # ≈ 6.2415×10²¹ MeV·s⁻¹·GW⁻¹
