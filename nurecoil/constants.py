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
# PDG: $G_F / (\hbar c)^3 = 1.1663788(6) \times 10^{-5}\ \mathrm{GeV}^{-2}$
G_F_GeV2: float = 1.1663788e-5   # GeV^{-2}  -- PDG reference value; do NOT use in formulas
G_F: float      = 1.1663788e-11  # MeV^{-2}  -- use this in all cross-section formulas
                                  # (= G_F_GeV2 * 1e-6, since 1 GeV^{-2} = 1e-6 MeV^{-2})

# ---------------------------------------------------------------------------
# Electroweak mixing
# ---------------------------------------------------------------------------
# Low-q² (on-shell / running) value used in CEvNS and EvES cross-section formulas.
# Sources: arXiv:2411.03122; Phys. Rev. D 110 (2024) 030001.
# Note: the MS-bar value at M_Z is 0.23122 (PDG 2024) — do NOT use that here.
sin2_theta_W: float = 0.23868  # $\sin^2\theta_W$ at $q^2 \to 0$

# ---------------------------------------------------------------------------
# Natural unit bridge
# ---------------------------------------------------------------------------
hbar_c: float = 197.3269804  # MeV * fm   ($\hbar c$; CODATA 2018)
# Useful derived quantity: $(\hbar c)^2$ in MeV^2 * fm^2
hbar_c_sq: float = hbar_c ** 2  # MeV^2 * fm^2

# ---------------------------------------------------------------------------
# Energy unit conversions  (internal energy unit: MeV)
# ---------------------------------------------------------------------------
keV_per_MeV: float = 1.0e3   # 1 MeV = 1 000 keV
eV_per_MeV: float  = 1.0e6   # 1 MeV = 1 000 000 eV

# ---------------------------------------------------------------------------
# Length unit conversions
# ---------------------------------------------------------------------------
fm_per_cm: float  = 1.0e13   # 1 cm = 1e13 fm   (used for $\sigma$ conversion via $(\hbar c)^2$)
cm_per_fm: float  = 1.0e-13  # 1 fm = 1e-13 cm
cm_per_m: float   = 1.0e2    # 1 m  = 100 cm    (baseline $L$: m -> cm at API boundary)

# ---------------------------------------------------------------------------
# Atomic mass unit
# ---------------------------------------------------------------------------
u_to_MeV: float = 931.494    # 1 u = 931.494 MeV/$c^2$  (CODATA 2018)

# ---------------------------------------------------------------------------
# Electron mass
# ---------------------------------------------------------------------------
m_e: float = 0.51099895  # MeV/$c^2$  (CODATA 2018; PDG 2024)
# Nuclear mass approximation: $M \approx A \times u\_to\_MeV$  (< 0.1 % error for A > 10)

# ---------------------------------------------------------------------------
# Reactor / flux conversions
# ---------------------------------------------------------------------------
MeV_per_fission: float = 205.0       # average energy released per fission [MeV]
                                      # (Huber & Mueller consensus value)

# 1 GW = 1e9 J/s;  1 MeV = 1.602176634e-13 J
GW_to_MeV_per_s: float = 1.0e9 / 1.602176634e-13  # ~6.2415e21  MeV * s^{-1} * GW^{-1}
