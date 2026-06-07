"""Tests for nurecoil.constants."""

import math
import nurecoil.constants as C


def test_G_F_unit_conversion():
    """G_F [MeV⁻²] = G_F_GeV2 × 1e-6."""
    assert math.isclose(C.G_F, C.G_F_GeV2 * 1e-6, rel_tol=1e-9)


def test_hbar_c_value():
    """ℏc ≈ 197.327 MeV·fm."""
    assert math.isclose(C.hbar_c, 197.3269804, rel_tol=1e-6)


def test_length_conversions():
    assert math.isclose(C.fm_per_cm * C.cm_per_fm, 1.0, rel_tol=1e-12)
    assert math.isclose(C.cm_per_m, 100.0, rel_tol=1e-12)


def test_energy_conversions():
    assert math.isclose(C.keV_per_MeV, 1e3, rel_tol=1e-12)
    assert math.isclose(C.eV_per_MeV,  1e6, rel_tol=1e-12)


def test_GW_to_MeV_per_s():
    """1 GW = 1e9 J/s; 1 MeV = 1.602176634e-13 J."""
    expected = 1.0e9 / 1.602176634e-13
    assert math.isclose(C.GW_to_MeV_per_s, expected, rel_tol=1e-9)
