"""Tests for nurecoil.nucleus."""

import pytest
from nurecoil.nucleus import Nucleus, ISOTOPE_TABLE, get_nucleus
import nurecoil.constants as C


def test_nucleus_properties():
    ge76 = Nucleus(Z=32, A=76)
    assert ge76.N == 44
    assert abs(ge76.M - 76 * C.u_to_MeV) < 1e-6


def test_nucleus_invalid():
    with pytest.raises(ValueError):
        Nucleus(Z=0, A=10)
    with pytest.raises(ValueError):
        Nucleus(Z=10, A=5)


def test_isotope_table_populated():
    assert "Ge76" in ISOTOPE_TABLE
    assert "Si28" in ISOTOPE_TABLE


def test_get_nucleus_found():
    n = get_nucleus("Ge76")
    assert n.Z == 32 and n.A == 76


def test_get_nucleus_not_found():
    with pytest.raises(KeyError, match="Unknown nucleus"):
        get_nucleus("Xx999")
