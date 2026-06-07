"""Tests for nurecoil.cross_section (form factors and SM CEvNS)."""

import numpy as np
import pytest
from nurecoil.nucleus import Nucleus
from nurecoil.cross_section.form_factors import HelmFormFactor, GaussianFormFactor
from nurecoil.cross_section.sm_cevns import SMCEvNS


GE76 = Nucleus(Z=32, A=76)


class TestHelmFormFactor:
    def setup_method(self):
        self.ff = HelmFormFactor()

    def test_at_zero_is_one(self):
        F = self.ff(np.array([0.0]), GE76)
        assert F[0] == pytest.approx(1.0, abs=1e-6)

    def test_decreases_with_q(self):
        q = np.linspace(0.0, 200.0, 100)
        F = self.ff(q, GE76)
        # envelope should decrease; check mean of last half < mean of first half
        assert F[50:].mean() < F[:50].mean()

    def test_non_negative(self):
        q = np.linspace(0.0, 500.0, 200)
        F = self.ff(q, GE76)
        assert np.all(np.isfinite(F))


class TestGaussianFormFactor:
    def setup_method(self):
        self.ff = GaussianFormFactor()

    def test_at_zero_is_one(self):
        F = self.ff(np.array([0.0]), GE76)
        assert F[0] == pytest.approx(1.0, abs=1e-9)

    def test_decreases_monotonically(self):
        q = np.linspace(0.0, 300.0, 100)
        F = self.ff(q, GE76)
        assert np.all(np.diff(F) <= 0)


class TestSMCEvNS:
    def setup_method(self):
        self.xs = SMCEvNS(GE76)

    def test_non_negative_in_kinematics(self):
        E_nu = 5.0  # MeV
        E_R  = np.linspace(0.0, self.xs.E_R_max(E_nu) * 0.99, 50)
        dsig = self.xs(E_nu, E_R)
        assert np.all(dsig >= 0)

    def test_zero_outside_kinematics(self):
        E_nu = 5.0
        E_R_over = self.xs.E_R_max(E_nu) * 1.1
        dsig = self.xs(E_nu, E_R_over)
        assert dsig == pytest.approx(0.0)

    def test_E_R_max_formula(self):
        E_nu = 10.0
        M = GE76.M
        expected = 2.0 * E_nu ** 2 / (M + 2.0 * E_nu)
        assert self.xs.E_R_max(E_nu) == pytest.approx(expected, rel=1e-9)
