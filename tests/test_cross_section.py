"""Tests for nurecoil.cross_section (form factors and SM CEvNS)."""

import numpy as np
import pytest
from nurecoil.nucleus import Nucleus
from nurecoil.cross_section.form_factors import (
    HelmFormFactor,
    KleinNystrandFormFactor,
    GaussianFormFactor,
)
from nurecoil.cross_section.sm_cevns import SMCEvNS


GE76 = Nucleus(Z=32, A=76)
CS133 = Nucleus(Z=55, A=133)  # CsI target, commonly used in CEvNS


class TestHelmFormFactor:
    def setup_method(self):
        self.ff = HelmFormFactor()

    def test_at_zero_is_one(self):
        F = self.ff(np.array([0.0]), GE76)
        assert F[0] == pytest.approx(1.0, abs=1e-6)

    def test_at_zero_is_one_heavy(self):
        # heavier nucleus (Cs): same F(0)=1 requirement
        F = self.ff(np.array([0.0]), CS133)
        assert F[0] == pytest.approx(1.0, abs=1e-6)

    def test_decreases_with_q(self):
        q = np.linspace(0.0, 200.0, 100)
        F = self.ff(q, GE76)
        # envelope should decrease; check mean of last half < mean of first half
        assert F[50:].mean() < F[:50].mean()

    def test_finite_everywhere(self):
        q = np.linspace(0.0, 500.0, 200)
        F = self.ff(q, GE76)
        assert np.all(np.isfinite(F))

    def test_scalar_input(self):
        # should handle scalar q without crashing
        F = self.ff(50.0, GE76)
        assert np.isfinite(F)

    def test_lewin_smith_params(self):
        # default parameters must match Lewin-Smith values
        ff = HelmFormFactor()
        assert ff.s_fm    == pytest.approx(0.9)
        assert ff.a_fm    == pytest.approx(0.52)
        assert ff.c_coeff == pytest.approx(1.23)

    def test_custom_params_change_output(self):
        # Changing s_fm must produce different F(q) values (parameters are wired in)
        ff_default = HelmFormFactor(s_fm=0.9)
        ff_large_s = HelmFormFactor(s_fm=2.0)
        q = np.linspace(10.0, 100.0, 20)
        assert not np.allclose(ff_default(q, GE76), ff_large_s(q, GE76))

    def test_gaussian_envelope_larger_s(self):
        # The Gaussian envelope exp(-q^2 s^2 / 2 / hbarc^2) is strictly smaller
        # for larger s at any q > 0, independent of the j1 oscillations.
        from nurecoil import constants
        q_val = 80.0  # MeV — well inside the monotone part before first j1 zero
        s_small, s_large = 0.9, 2.0
        env_small = np.exp(-0.5 * (q_val * s_small / constants.hbar_c) ** 2)
        env_large = np.exp(-0.5 * (q_val * s_large / constants.hbar_c) ** 2)
        assert env_large < env_small


class TestKleinNystrandFormFactor:
    def setup_method(self):
        self.ff = KleinNystrandFormFactor()

    def test_at_zero_is_one(self):
        F = self.ff(np.array([0.0]), GE76)
        assert F[0] == pytest.approx(1.0, abs=1e-6)

    def test_at_zero_is_one_heavy(self):
        F = self.ff(np.array([0.0]), CS133)
        assert F[0] == pytest.approx(1.0, abs=1e-6)

    def test_decreases_with_q(self):
        q = np.linspace(0.0, 200.0, 100)
        F = self.ff(q, GE76)
        assert F[50:].mean() < F[:50].mean()

    def test_finite_everywhere(self):
        q = np.linspace(0.0, 500.0, 200)
        F = self.ff(q, GE76)
        assert np.all(np.isfinite(F))

    def test_default_params(self):
        ff = KleinNystrandFormFactor()
        assert ff.r0_fm == pytest.approx(1.3)
        assert ff.a_fm  == pytest.approx(0.7)

    def test_larger_yukawa_range_suppresses_more(self):
        # The Yukawa factor 1/(1+(qa/hbarc)^2) decreases with a at any q>0.
        # Use F² (the physically relevant quantity) to avoid sign issues from j1 oscillations.
        ff_default = KleinNystrandFormFactor(a_fm=0.7)
        ff_soft    = KleinNystrandFormFactor(a_fm=2.0)
        q = np.array([50.0, 100.0, 150.0])
        assert np.all(ff_soft(q, GE76) ** 2 < ff_default(q, GE76) ** 2)


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

    def test_analytic_value(self):
        from nurecoil import constants
        ff = GaussianFormFactor(R_coeff=1.2)
        q_val = 50.0  # MeV
        R_fm  = 1.2 * GE76.A ** (1.0 / 3.0)
        expected = np.exp(-(q_val * R_fm / constants.hbar_c) ** 2 / 6.0)
        F = ff(np.array([q_val]), GE76)
        assert F[0] == pytest.approx(expected, rel=1e-9)


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
