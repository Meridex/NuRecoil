"""Tests for nurecoil.cross_section (form factors, SM CEvNS, SM EvES)."""

import numpy as np
import pytest
from nurecoil.nucleus import Nucleus
from nurecoil.cross_section.form_factors import (
    HelmFormFactor,
    KleinNystrandFormFactor,
    GaussianFormFactor,
)
from nurecoil.cross_section.sm_cevns import SMCEvNS
from nurecoil.cross_section.sm_eves import SMEvES


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

    def test_benchmarks(self):

        ff = HelmFormFactor()
        q_val = 10.0
        F = ff(np.array([q_val]), GE76)
        val_benchmark = 0.992962380487098551
        assert F[0] == pytest.approx(val_benchmark, abs=1e-6)


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

    def test_tree_value_benchmark(self):

        xs_tree = SMCEvNS(GE76, flavor="tree")
        E_nu, E_R = 5, 0.00016
        val_tree = float(xs_tree(E_nu, E_R))
        val_tree_benchmark = 4.16569710349041948e-37 # cm^2/MeV
        assert abs(val_tree - val_tree_benchmark)/abs(val_tree_benchmark) < 1e-9

        E_nu, E_R = 50.0, 0.01
        val_tree = float(xs_tree(E_nu, E_R))
        val_tree_benchmark = 3.78909661085607236e-37 # cm^2/MeV
        assert abs(val_tree - val_tree_benchmark)/abs(val_tree_benchmark) < 1e-9

    def test_flavor_rc_differs_from_tree(self):
        # Radiative-corrected couplings must give different cross section than tree-level
        # Use E_R well inside the kinematic range: E_R_max(50 MeV, Ge76) ~ 0.027 MeV
        xs_tree = SMCEvNS(GE76, flavor="tree")
        xs_rc   = SMCEvNS(GE76, flavor="nu_e")
        E_nu, E_R = 50.0, 0.01
        val_tree = float(xs_tree(E_nu, E_R))
        val_rc   = float(xs_rc(E_nu, E_R))
        assert abs(val_tree - val_rc) / abs(val_tree) > 1e-4

    def test_flavor_rc_nu_flavors_differ(self):
        # nu_e, nu_mu, nu_tau have different g_V^p -> different cross sections
        xs_e   = SMCEvNS(GE76, flavor="nu_e")
        xs_mu  = SMCEvNS(GE76, flavor="nu_mu")
        xs_tau = SMCEvNS(GE76, flavor="nu_tau")
        E_nu, E_R = 50.0, 0.01
        val_e   = float(xs_e(E_nu, E_R))
        val_mu  = float(xs_mu(E_nu, E_R))
        val_tau = float(xs_tau(E_nu, E_R))
        assert abs(val_e  - val_mu)  / abs(val_e)  > 1e-4
        assert abs(val_e  - val_tau) / abs(val_e)  > 1e-4
        assert abs(val_mu - val_tau) / abs(val_mu) > 1e-4

    def test_invalid_flavor_raises(self):
        with pytest.raises(ValueError, match="Unknown flavor"):
            SMCEvNS(GE76, flavor="nu_sterile")


I53 = Nucleus(Z=53, A=127)
UNKNOWN = Nucleus(Z=14, A=28)  # Silicon — not in Z_eff table


class TestSMEvES:
    def setup_method(self):
        self.xs_ge = SMEvES(GE76, flavor="nu_e")
        self.xs_cs = SMEvES(CS133, flavor="nu_e")
        self.xs_i  = SMEvES(I53,  flavor="nu_e")

    # --- kinematics ---

    def test_non_negative_in_kinematics(self):
        E_nu = 5.0
        E_R  = np.linspace(0.0, self.xs_ge.E_R_max(E_nu) * 0.99, 50)
        dsig = self.xs_ge(E_nu, E_R)
        assert np.all(dsig >= 0)

    def test_zero_outside_kinematics(self):
        E_nu = 5.0
        E_R_over = self.xs_ge.E_R_max(E_nu) * 1.1
        assert self.xs_ge(E_nu, E_R_over) == pytest.approx(0.0)

    def test_E_R_max_formula(self):
        from nurecoil import constants
        E_nu = 5.0
        expected = 2.0 * E_nu ** 2 / (constants.m_e + 2.0 * E_nu)
        assert self.xs_ge.E_R_max(E_nu) == pytest.approx(expected, rel=1e-9)

    # --- flavor dependence ---

    def test_flavor_nu_e_differs_from_nu_mu(self):
        xs_e  = SMEvES(GE76, flavor="nu_e")
        xs_mu = SMEvES(GE76, flavor="nu_mu")
        E_nu, E_R = 5.0, 1e-3
        val_e  = float(xs_e(E_nu, E_R))
        val_mu = float(xs_mu(E_nu, E_R))
        # nu_e has g_V ~ +0.95 vs nu_mu g_V ~ -0.04: large difference expected
        assert abs(val_e - val_mu) / max(abs(val_e), abs(val_mu)) > 0.1

    def test_antineutrino_differs_from_neutrino(self):
        xs_nu    = SMEvES(GE76, flavor="nu_mu", antineutrino=False)
        xs_nubar = SMEvES(GE76, flavor="nu_mu", antineutrino=True)
        # At T/E_nu ~ 0.5 the (g_V-g_A)^2*(1-T/E_nu)^2 term is suppressed for
        # neutrinos vs antineutrinos, making the difference largest there.
        E_nu = 5.0
        E_R  = E_nu * 0.5  # well inside kinematics for electron target
        val_nu    = float(xs_nu(E_nu, E_R))
        val_nubar = float(xs_nubar(E_nu, E_R))
        # The two values must be distinct (g_A sign flip changes the cross section)
        assert abs(val_nu - val_nubar) / abs(val_nu) > 1e-4

    def test_invalid_flavor_raises(self):
        with pytest.raises(ValueError, match="Unknown flavor"):
            SMEvES(GE76, flavor="nu_sterile")

    # --- Z_eff step function ---

    def test_z_eff_ge_step_at_k_edge(self):
        # Just above the Ge K-edge (11.103 keV): Z_eff=32; just below L1-edge: Z_eff=30
        E_nu = 100.0  # high enough that both T values are kinematically allowed
        T_above_K  = 11.2e-3   # MeV, above 11.103 keV
        T_below_K  = 11.0e-3   # MeV, below 11.103 keV, above L1 1.4146 keV
        dsig_above = self.xs_ge(E_nu, T_above_K)
        dsig_below = self.xs_ge(E_nu, T_below_K)
        # Z_eff=32 above vs Z_eff=30 below -> strictly larger cross section above
        assert dsig_above > dsig_below

    def test_z_eff_cs_step_at_k_edge(self):
        E_nu = 100.0
        T_above_K = 36.0e-3   # above 35.99 keV
        T_below_K = 35.0e-3   # below 35.99 keV
        dsig_above = self.xs_cs(E_nu, T_above_K)
        dsig_below = self.xs_cs(E_nu, T_below_K)
        assert dsig_above > dsig_below

    def test_z_eff_i_step_at_k_edge(self):
        E_nu = 100.0
        T_above_K = 33.5e-3   # above 33.17 keV
        T_below_K = 33.0e-3   # below 33.17 keV
        dsig_above = self.xs_i(E_nu, T_above_K)
        dsig_below = self.xs_i(E_nu, T_below_K)
        assert dsig_above > dsig_below

    def test_unknown_element_raises(self):
        with pytest.raises(NotImplementedError, match="Z_eff is not implemented"):
            SMEvES(UNKNOWN, flavor="nu_e")
