"""Tests for nurecoil.detector (quenching and resolution)."""

import numpy as np
import pytest
from nurecoil.nucleus import Nucleus
from nurecoil.detector.quenching import LindhardQuenching, ConstantQuenching
from nurecoil.detector.resolution import CDEXResolution, ConstantResolution
import nurecoil.constants as C


# --- Quenching ---

class TestLindhardQuenching:
    def setup_method(self):
        self.q = LindhardQuenching()
        self.ge = Nucleus(Z=32, A=76)

    def test_range(self):
        E_R = np.linspace(1e-4, 0.1, 50)  # MeV
        fq = self.q(E_R, self.ge)
        assert np.all(fq >= 0) and np.all(fq <= 1)

    def test_below_threshold_is_zero(self):
        fq = self.q(np.array([0.0]), self.ge)
        assert fq[0] == 0.0

    def test_increases_with_energy(self):
        E_R = np.array([0.001, 0.01, 0.1])  # MeV
        fq = self.q(E_R, self.ge)
        assert np.all(np.diff(fq) > 0)

    def test_default_k_semi_empirical(self):
        # Default k = 0.133 * Z^(2/3) * A^(-1/2)
        # Note: doc quotes k~0.157 for Ge using natural-abundance A~72.6;
        # Ge76 (A=76) gives k~0.1538. Both are in the expected range.
        q = LindhardQuenching()
        expected_k = 0.133 * 32 ** (2.0 / 3.0) * 76 ** (-0.5)
        assert expected_k == pytest.approx(0.1538, abs=0.001)
        # k=None should give same result as passing the computed k explicitly
        q_explicit = LindhardQuenching(k=expected_k)
        E_R = np.array([0.001, 0.01, 0.05])
        np.testing.assert_allclose(q(E_R, self.ge), q_explicit(E_R, self.ge), rtol=1e-9)

    def test_custom_k_changes_output(self):
        # Bonhomme et al. 2022: k=0.162 for Ge
        q_default = LindhardQuenching()
        q_bonhomme = LindhardQuenching(k=0.162)
        E_R = np.array([0.005, 0.01, 0.05])
        fq_default  = q_default(E_R, self.ge)
        fq_bonhomme = q_bonhomme(E_R, self.ge)
        # k=0.162 > 0.157 (default for Ge76) -> higher quenching
        assert np.all(fq_bonhomme > fq_default)

    def test_analytic_value_at_1keV(self):
        # Manual calculation at E_R = 1 keV = 1e-3 MeV for Ge76
        E_R_keV = 1.0
        Z, A = 32, 76
        epsilon = 11.5 * E_R_keV * Z ** (-7.0 / 3.0)
        g = 3.0 * epsilon ** 0.15 + 0.7 * epsilon ** 0.6 + epsilon
        k = 0.133 * Z ** (2.0 / 3.0) * A ** (-0.5)
        expected = k * g / (1.0 + k * g)
        fq = LindhardQuenching()(np.array([1e-3]), self.ge)
        assert fq[0] == pytest.approx(expected, rel=1e-9)


class TestConstantQuenching:
    def test_default_unity(self):
        q = ConstantQuenching()
        fq = q(np.array([0.1, 0.5, 1.0]), Nucleus(Z=14, A=28))
        np.testing.assert_array_equal(fq, 1.0)

    def test_custom_value(self):
        q = ConstantQuenching(value=0.3)
        fq = q(np.array([0.1]), Nucleus(Z=14, A=28))
        assert fq[0] == pytest.approx(0.3)


# --- Resolution ---

class TestCDEXResolution:
    def setup_method(self):
        self.res = CDEXResolution()

    def test_sigma_units_MeV(self):
        # At E_ee = 1 keV = 1e-3 MeV:
        # ΔE [eV] = 35.8 + 16.6×√1 = 52.4 eV  → 52.4e-6 MeV
        E_ee = 1e-3  # MeV = 1 keV
        sigma = self.res.sigma(np.array([E_ee]))[0]
        expected = (35.8 + 16.6 * 1.0) / C.eV_per_MeV
        assert sigma == pytest.approx(expected, rel=1e-6)

    def test_smear_normalisation(self):
        """Integral of f_res over E_det grid ≈ 1."""
        E_ee  = np.array([0.01])  # MeV
        E_det = np.linspace(0.0, 0.05, 5000)
        kernel = self.res.smear(E_ee, E_det)  # (1, 5000) MeV⁻¹
        integral = np.trapezoid(kernel[0], E_det)
        assert integral == pytest.approx(1.0, abs=1e-3)


class TestConstantResolution:
    def test_sigma_constant(self):
        res = ConstantResolution(sigma_MeV=0.001)
        sigma = res.sigma(np.array([0.01, 0.1, 1.0]))
        np.testing.assert_array_almost_equal(sigma, 0.001)
