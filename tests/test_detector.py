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
