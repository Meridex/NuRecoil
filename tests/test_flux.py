"""Tests for nurecoil.flux implementations."""

import numpy as np
import pytest
from nurecoil.flux.phenomenological import PhenomenologicalFlux
from nurecoil.flux.interpolated import InterpolatedFlux


E_NU = np.array([2.0, 4.0, 6.0])  # MeV
P    = 1.0   # GW
L    = 10.0  # m


class TestPhenomenologicalFlux:
    def setup_method(self):
        self.flux = PhenomenologicalFlux()

    def test_positive(self):
        phi = self.flux(E_NU, P, L)
        assert np.all(phi > 0)

    def test_shape(self):
        phi = self.flux(E_NU, P, L)
        assert phi.shape == E_NU.shape

    def test_scales_with_power(self):
        phi1 = self.flux(E_NU, 1.0, L)
        phi2 = self.flux(E_NU, 2.0, L)
        np.testing.assert_allclose(phi2, 2.0 * phi1, rtol=1e-9)

    def test_scales_with_distance(self):
        phi1 = self.flux(E_NU, P, 10.0)
        phi2 = self.flux(E_NU, P, 20.0)
        # 1/L² scaling: φ(2L) = φ(L)/4
        np.testing.assert_allclose(phi2, phi1 / 4.0, rtol=1e-9)

    def test_invalid_fractions(self):
        with pytest.raises(ValueError, match="sum to 1"):
            PhenomenologicalFlux({"U235": 0.5, "Pu239": 0.1})


class TestInterpolatedFlux:
    def setup_method(self):
        E = np.linspace(1.8, 8.0, 200)
        # Flat spectrum for easy verification
        S = np.ones_like(E) * 1.0  # # / MeV / fission
        self.flux = InterpolatedFlux(E, S)

    def test_positive_inside_range(self):
        phi = self.flux(np.array([3.0, 5.0]), P, L)
        assert np.all(phi > 0)

    def test_zero_outside_range(self):
        phi = self.flux(np.array([0.5, 9.0]), P, L)
        assert np.all(phi == 0.0)

    def test_invalid_table_unsorted(self):
        with pytest.raises(ValueError, match="strictly increasing"):
            InterpolatedFlux([3.0, 1.0, 2.0], [1.0, 1.0, 1.0])
