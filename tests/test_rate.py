"""Tests for final differential and binned event-rate calculations."""

import numpy as np
import pytest

from nurecoil import constants
from nurecoil.cross_section.base import CrossSectionBase
from nurecoil.cross_section.sm_cevns import SMCEvNS
from nurecoil.cross_section.sm_eves import SMEvES
from nurecoil.detector.quenching import ConstantQuenching, LindhardQuenching
from nurecoil.detector.resolution import ConstantResolution
from nurecoil.flux.base import FluxBase
from nurecoil.nucleus import Nucleus
from nurecoil.rate import (
    compute_binned_spectrum,
    compute_cevns_spectrum,
    compute_differential_rate,
    compute_eves_spectrum,
    compute_rate,
    target_count_from_mass,
)


GE76 = Nucleus(Z=32, A=76)


class ConstantFlux(FluxBase):
    E_min = 0.0
    E_max = 10.0

    def __init__(self, value=1.0):
        self.value = float(value)

    def _flux(self, E_nu, P, L):
        return np.full_like(E_nu, self.value * P / L**2)


class ConstantCrossSection(CrossSectionBase):
    def __init__(self, nucleus, value=1.0, recoil_max=1.0):
        super().__init__(nucleus)
        self.value = float(value)
        self.recoil_max = float(recoil_max)

    def E_R_max(self, E_nu):
        return self.recoil_max

    def __call__(self, E_nu, E_R):
        E_R_arr = np.asarray(E_R, dtype=float)
        return np.where(
            (E_R_arr >= 0.0) & (E_R_arr <= self.recoil_max),
            self.value,
            0.0,
        )


def _simple_inputs():
    return {
        "flux": ConstantFlux(value=2.0),
        "cross_section": ConstantCrossSection(GE76, value=3.0, recoil_max=0.004),
        "quenching": ConstantQuenching(1.0),
        "resolution": ConstantResolution(1e-3),
        "nucleus": GE76,
        "P": 1.0,
        "L": 1.0,
        "E_nu_range": (0.0, 0.01),
    }


def test_target_count_from_mass_uses_atomic_mass_number():
    assert target_count_from_mass(76.0, GE76) == pytest.approx(constants.N_A)


def test_target_count_from_mass_rejects_non_positive_mass():
    with pytest.raises(ValueError, match="target_mass_g"):
        target_count_from_mass(0.0, GE76)


def test_compute_rate_remains_alias_for_differential_rate():
    inputs = _simple_inputs()
    E_det = np.array([0.001, 0.002])

    direct = compute_differential_rate(E_det, N_T=5.0, **inputs)
    compat = compute_rate(E_det, N_T=5.0, **inputs)

    np.testing.assert_allclose(compat, direct)


def test_binned_spectrum_validates_bin_edges():
    inputs = _simple_inputs()

    with pytest.raises(ValueError, match="one-dimensional"):
        compute_binned_spectrum([[0.0, 0.001], [0.002, 0.003]], N_T=1.0, **inputs)

    with pytest.raises(ValueError, match="strictly increasing"):
        compute_binned_spectrum([0.0, 0.002, 0.001], N_T=1.0, **inputs)


def test_binned_spectrum_shape_non_negative_and_scales_with_N_T():
    inputs = _simple_inputs()
    bin_edges = np.array([0.0, 0.001, 0.002])

    spectrum = compute_binned_spectrum(bin_edges, N_T=2.0, **inputs)
    doubled = compute_binned_spectrum(bin_edges, N_T=4.0, **inputs)

    assert spectrum.shape == (2,)
    assert np.all(np.isfinite(spectrum))
    assert np.all(spectrum >= 0.0)
    np.testing.assert_allclose(doubled, 2.0 * spectrum, rtol=1e-10, atol=0.0)


def test_binned_spectrum_exposure_scales_counts():
    inputs = _simple_inputs()
    bin_edges = np.array([0.0, 0.001, 0.002])

    one_second = compute_binned_spectrum(bin_edges, N_T=2.0, exposure_s=1.0, **inputs)
    two_seconds = compute_binned_spectrum(bin_edges, N_T=2.0, exposure_s=2.0, **inputs)

    np.testing.assert_allclose(two_seconds, 2.0 * one_second, rtol=1e-10, atol=0.0)


def test_cevns_wrapper_defaults_to_lindhard_quenching_and_accepts_target_mass():
    flux = ConstantFlux(value=1.0)
    bin_edges = np.array([0.0, 1e-5, 2e-5])

    spectrum = compute_cevns_spectrum(
        bin_edges,
        flux,
        nucleus=GE76,
        target_mass_g=76.0,
        P=1.0,
        L=1.0,
        E_nu_range=(2.0, 2.01),
        resolution=ConstantResolution(1e-5),
    )

    assert isinstance(compute_cevns_spectrum.default_quenching(), LindhardQuenching)
    assert spectrum.shape == (2,)
    assert np.all(np.isfinite(spectrum))
    assert np.all(spectrum >= 0.0)


def test_eves_wrapper_defaults_to_unit_quenching_and_uses_atom_count_not_electron_count():
    flux = ConstantFlux(value=1.0)
    bin_edges = np.array([0.0, 1e-3, 2e-3])
    common = {
        "flux": flux,
        "nucleus": GE76,
        "P": 1.0,
        "L": 1.0,
        "E_nu_range": (2.0, 2.01),
        "resolution": ConstantResolution(1e-3),
        "cross_section": SMEvES(GE76),
    }

    from_mass = compute_eves_spectrum(bin_edges, target_mass_g=76.0, **common)
    from_count = compute_eves_spectrum(bin_edges, N_T=constants.N_A, **common)

    assert isinstance(compute_eves_spectrum.default_quenching(), ConstantQuenching)
    np.testing.assert_allclose(from_mass, from_count, rtol=1e-12)
    assert from_mass.shape == (2,)
    assert np.all(np.isfinite(from_mass))
    assert np.all(from_mass >= 0.0)


def test_wrappers_can_build_default_cross_sections():
    flux = ConstantFlux(value=1.0)

    cevns = compute_cevns_spectrum(
        [0.0, 1e-5],
        flux,
        nucleus=GE76,
        N_T=1.0,
        P=1.0,
        L=1.0,
        E_nu_range=(2.0, 2.01),
        resolution=ConstantResolution(1e-5),
    )
    eves = compute_eves_spectrum(
        [0.0, 1e-3],
        flux,
        nucleus=GE76,
        N_T=1.0,
        P=1.0,
        L=1.0,
        E_nu_range=(2.0, 2.01),
        resolution=ConstantResolution(1e-3),
    )

    assert isinstance(SMCEvNS(GE76), SMCEvNS)
    assert isinstance(SMEvES(GE76), SMEvES)
    assert np.all(np.isfinite(cevns))
    assert np.all(np.isfinite(eves))
