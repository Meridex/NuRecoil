"""Tests for nurecoil.flux implementations."""

import numpy as np
import pytest
from nurecoil.flux.phenomenological import (
    PhenomenologicalFlux,
    SPECTRUM_MODELS,
    SPECTRUM_HUBER_MUELLER,
)
from nurecoil.flux.reactor_mix import (
    ENERGY_PER_FISSION_MODELS,
    FISSION_FRACTION_PRESETS,
    ENERGY_PER_FISSION_MA_2013,
    FISSION_FRACTIONS_TYPICAL,
)
from nurecoil.flux.interpolated import InterpolatedFlux, IsotopeInterpolatedFlux
from nurecoil.flux import PhenoFlux


E_NU = np.array([2.0, 4.0, 6.0])  # MeV, all within the valid 1.8–8.0 MeV range
P    = 1.0   # GW
L    = 10.0  # m


class TestPhenomenologicalFlux:
    def setup_method(self):
        self.flux = PhenomenologicalFlux()

    # --- basic output sanity ---

    def test_positive(self):
        phi = self.flux(E_NU, P, L)
        assert np.all(phi > 0)

    def test_shape(self):
        phi = self.flux(E_NU, P, L)
        assert phi.shape == E_NU.shape

    def test_scalar_energy(self):
        phi = self.flux(np.array([3.0]), P, L)
        assert phi.shape == (1,)
        assert phi[0] > 0

    # --- physical scaling ---

    def test_scales_with_power(self):
        phi1 = self.flux(E_NU, 1.0, L)
        phi2 = self.flux(E_NU, 2.0, L)
        np.testing.assert_allclose(phi2, 2.0 * phi1, rtol=1e-9)

    def test_scales_with_distance(self):
        phi1 = self.flux(E_NU, P, 10.0)
        phi2 = self.flux(E_NU, P, 20.0)
        np.testing.assert_allclose(phi2, phi1 / 4.0, rtol=1e-9)

    def test_zero_outside_energy_range(self):
        phi = self.flux(np.array([0.5, 10.0]), P, L)
        assert np.all(phi == 0.0)

    # --- preset / model selection ---

    @pytest.mark.parametrize("preset", list(FISSION_FRACTION_PRESETS))
    def test_fission_fraction_presets(self, preset):
        flux = PhenomenologicalFlux(fission_fractions=preset)
        phi = flux(E_NU, P, L)
        assert np.all(phi > 0)

    @pytest.mark.parametrize("model", [m for m in SPECTRUM_MODELS if m != "huber_2011"])
    def test_spectrum_models(self, model):
        flux = PhenomenologicalFlux(spectrum_model=model)
        phi = flux(E_NU, P, L)
        assert np.all(phi > 0)

    @pytest.mark.parametrize("model", list(ENERGY_PER_FISSION_MODELS))
    def test_energy_per_fission_models(self, model):
        flux = PhenomenologicalFlux(energy_per_fission=model)
        phi = flux(E_NU, P, L)
        assert np.all(phi > 0)

    def test_huber_2011_no_u238_ok(self):
        # huber_2011 is fine when U238 fraction is zero
        flux = PhenomenologicalFlux(
            fission_fractions={"U235": 0.57, "U238": 0.0, "Pu239": 0.30, "Pu241": 0.13},
            spectrum_model="huber_2011",
        )
        phi = flux(E_NU, P, L)
        assert np.all(phi > 0)

    def test_huber_2011_with_u238_raises(self):
        with pytest.raises(ValueError, match="U238"):
            PhenomenologicalFlux(
                fission_fractions={"U235": 0.57, "U238": 0.07, "Pu239": 0.30, "Pu241": 0.06},
                spectrum_model="huber_2011",
            )

    def test_custom_fission_fractions(self):
        fracs = {"U235": 0.564, "U238": 0.076, "Pu239": 0.304, "Pu241": 0.056}
        flux = PhenomenologicalFlux(fission_fractions=fracs)
        phi = flux(E_NU, P, L)
        assert np.all(phi > 0)

    def test_custom_spectrum_coeffs(self):
        flux = PhenomenologicalFlux(spectrum_model=SPECTRUM_HUBER_MUELLER)
        phi = flux(E_NU, P, L)
        assert np.all(phi > 0)

    def test_custom_energy_per_fission(self):
        flux = PhenomenologicalFlux(energy_per_fission=ENERGY_PER_FISSION_MA_2013)
        phi = flux(E_NU, P, L)
        assert np.all(phi > 0)

    # --- e_bar sanity: weighted average should be in a physical range ---

    def test_e_bar_physical_range(self):
        flux = PhenomenologicalFlux()
        assert 200.0 < flux._e_bar < 215.0

    # --- validation errors ---

    def test_invalid_fractions_sum(self):
        with pytest.raises(ValueError, match="sum to 1"):
            PhenomenologicalFlux(fission_fractions={"U235": 0.5, "Pu239": 0.1})

    def test_invalid_fraction_preset_name(self):
        with pytest.raises(ValueError, match="Unknown fission fraction preset"):
            PhenomenologicalFlux(fission_fractions="nonexistent")

    def test_invalid_spectrum_model_name(self):
        with pytest.raises(ValueError, match="Unknown spectrum model"):
            PhenomenologicalFlux(spectrum_model="nonexistent")

    def test_invalid_energy_per_fission_name(self):
        with pytest.raises(ValueError, match="Unknown energy-per-fission model"):
            PhenomenologicalFlux(energy_per_fission="nonexistent")

    # --- alias ---

    def test_phenoflux_alias(self):
        assert PhenoFlux is PhenomenologicalFlux


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


# ---------------------------------------------------------------------------
# Helpers shared by IsotopeInterpolatedFlux tests
# ---------------------------------------------------------------------------

def _flat_spectra(isotopes=("U235", "U238", "Pu239", "Pu241"), E_lo=2.0, E_hi=8.0, val=1.0):
    """Return a spectra dict with flat unit spectra for the given isotopes."""
    E = np.linspace(E_lo, E_hi, 100)
    S = np.full_like(E, val)
    return {iso: (E.copy(), S.copy()) for iso in isotopes}


class TestIsotopeInterpolatedFlux:
    def setup_method(self):
        self.spectra = _flat_spectra()
        self.flux = IsotopeInterpolatedFlux(self.spectra)

    # --- basic sanity ---

    def test_positive_inside_range(self):
        phi = self.flux(np.array([3.0, 5.0]), P, L)
        assert np.all(phi > 0)

    def test_zero_outside_range(self):
        phi = self.flux(np.array([0.5, 9.0]), P, L)
        assert np.all(phi == 0.0)

    def test_scales_with_power(self):
        phi1 = self.flux(np.array([4.0]), 1.0, L)
        phi2 = self.flux(np.array([4.0]), 2.0, L)
        assert np.isclose(phi2 / phi1, 2.0, rtol=1e-6)

    def test_scales_with_distance(self):
        phi1 = self.flux(np.array([4.0]), P, 10.0)
        phi2 = self.flux(np.array([4.0]), P, 20.0)
        assert np.isclose(phi2 / phi1, 0.25, rtol=1e-6)

    # --- fission fraction presets ---

    @pytest.mark.parametrize("preset", ["typical", "ksnps", "conus", "daya_bay"])
    def test_all_presets(self, preset):
        f = IsotopeInterpolatedFlux(self.spectra, fission_fractions=preset)
        phi = f(np.array([4.0]), P, L)
        assert phi[0] > 0

    def test_custom_fission_fractions(self):
        fracs = {"U235": 0.6, "U238": 0.1, "Pu239": 0.2, "Pu241": 0.1}
        f = IsotopeInterpolatedFlux(self.spectra, fission_fractions=fracs)
        phi = f(np.array([4.0]), P, L)
        assert phi[0] > 0

    # --- energy per fission presets ---

    @pytest.mark.parametrize("epf", list(ENERGY_PER_FISSION_MODELS))
    def test_all_epf_models(self, epf):
        f = IsotopeInterpolatedFlux(self.spectra, energy_per_fission=epf)
        phi = f(np.array([4.0]), P, L)
        assert phi[0] > 0

    # --- partial isotope set ---

    def test_subset_of_isotopes(self):
        # Only U235 and Pu239 provided, zero U238/Pu241 fraction
        spectra = _flat_spectra(("U235", "Pu239"))
        fracs = {"U235": 0.6, "U238": 0.0, "Pu239": 0.4, "Pu241": 0.0}
        f = IsotopeInterpolatedFlux(spectra, fission_fractions=fracs)
        phi = f(np.array([4.0]), P, L)
        assert phi[0] > 0

    # --- E_min / E_max union ---

    def test_energy_range_is_union(self):
        spectra = {
            "U235":  (np.linspace(1.0, 6.0, 50), np.ones(50)),
            "Pu239": (np.linspace(2.0, 8.0, 50), np.ones(50)),
        }
        fracs = {"U235": 0.6, "U238": 0.0, "Pu239": 0.4, "Pu241": 0.0}
        f = IsotopeInterpolatedFlux(spectra, fission_fractions=fracs)
        assert f.E_min == pytest.approx(1.0)
        assert f.E_max == pytest.approx(8.0)

    # --- error paths ---

    def test_invalid_isotope_key(self):
        bad = {"Xe135": (np.linspace(2, 8, 10), np.ones(10))}
        with pytest.raises(ValueError, match="Unknown actinide"):
            IsotopeInterpolatedFlux(bad)

    def test_empty_spectra(self):
        with pytest.raises(ValueError, match="at least one"):
            IsotopeInterpolatedFlux({})

    def test_fractions_not_sum_to_one(self):
        fracs = {"U235": 0.5, "U238": 0.1, "Pu239": 0.1, "Pu241": 0.1}
        with pytest.raises(ValueError, match="sum to 1"):
            IsotopeInterpolatedFlux(self.spectra, fission_fractions=fracs)

    def test_invalid_fraction_preset(self):
        with pytest.raises(ValueError, match="Unknown fission fraction preset"):
            IsotopeInterpolatedFlux(self.spectra, fission_fractions="nonexistent")

    def test_invalid_epf_preset(self):
        with pytest.raises(ValueError, match="Unknown energy-per-fission model"):
            IsotopeInterpolatedFlux(self.spectra, energy_per_fission="nonexistent")

    def test_unsorted_table(self):
        bad = {"U235": ([3.0, 1.0, 2.0], [1.0, 1.0, 1.0])}
        fracs = {"U235": 1.0, "U238": 0.0, "Pu239": 0.0, "Pu241": 0.0}
        with pytest.raises(ValueError, match="strictly increasing"):
            IsotopeInterpolatedFlux(bad, fission_fractions=fracs)
