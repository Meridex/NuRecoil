"""Tests for nurecoil.flux implementations."""

import numpy as np
import pytest
from nurecoil.flux.data_loader import (
    load_spectrum,
    load_all_isotopes,
    list_sources,
    list_isotopes,
)
from nurecoil.flux.interpolated import make_flux
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


# ---------------------------------------------------------------------------
# TestDataLoader — unit tests for nurecoil.flux.data_loader
# ---------------------------------------------------------------------------

ISOSPEC_SOURCES = ["estienne2019", "mueller2011", "CEA2023", "vogel1989"]
COMPOSITE_SOURCES = ["kopeikin2012", "kopeikin1999"]
ALL_ISOTOPES = ["U235", "U238", "Pu239", "Pu241"]


class TestDataLoader:
    # --- list_sources / list_isotopes ---

    def test_list_sources_returns_list(self):
        sources = list_sources()
        assert isinstance(sources, list)
        assert len(sources) > 0

    def test_list_sources_contains_known(self):
        sources = list_sources()
        for s in ISOSPEC_SOURCES + COMPOSITE_SOURCES:
            assert s in sources, f"Expected source '{s}' in list_sources()"

    def test_list_isotopes_isospec(self):
        isos = list_isotopes("estienne2019")
        assert set(isos) == {"U235", "U238", "Pu239", "Pu241"}

    def test_list_isotopes_composite_returns_empty(self):
        isos = list_isotopes("kopeikin2012")
        assert isos == []

    # --- load_spectrum: isospec ---

    @pytest.mark.parametrize("source", ISOSPEC_SOURCES)
    @pytest.mark.parametrize("isotope", ALL_ISOTOPES)
    def test_load_spectrum_isospec_shape(self, source, isotope):
        E, S = load_spectrum(source, isotope)
        assert E.ndim == 1 and S.ndim == 1
        assert len(E) == len(S)
        assert len(E) > 1

    @pytest.mark.parametrize("source", ISOSPEC_SOURCES)
    @pytest.mark.parametrize("isotope", ALL_ISOTOPES)
    def test_load_spectrum_isospec_strictly_increasing(self, source, isotope):
        E, _ = load_spectrum(source, isotope)
        assert np.all(np.diff(E) > 0)

    @pytest.mark.parametrize("source", ISOSPEC_SOURCES)
    @pytest.mark.parametrize("isotope", ALL_ISOTOPES)
    def test_load_spectrum_isospec_positive(self, source, isotope):
        _, S = load_spectrum(source, isotope)
        assert np.all(S >= 0)

    # --- load_spectrum: isotope name aliases ---

    @pytest.mark.parametrize("alias", ["u235", "U235", "u-235", "U-235"])
    def test_load_spectrum_isotope_aliases(self, alias):
        E1, S1 = load_spectrum("estienne2019", "u235")
        E2, S2 = load_spectrum("estienne2019", alias)
        np.testing.assert_array_equal(E1, E2)
        np.testing.assert_array_equal(S1, S2)

    # --- load_spectrum: composite ---

    @pytest.mark.parametrize("source", COMPOSITE_SOURCES)
    def test_load_spectrum_composite(self, source):
        E, S = load_spectrum(source)
        assert E.ndim == 1 and S.ndim == 1
        assert len(E) > 1
        assert np.all(np.diff(E) > 0)
        assert np.all(S >= 0)

    # --- load_spectrum: error paths ---

    def test_load_spectrum_unknown_isotope(self):
        with pytest.raises(ValueError, match="Unknown isotope"):
            load_spectrum("estienne2019", "Xe135")

    def test_load_spectrum_unknown_source(self):
        with pytest.raises(ValueError):
            load_spectrum("nonexistent_source", "U235")

    def test_load_spectrum_composite_missing(self):
        # estienne2019 has no composite file
        with pytest.raises(ValueError, match="No composite file"):
            load_spectrum("estienne2019")

    # --- load_all_isotopes ---

    @pytest.mark.parametrize("source", ISOSPEC_SOURCES)
    def test_load_all_isotopes_keys(self, source):
        spectra = load_all_isotopes(source)
        assert set(spectra.keys()) == {"U235", "U238", "Pu239", "Pu241"}

    @pytest.mark.parametrize("source", ISOSPEC_SOURCES)
    def test_load_all_isotopes_shapes(self, source):
        spectra = load_all_isotopes(source)
        for iso, (E, S) in spectra.items():
            assert len(E) == len(S), f"{source}/{iso}: length mismatch"
            assert np.all(np.diff(E) > 0), f"{source}/{iso}: not strictly increasing"

    def test_load_all_isotopes_unknown_source(self):
        with pytest.raises(ValueError, match="No isospec files found"):
            load_all_isotopes("nonexistent_source")


# ---------------------------------------------------------------------------
# TestFromSource — integration tests for from_source classmethods
# ---------------------------------------------------------------------------

class TestFromSource:
    # --- InterpolatedFlux.from_composite ---

    @pytest.mark.parametrize("source", COMPOSITE_SOURCES)
    def test_interpolated_from_composite(self, source):
        from nurecoil.flux.interpolated import InterpolatedFlux
        flux = InterpolatedFlux.from_composite(source)
        E_mid = (flux.E_min + flux.E_max) / 2.0
        phi = flux(np.array([E_mid]), P=1.0, L=500.0)
        assert phi[0] > 0

    # --- IsotopeInterpolatedFlux.from_source ---

    @pytest.mark.parametrize("source", ISOSPEC_SOURCES)
    def test_isotope_interpolated_from_source(self, source):
        from nurecoil.flux.interpolated import IsotopeInterpolatedFlux
        flux = IsotopeInterpolatedFlux.from_source(source)
        E_mid = (flux.E_min + flux.E_max) / 2.0
        phi = flux(np.array([E_mid]), P=1.0, L=500.0)
        assert phi[0] > 0

    @pytest.mark.parametrize("preset", ["typical", "ksnps", "conus", "daya_bay"])
    def test_isotope_from_source_presets(self, preset):
        from nurecoil.flux.interpolated import IsotopeInterpolatedFlux
        flux = IsotopeInterpolatedFlux.from_source("estienne2019", fission_fractions=preset)
        E_mid = (flux.E_min + flux.E_max) / 2.0
        phi = flux(np.array([E_mid]), P=1.0, L=500.0)
        assert phi[0] > 0

    def test_isotope_from_source_unknown(self):
        from nurecoil.flux.interpolated import IsotopeInterpolatedFlux
        with pytest.raises(ValueError):
            IsotopeInterpolatedFlux.from_source("nonexistent")


# ---------------------------------------------------------------------------
# TestMakeFlux — unified constructor
# ---------------------------------------------------------------------------

class TestMakeFlux:
    @pytest.mark.parametrize("source", ISOSPEC_SOURCES)
    def test_isospec_returns_isotope_class(self, source):
        from nurecoil.flux.interpolated import IsotopeInterpolatedFlux
        flux = make_flux(source)
        assert isinstance(flux, IsotopeInterpolatedFlux)

    @pytest.mark.parametrize("source", COMPOSITE_SOURCES)
    def test_composite_returns_interpolated_class(self, source):
        from nurecoil.flux.interpolated import InterpolatedFlux
        flux = make_flux(source)
        assert isinstance(flux, InterpolatedFlux)

    @pytest.mark.parametrize("source", ISOSPEC_SOURCES + COMPOSITE_SOURCES)
    def test_make_flux_produces_positive_output(self, source):
        flux = make_flux(source)
        E_mid = np.array([(flux.E_min + flux.E_max) / 2.0])
        phi = flux(E_mid, P=1.0, L=500.0)
        assert phi[0] > 0

    def test_make_flux_fission_fraction_preset(self):
        flux = make_flux("estienne2019", fission_fractions="daya_bay")
        E_mid = np.array([4.0])
        phi = flux(E_mid, P=1.0, L=500.0)
        assert phi[0] > 0

    def test_make_flux_unknown_source_raises(self):
        with pytest.raises(ValueError):
            make_flux("nonexistent")


# ---------------------------------------------------------------------------
# TestIsotopeSpectrum — per-isotope access on IsotopeInterpolatedFlux
# ---------------------------------------------------------------------------

class TestIsotopeSpectrum:
    def setup_method(self):
        self.flux = IsotopeInterpolatedFlux.from_source("estienne2019")
        self.E = np.linspace(2.0, 8.0, 100)

    def test_isotopes_property(self):
        assert set(self.flux.isotopes) == {"U235", "U238", "Pu239", "Pu241"}

    @pytest.mark.parametrize("iso", ["U235", "U238", "Pu239", "Pu241"])
    def test_isotope_spectrum_positive(self, iso):
        S = self.flux.isotope_spectrum(iso, self.E)
        assert np.all(S >= 0)
        assert np.any(S > 0)

    @pytest.mark.parametrize("iso", ["U235", "U238", "Pu239", "Pu241"])
    def test_isotope_spectrum_zero_outside_range(self, iso):
        e_lo, e_hi = self.flux._e_ranges[iso]
        E_out = np.array([e_lo - 1.0, e_hi + 1.0])
        S = self.flux.isotope_spectrum(iso, E_out)
        assert np.all(S == 0.0)

    def test_isotope_spectrum_unknown_raises(self):
        with pytest.raises(ValueError, match="not in this instance"):
            self.flux.isotope_spectrum("Xe135", self.E)

    @pytest.mark.parametrize("iso", ["U235", "U238", "Pu239", "Pu241"])
    def test_isotope_flux_positive(self, iso):
        phi_iso = self.flux.isotope_flux(iso, self.E, P=1.0, L=500.0)
        assert np.all(phi_iso >= 0)
        assert np.any(phi_iso > 0)

    def test_isotope_flux_sums_to_total(self):
        # Sum of per-isotope fluxes should equal the total flux
        phi_total = self.flux(self.E, P=1.0, L=500.0)
        phi_sum = sum(
            self.flux.isotope_flux(iso, self.E, P=1.0, L=500.0)
            for iso in self.flux.isotopes
        )
        np.testing.assert_allclose(phi_sum, phi_total, rtol=1e-10)

    def test_isotope_flux_unknown_raises(self):
        with pytest.raises(ValueError, match="not in this instance"):
            self.flux.isotope_flux("Xe135", self.E, P=1.0, L=500.0)


# ---------------------------------------------------------------------------
# TestPhenoVsInterpolated — cross-class consistency checks
# ---------------------------------------------------------------------------

class TestPhenoVsInterpolated:
    """
    Sanity checks comparing PhenomenologicalFlux and IsotopeInterpolatedFlux.

    Both classes share the same fission-fraction weighting and e_bar
    bookkeeping.  Their absolute values differ (different spectral models),
    but they must agree on basic physical properties.
    """

    # Common evaluation grid — within the valid range of all models
    E = np.linspace(2.5, 7.5, 200)
    P = 1.0   # GW
    L = 10.0  # m

    def setup_method(self):
        from nurecoil.flux.phenomenological import SPECTRUM_MODELS
        self.pheno_models = list(SPECTRUM_MODELS)
        self.isospec_sources = [s for s in list_sources() if list_isotopes(s)]

    # --- both classes produce positive flux on the shared grid ---------------

    @pytest.mark.parametrize("model", ["huber_mueller", "mueller_2011", "vogel_1985"])
    def test_pheno_positive_on_shared_grid(self, model):
        from nurecoil.flux import PhenomenologicalFlux
        flux = PhenomenologicalFlux(spectrum_model=model)
        phi = flux(self.E, self.P, self.L)
        assert np.all(phi > 0)

    @pytest.mark.parametrize("source", ["estienne2019", "mueller2011", "CEA2023"])
    def test_interp_positive_on_shared_grid(self, source):
        flux = make_flux(source)
        phi = flux(self.E, self.P, self.L)
        assert np.all(phi > 0)

    # --- same fission fractions → same e_bar → same normalisation factor -----

    def test_same_ebar_same_preset(self):
        """Both classes built with the same preset must have equal e_bar."""
        from nurecoil.flux import PhenomenologicalFlux
        pheno = PhenomenologicalFlux(fission_fractions="typical")
        interp = make_flux("estienne2019", fission_fractions="typical")
        assert pheno._e_bar == pytest.approx(interp._e_bar, rel=1e-10)

    @pytest.mark.parametrize("preset", ["typical", "ksnps", "conus", "daya_bay"])
    def test_ebar_consistent_across_presets(self, preset):
        """e_bar must match between pheno and interpolated for every preset."""
        from nurecoil.flux import PhenomenologicalFlux
        pheno  = PhenomenologicalFlux(fission_fractions=preset)
        interp = make_flux("estienne2019", fission_fractions=preset)
        assert pheno._e_bar == pytest.approx(interp._e_bar, rel=1e-10)

    # --- ratio between pheno and interpolated is within physically plausible range

    @pytest.mark.parametrize("model", ["huber_mueller", "mueller_2011"])
    @pytest.mark.parametrize("source", ["estienne2019", "mueller2011"])
    def test_ratio_within_factor_of_two(self, model, source):
        """Pheno / interpolated ratio should stay within [0.5, 2.0] on 2.5–7.5 MeV."""
        from nurecoil.flux import PhenomenologicalFlux
        pheno  = PhenomenologicalFlux(spectrum_model=model)
        interp = make_flux(source)
        phi_p = pheno(self.E, self.P, self.L)
        phi_i = interp(self.E, self.P, self.L)
        ratio = phi_p / np.where(phi_i > 0, phi_i, np.nan)
        assert np.nanmin(ratio) > 0.5
        assert np.nanmax(ratio) < 2.0

    # --- changing preset shifts both classes in the same direction -----------

    def test_preset_shift_direction_consistent(self):
        """
        Switching from 'typical' to 'daya_bay' should change e_bar by the same
        amount in both pheno and interpolated (they share ReactorMix).
        """
        from nurecoil.flux import PhenomenologicalFlux
        p_typ = PhenomenologicalFlux(fission_fractions="typical")
        p_db  = PhenomenologicalFlux(fission_fractions="daya_bay")
        i_typ = make_flux("estienne2019", fission_fractions="typical")
        i_db  = make_flux("estienne2019", fission_fractions="daya_bay")
        assert (p_db._e_bar - p_typ._e_bar) == pytest.approx(
            i_db._e_bar - i_typ._e_bar, rel=1e-10
        )
