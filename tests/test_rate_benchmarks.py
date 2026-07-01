"""External benchmark comparisons for CEvNS/EvES binned rate spectra."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from nurecoil.cross_section.form_factors import HelmFormFactor
from nurecoil.cross_section.sm_cevns import SMCEvNS
from nurecoil.cross_section.sm_eves import SMEvES
from nurecoil.detector.quenching import ConstantQuenching, LindhardQuenching
from nurecoil.detector.resolution import CDEXResolution, ConstantResolution
from nurecoil.flux.interpolated import make_flux as make_interpolated_flux
from nurecoil.flux.phenomenological import PhenomenologicalFlux
from nurecoil.nucleus import Nucleus
from nurecoil.rate import compute_cevns_spectrum, compute_eves_spectrum


BENCHMARK_DIR = Path(__file__).parent / "data" / "rate_benchmarks"
BENCHMARK_FILES = sorted(BENCHMARK_DIR.glob("*.json"))
TEMPLATE_FILES = sorted(BENCHMARK_DIR.glob("*.json.template"))


def _load_case(path: Path) -> dict:
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def _make_flux(config: dict):
    model = config.get("model")
    if model == "phenomenological":
        return PhenomenologicalFlux(
            fission_fractions=config.get("fission_fractions"),
            spectrum_model=config.get("spectrum_model", "huber_mueller"),
            energy_per_fission=config.get("energy_per_fission", "ma_2013"),
        )
    if model == "interpolated":
        source = config.get("source")
        if not source:
            raise ValueError("interpolated flux requires a 'source' field.")
        return make_interpolated_flux(
            source=source,
            fission_fractions=config.get("fission_fractions"),
            energy_per_fission=config.get("energy_per_fission", "ma_2013"),
            extrapolate=config.get("extrapolate", False),
        )
    if model == "tabulated":
        from nurecoil.flux.tabulated import TabulatedFlux
        return TabulatedFlux(source=config.get("source"))
    if model == "conflux":
        from nurecoil.flux.conflux import ConfluxFlux
        return ConfluxFlux()
    raise ValueError(f"Unsupported flux model: {model!r}")


def _make_resolution(config: dict):
    model = config.get("model")
    if model == "cdex":
        return CDEXResolution()
    if model == "constant":
        return ConstantResolution(config["sigma_MeV"])
    raise ValueError(f"Unsupported resolution model: {model!r}")


def _make_quenching(config: dict):
    model = config.get("model")
    if model == "lindhard":
        return LindhardQuenching(k=config.get("k"))
    if model == "constant":
        return ConstantQuenching(value=config.get("value", 1.0))
    raise ValueError(f"Unsupported quenching model: {model!r}")


def _make_cross_section(observable: str, nucleus: Nucleus, config: dict):
    model = config.get("model")
    if observable == "cevns" and model == "sm_cevns":
        form_factor_name = config.get("form_factor", "helm")
        if form_factor_name != "helm":
            raise ValueError("Only Helm CEvNS form factor benchmarks are supported.")
        return SMCEvNS(
            nucleus,
            form_factor=HelmFormFactor(),
            flavor=config.get("flavor", "tree"),
        )
    if observable == "eves" and model == "sm_eves":
        return SMEvES(
            nucleus,
            flavor=config.get("flavor", "nu_e"),
            antineutrino=config.get("antineutrino", False),
        )
    raise ValueError(f"Unsupported cross-section model for {observable}: {model!r}")


def _reference_arrays(config: dict) -> tuple[np.ndarray, np.ndarray]:
    bins = config["reference"]["bins"]
    if not bins:
        raise ValueError("reference.bins must not be empty.")
    edges = [float(bins[0]["low_MeV"])]
    values = []
    for item in bins:
        if float(item["low_MeV"]) != pytest.approx(edges[-1]):
            raise ValueError("reference.bins must be contiguous and ordered.")
        edges.append(float(item["high_MeV"]))
        values.append(float(item["value"]))
    return np.asarray(edges, dtype=float), np.asarray(values, dtype=float)


@pytest.mark.parametrize("template_path", TEMPLATE_FILES, ids=lambda p: p.name)
def test_rate_benchmark_templates_are_valid(template_path: Path):
    case = _load_case(template_path)
    assert case["observable"] in {"cevns", "eves"}
    assert {"Z", "A"} <= set(case["target"])
    assert case["normalization"]["P_GW"] > 0.0
    assert case["normalization"]["L_m"] > 0.0
    assert case["reference"]["unit"] in {"counts_per_s_per_bin", "counts_per_bin"}
    bin_edges, values = _reference_arrays(case)
    assert len(bin_edges) == len(values) + 1


@pytest.mark.skipif(
    not BENCHMARK_FILES,
    reason=f"No external rate benchmark JSON files found in {BENCHMARK_DIR}.",
)
@pytest.mark.parametrize("benchmark_path", BENCHMARK_FILES, ids=lambda p: p.stem)
def test_rate_spectrum_matches_external_benchmark(benchmark_path: Path):
    case = _load_case(benchmark_path)
    observable = case["observable"]
    nucleus = Nucleus(**case["target"])
    normalization = case["normalization"]
    detector = case["detector"]
    reference = case["reference"]
    bin_edges, expected = _reference_arrays(case)

    flux = _make_flux(case["flux"])
    cross_section = _make_cross_section(observable, nucleus, case["cross_section"])
    resolution = _make_resolution(detector["resolution"])
    quenching = _make_quenching(detector["quenching"])
    E_nu_range = tuple(case["integration"]["E_nu_range_MeV"])

    common = {
        "flux": flux,
        "nucleus": nucleus,
        "P": normalization["P_GW"],
        "L": normalization["L_m"],
        "target_mass_g": normalization.get("target_mass_g"),
        "N_T": normalization.get("N_T"),
        "cross_section": cross_section,
        "quenching": quenching,
        "resolution": resolution,
        "E_nu_range": E_nu_range,
    }

    unit = reference["unit"]
    if unit == "counts_per_bin":
        if normalization.get("exposure_s") is None:
            raise ValueError("counts_per_bin benchmarks must set exposure_s.")
        common["exposure_s"] = normalization["exposure_s"]
    elif unit != "counts_per_s_per_bin":
        raise ValueError(f"Unsupported reference unit: {unit!r}")

    if observable == "cevns":
        actual = compute_cevns_spectrum(bin_edges, **common)
    elif observable == "eves":
        actual = compute_eves_spectrum(bin_edges, **common)
    else:
        raise ValueError(f"Unsupported observable: {observable!r}")

    np.testing.assert_allclose(
        actual,
        expected,
        rtol=reference.get("rtol", 0.02),
        atol=reference.get("atol", 0.0),
        err_msg=f"Benchmark mismatch: {benchmark_path.name}",
    )
