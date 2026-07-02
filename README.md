# NuRecoil

Python package for computing reactor-neutrino event rates, including
**CEvNS** (Coherent Elastic Neutrino-Nucleus Scattering) and
**EvES** (Elastic Neutrino-Electron Scattering),
targeting experiments such as CDEX and TEXONO.

## Features

- **Flux models** — phenomenological (Huber-Mueller etc.) and numerically
  interpolated tabulated spectra; unified factory via `load_flux()`
- **Cross sections** — SM CEvNS (`SMCEvNS`) and SM EvES (`SMEvES`)
- **Nuclear form factors** — Helm and Gaussian
- **Detector response** — Lindhard quenching factor, CDEX energy resolution
- **Rate integrators** — differential rate, binned spectrum, exposure-scaled counts
- **Bundled datasets** — six reactor antineutrino spectral datasets
  (Estienne 2019, Mueller 2011, Vogel 1989, CEA 2023, Kopeikin 1999/2012)

## Requirements

- Python ≥ 3.10
- `numpy >= 1.24`, `scipy >= 1.11`

## Install

```bash
# Editable install (recommended for development)
pip install -e .

# With dev tools (pytest, matplotlib)
pip install -e ".[dev]"

# With CONFLUX flux backend
pip install -e ".[conflux]"
```

> **Conda users** — create and activate the project environment first:
> ```bash
> conda env create -f environment.yml   # installs everything, including nurecoil
> conda activate NURECOIL
> ```

## Quick Start

### CEvNS event rate

```python
import numpy as np
from nurecoil import Nucleus, compute_rate
from nurecoil.flux import load_flux
from nurecoil.cross_section import SMCEvNS, HelmFormFactor
from nurecoil.detector import LindhardQuenching, CDEXResolution

nucleus = Nucleus(Z=32, A=76)                        # ⁷⁶Ge
flux    = load_flux("pheno:huber_mueller")            # default reactor flux
xs      = SMCEvNS(nucleus, HelmFormFactor())
quench  = LindhardQuenching()
res     = CDEXResolution()

E_det = np.linspace(1e-4, 5e-3, 200)                # detected energy grid [MeV]

rate = compute_rate(
    E_det, flux, xs, quench, res,
    nucleus=nucleus,
    N_T=1e27,   # number of target atoms
    P=3.0,      # thermal power [GW]
    L=30.0,     # baseline [m]
)
# rate shape: (200,)  units: counts / s / MeV
```

### Binned spectrum (with exposure)

```python
from nurecoil import compute_binned_spectrum

bin_edges = np.linspace(1e-4, 5e-3, 51)   # 50 bins [MeV]

counts = compute_binned_spectrum(
    bin_edges, flux, xs, quench, res,
    nucleus=nucleus,
    N_T=1e27,
    P=3.0, L=30.0,
    exposure_s=3.15e7,   # 1 year in seconds
)
# counts shape: (50,)  units: counts / bin
```

### Flux selection

```python
# Phenomenological spectrum models
flux = load_flux("pheno:huber_mueller")   # default — Huber+Mueller (2–8 MeV)
flux = load_flux("pheno:mueller_2011")
flux = load_flux("pheno:vogel_1985")

# Interpolated tabulated datasets
flux = load_flux("interp:estienne2019")   # 0.05–10.05 MeV
flux = load_flux("interp:CEA2023")        # 0.01–12.5 MeV
flux = load_flux("interp:kopeikin2012")   # 0.01–9 MeV (composite)

# Custom fission fractions / energy per fission
flux = load_flux("pheno:huber_mueller",
                 fission_fractions="daya_bay",
                 energy_per_fission="ma_2013")

# User-supplied CSV file
flux = load_flux("file:/path/to/spectrum.csv")
```

See `help(nurecoil.flux)` for the full reference of all `spec` values and kwargs.

## Module Overview

| Module | Contents |
|--------|----------|
| `nurecoil` | `Nucleus`, `compute_rate`, `compute_binned_spectrum`, `compute_differential_rate`, `compute_cevns_spectrum`, `compute_eves_spectrum`, `target_count_from_mass`, `constants` |
| `nurecoil.flux` | `load_flux`, `PhenomenologicalFlux`, `InterpolatedFlux`, `IsotopeInterpolatedFlux`, `TabulatedFlux`, `ReactorMix`, data loaders |
| `nurecoil.cross_section` | `SMCEvNS`, `SMEvES`, `HelmFormFactor`, `GaussianFormFactor` |
| `nurecoil.detector` | `LindhardQuenching`, `ConstantQuenching`, `CDEXResolution`, `ConstantResolution` |

## Units Convention

| Quantity | User-facing input | Internal |
|----------|-------------------|----------|
| Energy (neutrino / recoil) | MeV | MeV |
| Detected energy | MeV | MeV |
| Baseline | m | cm |
| Thermal power | GW | MeV/s |
| Cross section | — | cm² |
| Flux | — | # / MeV / cm² / s |
| Rate | — | counts / s / MeV |

## Examples

The `examples/` directory contains Jupytext-paired notebooks (tracked as `.py`
files; `.ipynb` generated locally):

| Script | Description |
|--------|-------------|
| `sm_cevns.py` | SM CEvNS cross section, event rate, binned spectrum |
| `sm_eves.py` | SM EvES cross section, CEvNS vs EvES comparison |
| `form_factors.py` | Helm and Gaussian form factor comparison |
| `phenomenological_flux.py` | Phenomenological flux model comparison |
| `interpolated_flux.py` | Interpolated tabulated flux datasets |
| `model_comparison.py` | Phenomenological vs. interpolated flux side-by-side |

Restore a notebook from its `.py` file:

```bash
jupytext --sync examples/<name>.py
```

## Run Tests

```bash
# With conda environment
conda run -n NURECOIL python -m pytest tests/ -v

# Or with activated environment
pytest tests/ -v
pytest tests/ -v --cov=nurecoil
```
