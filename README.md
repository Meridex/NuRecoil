# NuRecoil

Python package for computing CEvNS (Coherent Elastic Neutrino-Nucleus Scattering)
event rates for reactor-neutrino experiments such as CDEX and TEXONO.

## Install

```bash
# Standard install
pip install -e .

# With dev tools (pytest, matplotlib)
pip install -e ".[dev]"

# With CONFLUX flux backend
pip install -e ".[conflux]"
```

Requires Python ≥ 3.10.

## Quick start

```python
from nurecoil import Nucleus, compute_rate, constants
from nurecoil.flux import PhenomenologicalFlux
from nurecoil.cross_section import SMCEvNS, HelmFormFactor
from nurecoil.detector.quenching import LindhardQuenching
from nurecoil.detector.resolution import CDEXResolution
import numpy as np

nucleus = Nucleus(Z=32, A=76)          # ⁷⁶Ge
flux    = PhenomenologicalFlux()
xs      = SMCEvNS(nucleus, HelmFormFactor())
quench  = LindhardQuenching()
res     = CDEXResolution()

E_det   = np.linspace(1e-4, 5e-3, 100)  # MeV
N_T     = 1e27                           # number of target atoms

rate = compute_rate(
    E_det, flux, xs, quench, res,
    nucleus=nucleus, N_T=N_T,
    P=3.0,   # GW
    L=30.0,  # m
)
```

## Units convention

See `misc/plans/units_convention.md` for the full description.
All internal calculations use MeV (energy), fm/cm (length), and cm² (cross section).
User-facing inputs use the experimental convention: keV or MeV (energy), m (baseline), GW (power).

## Run tests

```bash
pytest --cov=nurecoil
```
