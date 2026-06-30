# Examples

All notebooks in this folder are tracked as **Jupytext percent-format `.py` files**.
The paired `.ipynb` files are generated locally and are listed in `.gitignore`.

## Workflow

Only `.py` files are committed to git. To work with a notebook:

```bash
# Restore the .ipynb from the .py (first time on a new machine)
jupytext --sync examples/<name>.py

# Or open the .py directly in JupyterLab/VS Code — Jupytext keeps both in sync on save.
```

To create a new example notebook, create the `.py` file with `# %%` cell markers
and the standard Jupytext header (copy from an existing file), then run:

```bash
jupytext --sync examples/<name>.py
```

Commit only the `.py` file.

## Notebooks

| File | Description |
|------|-------------|
| [sm_cevns.py](sm_cevns.py) | SM CEvNS differential cross section, radiative corrections / flavor comparison, total cross section, event rate, and binned detected-energy spectrum |
| [sm_eves.py](sm_eves.py) | SM EvES (elastic neutrino-electron scattering) differential cross section, flavor and antineutrino effects, $Z_{\rm eff}$ atomic ionization steps, CEvNS vs EvES comparison, and binned detected-energy spectrum |
| [form_factors.py](form_factors.py) | Nuclear form factor models (Helm, Klein-Nystrand, Gaussian) — comparison and parameter sensitivity |
| [phenomenological_flux.py](phenomenological_flux.py) | Phenomenological reactor antineutrino flux models |
| [interpolated_flux.py](interpolated_flux.py) | Interpolated (tabulated) reactor antineutrino flux models |
| [model_comparison.py](model_comparison.py) | Side-by-side comparison of phenomenological vs. interpolated flux models |

## Kernel

All notebooks use the `NURECOIL` conda environment as their kernel.
Activate it before launching Jupyter:

```bash
conda activate NURECOIL
jupyter lab
```
