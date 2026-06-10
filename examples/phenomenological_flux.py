# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: NURECOIL
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Phenomenological Reactor Antineutrino Flux
#
# This notebook demonstrates `PhenomenologicalFlux` (alias `PhenoFlux`) from `nurecoil`.
#
# The differential flux is:
#
# $$
# \frac{d\Phi}{dE_\nu} = \frac{1}{4\pi L^{2}} \frac{P_{\rm th}}{\bar{e}} \sum_{i} \frac{f_i}{F} S_i(E_\nu)
# $$
#
# where $\bar{e} = \sum_i (f_i/F)\,e_i$ and $S_i$ is the per-fission spectrum for isotope $i$.

# %%
import numpy as np
import matplotlib.pyplot as plt

from nurecoil.flux.phenomenological import (
    PhenomenologicalFlux,
    SPECTRUM_MODELS,
    ENERGY_PER_FISSION_MODELS,
    FISSION_FRACTION_PRESETS,
    SPECTRUM_HUBER_MUELLER,
    FISSION_FRACTIONS_TYPICAL,
    _eval_spectrum,
)

plt.rcParams.update({"figure.dpi": 120, "font.size": 11})

E = np.linspace(2.0, 8.0, 500)   # MeV
P = 1.0   # GW
L = 10.0  # m

ISOTOPE_COLORS = {"U235": "C0", "U238": "C1", "Pu239": "C2", "Pu241": "C3"}
ISOTOPE_LABELS = {
    "U235":  r"$^{235}$U",
    "U238":  r"$^{238}$U",
    "Pu239": r"$^{239}$Pu",
    "Pu241": r"$^{241}$Pu",
}

flux_default = PhenomenologicalFlux()  # huber_mueller + typical + ma_2013

# %% [markdown]
# ## 1. Per-Isotope Spectra
#
# Raw spectrum shape $S_i(E_\nu)$ [# / MeV / fission] for each isotope using Huber-Mueller coefficients.

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))

for iso, coeffs in SPECTRUM_HUBER_MUELLER.items():
    ax.semilogy(E, _eval_spectrum(coeffs, E),
                label=ISOTOPE_LABELS[iso], color=ISOTOPE_COLORS[iso])

ax.set_xlabel(r"$E_\nu$ [MeV]")
ax.set_ylabel(r"$S_i(E_\nu)$ [# / MeV / fission]")
ax.set_title("Per-isotope antineutrino spectra (Huber-Mueller)")
ax.legend()
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 2. Weighted Combined Spectrum vs Individual Contributions
#
# Dashed lines show each isotope's weighted contribution $(f_i/F)\,S_i$; solid black is the total. Using typical PWR fission fractions.

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))

for iso, coeffs in SPECTRUM_HUBER_MUELLER.items():
    frac = FISSION_FRACTIONS_TYPICAL[iso]
    ax.semilogy(E, frac * _eval_spectrum(coeffs, E), "--", alpha=0.7,
                label=rf"{ISOTOPE_LABELS[iso]} (f={frac:.0%})",
                color=ISOTOPE_COLORS[iso])

total = flux_default._spectrum_per_fission(E)
ax.semilogy(E, total, "k-", lw=2, label="Total")

ax.set_xlabel(r"$E_\nu$ [MeV]")
ax.set_ylabel(r"$\sum_i (f_i/F)\,S_i$ [# / MeV / fission]")
ax.set_title("Combined spectrum: isotope contributions (Typical PWR)")
ax.legend(fontsize=9)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 3. Comparison of Spectrum Models
#
# `huber_mueller` (default), `mueller_2011`, `schreckenbach_1985`.  
# Note: `huber_2011` lacks U238 and cannot be used with standard fission fractions.

# %%
MODEL_STYLES = {
    "huber_mueller":      ("C0", "-",   "Huber-Mueller (default)"),
    "mueller_2011":       ("C1", "--",  "Mueller 2011"),
    "schreckenbach_1985": ("C2", "-.", "Schreckenbach 1985"),
}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

for name, (color, ls, label) in MODEL_STYLES.items():
    phi = PhenomenologicalFlux(spectrum_model=name)(E, P, L)
    axes[0].semilogy(E, phi, ls=ls, color=color, label=label)

axes[0].set_xlabel(r"$E_\nu$ [MeV]")
axes[0].set_ylabel(r"$d\Phi/dE_\nu$ [# / MeV / cm$^2$ / s]")
axes[0].set_title("Spectrum model comparison")
axes[0].legend(fontsize=9)
axes[0].grid(True, which="both", alpha=0.3)

phi_ref = PhenomenologicalFlux(spectrum_model="huber_mueller")(E, P, L)
for name, (color, ls, label) in MODEL_STYLES.items():
    if name == "huber_mueller":
        continue
    phi = PhenomenologicalFlux(spectrum_model=name)(E, P, L)
    axes[1].plot(E, phi / phi_ref, ls=ls, color=color, label=label)

axes[1].axhline(1.0, color="C0", ls="-", label="Huber-Mueller (ref)")
axes[1].set_xlabel(r"$E_\nu$ [MeV]")
axes[1].set_ylabel("Ratio to Huber-Mueller")
axes[1].set_title("Spectrum model ratio")
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 4. Comparison of Fission Fraction Presets

# %%
PRESET_STYLES = {
    "typical":  ("C0", "-",   "Typical PWR"),
    "ksnps":    ("C1", "--",  "KSNPS"),
    "conus":    ("C2", "-.", "CONUS"),
    "daya_bay": ("C3", ":",   "Daya Bay"),
}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

for preset, (color, ls, label) in PRESET_STYLES.items():
    phi = PhenomenologicalFlux(fission_fractions=preset)(E, P, L)
    axes[0].semilogy(E, phi, ls=ls, color=color, label=label)

axes[0].set_xlabel(r"$E_\nu$ [MeV]")
axes[0].set_ylabel(r"$d\Phi/dE_\nu$ [# / MeV / cm$^2$ / s]")
axes[0].set_title("Fission fraction preset comparison")
axes[0].legend(fontsize=9)
axes[0].grid(True, which="both", alpha=0.3)

phi_ref = PhenomenologicalFlux(fission_fractions="typical")(E, P, L)
for preset, (color, ls, label) in PRESET_STYLES.items():
    if preset == "typical":
        continue
    phi = PhenomenologicalFlux(fission_fractions=preset)(E, P, L)
    axes[1].plot(E, phi / phi_ref, ls=ls, color=color, label=label)

axes[1].axhline(1.0, color="C0", ls="-", label="Typical (ref)")
axes[1].set_xlabel(r"$E_\nu$ [MeV]")
axes[1].set_ylabel("Ratio to Typical")
axes[1].set_title("Fission fraction ratio")
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 5. Effect of Energy-per-Fission Model on Normalization
#
# $\bar{e}$ enters only the overall normalization $P_{\rm th}/\bar{e}$, so different $e_i$ models produce a uniform ratio shift across all energies.

# %%
EPF_STYLES = {
    "ma_2013":         ("C0", "-",   "Ma 2013 (default)"),
    "bemporad_2002":   ("C1", "--",  "Bemporad 2002"),
    "meulenberg_1969": ("C2", "-.", "Meulenberg 1969"),
}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

for model, (color, ls, label) in EPF_STYLES.items():
    f = PhenomenologicalFlux(energy_per_fission=model)
    phi = f(E, P, L)
    axes[0].semilogy(E, phi, ls=ls, color=color,
                     label=rf"{label}  ($\bar{{e}}={f._e_bar:.2f}$ MeV)")

axes[0].set_xlabel(r"$E_\nu$ [MeV]")
axes[0].set_ylabel(r"$d\Phi/dE_\nu$ [# / MeV / cm$^2$ / s]")
axes[0].set_title(r"Effect of $e_i$ model on normalization")
axes[0].legend(fontsize=9)
axes[0].grid(True, which="both", alpha=0.3)

phi_ref = PhenomenologicalFlux(energy_per_fission="ma_2013")(E, P, L)
for model, (color, ls, label) in EPF_STYLES.items():
    if model == "ma_2013":
        continue
    phi = PhenomenologicalFlux(energy_per_fission=model)(E, P, L)
    axes[1].plot(E, phi / phi_ref, ls=ls, color=color, label=label)

axes[1].axhline(1.0, color="C0", ls="-", label="Ma 2013 (ref)")
axes[1].set_xlabel(r"$E_\nu$ [MeV]")
axes[1].set_ylabel(r"Ratio to Ma 2013")
axes[1].set_title(r"Normalization ratio ($\bar{e}$ effect)")
axes[1].legend(fontsize=9)
axes[1].grid(True, alpha=0.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 6. Distance and Power Scaling

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

for p_gw in [0.5, 1.0, 2.0, 5.0]:
    axes[0].semilogy(E, flux_default(E, p_gw, L), label=rf"$P={p_gw}$ GW")

axes[0].set_xlabel(r"$E_\nu$ [MeV]")
axes[0].set_ylabel(r"$d\Phi/dE_\nu$ [# / MeV / cm$^2$ / s]")
axes[0].set_title("Flux vs reactor power")
axes[0].legend()
axes[0].grid(True, which="both", alpha=0.3)

for l_m in [10, 25, 50, 100]:
    axes[1].semilogy(E, flux_default(E, P, l_m), label=rf"$L={l_m}$ m")

axes[1].set_xlabel(r"$E_\nu$ [MeV]")
axes[1].set_ylabel(r"$d\Phi/dE_\nu$ [# / MeV / cm$^2$ / s]")
axes[1].set_title(r"Flux vs baseline $L$")
axes[1].legend()
axes[1].grid(True, which="both", alpha=0.3)
fig.tight_layout()
plt.show()
